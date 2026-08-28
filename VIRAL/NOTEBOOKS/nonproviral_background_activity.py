"""Parallel PropagAtE-style activity calls with a non-proviral host background.

The expensive input is a set of sparse per-base depth files, one file per
sample/reference contig. Whole-host per-base coverage statistics are reconstructed from
the existing covstats tables.

For every tested interval, the background is the eligible bacterial genome
after removing the union of all explicit geNomad ``topology == 'Provirus'``
intervals.  If a tested interval extends outside that union, its additional
sequence is removed as well, so no tested base occurs in both foreground and
background.
"""

from __future__ import annotations

import argparse
import bisect
import math
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd


SAMPLES_36 = [
    "vFL", "vOL", "vF23",
    "vN1", "vN2", "vN3",
    "vZ1", "vZ2", "vZ3",
    "vF1A", "vF2A", "vF3A",
    "vF1B", "vF2B", "vF3B",
    "vF1C", "vF2C", "vF3C",
    "vO1A", "vO2A", "vO3A",
    "vO1B", "vO2B", "vO3B",
    "vO1C", "vO2C", "vO3C",
    "vS1B", "vS2B", "vS3B",
    "vS1C", "vS2C", "vS3C",
    "vNWd", "vSWd", "vSEd",
]

BACKGROUND_DEFINITION = (
    "eligible bacterial-genome sequence outside all explicit geNomad Provirus "
    "intervals; tested interval sequence outside those intervals was also excluded"
)


@dataclass(frozen=True)
class ActivityThresholds:
    min_mean_coverage: float = 1.0
    min_breadth: float = 0.50
    min_fold_change: float = 2.0
    min_cohen_d: float = 0.70
    contig_edge_mask_bp: int = 150


@dataclass
class ActivityRunResult:
    provirus_activity: pd.DataFrame
    provirus_summary: pd.DataFrame
    provirus_overview: pd.DataFrame
    votu_activity: pd.DataFrame
    interval_coverage_statistics: pd.DataFrame
    output_paths: dict[str, Path]

    @property
    def interval_moments(self) -> pd.DataFrame:
        """Backward-compatible alias for older notebook code."""
        return self.interval_coverage_statistics


_WORKER_CONTIG_SPECS: dict[str, dict[str, Any]] = {}
_WORKER_DEPTH_ROOT: Path | None = None


def _normalise_seq(values: pd.Series) -> pd.Series:
    return values.astype(str).str.strip().str.replace("-", "_", regex=False)


def _merge_intervals(intervals: Iterable[tuple[int, int]]) -> list[tuple[int, int]]:
    merged: list[list[int]] = []
    for start, end in sorted((int(a), int(b)) for a, b in intervals if int(b) >= int(a)):
        if not merged or start > merged[-1][1] + 1:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return [(start, end) for start, end in merged]


def _subtract_intervals(
    interval: tuple[int, int], masks: Iterable[tuple[int, int]]
) -> list[tuple[int, int]]:
    """Return inclusive parts of ``interval`` not covered by merged ``masks``."""
    start, end = interval
    if end < start:
        return []
    cursor = start
    pieces: list[tuple[int, int]] = []
    for mask_start, mask_end in _merge_intervals(masks):
        if mask_end < cursor:
            continue
        if mask_start > end:
            break
        if mask_start > cursor:
            pieces.append((cursor, min(end, mask_start - 1)))
        cursor = max(cursor, mask_end + 1)
        if cursor > end:
            break
    if cursor <= end:
        pieces.append((cursor, end))
    return pieces


def _interval_moments(
    positions: np.ndarray,
    depths: np.ndarray,
    intervals: Iterable[tuple[int, int]],
    min_coverage: float,
) -> tuple[int, float, float, int]:
    n_bases = 0
    depth_sum = 0.0
    depth_sum_squares = 0.0
    bases_at_min_coverage = 0
    for start, end in _merge_intervals(intervals):
        if end < start:
            continue
        left = int(np.searchsorted(positions, start, side="left"))
        right = int(np.searchsorted(positions, end, side="right"))
        interval_depths = depths[left:right]
        n_bases += end - start + 1
        depth_sum += float(np.sum(interval_depths))
        depth_sum_squares += float(np.dot(interval_depths, interval_depths))
        bases_at_min_coverage += int(np.count_nonzero(interval_depths >= min_coverage))
    return n_bases, depth_sum, depth_sum_squares, bases_at_min_coverage


def _read_sparse_depth(path: Path, valid_start: int, valid_end: int) -> tuple[np.ndarray, np.ndarray]:
    if not path.exists() or path.stat().st_size == 0 or valid_end < valid_start:
        return np.array([], dtype=np.int64), np.array([], dtype=float)

    # Read only the tested/masked intervals.  High-depth viromes can have tens
    # of megabytes per contig; streaming avoids loading irrelevant positions.
    selected_positions: list[int] = []
    selected_depths: list[float] = []
    with path.open() as handle:
        for line in handle:
            position_text, depth_text = line.split("\t", 1)
            position = int(position_text)
            if valid_start <= position <= valid_end:
                selected_positions.append(position)
                selected_depths.append(float(depth_text))
    positions = np.asarray(selected_positions, dtype=np.int64)
    depths = np.asarray(selected_depths, dtype=float)
    keep = (positions >= valid_start) & (positions <= valid_end)
    positions = positions[keep]
    depths = depths[keep]

    # The normal files contain one row per covered position.  Summing duplicates
    # makes the parser safe if depth chunks were ever concatenated.
    if positions.size > 1 and np.any(positions[1:] == positions[:-1]):
        unique_positions, inverse = np.unique(positions, return_inverse=True)
        unique_depths = np.zeros(unique_positions.size, dtype=float)
        np.add.at(unique_depths, inverse, depths)
        positions, depths = unique_positions, unique_depths
    return positions, depths


def _read_sparse_intervals(
    path: Path,
    intervals: list[tuple[int, int]],
    index_stride_bytes: int = 4 * 1024 * 1024,
) -> tuple[np.ndarray, np.ndarray]:
    """Read exact rows in inclusive intervals from a position-sorted text file.

    A small in-memory byte-offset index is constructed with one seek per 4 MiB.
    Queries then begin at the last checkpoint preceding each interval.  This
    avoids parsing the usually much larger parts of high-depth sparse files.
    """
    intervals = _merge_intervals(intervals)
    if not intervals or not path.exists() or path.stat().st_size == 0:
        return np.array([], dtype=np.int64), np.array([], dtype=float)

    file_size = path.stat().st_size
    checkpoint_positions = [0]
    checkpoint_offsets = [0]
    with path.open("rb") as handle:
        for approximate_offset in range(index_stride_bytes, file_size, index_stride_bytes):
            handle.seek(approximate_offset)
            handle.readline()  # discard the partial line at the seek position
            line_offset = handle.tell()
            line = handle.readline()
            if not line:
                break
            checkpoint_positions.append(int(line.split(b"\t", 1)[0]))
            checkpoint_offsets.append(line_offset)

        selected_positions: list[int] = []
        selected_depths: list[float] = []
        for start, end in intervals:
            checkpoint_index = max(
                0, bisect.bisect_right(checkpoint_positions, start) - 1
            )
            handle.seek(checkpoint_offsets[checkpoint_index])
            for line in handle:
                position_text, depth_text = line.split(b"\t", 1)
                position = int(position_text)
                if position < start:
                    continue
                if position > end:
                    break
                selected_positions.append(position)
                selected_depths.append(float(depth_text))

    positions = np.asarray(selected_positions, dtype=np.int64)
    depths = np.asarray(selected_depths, dtype=float)
    if positions.size > 1 and np.any(positions[1:] == positions[:-1]):
        unique_positions, inverse = np.unique(positions, return_inverse=True)
        unique_depths = np.zeros(unique_positions.size, dtype=float)
        np.add.at(unique_depths, inverse, depths)
        positions, depths = unique_positions, unique_depths
    return positions, depths


def _init_worker(contig_specs: dict[str, dict[str, Any]], depth_root: str) -> None:
    global _WORKER_CONTIG_SPECS, _WORKER_DEPTH_ROOT
    _WORKER_CONTIG_SPECS = contig_specs
    _WORKER_DEPTH_ROOT = Path(depth_root)


def _process_sample_depths(
    sample_and_settings: tuple[str, float, int]
) -> tuple[
    str, list[dict[str, Any]], list[dict[str, Any]]
]:
    sample, min_coverage, edge_mask_bp = sample_and_settings
    if _WORKER_DEPTH_ROOT is None:
        raise RuntimeError("Depth worker was not initialised.")

    region_rows: list[dict[str, Any]] = []
    union_rows: list[dict[str, Any]] = []
    for target_seq, spec in _WORKER_CONTIG_SPECS.items():
        seq_length = int(spec["seq_length"])
        valid_start = min(edge_mask_bp, max((seq_length - 1) // 2, 0))
        valid_end = seq_length - valid_start - 1
        clipped_proviruses = [
            (max(start, valid_start), min(end, valid_end))
            for start, end in spec["provirus_union"]
            if min(end, valid_end) >= max(start, valid_start)
        ]
        tested_intervals = [
            (
                max(int(region["start0"]), valid_start),
                min(int(region["end0"]), valid_end),
            )
            for region in spec["test_regions"]
            if min(int(region["end0"]), valid_end)
            >= max(int(region["start0"]), valid_start)
        ]
        intervals_needed = _merge_intervals(
            clipped_proviruses + tested_intervals
        )
        positions, depths = _read_sparse_intervals(
            _WORKER_DEPTH_ROOT / sample / f"{target_seq}.tsv",
            intervals_needed,
        )
        union_n, union_sum, union_sumsq, _ = _interval_moments(
            positions, depths, clipped_proviruses, min_coverage
        )
        union_rows.append(
            {
                "sample": sample,
                "genome": spec["genome"],
                "target_seq": target_seq,
                "provirus_union_n_bases": union_n,
                "provirus_union_depth_sum": union_sum,
                "provirus_union_depth_sum_squares": union_sumsq,
            }
        )

        for region in spec["test_regions"]:
            start = max(int(region["start0"]), valid_start)
            end = min(int(region["end0"]), valid_end)
            inside = _interval_moments(
                positions, depths, [(start, end)], min_coverage
            )
            outside_provirus = _subtract_intervals(
                (start, end), clipped_proviruses
            )
            extra = _interval_moments(
                positions, depths, outside_provirus, min_coverage
            )
            region_rows.append(
                {
                    "sample": sample,
                    "region_key": region["region_key"],
                    "region_type": region["region_type"],
                    "genome": spec["genome"],
                    "target_seq": target_seq,
                    "region_n_bases": inside[0],
                    "region_depth_sum": inside[1],
                    "region_depth_sum_squares": inside[2],
                    "region_bases_at_min_coverage": inside[3],
                    "test_interval_outside_provirus_n_bases": extra[0],
                    "test_interval_outside_provirus_depth_sum": extra[1],
                    "test_interval_outside_provirus_depth_sum_squares": extra[2],
                }
            )

    return sample, region_rows, union_rows


def _cohen_d_from_moments(
    inside_n: int,
    inside_sum: float,
    inside_sumsq: float,
    outside_n: int,
    outside_sum: float,
    outside_sumsq: float,
) -> float:
    if inside_n < 2 or outside_n < 2:
        return np.nan
    inside_variance = max(
        (inside_sumsq - inside_sum**2 / inside_n) / (inside_n - 1), 0.0
    )
    outside_variance = max(
        (outside_sumsq - outside_sum**2 / outside_n) / (outside_n - 1), 0.0
    )
    pooled_variance = (
        (inside_n - 1) * inside_variance
        + (outside_n - 1) * outside_variance
    ) / (inside_n + outside_n - 2)
    if pooled_variance == 0:
        inside_mean = inside_sum / inside_n
        outside_mean = outside_sum / outside_n
        if inside_mean == outside_mean:
            return 0.0
        return np.inf if inside_mean > outside_mean else -np.inf
    return float(
        (inside_sum / inside_n - outside_sum / outside_n)
        / math.sqrt(pooled_variance)
    )


def _covstats_path(
    references_dir: Path,
    sample: str,
    mapping_reference_prefix: str,
) -> Path:
    return references_dir / f"bowtie2_{mapping_reference_prefix}_{sample}_tot_covstats.txt"


def _read_host_moments(
    references_dir: Path,
    eligibility: pd.DataFrame,
    samples: list[str],
    edge_mask_bp: int,
    exact_host_moments_path: Path | None = None,
    mapping_reference_prefix: str = "strains_in_microbial_fraction_2018",
) -> pd.DataFrame:
    eligible = eligibility.loc[eligibility["eligible_for_votu_blocks"]].copy()
    eligible["target_seq"] = _normalise_seq(eligible["target_seq"])
    eligible = eligible[["target_seq", "genome", "sequence_length"]].drop_duplicates("target_seq")

    rows: list[pd.DataFrame] = []
    for sample in samples:
        path = _covstats_path(references_dir, sample, mapping_reference_prefix)
        if not path.exists():
            raise FileNotFoundError(path)
        covstats = pd.read_csv(path, sep="\t")
        if covstats.shape[1] != 8:
            raise ValueError(f"Expected eight covstats columns in {path}; found {covstats.shape[1]}")
        covstats.columns = [
            "target_seq", "mean_depth", "length", "covered_bases",
            "read_count", "variance", "trimmed_mean", "rpkm",
        ]
        covstats["target_seq"] = _normalise_seq(covstats["target_seq"])
        covstats = covstats.merge(eligible, on="target_seq", how="inner", validate="one_to_one")
        if covstats.empty:
            raise ValueError(f"No eligible contigs matched covstats for {sample}")

        # The 150-bp edge mask is represented by scaling each contig's covstats
        # coverage statistics to its effective length. Exact interval values still come
        # from per-base depth.  On >=100-kb contigs this approximation affects
        # at most 0.3% of a contig and avoids reading every background contig.
        covstats["effective_length"] = (
            covstats["length"].astype(int) - 2 * edge_mask_bp
        ).clip(lower=0)
        covstats["host_depth_sum"] = (
            covstats["mean_depth"].astype(float) * covstats["effective_length"]
        )
        covstats["host_depth_sum_squares"] = (
            covstats["variance"].astype(float)
            + covstats["mean_depth"].astype(float) ** 2
        ) * covstats["effective_length"]
        grouped = (
            covstats.groupby("genome", as_index=False)
            .agg(
                host_n_bases=("effective_length", "sum"),
                host_depth_sum=("host_depth_sum", "sum"),
                host_depth_sum_squares=("host_depth_sum_squares", "sum"),
                n_eligible_host_contigs=("target_seq", "nunique"),
            )
            .assign(sample=sample)
        )
        rows.append(grouped)
    host_moments = pd.concat(rows, ignore_index=True)
    if exact_host_moments_path is not None and exact_host_moments_path.exists():
        exact = pd.read_csv(exact_host_moments_path)
        exact = exact.loc[exact["sample"].isin(samples)].rename(
            columns={
                "masked_genome_n_bases": "host_n_bases",
                "masked_genome_depth_sum": "host_depth_sum",
                "masked_genome_depth_sum_squares": "host_depth_sum_squares",
            }
        )
        exact_columns = [
            "sample", "genome", "host_n_bases", "host_depth_sum",
            "host_depth_sum_squares",
        ]
        exact = exact[exact_columns].drop_duplicates(["sample", "genome"])
        host_moments = host_moments.merge(
            exact,
            on=["sample", "genome"],
            how="left",
            suffixes=("", "_exact"),
            validate="one_to_one",
        )
        for column in ["host_n_bases", "host_depth_sum", "host_depth_sum_squares"]:
            host_moments[column] = host_moments[f"{column}_exact"].fillna(
                host_moments[column]
            )
            host_moments = host_moments.drop(columns=f"{column}_exact")
    return host_moments


def _prepare_proviruses(
    prophage_regions_path: Path,
    eligibility: pd.DataFrame,
    detected_genomes_path: Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    regions = pd.read_csv(prophage_regions_path, sep="\t")
    regions["target_seq"] = _normalise_seq(regions["target_seq"])
    explicit = regions.loc[regions["topology"].eq("Provirus")].copy()

    eligible = eligibility.loc[eligibility["eligible_for_votu_blocks"]].copy()
    eligible["target_seq"] = _normalise_seq(eligible["target_seq"])
    eligible = eligible[["target_seq", "genome", "sequence_length"]].drop_duplicates("target_seq")
    explicit = explicit.merge(eligible, on="target_seq", how="inner", validate="many_to_one")
    explicit["start0"] = explicit["start_1based"].astype(int) - 1
    explicit["end0"] = explicit["end_1based_inclusive"].astype(int) - 1
    explicit["region_key"] = (
        "provirus::" + explicit["genome"].astype(str)
        + "::" + explicit["target_seq"].astype(str)
        + "::" + explicit["start_1based"].astype(int).astype(str)
        + "-" + explicit["end_1based_inclusive"].astype(int).astype(str)
    )
    explicit["region_type"] = "geNomad_provirus"
    if explicit["region_key"].duplicated().any():
        raise ValueError("Explicit geNomad provirus identifiers are not unique.")

    detected = pd.read_csv(detected_genomes_path, usecols=["genome"])["genome"].astype(str)
    tested = explicit.loc[explicit["genome"].isin(set(detected))].copy()
    return explicit.reset_index(drop=True), tested.reset_index(drop=True)


def _prepare_votu_regions(votu_activity_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    old_activity = pd.read_csv(votu_activity_path)
    required = {
        "genome", "sample", "target_seq", "block_start_1based",
        "block_end_1based", "source", "is_active_votu_from_list",
    }
    missing = required.difference(old_activity.columns)
    if missing:
        raise ValueError(f"vOTU activity table is missing columns: {sorted(missing)}")
    old_activity["target_seq"] = _normalise_seq(old_activity["target_seq"])

    identity = [
        "genome", "target_seq", "block_start_1based", "block_end_1based",
        "source", "is_active_votu_from_list",
    ]
    regions = old_activity[identity].drop_duplicates().copy()
    regions["start0"] = regions["block_start_1based"].astype(int) - 1
    regions["end0"] = regions["block_end_1based"].astype(int) - 1
    active_label = np.where(
        regions["is_active_votu_from_list"].astype(bool), "active_list", "not_active_list"
    )
    regions["region_key"] = (
        "votu::" + regions["genome"].astype(str)
        + "::" + regions["target_seq"].astype(str)
        + "::" + regions["block_start_1based"].astype(int).astype(str)
        + "-" + regions["block_end_1based"].astype(int).astype(str)
        + "::" + regions["source"].astype(str)
        + "::" + active_label
    )
    regions["region_type"] = "vOTU_block"
    if regions["region_key"].duplicated().any():
        raise ValueError("vOTU block identifiers are not unique.")
    return regions.reset_index(drop=True), old_activity


def _build_contig_specs(
    all_proviruses: pd.DataFrame,
    tested_proviruses: pd.DataFrame,
    votu_regions: pd.DataFrame,
    eligibility: pd.DataFrame,
) -> dict[str, dict[str, Any]]:
    test_regions = pd.concat(
        [
            tested_proviruses[["region_key", "region_type", "genome", "target_seq", "start0", "end0"]],
            votu_regions[["region_key", "region_type", "genome", "target_seq", "start0", "end0"]],
        ],
        ignore_index=True,
    )
    relevant_genomes = set(test_regions["genome"].astype(str))
    mask_regions = all_proviruses.loc[all_proviruses["genome"].isin(relevant_genomes)].copy()
    eligible = eligibility.loc[eligibility["eligible_for_votu_blocks"]].copy()
    eligible["target_seq"] = _normalise_seq(eligible["target_seq"])
    length_lookup = eligible.drop_duplicates("target_seq").set_index("target_seq")["sequence_length"].to_dict()
    genome_lookup = eligible.drop_duplicates("target_seq").set_index("target_seq")["genome"].to_dict()
    required_contigs = sorted(set(test_regions["target_seq"]) | set(mask_regions["target_seq"]))

    missing_lengths = sorted(set(required_contigs) - set(length_lookup))
    if missing_lengths:
        raise ValueError(f"Test/mask contigs are not eligible host contigs: {missing_lengths[:5]}")

    specs: dict[str, dict[str, Any]] = {}
    for target_seq in required_contigs:
        tests = test_regions.loc[test_regions["target_seq"].eq(target_seq)]
        masks = mask_regions.loc[mask_regions["target_seq"].eq(target_seq)]
        provirus_intervals = _merge_intervals(zip(masks["start0"], masks["end0"]))
        test_records = tests[
            ["region_key", "region_type", "start0", "end0"]
        ].to_dict("records")
        specs[target_seq] = {
            "genome": str(genome_lookup[target_seq]),
            "seq_length": int(length_lookup[target_seq]),
            "provirus_union": provirus_intervals,
            "test_regions": test_records,
        }
    return specs


def _calculate_activity(
    region_moments: pd.DataFrame,
    union_moments: pd.DataFrame,
    host_moments: pd.DataFrame,
    thresholds: ActivityThresholds,
) -> pd.DataFrame:
    union_by_genome = (
        union_moments.groupby(["sample", "genome"], as_index=False)
        .agg(
            provirus_union_n_bases=("provirus_union_n_bases", "sum"),
            provirus_union_depth_sum=("provirus_union_depth_sum", "sum"),
            provirus_union_depth_sum_squares=("provirus_union_depth_sum_squares", "sum"),
        )
    )
    activity = (
        region_moments
        .merge(union_by_genome, on=["sample", "genome"], how="left", validate="many_to_one")
        .merge(host_moments, on=["sample", "genome"], how="left", validate="many_to_one")
    )
    if activity[["host_n_bases", "host_depth_sum", "host_depth_sum_squares"]].isna().any().any():
        missing = activity.loc[activity["host_n_bases"].isna(), ["sample", "genome"]].drop_duplicates()
        raise ValueError(
            "Missing host per-base coverage statistics for interval rows: "
            f"{missing.head().to_dict('records')}"
        )

    activity["background_n_bases"] = (
        activity["host_n_bases"]
        - activity["provirus_union_n_bases"]
        - activity["test_interval_outside_provirus_n_bases"]
    )
    activity["background_depth_sum"] = (
        activity["host_depth_sum"]
        - activity["provirus_union_depth_sum"]
        - activity["test_interval_outside_provirus_depth_sum"]
    )
    activity["background_depth_sum_squares"] = (
        activity["host_depth_sum_squares"]
        - activity["provirus_union_depth_sum_squares"]
        - activity["test_interval_outside_provirus_depth_sum_squares"]
    )
    if (activity["background_n_bases"] <= 1).any():
        bad = activity.loc[activity["background_n_bases"] <= 1, ["genome", "region_key"]]
        raise ValueError(f"No usable non-proviral background: {bad.head().to_dict('records')}")

    # Covstats values are stored with finite decimal precision, whereas interval
    # coverage statistics are summed from integer per-base depths. In a low-depth genome
    # whose reads are almost entirely proviral, subtracting the exact interval
    # from the rounded whole-host statistic can produce a very small negative
    # remainder. Constrain those values to their mathematical lower bounds.
    activity["background_depth_sum"] = activity["background_depth_sum"].clip(lower=0.0)
    minimum_background_sumsq = (
        activity["background_depth_sum"] ** 2 / activity["background_n_bases"]
    )
    activity["background_depth_sum_squares"] = np.maximum(
        activity["background_depth_sum_squares"], minimum_background_sumsq
    )

    activity["region_mean_coverage"] = (
        activity["region_depth_sum"] / activity["region_n_bases"]
    )
    activity["region_coverage_breadth"] = (
        activity["region_bases_at_min_coverage"] / activity["region_n_bases"]
    )
    activity["nonproviral_host_mean_coverage"] = (
        activity["background_depth_sum"] / activity["background_n_bases"]
    )
    activity["fold_change_inside_vs_background"] = (
        activity["region_mean_coverage"] + 1e-12
    ) / (activity["nonproviral_host_mean_coverage"] + 1e-12)
    activity["log2_fold_change_inside_vs_background"] = np.log2(
        activity["fold_change_inside_vs_background"]
    )
    activity["cohen_d_inside_vs_background"] = [
        _cohen_d_from_moments(
            int(row.region_n_bases), float(row.region_depth_sum),
            float(row.region_depth_sum_squares), int(row.background_n_bases),
            float(row.background_depth_sum), float(row.background_depth_sum_squares),
        )
        for row in activity.itertuples(index=False)
    ]

    activity["passes_min_mean_coverage"] = activity["region_mean_coverage"].ge(
        thresholds.min_mean_coverage
    )
    activity["passes_min_breadth"] = activity["region_coverage_breadth"].ge(
        thresholds.min_breadth
    )
    activity["passes_coverage_ratio"] = activity["fold_change_inside_vs_background"].ge(
        thresholds.min_fold_change
    )
    activity["passes_cohen_d"] = activity["cohen_d_inside_vs_background"].ge(
        thresholds.min_cohen_d
    )
    present = activity["region_depth_sum"].gt(0)
    active = (
        present
        & activity["passes_min_mean_coverage"]
        & activity["passes_min_breadth"]
        & activity["passes_coverage_ratio"]
        & activity["passes_cohen_d"]
    )
    ambiguous = (
        present
        & ~active
        & activity["passes_coverage_ratio"]
        & activity["passes_cohen_d"]
    )
    activity["activity_status"] = np.select(
        [~present, active, ambiguous],
        ["NOT_PRESENT", "ACTIVE", "AMBIGUOUS"],
        default="DORMANT",
    )
    activity["ACTIVE"] = active
    activity["background_definition"] = BACKGROUND_DEFINITION
    return activity


def _merge_votu_activity(old_activity: pd.DataFrame, recalculated: pd.DataFrame) -> pd.DataFrame:
    rows = recalculated.loc[recalculated["region_type"].eq("vOTU_block")].copy()
    identity = [
        "genome", "target_seq", "block_start_1based", "block_end_1based",
        "source", "is_active_votu_from_list",
    ]
    keys = old_activity[identity].drop_duplicates().copy()
    keys["region_key"] = (
        "votu::" + keys["genome"].astype(str)
        + "::" + keys["target_seq"].astype(str)
        + "::" + keys["block_start_1based"].astype(int).astype(str)
        + "-" + keys["block_end_1based"].astype(int).astype(str)
        + "::" + keys["source"].astype(str)
        + "::" + np.where(keys["is_active_votu_from_list"].astype(bool), "active_list", "not_active_list")
    )
    old = old_activity.merge(keys, on=identity, how="left", validate="many_to_one")
    update_columns = [
        "region_key", "sample", "region_n_bases", "region_bases_at_min_coverage",
        "region_mean_coverage", "nonproviral_host_mean_coverage",
        "region_coverage_breadth", "fold_change_inside_vs_background",
        "log2_fold_change_inside_vs_background", "cohen_d_inside_vs_background",
        "passes_min_mean_coverage", "passes_min_breadth", "passes_coverage_ratio",
        "passes_cohen_d", "activity_status", "ACTIVE", "background_n_bases",
        "background_definition",
    ]
    updated = old.drop(
        columns=[
            column for column in update_columns
            if column in old.columns and column not in {"region_key", "sample"}
        ],
        errors="ignore",
    ).merge(rows[update_columns], on=["region_key", "sample"], how="left", validate="one_to_one")
    if updated["activity_status"].isna().any():
        missing = updated.loc[updated["activity_status"].isna(), ["region_key", "sample"]]
        raise ValueError(f"Missing recalculated vOTU rows: {missing.head().to_dict('records')}")
    updated = updated.rename(
        columns={
            "region_n_bases": "block_effective_length_bp",
            "region_bases_at_min_coverage": "block_bases_at_min_coverage",
            "region_mean_coverage": "block_mean_coverage",
            "nonproviral_host_mean_coverage": "host_mean_coverage",
            "region_coverage_breadth": "block_coverage_breadth",
        }
    )
    return updated.drop(columns="region_key").sort_values(
        ["genome", "source", "block_index", "sample"]
    ).reset_index(drop=True)


def _summarise_proviruses(activity: pd.DataFrame, metadata: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    active_samples = (
        activity.loc[activity["ACTIVE"]]
        .groupby("region_key")["sample"]
        .agg(lambda values: ", ".join(sorted(values.astype(str))))
    )
    summary = (
        activity.groupby("region_key", as_index=False)
        .agg(
            n_samples_evaluated=("sample", "nunique"),
            n_samples_ACTIVE=("ACTIVE", "sum"),
            n_samples_AMBIGUOUS=("activity_status", lambda x: x.eq("AMBIGUOUS").sum()),
            n_samples_DORMANT=("activity_status", lambda x: x.eq("DORMANT").sum()),
            n_samples_NOT_PRESENT=("activity_status", lambda x: x.eq("NOT_PRESENT").sum()),
            max_region_mean_coverage=("region_mean_coverage", "max"),
            max_region_coverage_breadth=("region_coverage_breadth", "max"),
            max_fold_change_inside_vs_background=("fold_change_inside_vs_background", "max"),
            max_cohen_d_inside_vs_background=("cohen_d_inside_vs_background", "max"),
        )
        .merge(active_samples.rename("ACTIVE_samples"), on="region_key", how="left")
    )
    summary["ACTIVE_samples"] = summary["ACTIVE_samples"].fillna("")
    summary["propagate_style_active"] = summary["n_samples_ACTIVE"].gt(0)
    meta_columns = [
        "region_key", "genome", "target_seq", "start_1based",
        "end_1based_inclusive", "predicted_region_length", "virus_score",
        "n_hallmarks", "taxonomy",
    ]
    summary = metadata[meta_columns].merge(summary, on="region_key", how="left", validate="one_to_one")

    present_calls = activity["activity_status"].ne("NOT_PRESENT")
    present_intervals = summary["n_samples_NOT_PRESENT"].lt(summary["n_samples_evaluated"])
    overview = pd.DataFrame(
        [
            {"metric": "Eligible explicit geNomad Provirus intervals tested", "value": len(summary)},
            {"metric": "Bacterial reference genomes tested", "value": summary["genome"].nunique()},
            {"metric": "Provirus intervals detected in >=1 virome", "value": int(present_intervals.sum())},
            {"metric": "Provirus intervals ACTIVE in >=1 virome", "value": int(summary["propagate_style_active"].sum())},
            {
                "metric": "Detected provirus intervals never ACTIVE",
                "value": int((present_intervals & ~summary["propagate_style_active"]).sum()),
            },
            {
                "metric": "Detected provirus intervals ACTIVE in >=1 virome (%)",
                "value": 100 * summary.loc[present_intervals, "propagate_style_active"].mean(),
            },
            {
                "metric": "Detected provirus intervals never ACTIVE (%)",
                "value": 100 * (~summary.loc[present_intervals, "propagate_style_active"]).mean(),
            },
            {"metric": "Detected provirus-virome combinations", "value": int(present_calls.sum())},
            {"metric": "ACTIVE provirus-virome combinations", "value": int(activity["ACTIVE"].sum())},
            {
                "metric": "Detected provirus-virome combinations ACTIVE (%)",
                "value": 100 * activity["ACTIVE"].sum() / present_calls.sum(),
            },
        ]
    )
    return summary, overview


def run_parallel_activity_analysis(
    *,
    viral_fraction_dir: Path,
    bacterial_fraction_dir: Path,
    output_dir: Path,
    votu_activity_path: Path,
    host_eligibility_path: Path,
    detected_genomes_path: Path,
    exact_host_coverage_statistics_path: Path | None = None,
    mapping_reference_prefix: str = "strains_in_microbial_fraction_2018",
    depth_root: Path | None = None,
    prophage_regions_path: Path | None = None,
    samples: list[str] | None = None,
    workers: int | None = None,
    thresholds: ActivityThresholds | None = None,
    include_votu_blocks: bool = True,
    exact_host_moments_path: Path | None = None,
) -> ActivityRunResult:
    samples = list(samples or SAMPLES_36)
    thresholds = thresholds or ActivityThresholds()
    if (
        exact_host_coverage_statistics_path is not None
        and exact_host_moments_path is not None
    ):
        raise ValueError(
            "Provide exact_host_coverage_statistics_path only; "
            "exact_host_moments_path is retained only for older calls."
        )
    exact_host_coverage_statistics_path = (
        exact_host_coverage_statistics_path
        if exact_host_coverage_statistics_path is not None
        else exact_host_moments_path
    )
    workers = workers or min(8, len(samples), os.cpu_count() or 1)
    output_dir.mkdir(parents=True, exist_ok=True)

    references_dir = viral_fraction_dir / "06_MAPPING" / "REFERENCES"
    depth_root = Path(depth_root) if depth_root is not None else (
        references_dir / "STRAINS_IN_MICROBIAL_SPLIT"
    )
    prophage_regions_path = Path(prophage_regions_path) if prophage_regions_path is not None else (
        bacterial_fraction_dir / "MICROBIAL_GENOMES_PHYLLOVIR"
        / "strains_in_microbial_fraction_2018_geNomad_virus_masked_regions.tsv"
    )
    eligibility = pd.read_csv(host_eligibility_path)
    all_proviruses, tested_proviruses = _prepare_proviruses(
        prophage_regions_path, eligibility, detected_genomes_path
    )
    all_votu_regions, old_votu_activity = _prepare_votu_regions(votu_activity_path)
    old_votu_activity = old_votu_activity.loc[
        old_votu_activity["sample"].isin(samples)
    ].copy()
    votu_regions = (
        all_votu_regions
        if include_votu_blocks
        else all_votu_regions.iloc[0:0].copy()
    )
    contig_specs = _build_contig_specs(
        all_proviruses, tested_proviruses, votu_regions, eligibility
    )

    print(
        f"Testing {len(tested_proviruses)} explicit geNomad proviruses from "
        f"{tested_proviruses['genome'].nunique()} detected genomes and "
        f"{len(votu_regions)} vOTU blocks across {len(samples)} viromes."
    )
    print(
        f"Reading sparse depth for {len(contig_specs)} interval-bearing contigs "
        f"with {workers} sample workers."
    )
    sample_jobs = [
        (
            sample,
            thresholds.min_mean_coverage,
            thresholds.contig_edge_mask_bp,
        )
        for sample in samples
    ]
    region_rows: list[dict[str, Any]] = []
    union_rows: list[dict[str, Any]] = []
    with ProcessPoolExecutor(
        max_workers=workers,
        initializer=_init_worker,
        initargs=(contig_specs, str(depth_root)),
    ) as executor:
        futures = {executor.submit(_process_sample_depths, job): job[0] for job in sample_jobs}
        for completed, future in enumerate(as_completed(futures), start=1):
            sample, sample_region_rows, sample_union_rows = future.result()
            region_rows.extend(sample_region_rows)
            union_rows.extend(sample_union_rows)
            print(
                f"Completed per-base coverage statistics for {sample} "
                f"({completed}/{len(samples)})."
            )

    region_moments = pd.DataFrame(region_rows)
    union_moments = pd.DataFrame(union_rows)
    relevant_genomes = set(region_moments["genome"].astype(str))
    eligibility_for_hosts = eligibility.loc[eligibility["genome"].astype(str).isin(relevant_genomes)]
    host_moments = _read_host_moments(
        references_dir,
        eligibility_for_hosts,
        samples,
        thresholds.contig_edge_mask_bp,
        exact_host_coverage_statistics_path,
        mapping_reference_prefix,
    )
    activity = _calculate_activity(region_moments, union_moments, host_moments, thresholds)

    provirus_activity = activity.loc[activity["region_type"].eq("geNomad_provirus")].copy()
    provirus_activity = provirus_activity.merge(
        tested_proviruses[
            [
                "region_key", "start_1based", "end_1based_inclusive",
                "predicted_region_length", "virus_score", "n_hallmarks", "taxonomy",
            ]
        ],
        on="region_key",
        how="left",
        validate="many_to_one",
    )
    provirus_summary, provirus_overview = _summarise_proviruses(
        provirus_activity, tested_proviruses
    )
    votu_activity = (
        _merge_votu_activity(old_votu_activity, activity)
        if include_votu_blocks
        else old_votu_activity.copy()
    )
    prefix = "10_nonproviral_background"
    paths = {
        "interval_coverage_statistics": (
            output_dir / f"{prefix}_interval_per_base_coverage_statistics.csv"
        ),
        "provirus_activity": output_dir / f"{prefix}_genomad_provirus_activity_by_sample.csv",
        "provirus_summary": output_dir / f"{prefix}_genomad_provirus_activity_summary.csv",
        "provirus_overview": output_dir / f"{prefix}_genomad_provirus_activity_overview.csv",
    }
    if include_votu_blocks:
        paths["votu_activity"] = output_dir / f"{prefix}_votu_block_sample_activity.csv"
    activity.to_csv(paths["interval_coverage_statistics"], index=False)
    provirus_activity.to_csv(paths["provirus_activity"], index=False)
    provirus_summary.to_csv(paths["provirus_summary"], index=False)
    provirus_overview.to_csv(paths["provirus_overview"], index=False)
    if include_votu_blocks:
        votu_activity.to_csv(paths["votu_activity"], index=False)
    return ActivityRunResult(
        provirus_activity=provirus_activity,
        provirus_summary=provirus_summary,
        provirus_overview=provirus_overview,
        votu_activity=votu_activity,
        interval_coverage_statistics=activity,
        output_paths=paths,
    )


def _default_paths() -> dict[str, Path]:
    viral_fraction = Path("/home/lmf/PhylloVir/VIRAL_WORLD/VIRAL_FRACTION")
    bacterial_fraction = viral_fraction.parent / "BACTERIAL_FRACTION"
    supplementary = viral_fraction / "FIGURES_AND_TABLES" / "SUPPLEMENTARY"
    output = supplementary / "10_vOTU_block_activity" / "strain_depth_profiles"
    return {
        "viral_fraction_dir": viral_fraction,
        "bacterial_fraction_dir": bacterial_fraction,
        "output_dir": output,
        "votu_activity_path": output / "10_supplemented_votu_blocks_top0_rpkm_plus_votu_block_genomes_votu_block_sample_activity.csv",
        "host_eligibility_path": supplementary / "09_vOTU_block_test_set" / "09_host_contig_eligibility.csv",
        "detected_genomes_path": supplementary / "08_microbial_signal_in_viromes" / "08_detected_bacterial_reference_genomes_breadth30.csv",
        "exact_host_coverage_statistics_path": output / "10_supplemented_votu_blocks_top0_rpkm_plus_votu_block_genomes_split_file_summary.csv",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=min(8, os.cpu_count() or 1))
    args = parser.parse_args()
    result = run_parallel_activity_analysis(**_default_paths(), workers=args.workers)
    print(result.provirus_overview.to_string(index=False))
    print("Saved outputs:")
    for label, path in result.output_paths.items():
        print(f"  {label}: {path}")


if __name__ == "__main__":
    main()
