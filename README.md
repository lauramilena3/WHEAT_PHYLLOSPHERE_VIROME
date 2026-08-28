# Wheat Phyllosphere Virome

This repository contains the analysis code supporting the manuscript:

**Active Temperate Phages Linked to Core Bacterial Taxa Dominate the Wheat Phyllosphere Virome**

The study combines paired viral-fraction and microbial metagenomic sequencing of the wheat phyllosphere to characterize viral community composition, active viral operational taxonomic units (avOTUs), phage lifestyle, viral and microbial prevalence, and relationships between phage and bacterial communities.

## Analysis overview

The computational analysis consisted of two main components:

1. **Upstream metagenomic processing with MOSAIC**, including read quality control, assembly, viral identification, vOTU generation and filtering, read mapping, and abundance estimation.
2. **Downstream analyses in Jupyter notebooks**, including identification of active vOTUs, abundance normalization, prevalence analyses, phage lifestyle analyses, microbial community analyses, host-prediction summaries, statistical analyses, and generation of manuscript figures and tables.

The MOSAIC workflow is available at:

https://github.com/lauramilena3/MOSAIC

## Repository structure

```text
WHEAT_PHYLLOSPHERE_VIROME/
├── VIRAL/
│   └── NOTEBOOKS/
│       ├── 01_MOSAIC_01_QC.tot.ipynb
│       ├── 02_MOSAIC_03_assembly_short.tot.ipynb
│       ├── 03_MOSAIC_04_viral_ID_tot.ipynb
│       ├── 04_MOSAIC_05_vOTU_representative.tot.ipynb
│       ├── 05_MOSAIC_05_vOTU_filtering.tot.ipynb
│       ├── 06_MOSAIC_07_mapping_statistics_tot.ipynb
│       ├── 07_MOSAIC_07_subsampling.ipynb
│       ├── 08_microbial_signal_in_viromes.ipynb
│       ├── 09_vOTU_block_test_set.ipynb
│       ├── 10_vOTU_block_activity.ipynb
│       ├── 11_rank_abundance_threshold_selection.ipynb
│       ├── 12_MOSAIC_07_Normalise_active_vOTUs.tot.ipynb
│       ├── 13_primary_results_lifestyle_abundance.ipynb
│       ├── 14_rank_abundance_threshold_robustness.ipynb
│       └── 15_vOTU_quality_robustness.ipynb
│
├── MICROBIAL/
│   └── NOTEBOOKS/
│       ├── MOSAIC_01_QC.tot.ipynb
│       ├── MOSAIC_03_assembly_short.tot.ipynb
│       ├── MOSAIC_04_viral_ID_tot.ipynb
│       ├── MOSAIC_05-1_vOTU_representative.tot.ipynb
│       ├── MOSAIC_05-2_vOTU_filtering.tot.ipynb
│       └── MOSAIC_07_microbial_abundance.ipynb
│
├── environment.yaml
└── README.md
```

### Viral analyses

The notebooks under `VIRAL/NOTEBOOKS/` document the processing and downstream analysis of the viral-fraction metagenomes. These analyses include:

* sequencing and assembly quality summaries
* viral sequence identification
* generation and filtering of representative vOTUs
* read-mapping statistics
* subsampling analyses
* evaluation of microbial signal in the viral fraction
* selection of active vOTUs
* rank-abundance threshold selection
* abundance normalization
* phage lifestyle and abundance analyses
* robustness analyses for activity thresholds and vOTU quality
* generation of manuscript figures, statistics, and supplementary outputs

### Microbial analyses

The notebooks under `MICROBIAL/NOTEBOOKS/` contain analyses associated with the paired microbial metagenomes, including:

* sequencing and assembly quality summaries
* sequence classification and filtering
* microbial abundance estimation
* taxonomic summaries
* microbial community composition and prevalence analyses
* data used for comparisons between bacterial and viral communities

## Upstream MOSAIC processing

The viral-fraction sequencing data were processed with the MOSAIC Snakemake workflow.

### Workflow run recorded on 15 July 2025

The following MOSAIC/Snakemake configuration was recorded on **15 July 2025**:

```bash
snakemake \
    --use-conda \
    -p \
    runWorkflow \
    --config \
        input_dir=/home/lmf/PhylloVir/VIRAL_WORLD/VIRAL_FRACTION/00_RAW_DATA \
        ecc_memory=16000 \        additional_reference_contigs=/home/lmf/PhylloVir/VIRAL_WORLD/VIRAL_FRACTION/erwinia_pseudomonas_and_reference_phages.fasta \
        assembly_stats=True \
        subassembly=True \
        min_votu_length=10000 \
    -j 144 \
    -k \
    -n
```

This configuration specified:

* `input_dir`: directory containing the raw viral-fraction sequencing reads
* `ecc_memory=16000`: memory setting used by the workflow
* `subassembly=True`: enabled the subassembly procedure
* `additional_reference_contigs`: additional *Erwinia*, *Pseudomonas*, and reference-phage sequences supplied to the workflow
* `assembly_stats=True`: enabled assembly-statistics calculations
* `min_votu_length=10000`: minimum vOTU sequence length of 10 kb
* `-j 144`: allowed up to 144 concurrent jobs/cores according to the Snakemake execution environment
* `-k`: instructed Snakemake to continue with independent jobs when possible after an error
* `--use-conda`: used the Conda environments specified by the workflow

The absolute paths above correspond to the original computational environment and must be changed when reproducing the workflow on another system.

## Computational environment

A Conda environment containing the software used for analysis is provided as:

```text
environment.yaml
```

The environment currently specifies, among other dependencies:

* Snakemake 7.18.2
* Python 3.11
* pandas
* NumPy
* Matplotlib
* seaborn
* Biopython
* scikit-bio
* Jupyter

It can be created, for example, with:

```bash
mamba env create -f environment.yaml -n wheat_phyllosphere_virome
conda activate wheat_phyllosphere_virome
```

## Data availability

Sequencing data generated and analyzed in this study are available through NCBI under BioProject:

**PRJNA1428567**

Large intermediate files and raw sequencing data are not stored in this repository.

## Manuscript outputs

The notebooks generate processed data, summary statistics, supplementary tables, and figure source data used in the manuscript, including analyses of:

* active vOTU abundance and prevalence
* temperate and virulent phage lifestyles
* viral community composition
* microbial genus-level abundance and prevalence
* phage-host associations
* bacterial-viral community concordance
* robustness of active-vOTU definitions
* main and supplementary manuscript figures

## Reproducibility

For reproducibility, the repository provides:

* the analysis notebooks used for viral and microbial analyses
* the MOSAIC workflow configuration used for upstream processing
* the date of the recorded upstream workflow run
* a Conda environment specifying the main analysis dependencies
* the NCBI BioProject accession for the sequencing data

Because database contents and software repositories can change over time, reproducing the original upstream analysis should use the MOSAIC version and reference databases corresponding as closely as possible to those used in July 2025.

## Citation

If using this repository or associated data, please cite the corresponding manuscript:

**Active Temperate Phages Linked to Core Bacterial Taxa Dominate the Wheat Phyllosphere Virome**

Full citation information will be added upon publication.
