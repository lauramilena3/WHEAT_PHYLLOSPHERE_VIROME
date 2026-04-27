# Wheat Phyllosphere Virome

This repository contains the downstream analysis code supporting the manuscript:

**Active Temperate Phages Linked to Core Bacterial Taxa Dominate the Wheat Phyllosphere Virome**

The analyses use paired wheat phyllosphere virome and microbial metagenome data to characterize active viral operational taxonomic units (avOTUs), phage lifestyle distributions, vOTU prevalence, microbial genus prevalence, and phage-host community relationships.

## Repository structure

```text
VIRAL/NOTEBOOKS/      Downstream viral analyses, including vOTU filtering, active-vOTU selection, abundance normalization, mapping summaries, prevalence analyses, and viral figure/table generation.

MICROBIAL/NOTEBOOKS/  Microbial abundance, taxonomy, prevalence, and microbial figure/table generation.
```

## Data availability

Sequencing data generated and analyzed in this study are available under NCBI BioProject:

**PRJNA1428567**

The upstream quality-control, assembly, vOTU clustering, annotation, and viral abundance workflow was implemented with the MOSAIC pipeline:

<https://github.com/lauramilena3/MOSAIC>

## Manuscript outputs

The notebooks generate the processed data, summary statistics, supplementary tables, and figure inputs used in the manuscript, including:

- active vOTU abundance and prevalence analyses
- lifestyle comparisons between temperate and virulent avOTUs
- microbial genus-level abundance and prevalence summaries
- host-prediction summaries
- bacterial-viral community concordance analyses
- main and supplementary figure source data

## Citation

If using this repository, please cite the associated manuscript.
