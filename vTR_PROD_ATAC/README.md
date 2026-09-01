# PROD-ATAC analysis for the vTR study

This directory contains the core ArchR-based R code used to generate the peak-accessibility, motif-enrichment, and co-accessibility datasets for the associated vTR manuscript. It represents a manuscript-relevant subset of the broader PROD-ATAC analysis.

## Workflow

The script performs cell quality control, vTR assignment, LSI/UMAP analysis, peak calling, differential-accessibility testing, motif enrichment, and co-accessibility analysis using the human `hg38` reference genome.

## Requirements

- R with `ArchR`, `BSgenome.Hsapiens.UCSC.hg38`, `Biostrings`, `ggplot2`, `dplyr`, `data.table`, `stringr`, and `chromVARmotifs`
- MACS3
- Eight ATAC-seq fragment files (`.tsv.gz`)
- Eight sample-specific association tables containing `CBC` and `vTR` columns

## Usage

Update the fragment-file, association-table, and MACS3 paths in `PROD_ATAC_vTR_analysis.R`, then run:

```bash
Rscript PROD_ATAC_vTR_analysis.R
```

The supplied script retains cells with TSS enrichment ≥7 and ≥40,000 fragments. Doublet scores are calculated for ArchR compatibility but are not used for filtering. Cells without a unique vTR assignment are excluded from vTR-level analyses.
