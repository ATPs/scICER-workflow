# scICER Workflow: Example Workflows and Manuscript Figure Scripts for Integrated Single-Cell Analyses

## **Introduction**
**scICER** is an R package for evaluating cluster consistency in single-cell RNA-seq data. It reimplements the scICE methodology to assess the stability of cluster labels across stochastic runs. Fully integrated with Seurat, scICER works on embeddings corrected by **Harmony** or **scVI** enabling consistent clustering analysis in complex, multi-sample datasets.

## Data sources, preprocessing, and conversion
The repository includes source tables, dataset preparation notes, and preprocessing scripts for the benchmark analyses, pancreas integration examples, the large pancreas atlas from **GSE221156**, and the runtime comparison panels.

The corresponding materials are organized under `data/`, including:
- `data/README.md`
- `data/preprocessing/`
- `data/tables/`

## Codes examples to run scICER
The `codes/` directory contains example R scripts demonstrating the use of **scICER** in several common analysis settings.

The workflow examples include:
- `codes/run_scICER.R`
- `codes/harmony_integration_scICER.R`
- `codes/scvi_integration_scICER.R`

These scripts are written as reusable examples and use generic inputs rather than manuscript-specific output paths.

## Manuscript figure scripts
The `codes/figures/` directory contains the R scripts used to reproduce the manuscript figures. Scripts are organized by figure, and each `.R` file corresponds to one figure block.

The figure scripts are accompanied by:
- `codes/figures/README.md`
- source tables in `data/tables/`
