# Data for scICER Workflow

This folder provides the source tables and dataset preparation instructions used by the example workflows and manuscript figure scripts in this repository.

For the main analysis scripts, see `../codes/`.

## 1. Benchmark summary tables for the benchmark comparison panels

### Source
These tables were prepared from the benchmark result summaries generated from the 48 scICE benchmark datasets and the extended workflow-plus-graph comparison analyses.

### Included files

- `tables/benchmark/benchmark_summary_primary.tsv`
  Summary table used by:
  - `codes/figures/Figure1/Fig1_AB_scice_vs_scicer_ratio.R`
  - `codes/figures/FigureS1/FigS1_consistent_cluster_ranges.R`
  - `codes/figures/FigureS2/FigS2_AB_paper_vs_scice_ranges.R`
  - `codes/figures/FigureS2/FigS2_CD_julia_vs_r_umap_ranges.R`
  - `codes/figures/FigureS2/FigS2_EF_julia_vs_r_umap_ratio.R`

- `tables/benchmark/benchmark_summary_workflow_graph.tsv`
  Summary table used by:
  - `codes/figures/Figure1/Fig1_CD_workflow_graph_ratio.R`
  - `codes/figures/FigureS3/FigS3_workflow_graph_ranges.R`

### Table contents
Both tables contain the figure-level summary fields used in the plotting scripts:

- `sample`
- `sample_type`
- `method`
- `graph`
- `Optimal_Cluster_by_ECS`
- `consistent_cluster`
- `Consistent_Clustering_Ratio`

## 2. Iterative analysis datasets

### Source
These datasets were used for the iterative scICER analyses shown in the iterative clustering panels.

### Expected local files

- `data/generated/iterative/ZhengMix_4757/seurat.scICER.qs`
- `data/generated/iterative/sim_tcell_2/seurat.scICER.qs`

### Used by

- `codes/figures/FigureS4/FigS4_AH_zhengmix4757_ic.R`
- `codes/figures/FigureS4/FigS4_BC_zhengmix4757_labels.R`
- `codes/figures/FigureS4/FigS4_GI_zhengmix4757_iterative_graph.R`
- `codes/figures/FigureS4/FigS4_DEFJKL_simTcell2_iterative.R`

## 3. Runtime and memory source tables

### Source
These tables were prepared from the runtime benchmark summaries and the memory monitoring outputs used for the pancreas runtime comparison.

### Included files

- `tables/runtime/scaling/runtime_summary_scicer.tsv`
  Runtime and peak-memory summary across datasets and thread counts for the Figure S7A-B comparison.

- `tables/runtime/scaling/runtime_summary_scice_reference.tsv`
  Reference scICE timing table used in the total-time panel of Figure S7A-B.

- `tables/runtime/pancreas/full_usage_memory.tsv`
  Memory trace for the full-object pancreas run.

- `tables/runtime/pancreas/full_window.tsv`
  Start and end boundaries used to isolate the clustering step from the full-object trace.

- `tables/runtime/pancreas/diet_usage_memory.tsv`
  Memory trace for the DietSeurat-pruned pancreas run.

- `tables/runtime/pancreas/diet_window.tsv`
  Start and end boundaries used to isolate the clustering step from the DietSeurat trace.

### Used by

- `codes/figures/FigureS7/FigS7_AB_runtime_scaling.R`
- `codes/figures/FigureS7/FigS7_C_dietseurat_memory.R`

## 4. scICE benchmark datasets

### Source
These datasets are from the article:
**"scICE: enhancing clustering reliability and efficiency of scRNA-seq data with multi-cluster label consistency evaluation"**

Code and data availability:
- https://www.nature.com/articles/s41467-025-60702-8#code-availability

Download from Zenodo:
- https://zenodo.org/records/15113898

### Notes
The public figure scripts in this repository use the prepared benchmark summary tables listed above.

## 5. Human pancreas integration benchmark dataset

### Source
These data are from the article:
**"Benchmarking atlas-level data integration in single-cell genomics"**

Data availability:
- https://www.nature.com/articles/s41592-021-01336-8#data-availability

Download from Figshare:
- https://figshare.com/articles/dataset/Benchmarking_atlas-level_data_integration_in_single-cell_genomics_-_integration_task_datasets_Immune_and_pancreas_/12420968

### Expected local files

- `data/generated/human_pancreas_norm_complexBatch/seurat.sclens.qs`
- `data/generated/human_pancreas_norm_complexBatch/seurat.scvi.HVG.qs`
- `data/generated/no_integration_pancreas/scICE_batch_data.rds`
- `data/generated/no_integration_pancreas/scICE_celltype_data.rds`

### Used by

- `codes/figures/Figure1/Fig1_EFGJK_harmony_pancreas.R`
- `codes/figures/Figure1/Fig1_HILM_scvi_pancreas.R`
- `codes/figures/FigureS5/FigS5_original_vs_sclens_pancreas.R`

## 6. Large pancreas atlas from GSE221156

### Source
Dataset accession:
- `GSE221156`

Download URL:
- https://datasets.cellxgene.cziscience.com/74a82e6d-2c1a-4328-9011-c6a437747854.h5ad

### Expected local files

- `data/raw/pancreas_gse221156/74a82e6d-2c1a-4328-9011-c6a437747854.h5ad`
- `data/generated/pancreas_gse221156/pancreas_gse221156_raw.rds`
- `data/generated/pancreas_gse221156/pancreatic_harmony.qs`

### Used by

- `codes/figures/FigureS6/FigS6_gse221156_pancreas_atlas.R`
- `codes/figures/FigureS7/FigS7_C_dietseurat_memory.R`

### Preparation
The pancreas atlas preprocessing scripts are provided in `data/preprocessing/`:

- `preprocessing/download_pancreatic_gse221156.sh`
- `preprocessing/convert_h5ad_to_seurat.R`
- `preprocessing/prepare_pancreatic_harmony.R`

These scripts download the public `.h5ad` file, convert it to a Seurat object, and prepare the Harmony-based object used by the example scripts.
