# Manuscript Figure Scripts

Run these scripts from the repository root so that the relative `data/...` paths resolve correctly.

## Figure 1

- `Figure1/Fig1_AB_scice_vs_scicer_ratio.R`: consistent-clustering ratio comparison between scICE and scICER.
- `Figure1/Fig1_CD_workflow_graph_ratio.R`: workflow-plus-graph comparison across the benchmark collection.
- `Figure1/Fig1_EFGJK_harmony_pancreas.R`: Harmony integration example on the pancreas benchmark.
- `Figure1/Fig1_HILM_scvi_pancreas.R`: scVI integration example on the pancreas benchmark.

## Figure S1

- `FigureS1/FigS1_consistent_cluster_ranges.R`: consistent-cluster range comparison for the main benchmark table.

## Figure S2

- `FigureS2/FigS2_AB_paper_vs_scice_ranges.R`: original scICE ranges versus the reproduced scICE ranges.
- `FigureS2/FigS2_CD_julia_vs_r_umap_ranges.R`: cluster-range comparison between the two UMAP workflows.
- `FigureS2/FigS2_EF_julia_vs_r_umap_ratio.R`: ratio summary for the same UMAP comparison.

## Figure S3

- `FigureS3/FigS3_workflow_graph_ranges.R`: consistent-cluster ranges for the workflow-plus-graph benchmark panel.

## Figure S4

- `FigureS4/FigS4_AH_zhengmix4757_ic.R`: IC distributions for the parent and iterative ZhengMix_4757 analyses.
- `FigureS4/FigS4_BC_zhengmix4757_labels.R`: cell-type and cluster labels for ZhengMix_4757.
- `FigureS4/FigS4_GI_zhengmix4757_iterative_graph.R`: iterative subgraph visualization for ZhengMix_4757.
- `FigureS4/FigS4_DEFJKL_simTcell2_iterative.R`: iterative analysis panels for `sim_Tcell_2`.

## Figure S5

- `FigureS5/FigS5_original_vs_sclens_pancreas.R`: original versus integrated pancreas embeddings.

## Figure S6

- `FigureS6/FigS6_gse221156_pancreas_atlas.R`: large pancreas atlas analysis based on GSE221156.

## Figure S7

- `FigureS7/FigS7_AB_runtime_scaling.R`: runtime and peak-memory scaling across thread counts.
- `FigureS7/FigS7_C_dietseurat_memory.R`: memory comparison between the full and DietSeurat-pruned pancreas runs.

## Shared helper

- `utils.R`
