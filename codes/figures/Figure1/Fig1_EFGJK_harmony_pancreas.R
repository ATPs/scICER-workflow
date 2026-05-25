#!/usr/bin/env Rscript

library(Seurat)
library(ggplot2)
library(qs)
library(ClustAssess)
library(scICER)
library(harmony)
library(patchwork)
source("codes/figures/utils.R")

input_file <- "data/generated/human_pancreas_norm_complexBatch/seurat.sclens.qs"
if (!file.exists(input_file)) {
  stop("Missing input file: ", input_file, "
See data/README.md for preparation details.")
}

sample_obj <- qs::qread(input_file)
sample_obj <- RunHarmony(
  object = sample_obj,
  group.by.vars = "tech",
  reduction = "pca",
  dims.use = seq_len(ncol(Embeddings(sample_obj, "pca"))),
  assay.use = "RNA",
  project.dim = FALSE,
  verbose = TRUE,
  reduction.save = "harmony_pca"
)
sample_obj <- RunUMAP(
  object = sample_obj,
  dims = seq_len(ncol(Embeddings(sample_obj, "harmony_pca"))),
  reduction = "harmony_pca",
  reduction.name = "harmony_umap",
  verbose = TRUE
)
sample_obj <- FindNeighbors(
  object = sample_obj,
  reduction = "harmony_pca",
  dims = seq_len(ncol(Embeddings(sample_obj, "harmony_pca"))),
  graph.name = c("harmony_nn", "harmony_snn")
)

scice_results <- scICE_clustering(
  object = sample_obj,
  cluster_range = 2:20,
  remove_threshold = Inf,
  n_workers = 8,
  n_trials = 30,
  n_bootstrap = 200,
  seed = 123,
  verbose = TRUE,
  graph_name = "harmony_snn"
)

score_table <- data.frame(cluster_number = scice_results$n_cluster, ic_score = scice_results$ic, is_consistent = scice_results$ic <= 1.005)
label_df <- get_robust_labels(scice_results, return_seurat = FALSE, threshold = Inf)
label_df$celltype <- sample_obj@meta.data$celltype
for (k in 2:20) {
  cluster_col <- paste0("clusters_", k)
  if (cluster_col %in% names(label_df)) {
    score_table[score_table$cluster_number == k, "ECS_score"] <- element_sim(label_df$celltype, label_df[[cluster_col]])
  }
}
valid_clusters <- subset(score_table, is_consistent)
optimal_cluster <- valid_clusters$cluster_number[which.max(valid_clusters$ECS_score)]

sample_obj <- get_robust_labels(scice_results, return_seurat = TRUE, threshold = Inf)
sample_obj$optimal_cluster <- sample_obj@meta.data[[paste0("clusters_", optimal_cluster)]]

base_theme <- theme(axis.ticks.length = unit(-0.15, "cm"), axis.text = element_text(margin = margin(3, 3, 3, 3)), panel.border = element_rect(color = "black", fill = NA, linewidth = 1), legend.position = "bottom", legend.direction = "vertical", legend.box = "vertical")

p1 <- DimPlot(sample_obj, reduction = "umap", group.by = "tech", pt.size = 3, alpha = 1, raster = TRUE, raster.dpi = c(1200, 1200)) + base_theme
p2 <- DimPlot(sample_obj, reduction = "harmony_umap", group.by = "tech", pt.size = 3, alpha = 1, raster = TRUE, raster.dpi = c(1200, 1200)) + base_theme
p3 <- DimPlot(sample_obj, reduction = "harmony_umap", group.by = "celltype", pt.size = 3, alpha = 1, raster = TRUE, raster.dpi = c(1200, 1200)) + base_theme
p4 <- DimPlot(sample_obj, reduction = "harmony_umap", group.by = "optimal_cluster", pt.size = 3, alpha = 1, raster = TRUE, raster.dpi = c(1200, 1200)) + base_theme
p_ic <- format_ic_plot(scice_results, threshold = 1.005)

(p1 + p2) / (p3 + p4)
p_ic
