#!/usr/bin/env Rscript

library(Seurat)
library(ggplot2)
library(qs)
library(scICER)
library(ClustAssess)
library(uwot)
library(patchwork)
source("codes/figures/utils.R")

input_file <- "data/generated/human_pancreas_norm_complexBatch/seurat.scvi.HVG.qs"
if (!file.exists(input_file)) {
  stop("Missing input file: ", input_file, "
See data/README.md for preparation details.")
}

sample_obj <- qs::qread(input_file)
n_scvi_dims <- ncol(sample_obj@reductions$scvi)

sample_obj <- FindNeighbors(sample_obj, reduction = "scvi", dims = seq_len(n_scvi_dims), graph.name = c("scvi_nn", "scvi_snn"))

umap_result <- uwot::umap(
  Embeddings(sample_obj, reduction = "scvi"),
  ret_model = TRUE,
  ret_extra = c("fgraph"),
  min_dist = 0.1,
  metric = "cosine",
  n_neighbors = 15L,
  n_components = 2L
)
sample_obj[["umap_scvi_umapR"]] <- CreateDimReducObject(embeddings = umap_result$embedding, key = "UMAP_", assay = DefaultAssay(sample_obj))
rownames(umap_result$fgraph) <- rownames(sample_obj@meta.data)
colnames(umap_result$fgraph) <- rownames(sample_obj@meta.data)
sample_obj[["scvi_umap_graph"]] <- as.Graph(umap_result$fgraph)

scice_results <- scICE_clustering(
  object = sample_obj,
  cluster_range = 2:20,
  remove_threshold = 1.005,
  n_workers = 8,
  n_trials = 30,
  n_bootstrap = 200,
  seed = 123,
  verbose = TRUE,
  graph_name = "scvi_umap_graph"
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

p1 <- DimPlot(sample_obj, reduction = "umap", group.by = "tech", pt.size = 3, alpha = 1, raster = TRUE) + base_theme
p2 <- DimPlot(sample_obj, reduction = "umap_scvi_umapR", group.by = "tech", pt.size = 3, alpha = 1, raster = TRUE) + base_theme
p3 <- DimPlot(sample_obj, reduction = "umap_scvi_umapR", group.by = "celltype", pt.size = 3, alpha = 1, raster = TRUE) + base_theme
p4 <- DimPlot(sample_obj, reduction = "umap_scvi_umapR", group.by = "optimal_cluster", pt.size = 3, alpha = 1, raster = TRUE) + base_theme
p_ic <- format_ic_plot(scice_results, threshold = 1.005, y_limits = c(NA, NA))

(p1 + p2) / (p3 + p4)
p_ic
