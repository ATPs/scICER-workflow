#!/usr/bin/env Rscript

library(Seurat)
library(ggplot2)
library(qs)
library(scICER)
library(uwot)
library(patchwork)
source("codes/figures/utils.R")

input_file <- "data/generated/iterative/ZhengMix_4757/seurat.scICER.qs"
if (!file.exists(input_file)) {
  stop("Missing input file: ", input_file, "
See data/README.md for preparation details.")
}

sample_obj <- qs::qread(input_file)
sub_obj <- subset(sample_obj, subset = clusters_6 == 4)
sub_obj@meta.data <- sub_obj@meta.data[, !grepl("^clusters_", colnames(sub_obj@meta.data)), drop = FALSE]
if (!"pca" %in% Reductions(sub_obj)) {
  sub_obj <- RunPCA(sub_obj, npcs = 30, verbose = FALSE)
}

umap_result <- uwot::umap(Embeddings(sub_obj, reduction = "pca"), ret_model = TRUE, ret_extra = c("fgraph"), min_dist = 0.1, metric = "cosine", n_neighbors = 15L, n_components = 2L)
sub_obj[["umap2"]] <- CreateDimReducObject(embeddings = umap_result$embedding, key = "UMAP_", assay = DefaultAssay(sub_obj))
rownames(umap_result$fgraph) <- rownames(sub_obj@meta.data)
colnames(umap_result$fgraph) <- rownames(sub_obj@meta.data)
sub_obj[["umap_graph"]] <- as.Graph(umap_result$fgraph)

scice_sub <- scICE_clustering(object = sub_obj, cluster_range = 2:20, remove_threshold = Inf, n_workers = 8, n_trials = 30, n_bootstrap = 200, seed = 123, verbose = TRUE, graph_name = "umap_graph")
scice_main <- scICE_clustering(object = sample_obj, cluster_range = 2:20, remove_threshold = 1.005, n_workers = 8, n_trials = 30, n_bootstrap = 200, seed = 123, graph_name = "umap_graph", verbose = TRUE)

format_ic_plot(scice_main, threshold = 1.005, y_limits = c(NA, NA)) /
  format_ic_plot(scice_sub, threshold = 1.011, y_limits = c(NA, NA))
