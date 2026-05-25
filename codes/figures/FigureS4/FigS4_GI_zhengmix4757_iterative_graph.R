#!/usr/bin/env Rscript

library(Seurat)
library(ggplot2)
library(qs)
library(scICER)
library(uwot)
library(igraph)
library(ggraph)
library(ggrastr)
library(patchwork)

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
sub_obj <- get_robust_labels(scice_sub, return_seurat = TRUE, threshold = Inf)

p_umap <- DimPlot(sub_obj, reduction = "umap2", group.by = "clusters_2", pt.size = 3, alpha = 1, raster = TRUE) + scale_color_manual(values = c("1" = "#AD002AFF", "2" = "#FDAF91FF")) + theme(axis.ticks.length = unit(-0.15, "cm"), axis.text = element_text(margin = margin(3, 3, 3, 3)), panel.border = element_rect(color = "black", fill = NA, linewidth = 1))

graph_obj <- graph_from_adjacency_matrix(sub_obj@graphs$umap_graph, mode = "undirected", weighted = TRUE)
layout_df <- as.data.frame(Embeddings(sub_obj, "umap2"))
layout_df$cell <- rownames(layout_df)
V(graph_obj)$umap_1 <- layout_df[match(V(graph_obj)$name, layout_df$cell), 1]
V(graph_obj)$umap_2 <- layout_df[match(V(graph_obj)$name, layout_df$cell), 2]

p_graph <- ggraph(graph_obj, layout = "manual", x = V(graph_obj)$umap_1, y = V(graph_obj)$umap_2) +
  ggrastr::rasterise(geom_edge_link(alpha = 0.05, color = "gray50", linewidth = 0.2), dpi = 300) +
  ggrastr::rasterise(geom_node_point(size = 0.01, color = "black"), dpi = 300) +
  theme_void() +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1), plot.margin = margin(5, 5, 5, 5))

p_umap / p_graph
