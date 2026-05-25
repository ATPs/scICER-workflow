#!/usr/bin/env Rscript

library(Seurat)
library(ggplot2)
library(qs)
library(scICER)
library(patchwork)
source("codes/figures/utils.R")

input_file <- "data/generated/pancreas_gse221156/pancreatic_harmony.qs"
if (!file.exists(input_file)) {
  stop("Missing input file: ", input_file, "
Prepare it with the scripts listed in data/README.md.")
}

data_obj <- qs::qread(input_file)
scice_results <- scICE_clustering(
  object = data_obj,
  cluster_range = 2:20,
  remove_threshold = Inf,
  n_workers = 8,
  n_trials = 15,
  n_bootstrap = 100,
  seed = 123,
  verbose = TRUE,
  graph_name = "harmony_nn",
  objective_function = "modularity"
)

data_obj <- get_robust_labels(scice_results, return_seurat = TRUE, threshold = Inf)
base_theme <- theme(axis.ticks.length = unit(-0.15, "cm"), axis.text = element_text(margin = margin(3, 3, 3, 3)), panel.border = element_rect(color = "black", fill = NA, linewidth = 1), legend.position = "bottom", legend.direction = "vertical", legend.box = "vertical")

p1 <- DimPlot(data_obj, reduction = "harmony_umap", group.by = "cell_type", pt.size = 3, alpha = 1, raster = TRUE, raster.dpi = c(1200, 1200)) + base_theme
p2 <- format_ic_plot(scice_results, threshold = 1.005)
p3 <- DimPlot(data_obj, reduction = "harmony_umap", group.by = "clusters_12", pt.size = 3, alpha = 1, raster = TRUE, raster.dpi = c(1200, 1200)) + base_theme

p1
p2
p3
