#!/usr/bin/env Rscript

library(Seurat)
library(ggplot2)
library(qs)
library(patchwork)

input_file <- "data/generated/iterative/ZhengMix_4757/seurat.scICER.qs"
if (!file.exists(input_file)) {
  stop("Missing input file: ", input_file, "
See data/README.md for preparation details.")
}

sample_obj <- qs::qread(input_file)
celltype_colors <- c("b-cells" = "#00468BFF", "cd14-monocytes" = "#ED0000FF", "cd4-t-helper" = "#EFC000FF", "cd56-nk" = "#0099B4FF", "memory-t" = "#925E9FFF", "naive-cytotoxic" = "#7E6148FF", "naive-t" = "#AD002AFF", "regulatory-t" = "#FDAF91FF")
cluster_colors <- c("1" = "#00468BFF", "2" = "#ED0000FF", "3" = "#925E9FFF", "4" = "#E18727FF", "5" = "#0099B4FF", "6" = "#7E6148FF")
base_theme <- theme(axis.ticks.length = unit(-0.15, "cm"), axis.text = element_text(margin = margin(3, 3, 3, 3)), panel.border = element_rect(color = "black", fill = NA, linewidth = 1))

p1 <- DimPlot(sample_obj, reduction = "umap", group.by = "celltype", pt.size = 2, alpha = 1, raster = TRUE) + scale_color_manual(values = celltype_colors) + base_theme
p2 <- DimPlot(sample_obj, reduction = "umap", group.by = "clusters_6", pt.size = 2, alpha = 1, raster = TRUE) + scale_color_manual(values = cluster_colors) + base_theme

p1 / p2
