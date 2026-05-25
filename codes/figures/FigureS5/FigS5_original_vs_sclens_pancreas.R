#!/usr/bin/env Rscript

library(Seurat)
library(ggplot2)
library(qs)
library(patchwork)

batch_rds <- "data/generated/no_integration_pancreas/scICE_batch_data.rds"
celltype_rds <- "data/generated/no_integration_pancreas/scICE_celltype_data.rds"
reference_qs <- "data/generated/human_pancreas_norm_complexBatch/seurat.sclens.qs"

required_files <- c(batch_rds, celltype_rds, reference_qs)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing input files:
", paste(missing_files, collapse = "
"), "
See data/README.md and data/preprocessing/convert_h5ad_to_seurat.R.")
}

batch_obj <- readRDS(batch_rds)
celltype_obj <- readRDS(celltype_rds)
reference_obj <- qs::qread(reference_qs)

tech_colors <- c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF", "#91D1C2FF", "#DC0000FF", "#7E6148FF")
celltype_colors <- c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF", "#91D1C2FF", "#DC0000FF", "#7E6148FF", "#925E9FFF", "#95CC5EFF", "#AD002AFF", "#ADB6B6FF", "#F7C530FF")
base_theme <- theme(axis.ticks.length = unit(-0.15, "cm"), axis.text = element_text(margin = margin(3, 3, 3, 3)), panel.border = element_rect(color = "black", fill = NA, linewidth = 1))

p1 <- DimPlot(batch_obj, reduction = "umap", group.by = "cell_id", pt.size = 3, alpha = 1, raster = TRUE, raster.dpi = c(1200, 1200)) + scale_color_manual(values = tech_colors) + base_theme
p2 <- DimPlot(celltype_obj, reduction = "umap", group.by = "cell_id", pt.size = 3, alpha = 1, raster = TRUE, raster.dpi = c(1200, 1200)) + scale_color_manual(values = celltype_colors) + base_theme
p3 <- DimPlot(reference_obj, reduction = "umap", group.by = "celltype", pt.size = 3, alpha = 1, raster = TRUE, raster.dpi = c(1200, 1200)) + scale_color_manual(values = celltype_colors) + base_theme
p4 <- DimPlot(reference_obj, reduction = "umap", group.by = "tech", pt.size = 3, alpha = 1, raster = TRUE, raster.dpi = c(1200, 1200)) + scale_color_manual(values = tech_colors) + base_theme

(p1 / p2) / (p3 / p4)
