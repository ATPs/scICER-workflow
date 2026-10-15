#!/usr/bin/env Rscript

library(Seurat)
library(ggplot2)
library(qs)
library(scICER)
library(patchwork)

harmony_file <- "data/generated/human_pancreas_norm_complexBatch/seurat.sclens.qs"
scvi_file <- "data/generated/human_pancreas_norm_complexBatch/seurat.scvi.HVG.qs"

harmony_results_file <- "data/generated/human_pancreas_norm_complexBatch/harmony_scice_results.qs"
scvi_results_file <- "data/generated/human_pancreas_norm_complexBatch/scvi_scice_results.qs"

# Harmony
sample_obj <- qs::qread(harmony_file)
scice_results <- qs::qread(harmony_results_file)
sample_obj <- get_robust_labels(scice_results, return_seurat = TRUE, object = sample_obj, threshold = Inf)

tab <- table(sample_obj$celltype, sample_obj$clusters_9)
dat <- as.data.frame(prop.table(tab, 1) * 100)
names(dat) <- c("celltype", "cluster", "pct")

n <- rowSums(tab)
dat$celltype <- factor(dat$celltype, levels = rev(rownames(tab)),
                       labels = paste0(rev(rownames(tab)), " (n=", rev(n), ")"))
dat$cluster <- factor(paste0("C", dat$cluster),
                      levels = paste0("C", colnames(tab)))

p1 <- ggplot(dat, aes(cluster, celltype, fill = pct)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = ifelse(pct >= 5, sprintf("%.1f%%", pct), "")),
            size = 12 / ggplot2::.pt, color = "black") +
  scale_fill_gradient(low = "#F5F7F4", high = "#167C79",
                      limits = c(0, 100), name = "Proportion (%)") +
  labs(x = "scICER cluster (K=9)", y = "Annotated cell type",
       title = "Harmony: K=9 cell-type correspondence") +
  theme_minimal(base_size = 12) +
  theme(
    text = element_text(size = 12, color = "black"),
    axis.text = element_text(size = 12, color = "black"),
    legend.text = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 12, color = "black"),
    plot.title = element_text(size = 12, color = "black"),
    panel.grid = element_blank()
  )

# scVI
sample_obj <- qs::qread(scvi_file)
scice_results <- qs::qread(scvi_results_file)
sample_obj <- get_robust_labels(scice_results, return_seurat = TRUE, object = sample_obj, threshold = Inf)

tab <- table(sample_obj$celltype, sample_obj$clusters_9)
dat <- as.data.frame(prop.table(tab, 1) * 100)
names(dat) <- c("celltype", "cluster", "pct")

n <- rowSums(tab)
dat$celltype <- factor(dat$celltype, levels = rev(rownames(tab)),
                       labels = paste0(rev(rownames(tab)), " (n=", rev(n), ")"))
dat$cluster <- factor(paste0("C", dat$cluster),
                      levels = paste0("C", colnames(tab)))

p2 <- ggplot(dat, aes(cluster, celltype, fill = pct)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = ifelse(pct >= 5, sprintf("%.1f%%", pct), "")),
            size = 12 / ggplot2::.pt, color = "black") +
  scale_fill_gradient(low = "#F5F7F4", high = "#167C79",
                      limits = c(0, 100), name = "Proportion (%)") +
  labs(x = "scICER cluster (K=9)", y = "Annotated cell type",
       title = "scVI: K=9 cell-type correspondence") +
  theme_minimal(base_size = 12) +
  theme(
    text = element_text(size = 12, color = "black"),
    axis.text = element_text(size = 12, color = "black"),
    legend.text = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 12, color = "black"),
    plot.title = element_text(size = 12, color = "black"),
    panel.grid = element_blank()
  )

p1
p2