#!/usr/bin/env Rscript

library(Seurat)
library(qs)
library(scICER)
library(ClustAssess)
library(uwot)
library(patchwork)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript codes/scvi_integration_scICER.R <input.qs> [batch_var] [celltype_var]")
}

input_file <- args[1]
batch_var <- if (length(args) >= 2) args[2] else "tech"
celltype_var <- if (length(args) >= 3) args[3] else "celltype"

score_scice <- function(scice_results, metadata, celltype_var, threshold = 1.005) {
  score_table <- data.frame(
    cluster_number = scice_results$n_cluster,
    ic_score = scice_results$ic,
    is_consistent = scice_results$ic <= threshold
  )

  label_df <- get_robust_labels(scice_results, return_seurat = FALSE, threshold = Inf)
  if (celltype_var %in% colnames(metadata)) {
    label_df$celltype_label <- metadata[[celltype_var]]
    for (k in score_table$cluster_number) {
      cluster_col <- paste0("clusters_", k)
      if (cluster_col %in% names(label_df)) {
        score_table[score_table$cluster_number == k, "ECS_score"] <- element_sim(
          label_df$celltype_label,
          label_df[[cluster_col]]
        )
      }
    }
  }

  consistent <- subset(score_table, is_consistent)
  if (nrow(consistent) == 0) {
    optimal_cluster <- score_table$cluster_number[which.min(score_table$ic_score)]
  } else if ("ECS_score" %in% colnames(consistent) && any(!is.na(consistent$ECS_score))) {
    optimal_cluster <- consistent$cluster_number[which.max(consistent$ECS_score)]
  } else {
    optimal_cluster <- consistent$cluster_number[which.min(consistent$ic_score)]
  }

  list(score_table = score_table, optimal_cluster = optimal_cluster)
}

sample_obj <- qs::qread(input_file)
if (!"scvi" %in% Reductions(sample_obj)) {
  stop("The input Seurat object must contain a reduction named 'scvi'.")
}
if (!batch_var %in% colnames(sample_obj@meta.data)) {
  stop("Missing batch variable in metadata: ", batch_var)
}

n_scvi_dims <- ncol(Embeddings(sample_obj, "scvi"))
sample_obj <- FindNeighbors(
  sample_obj,
  reduction = "scvi",
  dims = seq_len(n_scvi_dims),
  graph.name = c("scvi_nn", "scvi_snn")
)

umap_result <- uwot::umap(
  Embeddings(sample_obj, reduction = "scvi"),
  ret_model = TRUE,
  ret_extra = c("fgraph"),
  min_dist = 0.1,
  metric = "cosine",
  n_neighbors = 15L,
  n_components = 2L
)

sample_obj[["umap_scvi"]] <- CreateDimReducObject(
  embeddings = umap_result$embedding,
  key = "UMAP_",
  assay = DefaultAssay(sample_obj)
)
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

scored <- score_scice(scice_results, sample_obj@meta.data, celltype_var)
sample_obj <- get_robust_labels(scice_results, return_seurat = TRUE, threshold = Inf)
sample_obj$optimal_cluster <- sample_obj@meta.data[[paste0("clusters_", scored$optimal_cluster)]]

base_theme <- theme(
  axis.ticks.length = unit(-0.15, "cm"),
  axis.text = element_text(margin = margin(3, 3, 3, 3)),
  panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
  legend.position = "bottom",
  legend.direction = "vertical",
  legend.box = "vertical",
  panel.grid = element_blank()
)

plots <- list(
  batch_after = DimPlot(sample_obj, reduction = "umap_scvi", group.by = batch_var, pt.size = 3, raster = TRUE) + base_theme,
  cluster_plot = DimPlot(sample_obj, reduction = "umap_scvi", group.by = "optimal_cluster", pt.size = 3, raster = TRUE) + base_theme,
  ic_plot = plot_ic(scice_results, threshold = 1.005)
)

if ("umap" %in% Reductions(sample_obj)) {
  plots$batch_before <- DimPlot(
    sample_obj,
    reduction = "umap",
    group.by = batch_var,
    pt.size = 3,
    raster = TRUE
  ) + base_theme
}

if (celltype_var %in% colnames(sample_obj@meta.data)) {
  plots$celltype_plot <- DimPlot(
    sample_obj,
    reduction = "umap_scvi",
    group.by = celltype_var,
    pt.size = 3,
    raster = TRUE
  ) + base_theme
}

message("Optimal cluster: ", scored$optimal_cluster)
print(scored$score_table)

result <- list(
  seurat_object = sample_obj,
  scice_results = scice_results,
  score_table = scored$score_table,
  optimal_cluster = scored$optimal_cluster,
  plots = plots
)

invisible(result)
