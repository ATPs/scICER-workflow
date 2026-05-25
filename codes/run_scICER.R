#!/usr/bin/env Rscript

library(Seurat)
library(qs)
library(scICER)
library(ClustAssess)
library(uwot)
library(SeuratDisk)
library(reticulate)
library(patchwork)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript codes/run_scICER.R <input_file> <RNA|SCT|SCLENS> [graph_name]")
}

input_file <- args[1]
method <- toupper(args[2])
graph_name <- if (length(args) >= 3) args[3] else NULL

build_umap_graph <- function(object, reduction = "pca") {
  embedding <- Embeddings(object, reduction = reduction)
  umap_result <- uwot::umap(
    embedding,
    ret_model = TRUE,
    ret_extra = c("fgraph"),
    min_dist = 0.1,
    metric = "cosine",
    n_neighbors = 15L,
    n_components = 2L
  )

  object[["umap"]] <- CreateDimReducObject(
    embeddings = umap_result$embedding,
    key = "UMAP_",
    assay = DefaultAssay(object)
  )

  rownames(umap_result$fgraph) <- rownames(object@meta.data)
  colnames(umap_result$fgraph) <- rownames(object@meta.data)
  object[["umap_graph"]] <- as.Graph(umap_result$fgraph)
  object
}

score_scice <- function(scice_results, metadata, threshold = 1.005) {
  score_table <- data.frame(
    cluster_number = scice_results$n_cluster,
    ic_score = scice_results$ic,
    is_consistent = scice_results$ic <= threshold
  )

  label_df <- get_robust_labels(scice_results, return_seurat = FALSE, threshold = Inf)
  if ("celltype" %in% colnames(metadata)) {
    label_df$celltype <- metadata$celltype
    for (k in score_table$cluster_number) {
      cluster_col <- paste0("clusters_", k)
      if (cluster_col %in% names(label_df)) {
        score_table[score_table$cluster_number == k, "ECS_score"] <- element_sim(
          label_df$celltype,
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

if (method %in% c("RNA", "SCT")) {
  if (grepl("\\.qs$", input_file, ignore.case = TRUE)) {
    sample_obj <- qs::qread(input_file)
  } else if (grepl("\\.csv(\\.gz)?$", input_file, ignore.case = TRUE)) {
    expr_matrix <- read.csv(input_file, row.names = 1, check.names = FALSE)
    sample_obj <- CreateSeuratObject(counts = expr_matrix, project = tools::file_path_sans_ext(basename(input_file)))
  } else {
    stop("RNA and SCT examples expect a .qs, .csv, or .csv.gz input.")
  }

  if (method == "RNA") {
    sample_obj <- NormalizeData(sample_obj)
    sample_obj <- FindVariableFeatures(sample_obj)
    sample_obj <- ScaleData(sample_obj)
  } else {
    sample_obj <- PercentageFeatureSet(sample_obj, pattern = "^MT-", col.name = "percent.mt")
    sample_obj <- SCTransform(sample_obj, vars.to.regress = "percent.mt", verbose = FALSE)
  }

  sample_obj <- RunPCA(sample_obj)
  dims_to_use <- seq_len(min(30, ncol(Embeddings(sample_obj, "pca"))))
  sample_obj <- RunUMAP(sample_obj, dims = dims_to_use)
  sample_obj <- FindNeighbors(sample_obj, dims = dims_to_use)
  sample_obj <- build_umap_graph(sample_obj, reduction = "pca")
  default_graph <- if (method == "RNA") "RNA_snn" else "SCT_snn"
} else if (method == "SCLENS") {
  if (!grepl("\\.h5ad$", input_file, ignore.case = TRUE)) {
    stop("The SCLENS example expects an .h5ad input file.")
  }

  python_bin <- Sys.getenv("SCLENS_PYTHON", unset = "")
  if (nzchar(python_bin)) {
    reticulate::use_python(python_bin, required = TRUE)
  }

  SeuratDisk::Convert(input_file, dest = "h5seurat", overwrite = TRUE)
  h5seurat_file <- sub("\\.h5ad$", ".h5seurat", input_file, ignore.case = TRUE)
  sample_obj <- LoadH5Seurat(h5seurat_file, meta.data = FALSE, misc = TRUE)

  scanpy <- reticulate::import("scanpy")
  adata <- scanpy$read_h5ad(input_file)
  sample_obj <- AddMetaData(sample_obj, metadata = py_to_r(adata$obs))

  dims_to_use <- seq_len(ncol(Embeddings(sample_obj, "pca")))
  sample_obj <- FindNeighbors(sample_obj, reduction = "pca", dims = dims_to_use)
  sample_obj <- build_umap_graph(sample_obj, reduction = "pca")
  default_graph <- "RNA_snn"
} else {
  stop("Method must be one of RNA, SCT, or SCLENS.")
}

graph_name <- if (is.null(graph_name)) default_graph else graph_name
if (!graph_name %in% names(sample_obj@graphs)) {
  stop("Requested graph not found. Available graphs: ", paste(names(sample_obj@graphs), collapse = ", "))
}

scice_results <- scICE_clustering(
  object = sample_obj,
  cluster_range = 2:20,
  remove_threshold = 1.005,
  n_workers = 8,
  n_trials = 15,
  n_bootstrap = 100,
  seed = 123,
  verbose = TRUE,
  graph_name = graph_name
)

scored <- score_scice(scice_results, sample_obj@meta.data)
sample_obj <- get_robust_labels(scice_results, return_seurat = TRUE, threshold = Inf)
sample_obj$optimal_cluster <- sample_obj@meta.data[[paste0("clusters_", scored$optimal_cluster)]]

plots <- list(
  ic_plot = plot_ic(scice_results, threshold = 1.005),
  cluster_plot = DimPlot(sample_obj, reduction = "umap", group.by = "optimal_cluster", label = TRUE)
)

if ("celltype" %in% colnames(sample_obj@meta.data)) {
  plots$comparison_plot <- DimPlot(sample_obj, reduction = "umap", group.by = "celltype", label = TRUE) +
    DimPlot(sample_obj, reduction = "umap", group.by = "optimal_cluster", label = TRUE)
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
