#!/usr/bin/env Rscript

library(Seurat)
library(qs)
library(harmony)
library(biomaRt)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript data/preprocessing/prepare_pancreatic_harmony.R <input.rds> <output.qs>")
}

input_file <- args[1]
output_file <- args[2]

pancreas_obj <- readRDS(input_file)

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_annotations <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  mart = ensembl
)
ensembl_to_symbol <- setNames(gene_annotations$hgnc_symbol, gene_annotations$ensembl_gene_id)

current_ids <- rownames(pancreas_obj)
new_symbols <- ensembl_to_symbol[current_ids]
unmatched <- is.na(new_symbols) | new_symbols == ""
new_symbols[unmatched] <- current_ids[unmatched]

count_data <- LayerData(pancreas_obj, layer = "counts")
rownames(count_data) <- new_symbols
pancreas_obj[["RNA"]] <- CreateAssayObject(counts = count_data, min.cells = 0, min.features = 0)

pancreas_obj <- NormalizeData(pancreas_obj)
pancreas_obj <- FindVariableFeatures(pancreas_obj, selection.method = "vst", nfeatures = 3000)
pancreas_obj <- ScaleData(pancreas_obj)
pancreas_obj <- RunPCA(pancreas_obj)
pancreas_obj <- RunHarmony(
  object = pancreas_obj,
  group.by.vars = c("assay", "donor_id"),
  reduction = "pca",
  dims.use = seq_len(ncol(Embeddings(pancreas_obj, "pca"))),
  assay.use = "RNA",
  project.dim = FALSE,
  verbose = TRUE,
  reduction.save = "harmony_pca"
)
pancreas_obj <- RunUMAP(
  object = pancreas_obj,
  dims = seq_len(ncol(Embeddings(pancreas_obj, "harmony_pca"))),
  reduction = "harmony_pca",
  reduction.name = "harmony_umap",
  verbose = TRUE
)
pancreas_obj <- FindNeighbors(
  object = pancreas_obj,
  reduction = "harmony_pca",
  dims = seq_len(ncol(Embeddings(pancreas_obj, "harmony_pca"))),
  graph.name = c("harmony_nn", "harmony_snn")
)

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
qs::qsave(pancreas_obj, output_file)
