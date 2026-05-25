#!/usr/bin/env Rscript

library(sceasy)
library(reticulate)
library(BiocParallel)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript data/preprocessing/convert_h5ad_to_seurat.R <input.h5ad> <output.rds> [conda_env]")
}

input_file <- args[1]
output_file <- args[2]
conda_env <- if (length(args) >= 3) args[3] else Sys.getenv("SCEASY_CONDA_ENV", unset = "")

register(MulticoreParam(workers = 4, progressbar = TRUE))
if (nzchar(conda_env)) {
  reticulate::use_condaenv(conda_env, required = TRUE)
}

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
sceasy::convertFormat(
  obj = input_file,
  from = "anndata",
  to = "seurat",
  outFile = output_file
)
