#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1, unname(tools::md5sum(args)) == "9945e1b94d62375d4b845ac8cc2bbae9")
counts <- Seurat::Read10X_h5(args[[1]], use.names = TRUE, unique.features = TRUE)
if (is.list(counts)) counts <- counts[["Gene Expression"]]
set.seed(6152026)
indices <- sort(sample(which(Matrix::colSums(counts) > 0), 256))
subset <- counts[, indices, drop = FALSE]
totals <- Matrix::colSums(subset)
scaled <- subset %*% Matrix::Diagonal(x = 10000 / totals)
scaled@x <- log1p(scaled@x)
means <- Matrix::rowMeans(scaled)
variance <- Matrix::rowMeans(scaled ^ 2) - means ^ 2
markers <- c("MS4A1", "CD79A", "CD3D", "CD3E", "NKG7", "GNLY", "LST1", "S100A9")
stopifnot(all(markers %in% rownames(subset)))
genes <- head(unique(c(markers, rownames(subset)[order(variance, decreasing = TRUE)])), 160)
payload <- list(
  genes = genes,
  markers = markers,
  totals = unname(totals),
  cells = unname(t(as.matrix(subset[genes, , drop = FALSE])))
)
writeLines(paste0("const PBMC = ", jsonlite::toJSON(payload, matrix = "rowmajor"), ";"), "dist/pbmc.js")
