suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(scuttle)
  library(matrixStats)
})

bundle <- readRDS("data-raw/singler_results.rds")
query <- bundle$query
if (ncol(query) < 3L) stop("计算二维 PCA 至少需要三个细胞。")
normalized <- logNormCounts(query)
logcounts <- as.matrix(assay(normalized, "logcounts"))
variances <- rowVars(logcounts)
eligible <- which(is.finite(variances) & variances > 0)
if (length(eligible) < 2L) stop("可用于 PCA 的变异基因不足。")
ordered <- eligible[order(-variances[eligible], rownames(logcounts)[eligible])]
features <- head(ordered, 500L)
pca <- prcomp(t(logcounts[features, , drop = FALSE]), center = TRUE, scale. = FALSE, rank. = 2L)
varianceFraction <- pca$sdev^2 / sum(variances[features])
coordinates <- data.frame(
  cell_id = colnames(query),
  pc1 = pca$x[, 1],
  pc2 = pca$x[, 2],
  stringsAsFactors = FALSE
)
if (any(!is.finite(coordinates$pc1)) || any(!is.finite(coordinates$pc2))) stop("PCA 坐标不完整。")
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
write.table(coordinates, "data/processed/pca.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
writeLines(c(
  "method=prcomp; query only; centered log-normalized counts; no variance scaling",
  "features=up to 500 genes with highest variance; ties broken by gene identifier",
  paste0("selected_features=", length(features)),
  paste0("pc1_variance_fraction=", format(varianceFraction[1], digits = 12)),
  paste0("pc2_variance_fraction=", format(varianceFraction[2], digits = 12)),
  "purpose=cell overview only; not a SingleR score, model step, or independent validation"
), "data/processed/pca_provenance.txt")
cat("已导出 query PCA 坐标。\n")
