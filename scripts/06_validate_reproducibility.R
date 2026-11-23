required <- c("data/processed/cells.tsv", "data/processed/scores.tsv", "data/processed/markers.tsv", "data-raw/singler_results.rds", "data/processed/marker_expression.tsv", "data/processed/reference_expression.tsv")
missing <- required[!file.exists(required)]
if (length(missing) > 0L) stop("缺少导出文件：", paste(missing, collapse = ", "))
cells <- read.delim(required[1], check.names = FALSE, na.strings = c("", "NA"))
scores <- read.delim(required[2], check.names = FALSE)
markers <- read.delim(required[3], check.names = FALSE)
expression <- read.delim(required[5], check.names = FALSE)
referenceExpression <- read.delim(required[6], check.names = FALSE)
bundle <- readRDS(required[4])
prediction <- bundle$prediction
query <- bundle$query
score_matrix <- as.matrix(prediction$scores)
if (anyDuplicated(cells$cell_id)) stop("cells.tsv 中 cell_id 不唯一。")
if (nrow(scores) != nrow(cells) * ncol(score_matrix)) stop("scores.tsv 行数不符合长表预期。")
if (anyDuplicated(scores[c("cell_id", "label")])) stop("scores.tsv 中 cell-label 组合重复。")
if (!all(c("delta_next", "pruned_label", "singleR_label") %in% names(cells))) stop("cells.tsv 缺少 SingleR 字段。")
if (!identical(cells$cell_id, colnames(query))) stop("cells.tsv 的 query cell 顺序与 SingleR 输入不一致。")
expected_cells <- rep(colnames(query), times = ncol(score_matrix))
expected_labels <- rep(colnames(score_matrix), each = nrow(score_matrix))
if (!identical(scores$cell_id, expected_cells)) stop("scores.tsv 的 cell 顺序不匹配 score matrix。")
if (!identical(scores$label, expected_labels)) stop("scores.tsv 的 label 顺序不匹配 score matrix。")
if (any(!is.finite(scores$score))) stop("scores.tsv 含非有限 score。")
expected_scores <- as.vector(score_matrix)
if (!isTRUE(all.equal(scores$score, expected_scores, tolerance = 1e-12))) stop("scores.tsv 的 score 与 SingleR 原始矩阵不一致。")
for (cell_id in unique(scores$cell_id)) {
  rows <- scores$cell_id == cell_id
  if (!identical(sort(scores$rank[rows]), seq_len(ncol(score_matrix)))) stop("rank 不完整：", cell_id)
  if (!identical(scores$rank[rows], rank(-scores$score[rows], ties.method = "first"))) stop("rank 与 score 降序不一致：", cell_id)
}
if (!isTRUE(all.equal(cells$singleR_label, as.character(prediction$labels)))) stop("cells.tsv 的 labels 与 SingleR 原始结果不一致。")
expected_pruned <- as.character(prediction$pruned.labels)
actual_pruned <- as.character(cells$pruned_label)
same_pruned <- (is.na(actual_pruned) & is.na(expected_pruned)) |
  (!is.na(actual_pruned) & !is.na(expected_pruned) & actual_pruned == expected_pruned)
if (!all(same_pruned)) stop("cells.tsv 的 pruned labels 与 SingleR 原始结果不一致。")
if (!isTRUE(all.equal(cells$delta_next, as.numeric(prediction$delta.next), tolerance = 1e-12))) stop("cells.tsv 的 delta_next 与 SingleR 原始结果不一致。")
if (!isTRUE(all.equal(cells$delta_median, as.numeric(SingleR::getDeltaFromMedian(prediction)), tolerance = 1e-12))) stop("cells.tsv 的 delta_median 与 SingleR 结果不一致。")
if (!isTRUE(all.equal(cells$is_pruned, is.na(prediction$pruned.labels)))) stop("cells.tsv 的 pruning 标志不一致。")
if (nrow(markers) > 0L && anyDuplicated(markers[c("label", "comparison", "gene_id", "marker_rank")])) stop("markers.tsv 存在重复 marker 记录。")
if (!all(c("label", "comparison", "gene_id", "marker_rank") %in% names(markers))) stop("markers.tsv 缺少成对 marker 字段。")
marker_data <- S4Vectors::metadata(prediction)$de.genes
expected_markers <- list()
for (label in names(marker_data)) {
  for (comparison in names(marker_data[[label]])) {
    genes <- marker_data[[label]][[comparison]]
    if (length(genes) > 0L) expected_markers[[length(expected_markers) + 1L]] <- data.frame(
      label = label, comparison = comparison, gene_id = genes, marker_rank = seq_along(genes)
    )
  }
}
if (length(expected_markers) == 0L) stop("SingleR 没有可验证的 marker。")
expected_markers <- do.call(rbind, expected_markers)
rownames(expected_markers) <- NULL
if (!isTRUE(all.equal(markers, expected_markers, check.attributes = FALSE))) stop("markers.tsv 与 SingleR marker metadata 不一致。")
if (file.exists("data/processed/pca.tsv")) {
  pca <- read.delim("data/processed/pca.tsv", check.names = FALSE)
  if (!identical(pca$cell_id, colnames(query)) || any(!is.finite(pca$pc1)) || any(!is.finite(pca$pc2))) {
    stop("PCA 坐标与 query 细胞不匹配或包含非有限值。")
  }
}
genes <- intersect(unique(markers$gene_id), rownames(query))
expected_expression <- as.matrix(SummarizedExperiment::assay(query, "counts")[genes, , drop = FALSE])
if (!identical(expression$cell_id, rep(colnames(query), each = length(genes))) ||
    !identical(expression$gene_id, rep(genes, times = ncol(query))) ||
    !isTRUE(all.equal(expression$count, as.vector(expected_expression), check.attributes = FALSE))) {
  stop("marker_expression.tsv 与 query counts 不一致。")
}
reference <- bundle$reference
referenceMatrix <- as.matrix(SummarizedExperiment::assay(reference, "logcounts")[genes, , drop = FALSE])
referenceLabels <- as.character(reference$cell_type)
expectedReference <- do.call(rbind, lapply(sort(unique(referenceLabels)), function(label) {
  indices <- which(referenceLabels == label)
  data.frame(label = label, gene_id = genes, mean_log_expression = rowMeans(referenceMatrix[, indices, drop = FALSE]), reference_cells = length(indices))
}))
rownames(expectedReference) <- NULL
if (!isTRUE(all.equal(referenceExpression, expectedReference, tolerance = 1e-12, check.attributes = FALSE))) {
  stop("reference_expression.tsv 与参考表达矩阵不一致。")
}
dir.create("validation", showWarnings = FALSE)
writeLines(c(
  "export_against_saved_run=PASS",
  paste0("query_cells=", nrow(cells)),
  paste0("labels=", ncol(score_matrix)),
  paste0("scores_rows=", nrow(scores)),
  paste0("marker_rows=", nrow(markers)),
  paste0("pruned_count=", sum(cells$is_pruned)),
  "score_stage=pre_finetuning_initial_scores",
  "marker_semantics=pairwise_label_comparison_de_genes"
), "validation/reproducibility_validation.txt")
cat("导出与已保存的 SingleR 结果一致；此检查未重新运行 SingleR。\n")
