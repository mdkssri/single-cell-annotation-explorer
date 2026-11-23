suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(SingleR)
})

bundle <- readRDS("data-raw/singler_results.rds")
prediction <- bundle$prediction
query <- bundle$query
split <- read.delim("data-raw/split_cells.tsv", check.names = FALSE)
query_split <- split[split$role == "query", ]
query_split <- query_split[match(colnames(query), query_split$cell_id), ]
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

cells <- data.frame(
  cell_id = colnames(query),
  patient = query_split$patient,
  sample = query_split$sample,
  provided_annotation = query_split$provided_annotation,
  singleR_label = as.character(prediction$labels),
  pruned_label = as.character(prediction$pruned.labels),
  delta_next = as.numeric(prediction$delta.next),
  delta_median = as.numeric(getDeltaFromMedian(prediction)),
  is_pruned = is.na(prediction$pruned.labels)
)
write.table(cells, "data/processed/cells.tsv", sep = "\t", row.names = FALSE, quote = FALSE, na = "")

score_matrix <- as.matrix(prediction$scores)
scores <- data.frame(
  cell_id = rep(colnames(query), times = ncol(score_matrix)),
  label = rep(colnames(score_matrix), each = nrow(score_matrix)),
  score = as.vector(score_matrix),
  stringsAsFactors = FALSE
)
scores$rank <- ave(-scores$score, scores$cell_id, FUN = function(x) rank(x, ties.method = "first"))
scores <- scores[, c("cell_id", "label", "rank", "score")]
write.table(scores, "data/processed/scores.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

marker_data <- metadata(prediction)$de.genes
marker_rows <- list()
for (label in names(marker_data)) {
  comparisons <- marker_data[[label]]
  for (comparison in names(comparisons)) {
    genes <- comparisons[[comparison]]
    if (length(genes) > 0L) {
      marker_rows[[length(marker_rows) + 1L]] <- data.frame(
        label = label,
        comparison = comparison,
        gene_id = genes,
        marker_rank = seq_along(genes)
      )
    }
  }
}
markers <- if (length(marker_rows) > 0L) do.call(rbind, marker_rows) else
  data.frame(label = character(), comparison = character(), gene_id = character(), marker_rank = integer())
write.table(markers, "data/processed/markers.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

marker_genes <- intersect(unique(markers$gene_id), rownames(query))
expression_matrix <- as.matrix(assay(query, "counts")[marker_genes, , drop = FALSE])
expression <- data.frame(
  cell_id = rep(colnames(query), each = length(marker_genes)),
  gene_id = rep(marker_genes, times = ncol(query)),
  count = as.vector(expression_matrix),
  stringsAsFactors = FALSE
)
write.table(expression, "data/processed/marker_expression.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

reference <- bundle$reference
reference_matrix <- as.matrix(assay(reference, "logcounts")[marker_genes, , drop = FALSE])
reference_labels <- as.character(reference$cell_type)
reference_summary <- do.call(rbind, lapply(sort(unique(reference_labels)), function(label) {
  indices <- which(reference_labels == label)
  data.frame(label = label, gene_id = marker_genes,
             mean_log_expression = rowMeans(reference_matrix[, indices, drop = FALSE]),
             reference_cells = length(indices))
}))
write.table(reference_summary, "data/processed/reference_expression.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

cat("预计算结果已导出到 data/processed。\n")
