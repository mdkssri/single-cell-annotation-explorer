runTenX <-
function (project, input, out)
{
    suppressPackageStartupMessages({
        library(SingleCellExperiment)
        library(SingleR)
        library(scuttle)
    })
    stopifnot(!dir.exists(out))
    dir.create(out)
    tmp <- tempfile()
    dir.create(tmp)
    on.exit(unlink(tmp, recursive = TRUE))
    csv <- "cell_types/cell_types.csv"
    untar(file.path(input, "donor1-cell-types.tar.gz"), files = csv,
        exdir = tmp, tar = "internal")
    ann <- read.csv(file.path(tmp, csv))
    readCounts <- function(file) {
        x <- Seurat::Read10X_h5(file.path(input, file), use.names = FALSE,
            unique.features = FALSE)
        if (is.list(x))
            x <- x[["Gene Expression"]]
        stopifnot(inherits(x, "sparseMatrix"), !anyDuplicated(rownames(x)),
            !anyDuplicated(colnames(x)), Matrix::nnzero(x) >
                0)
        x
    }
    ref <- readCounts("donor1-counts.h5")
    test <- readCounts("donor2-counts.h5")
    shapes <- c(paste(dim(ref), collapse = " x "), paste(dim(test),
        collapse = " x "))
    stopifnot(!anyDuplicated(ann[["barcode"]]), all(ann[["barcode"]] %in%
        colnames(ref)))
    matched <- match(colnames(ref), ann[["barcode"]])
    labs <- trimws(ann[["coarse_cell_type"]][matched])
    group <- split(which(!is.na(labs) & nzchar(labs) & tolower(labs) !=
        "cell"), labs[!is.na(labs) & nzchar(labs) & tolower(labs) !=
        "cell"])
    group <- group[lengths(group) >= 30L]
    stopifnot(length(group) >= 3L)
    set.seed(6152026, kind = "Mersenne-Twister", normal.kind = "Inversion",
        sample.kind = "Rejection")
    ri <- sort(unlist(lapply(group, function(i) if (length(i) >
        180L)
        sample(i, 180L)
    else i), use.names = FALSE))
    qi <- sort(sample(seq_len(ncol(test)), 900L))
    genes <- intersect(rownames(ref), rownames(test))
    stopifnot(length(genes) > 1000L)
    rawIds <- genes
    features <- rhdf5::h5read(file.path(input, "donor1-counts.h5"),
        "matrix/features/id")
    symbols <- rhdf5::h5read(file.path(input, "donor1-counts.h5"),
        "matrix/features/name")[match(genes, features)]
    namesOk <- !is.na(symbols) & nzchar(symbols) & !duplicated(symbols) &
        !duplicated(symbols, fromLast = TRUE) & !symbols %in%
        genes
    genes[namesOk] <- symbols[namesOk]
    stopifnot(!anyDuplicated(genes))
    ref <- ref[rawIds, ri, drop = FALSE]
    test <- test[rawIds, qi, drop = FALSE]
    rownames(ref) <- rownames(test) <- genes
    colnames(ref) <- sprintf("R%04d", seq_along(ri))
    colnames(test) <- sprintf("Q%04d", seq_along(qi))
    reference <- logNormCounts(SingleCellExperiment(list(counts = ref),
        colData = S4Vectors::DataFrame(cell_type = labs[ri])))
    query <- SingleCellExperiment(list(counts = test))
    classify <- function() SingleR(query, reference, labels = reference[["cell_type"]],
        assay.type.test = "counts", assay.type.ref = "logcounts",
        de.method = "wilcox", quantile = 0.80000000000000004,
        fine.tune = TRUE, prune = TRUE)
    prediction <- classify()
    again <- classify()
    fields <- c("labels", "pruned.labels", "delta.next", "scores")
    stopifnot(all(vapply(fields, function(n) isTRUE(all.equal(prediction[[n]],
        again[[n]], tolerance = 9.9999999999999998e-13)), logical(1))),
        identical(S4Vectors::metadata(prediction)[["de.genes"]],
            S4Vectors::metadata(again)[["de.genes"]]), all(is.finite(prediction[["scores"]])),
        sd(as.numeric(prediction[["scores"]])) > 0)
    old <- getwd()
    on.exit(setwd(old), add = TRUE)
    setwd(out)
    dir.create("data-raw")
    write.table(data.frame(cell_id = colnames(query), role = "query",
        patient = "Donor 2", sample = "10x GEM-X", provided_annotation = NA_character_),
        "data-raw/split_cells.tsv", sep = intToUtf8(9), quote = FALSE,
        row.names = FALSE)
    saveRDS(list(prediction = prediction, reference = reference,
        query = query, genes = genes), "data-raw/singler_results.rds")
    for (file in c("05_export_precomputed.R", "09_export_pca.R",
        "06_validate_reproducibility.R")) source(file.path(project,
        "scripts", file), local = new.env())
    report <- c(paste("Original gene x cell dimensions:", paste(shapes,
        collapse = "; ")), paste("Annotated barcodes matched:",
        sum(!is.na(matched)), "of", length(matched)), paste("Reference cells:",
        length(ri)), paste("Query cells:", length(qi)), paste("Shared genes:",
        length(genes)), paste("Reference labels:", paste(names(table(labs[ri])),
        table(labs[ri]), collapse = "; ")), paste("Pruned:",
        sum(is.na(prediction[["pruned.labels"]]))), paste("Minimum delta.next:",
        min(prediction[["delta.next"]])), paste("Highest initial label differs from final:",
        sum(colnames(prediction[["scores"]])[max.col(prediction[["scores"]],
            ties.method = "first")] != prediction[["labels"]])),
        paste("SingleR version:", packageVersion("SingleR")),
        "Two computations and export validation: PASS")
    writeLines(report, "validation/10x-run.txt")
    writeLines(capture.output(sessionInfo()), "validation/10x-session.txt")
    cat(paste(report, collapse = intToUtf8(10)), intToUtf8(10))
}
args <- commandArgs(TRUE); stopifnot(length(args) == 3L); do.call(runTenX, list(normalizePath(args[[1]], mustWork = TRUE), normalizePath(args[[2]], mustWork = TRUE), args[[3]]))
