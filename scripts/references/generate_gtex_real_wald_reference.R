#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("usage: Rscript generate_gtex_real_wald_reference.R FIXTURE_DIR")
}
fixture_dir <- normalizePath(args[[1]], mustWork = TRUE)

if (!requireNamespace("DESeq2", quietly = TRUE) ||
    !requireNamespace("S4Vectors", quietly = TRUE) ||
    !requireNamespace("jsonlite", quietly = TRUE) ||
    !requireNamespace("BiocManager", quietly = TRUE)) {
  stop("DESeq2, S4Vectors, jsonlite, and BiocManager are required; use the pinned DESeq2 R 4.6.1 environment")
}
if (!identical(as.character(utils::packageVersion("DESeq2")), "1.52.0")) {
  stop("reference generation requires DESeq2 1.52.0")
}

read_tsv <- function(path) {
  utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE,
                    na.strings = character(), quote = "", comment.char = "")
}
write_tsv <- function(x, path) {
  utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE,
                     col.names = TRUE, na = "NA")
}

counts_table <- read_tsv(file.path(fixture_dir, "counts.tsv"))
design_table <- read_tsv(file.path(fixture_dir, "design_matrix.tsv"))
if (!identical(counts_table[[1]], counts_table[[1]][!duplicated(counts_table[[1]])])) {
  stop("gene IDs in counts.tsv must be unique")
}
samples <- names(counts_table)[-1L]
if (!identical(design_table$sample, samples)) {
  stop("design rows and count columns are not in the same sample order")
}

count_matrix <- as.matrix(counts_table[, -1L, drop = FALSE])
storage.mode(count_matrix) <- "integer"
rownames(count_matrix) <- counts_table[[1]]
col_data <- S4Vectors::DataFrame(row.names = samples)
design_matrix <- as.matrix(design_table[, -1L, drop = FALSE])
storage.mode(design_matrix) <- "double"
rownames(design_matrix) <- samples
colnames(design_matrix) <- names(design_table)[-1L]

dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = col_data,
  design = design_matrix
)
dds <- DESeq2::DESeq(
  dds,
  fitType = "parametric",
  betaPrior = FALSE,
  minReplicatesForReplace = Inf,
  quiet = TRUE
)
result_names <- DESeq2::resultsNames(dds)
coefficient <- result_names[[length(result_names)]]
result <- DESeq2::results(
  dds,
  name = coefficient,
  cooksCutoff = FALSE,
  independentFiltering = FALSE
)

reference <- data.frame(
  gene = rownames(result),
  dispersion = DESeq2::dispersions(dds),
  baseMean = result$baseMean,
  log2FoldChange = result$log2FoldChange,
  lfcSE = result$lfcSE,
  stat = result$stat,
  pvalue = result$pvalue,
  check.names = FALSE
)
write_tsv(reference, file.path(fixture_dir, "reference_results.tsv"))
write_tsv(
  data.frame(sample = samples, sizeFactor = DESeq2::sizeFactors(dds)),
  file.path(fixture_dir, "reference_size_factors.tsv")
)

source <- jsonlite::fromJSON(file.path(fixture_dir, "source.json"))
source$reference <- list(
  R = R.version.string,
  DESeq2 = as.character(utils::packageVersion("DESeq2")),
  Bioconductor = as.character(BiocManager::version()),
  coefficient = coefficient,
  fitType = "parametric",
  betaPrior = FALSE,
  minReplicatesForReplace = "Inf",
  rows_with_finite_dispersion = sum(is.finite(reference$dispersion))
)
jsonlite::write_json(source, file.path(fixture_dir, "source.json"), auto_unbox = TRUE,
                     pretty = TRUE, digits = 17)

cat(sprintf("Wrote DESeq2 reference for %s genes x %s samples (%s)\n",
            nrow(reference), ncol(count_matrix), coefficient))
