# DESeq(betaPrior = TRUE) references for `Analysis::beta_prior` (hardening 4.3d).
#
# Reuses the counts of two breadth cases (tests/data/breadth) and rebuilds
# their sample tables from the committed design matrices:
#
#   three_vs_three  ~ condition           (A, B)
#   many_samples    ~ batch + condition   (b1..b5; A, B)
#
# For each case runs DESeq(dds, betaPrior = TRUE, minReplicatesForReplace = Inf)
# (no Cook's replacement, which `run` does not combine with the beta prior)
# with the defaults otherwise and results(dds) (DESeq2 picks the last condition
# level against the first through an expanded-model contrast), and writes
# <case>/reference_results_beta_prior.tsv and
# <case>/reference_beta_prior_variance.tsv.
#
# Usage (from the repository root):
#   Rscript --vanilla scripts/references/generate_beta_prior_references.R \
#     crates/rsdeseq2/tests/data/breadth
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
root <- args[[1]]
suppressPackageStartupMessages(library(DESeq2))

read_tsv <- function(path) {
  read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
}
fmt <- function(v) ifelse(is.na(v), NA, sprintf("%.17g", v))

sample_table <- function(design) {
  condition <- factor(ifelse(design$conditionB == 1, "B", "A"), levels = c("A", "B"))
  table <- data.frame(condition = condition, row.names = design$sample)
  batch_columns <- grep("^batch", names(design), value = TRUE)
  if (length(batch_columns) > 0) {
    batch <- rep("b1", nrow(design))
    for (column in batch_columns) {
      batch[design[[column]] == 1] <- sub("^batch", "", column)
    }
    table$batch <- factor(batch, levels = c("b1", sub("^batch", "", batch_columns)))
  }
  table
}

run_case <- function(name, formula) {
  dir <- file.path(root, name)
  counts <- read_tsv(file.path(dir, "counts.tsv"))
  matrix <- as.matrix(counts[, -1])
  rownames(matrix) <- counts$gene
  storage.mode(matrix) <- "integer"
  col_data <- sample_table(read_tsv(file.path(dir, "design_matrix.tsv")))
  dds <- DESeqDataSetFromMatrix(matrix, col_data, design = formula)
  dds <- DESeq(dds, betaPrior = TRUE, minReplicatesForReplace = Inf, quiet = TRUE)
  res <- results(dds)
  write.table(data.frame(
    gene = rownames(res),
    baseMean = fmt(res$baseMean), log2FoldChange = fmt(res$log2FoldChange),
    lfcSE = fmt(res$lfcSE), stat = fmt(res$stat), pvalue = fmt(res$pvalue),
    padj = fmt(res$padj), dispersion = fmt(mcols(dds)$dispersion)
  ), file.path(dir, "reference_results_beta_prior.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  prior <- attr(dds, "betaPriorVar")
  write.table(data.frame(coefficient = names(prior), variance = fmt(prior)),
              file.path(dir, "reference_beta_prior_variance.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  message(name, ": ", paste(resultsNames(dds), collapse = ", "),
          "; betaPriorVar = ", paste(signif(attr(dds, "betaPriorVar"), 6), collapse = ", "))
}

run_case("three_vs_three", ~ condition)
run_case("many_samples", ~ batch + condition)
