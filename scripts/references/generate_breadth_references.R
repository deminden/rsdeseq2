# Reference breadth cases (hardening 6.6): DESeq() + results() defaults on
# three simulated shapes that the GTEx panels do not cover.
#
#   three_vs_three  2,000 genes x 6 samples (3 vs 3), the smallest common design
#   sparse          1,000 genes x 48 cells, single-cell-like (mostly zeros),
#                   sfType = "poscounts" because every gene has a zero
#   many_samples    400 genes x 500 samples, two groups plus a batch
#
# Usage (from the repository root):
#   Rscript --vanilla scripts/references/generate_breadth_references.R OUTDIR
# Writes OUTDIR/<case>/{counts,design_matrix,reference_results,
# reference_size_factors}.tsv and provenance.json.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
out_root <- args[[1]]
suppressPackageStartupMessages(library(DESeq2))

write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}
fmt <- function(v) ifelse(is.na(v), NA, sprintf("%.17g", v))

simulate <- function(seed, genes, samples, group, mean_log, dispersion, effect, zero_inflation = 0) {
  set.seed(seed)
  base <- exp(rnorm(genes, mean_log, 1.5))
  disp <- dispersion * (1 + 2 / sqrt(base))
  lfc <- ifelse(runif(genes) < 0.1, rnorm(genes, 0, effect), 0)
  sf <- exp(rnorm(samples, 0, 0.25))
  mu <- outer(base, sf) * 2^(outer(lfc, as.numeric(group == levels(group)[2])))
  counts <- matrix(rnbinom(genes * samples, mu = mu, size = 1 / disp), genes, samples)
  if (zero_inflation > 0) {
    counts[runif(length(counts)) < zero_inflation] <- 0L
  }
  storage.mode(counts) <- "integer"
  rownames(counts) <- sprintf("gene%04d", seq_len(genes))
  colnames(counts) <- sprintf("s%03d", seq_len(samples))
  counts
}

run_case <- function(name, counts, col_data, design, sf_type = "ratio") {
  dir <- file.path(out_root, name)
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  keep <- rowSums(counts) > 0
  counts <- counts[keep, , drop = FALSE]
  dds <- DESeqDataSetFromMatrix(counts, col_data, design = design)
  dds <- DESeq(dds, sfType = sf_type, quiet = TRUE)
  res <- results(dds)
  # The replace column exists only when some sample was replaceable.
  replaced <- if (is.null(mcols(dds)$replace)) rep(FALSE, nrow(dds)) else mcols(dds)$replace %in% TRUE
  model <- model.matrix(design, col_data)
  write_tsv(data.frame(gene = rownames(counts), counts, check.names = FALSE),
            file.path(dir, "counts.tsv"))
  write_tsv(data.frame(sample = colnames(counts), model, check.names = FALSE),
            file.path(dir, "design_matrix.tsv"))
  write_tsv(data.frame(
    gene = rownames(res),
    baseMean = fmt(res$baseMean), log2FoldChange = fmt(res$log2FoldChange),
    lfcSE = fmt(res$lfcSE), stat = fmt(res$stat), pvalue = fmt(res$pvalue),
    padj = fmt(res$padj), dispersion = fmt(mcols(dds)$dispersion),
    replaced = as.integer(replaced)
  ), file.path(dir, "reference_results.tsv"))
  write_tsv(data.frame(sample = colnames(dds), size_factor = fmt(sizeFactors(dds))),
            file.path(dir, "reference_size_factors.tsv"))
  writeLines(jsonlite::toJSON(list(
    generator = "scripts/references/generate_breadth_references.R",
    case = name,
    r_version = R.version.string,
    deseq2_version = as.character(packageVersion("DESeq2")),
    platform = R.version$platform,
    genes = nrow(counts), samples = ncol(counts),
    design = paste(deparse(design), collapse = ""),
    result = resultsNames(dds)[length(resultsNames(dds))],
    size_factor_type = sf_type,
    dispersion_fit_type = attr(dispersionFunction(dds), "fitType"),
    replaced_genes = sum(replaced),
    significant_0.1 = sum(res$padj < 0.1, na.rm = TRUE),
    significant_0.05 = sum(res$padj < 0.05, na.rm = TRUE)
  ), auto_unbox = TRUE, pretty = TRUE), file.path(dir, "provenance.json"))
  cat(sprintf("%s: %d genes x %d samples, fit %s, %d replaced, %d significant at 0.1\n",
              name, nrow(counts), ncol(counts), attr(dispersionFunction(dds), "fitType"),
              sum(replaced), sum(res$padj < 0.1, na.rm = TRUE)))
}

group <- factor(rep(c("A", "B"), each = 3))
run_case("three_vs_three",
         simulate(6601, 2000, 6, group, 4, 0.08, 2),
         data.frame(condition = group, row.names = sprintf("s%03d", 1:6)),
         ~ condition)

group <- factor(rep(c("A", "B"), each = 24))
run_case("sparse",
         simulate(6602, 1000, 48, group, 0, 0.5, 2, zero_inflation = 0.4),
         data.frame(condition = group, row.names = sprintf("s%03d", 1:48)),
         ~ condition, sf_type = "poscounts")

group <- factor(rep(c("A", "B"), times = 250))
batch <- factor(rep(c("b1", "b2", "b3", "b4", "b5"), each = 100))
run_case("many_samples",
         simulate(6603, 400, 500, group, 3, 0.15, 1),
         data.frame(batch = batch, condition = group, row.names = sprintf("s%03d", 1:500)),
         ~ batch + condition)
