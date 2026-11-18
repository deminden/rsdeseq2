# End-to-end DESeq2 reference on a GTEx kidney-cortex panel (hardening 6.5).
#
# Unlike the compact Wald panel, nothing is supplied: size factors,
# dispersions (parametric trend), Cook's outlier replacement and independent
# filtering all run with DESeq() and results() defaults.
#
# Usage (from the repository root):
#   Rscript --vanilla scripts/references/generate_gtex_end_to_end_reference.R OUTDIR
# Writes OUTDIR/counts.tsv, OUTDIR/design_matrix.tsv, OUTDIR/reference_results.tsv,
# OUTDIR/reference_size_factors.tsv and OUTDIR/provenance.json.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
out <- args[[1]]
dir.create(out, showWarnings = FALSE, recursive = TRUE)
suppressPackageStartupMessages(library(DESeq2))
src <- "crates/rsdeseq2/tests/data/gtex_v8_real_wald"
counts <- as.matrix(read.delim(file.path(src, "counts.tsv"), row.names = 1, check.names = FALSE))
meta <- read.delim(file.path(src, "sample_metadata.tsv"), row.names = 1, check.names = FALSE)
stopifnot(identical(rownames(meta), colnames(counts)))
# Deterministic panel: every 13th gene with at least one count, original order.
expressed <- which(rowSums(counts) > 0)
panel <- expressed[seq(1L, length(expressed), by = 13L)]
counts <- counts[panel, , drop = FALSE]
storage.mode(counts) <- "integer"
# Sex has a strong real signal (sex-chromosome genes), so the panel has
# significant genes whose calls can flip, unlike the random condition split.
meta$sex <- factor(meta$sex, levels = c("female", "male"))
dds <- DESeqDataSetFromMatrix(counts, meta, design = ~ sex)
dds <- DESeq(dds, quiet = TRUE)
res <- results(dds, name = "sex_male_vs_female")
write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}
write_tsv(data.frame(gene = rownames(counts), counts, check.names = FALSE),
          file.path(out, "counts.tsv"))
fmt <- function(v) ifelse(is.na(v), NA, sprintf("%.17g", v))
write_tsv(data.frame(
  gene = rownames(res),
  baseMean = fmt(res$baseMean), log2FoldChange = fmt(res$log2FoldChange),
  lfcSE = fmt(res$lfcSE), stat = fmt(res$stat), pvalue = fmt(res$pvalue),
  padj = fmt(res$padj), dispersion = fmt(mcols(dds)$dispersion),
  replaced = as.integer(mcols(dds)$replace %in% TRUE)
), file.path(out, "reference_results.tsv"))
write_tsv(data.frame(sample = colnames(dds), size_factor = fmt(sizeFactors(dds))),
          file.path(out, "reference_size_factors.tsv"))
write_tsv(data.frame(sample = colnames(dds), Intercept = 1,
                     sexmale = as.integer(meta$sex == "male")),
          file.path(out, "design_matrix.tsv"))
writeLines(jsonlite::toJSON(list(
  generator = "scripts/references/generate_gtex_end_to_end_reference.R",
  r_version = R.version.string,
  deseq2_version = as.character(packageVersion("DESeq2")),
  platform = R.version$platform,
  genes = nrow(counts), samples = ncol(counts),
  selection = "every 13th gene with a nonzero row sum, original order",
  design = "~ sex", result = "sex_male_vs_female",
  defaults = "DESeq() and results() defaults: parametric fit, minReplicatesForReplace = 7, Cook's cutoff and independent filtering on, alpha = 0.1",
  replaced_genes = sum(mcols(dds)$replace %in% TRUE),
  dispersion_fit_type = attr(dispersionFunction(dds), "fitType"),
  significant_0.1 = sum(res$padj < 0.1, na.rm = TRUE),
  significant_0.05 = sum(res$padj < 0.05, na.rm = TRUE)
), auto_unbox = TRUE, pretty = TRUE), file.path(out, "provenance.json"))
