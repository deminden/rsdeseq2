#!/usr/bin/env Rscript
# Wald per-gene inputs/results retain full-cohort dispersions and normalization.
# BH must be recomputed on the selected family of hypotheses.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==1L)
path <- file.path(args[[1]], 'reference_results.tsv')
x <- read.delim(path, check.names=FALSE)
x$padj <- stats::p.adjust(x$pvalue, method='BH')
write.table(x, path, sep='\t', quote=FALSE, row.names=FALSE, na='NA')
