#!/usr/bin/env Rscript
# Trace original DESeq2 execution; no candidate code or numerical output used.
args <- commandArgs(TRUE)
if (length(args) != 1L) stop("usage: FRESH_OUTPUT_JSON")
stopifnot(as.character(getRversion()) == "4.6.2",
  as.character(packageVersion("DESeq2")) == "1.52.0", !file.exists(args[[1L]]))
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[[1L]])
root <- normalizePath(file.path(dirname(gsub("~+~", " ", script, fixed=TRUE)), "..", ".."))
source(file.path(root, "r-pkg/rsdeseq2/tests/testthat/helper-fit-data.R"))
source(file.path(root, "r-pkg/rsdeseq2/tests/testthat/helper-formula.R"))
bits <- function(v) lapply(as.numeric(v), function(x) paste0("bits:", hex(x)))
rows <- function(m) unname(as.vector(t(m)))
counts <- .rsdeseq2_test_counts()
counts[40L, ] <- 0
counts[1L, 8L] <- 2000
counts[2L, ] <- 0
counts[2L, 8L] <- 1000
size <- exp(seq(log(.8), log(1.2), length.out=8L))
factors <- outer(seq(.9, 1.1, length.out=40L), size)
weights <- matrix(rep(c(.3, 1, .7, 1, 1, .4, 1, .8), each=40L), 40L, 8L)
design <- cbind(Intercept=1, conditionB=rep(0:1, each=4L))
col_data <- data.frame(condition=factor(rep(c("A", "B"), each=4L)), row.names=colnames(counts))
.replacement_execution_capture <- new.env()
ns <- asNamespace("DESeq2")
suppressMessages(trace("refitWithoutOutliers", where=ns, print=FALSE,
  tracer=quote({.replacement_execution_capture$refits <- .replacement_execution_capture$refits + 1L})))
suppressMessages(trace("estimateDispersionsGeneEst", where=ns, print=FALSE,
  tracer=quote({.replacement_execution_capture$gene_fits <- .replacement_execution_capture$gene_fits + 1L})))
run <- function(test, use_factors, use_weights, minimum) {
  dds <- DESeq2::DESeqDataSetFromMatrix(counts, col_data, ~condition)
  DESeq2::sizeFactors(dds) <- size
  if (use_factors) DESeq2::normalizationFactors(dds) <- factors
  if (use_weights) SummarizedExperiment::assay(dds, "weights", withDimnames=FALSE) <- weights
  .replacement_execution_capture$refits <- 0L
  .replacement_execution_capture$gene_fits <- 0L
  parameters <- list(object=dds, test=test, fitType="mean", minReplicatesForReplace=minimum, quiet=TRUE)
  if (test == "LRT") parameters$reduced <- ~1
  fitted <- do.call(DESeq2::DESeq, parameters)
  table <- DESeq2::results(fitted, cooksCutoff=FALSE, independentFiltering=FALSE)
  list(table=as.data.frame(table), refits=.replacement_execution_capture$refits,
    gene_fits=.replacement_execution_capture$gene_fits,
    replacements=sum(S4Vectors::mcols(fitted)$replace, na.rm=TRUE))
}
cases <- list()
for (test in c("Wald", "LRT")) for (use_factors in c(FALSE, TRUE)) for (use_weights in c(FALSE, TRUE)) {
  disabled <- run(test, use_factors, use_weights, Inf)
  ineligible <- run(test, use_factors, use_weights, 7L)
  stopifnot(disabled$refits == 0L, ineligible$refits == 0L,
    disabled$gene_fits == 1L, ineligible$gene_fits == 1L,
    ineligible$replacements == 0L, identical(disabled$table, ineligible$table))
  cases[[length(cases)+1L]] <- list(name=paste(test, use_factors, use_weights, sep="_"),
    test=test, normalization_factors=use_factors, weights=use_weights,
    refit_calls=ineligible$refits, gene_fit_calls=ineligible$gene_fits,
    replacements=ineligible$replacements, complete_results_identical_to_disabled=TRUE)
}
positive <- lapply(c("Wald", "LRT"), function(test) {
  value <- run(test, FALSE, FALSE, 3L)
  stopifnot(value$refits == 1L, value$gene_fits == 2L, value$replacements == 2L)
  list(test=test, refit_calls=value$refits, gene_fit_calls=value$gene_fits, replacements=value$replacements)
})
suppressMessages(untrace("refitWithoutOutliers", where=ns))
suppressMessages(untrace("estimateDispersionsGeneEst", where=ns))
jsonlite::write_json(list(source=list(R=as.character(getRversion()), DESeq2=as.character(packageVersion("DESeq2"))),
  genes=40L, samples=8L, counts=as.list(as.integer(rows(counts))), design=bits(rows(design)),
  size_factors=bits(size), normalization_factors=bits(rows(factors)), weights=bits(rows(weights)),
  cases=cases, positive_controls=positive), args[[1L]], auto_unbox=TRUE, pretty=TRUE, null="null", digits=NA)
