#!/usr/bin/env Rscript
# Original DESeq2 fits only. Candidate fitting/metadata never enters this oracle.
args <- commandArgs(TRUE)
if (length(args) != 1L) stop("usage: FRESH_OUTPUT_JSON")
stopifnot(as.character(getRversion()) == "4.6.2",
          as.character(packageVersion("DESeq2")) == "1.52.0", !file.exists(args[[1L]]))
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[[1L]])
root <- normalizePath(file.path(dirname(gsub("~+~", " ", script, fixed=TRUE)), "..", ".."))
source(file.path(root, "r-pkg/rsdeseq2/tests/testthat/helper-fit-data.R"))
source(file.path(root, "r-pkg/rsdeseq2/tests/testthat/helper-formula.R"))
bits <- function(x) lapply(as.numeric(x), function(v) paste0("bits:", hex(v)))
rows <- function(x) unname(as.vector(t(x)))
counts <- .rsdeseq2_test_counts()
counts[40L, ] <- 0
weights <- matrix(rep(c(.3, 1, .7, 1, 1, .4, 1, .8), each=40L), 40L, 8L)
weights[39L, ] <- c(1, rep(0, 7L))
size <- exp(seq(log(.8), log(1.2), length.out=8L))
factors <- outer(seq(.9, 1.1, length.out=40L), size)
design <- cbind(Intercept=1, conditionB=rep(0:1, each=4L))
col_data <- data.frame(condition=factor(rep(c("A", "B"), each=4L)), row.names=colnames(counts))
cases <- list()
for (test in c("Wald", "LRT")) for (kind in c("disabled", "no_refit", "replacement", "new_zero_only")) {
  input <- counts
  if (kind == "replacement") input[1L, 8L] <- 2000
  if (kind %in% c("replacement", "new_zero_only")) {
    input[2L, ] <- 0
    input[2L, 8L] <- 1000
  }
  dds <- DESeq2::DESeqDataSetFromMatrix(input, col_data, ~condition)
  DESeq2::sizeFactors(dds) <- size
  if (test == "LRT") DESeq2::normalizationFactors(dds) <- factors
  SummarizedExperiment::assay(dds, "weights", withDimnames=FALSE) <- weights
  parameters <- list(object=dds, test=test, fitType="mean",
    minReplicatesForReplace=if (kind == "disabled") Inf else 3L, quiet=TRUE)
  if (test == "LRT") parameters$reduced <- ~1
  fitted <- suppressWarnings(do.call(DESeq2::DESeq, parameters))
  metadata <- S4Vectors::mcols(fitted)
  result <- DESeq2::results(fitted, cooksCutoff=FALSE, independentFiltering=FALSE)
  convergence <- metadata[[if (test == "Wald") "betaConv" else "fullBetaConv"]]
  n_refit <- sum(metadata$replace, na.rm=TRUE)
  stopifnot(nrow(result) == nrow(input), length(result$stat) == nrow(input),
            length(result$pvalue) == nrow(input), length(convergence) == nrow(input),
            isTRUE(metadata$weightsFail[[39L]]),
            isTRUE(metadata$allZero[[39L]]) == (n_refit == 0L),
            is.na(result$pvalue[[39L]]), is.na(result$stat[[39L]]),
            is.na(convergence[[39L]]),
            n_refit == switch(kind, disabled=0L, no_refit=0L, replacement=2L, new_zero_only=1L))
  cases[[length(cases)+1L]] <- list(name=paste(test, kind, sep="_"), test=test, kind=kind,
    counts=as.list(as.integer(rows(input))), normalization_factors=test == "LRT",
    all_zero=as.list(unname(metadata$allZero)), weights_fail=as.list(unname(metadata$weightsFail)),
    n_refit=n_refit, failed_mean=bits(metadata$baseMean[[39L]])[[1L]],
    failed_variance=bits(metadata$baseVar[[39L]])[[1L]],
    failed_stat_missing=is.na(result$stat[[39L]]), failed_pvalue_missing=is.na(result$pvalue[[39L]]),
    failed_convergence_missing=is.na(convergence[[39L]]))
}
jsonlite::write_json(list(source=list(R=as.character(getRversion()), DESeq2=as.character(packageVersion("DESeq2"))),
  genes=40L, samples=8L, failed_row=38L, design=bits(rows(design)), size_factors=bits(size),
  factors=bits(rows(factors)), weights=bits(rows(weights)), cases=cases),
  args[[1L]], auto_unbox=TRUE, pretty=TRUE, null="null", digits=NA)
