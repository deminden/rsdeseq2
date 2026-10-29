#!/usr/bin/env Rscript
# Deterministic workflow/transform profiling; outputs are precision evidence too.
args <- commandArgs(TRUE)
if (length(args) != 6L) stop("usage: LIB WORKFLOW GENES SAMPLES OUTPUT TIMINGS")
lib <- normalizePath(args[[1L]])
workflow <- args[[2L]]; genes <- as.integer(args[[3L]]); samples <- as.integer(args[[4L]])
suppressPackageStartupMessages({library(rsdeseq2, lib.loc = lib); library(DESeq2); library(S4Vectors)})
stopifnot(getRversion() == "4.6.2", samples %% 2L == 0L, genes >= 1L)
stages <- data.frame(stage = character(), elapsed_s = numeric())
timed <- function(name, expression) {
  start <- proc.time()[["elapsed"]]
  value <- force(expression)
  stages[nrow(stages) + 1L, ] <<- list(name, proc.time()[["elapsed"]] - start)
  write.table(stages, args[[6L]], sep = "\t", quote = FALSE, row.names = FALSE)
  value
}
set.seed(20261029)
condition <- factor(rep(c("A", "B"), each = samples / 2L))
batch <- factor(rep(seq_len(4L), length.out = samples))
base <- exp(seq(log(10), log(1000), length.out = genes))
dispersions <- seq(.03, .3, length.out = genes)
sizeFactors <- exp(seq(log(.8), log(1.2), length.out = samples))
means <- outer(base, sizeFactors)
means[, condition == "B"] <- means[, condition == "B"] * rep(c(1, 1.5, .75, 2), length.out = genes)
counts <- matrix(rnbinom(genes * samples, mu = as.vector(means), size = rep(1/dispersions, samples)), nrow = genes)
rownames(counts) <- paste0("gene", seq_len(genes)); colnames(counts) <- paste0("sample", seq_len(samples))
if (workflow == "replacement") {
  script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]])
  script <- gsub("~+~", " ", script, fixed = TRUE)
  source(file.path(dirname(normalizePath(script)), "..", "r-pkg/rsdeseq2/tests/testthat/helper-fit-data.R"))
  counts <- .rsdeseq2_test_counts(); counts[1L, 8L] <- 2000L
  counts[2L, ] <- 0L; counts[2L, 8L] <- 1000L
  genes <- nrow(counts); samples <- ncol(counts)
  condition <- factor(rep(c("A", "B"), each = 4L))
  batch <- factor(rep(1:4, 2)); sizeFactors <- NULL
}
cd <- data.frame(condition = condition, batch = batch, row.names = colnames(counts))
formula <- if (workflow %in% c("rlog", "vst", "replacement")) ~ condition else ~ batch + condition
dds <- DESeq2::DESeqDataSetFromMatrix(counts, cd, formula)
fitArgs <- list(object = dds, fitType = "mean", sizeFactors = sizeFactors, nproc = 1L,
                minReplicatesForReplace = if (workflow == "replacement") 3L else if (workflow == "replacement_none") 7L else Inf)
if (workflow == "lrt") { fitArgs$test <- "LRT"; fitArgs$reduced <- ~ batch }
fit <- timed("fit", do.call(rsdeseq2::DESeq, fitArgs))
state <- S4Vectors::metadata(fit)$rsdeseq2RustFit
if (workflow == "replacement") stopifnot(state$nRefit > 0L)
if (workflow == "replacement_none") stopifnot(state$nRefit == 0L)
value <- switch(workflow,
  rlog = timed("rlog", {out <- rsdeseq2::rlog(fit, blind = FALSE, fitType = "mean", nproc = 1L); list(values = SummarizedExperiment::assay(out), metadata = S4Vectors::metadata(out))}),
  vst = timed("vst", {out <- rsdeseq2::varianceStabilizingTransformation(fit, blind = FALSE, fitType = "mean", nproc = 1L); SummarizedExperiment::assay(out)}),
  multi_results = timed("ten_results", lapply(seq_len(10L), function(i) as.data.frame(rsdeseq2::results(fit, name = rsdeseq2::resultsNames(fit)[[1L + i %% ncol(state$beta)]])))),
  shrink_normal = timed("shrink_normal", as.data.frame(rsdeseq2::lfcShrink(fit, coef = tail(rsdeseq2::resultsNames(fit), 1L), type = "normal"))),
  shrink_ashr = timed("shrink_ashr", as.data.frame(rsdeseq2::lfcShrink(fit, coef = tail(rsdeseq2::resultsNames(fit), 1L), type = "ashr"))),
  shrink_apeglm = timed("shrink_apeglm", as.data.frame(rsdeseq2::lfcShrink(fit, coef = tail(rsdeseq2::resultsNames(fit), 1L), type = "apeglm"))),
  timed("results", as.data.frame(rsdeseq2::results(fit))))
saveRDS(list(counts = counts, fitted = state, output = value), args[[5L]], compress = FALSE)
