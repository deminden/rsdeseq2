#!/usr/bin/env Rscript
# Independent original-DESeq2 reference. Never load or call rsdeseq2 here.
args <- commandArgs(TRUE)
if (length(args) != 1L) stop("usage: generate_frozen_local_vst_reference.R OUTPUT_JSON")
suppressPackageStartupMessages(library(DESeq2))
stopifnot(as.character(getRversion()) == "4.6.2",
          as.character(packageVersion("DESeq2")) == "1.52.0",
          as.character(packageVersion("locfit")) == "1.5.9.12")
numbers <- function(x) as.list(sprintf("%.17g", as.numeric(x)))
cases <- list()
for (profile in c("smooth", "irregular")) {
  n <- if (profile == "smooth") 60L else 120L
  position <- seq(0, 1, length.out = n)
  if (profile == "irregular") position <- position ^ 2.4
  means <- exp(log(.3) + (log(5000) - log(.3)) * position)
  if (profile == "irregular") means[20L] <- means[19L]
  disps <- if (profile == "smooth") {
    .03 + 1.5 / means + .015 * (1 + sin(log(means)))
  } else {
    .025 + 2 / means + .012 * (1 + cos(log(means) * 2.3))
  }
  trend <- DESeq2:::localDispersionFit(means, disps, 1e-8)
  attr(trend, "fitType") <- "local"
  path <- tempfile(fileext = ".rds")
  saveRDS(trend, path)
  reloaded <- readRDS(path)
  unlink(path)
  fit <- get("fit", environment(trend), inherits = FALSE)
  log_queries <- sort(unique(c(seq(-8, 15, length.out = 79L), fit$eva$xev)))
  queries <- exp(log_queries)
  predictions <- trend(queries)
  stopifnot(identical(predictions, reloaded(queries)))
  for (normalization in c("size_factors", "normalization_factors")) {
    counts <- matrix(as.integer(round(exp(seq(0, log(1e6), length.out = 48L)))), 12L, 4L)
    counts[1L, 1L] <- 0L
    rownames(counts) <- paste0("gene", seq_len(nrow(counts)))
    colnames(counts) <- paste0("sample", seq_len(ncol(counts)))
    dds <- DESeqDataSetFromMatrix(counts,
      data.frame(condition = factor(c("A", "A", "B", "B")), row.names = colnames(counts)),
      ~ condition)
    sf <- c(.7, 1.1, 1.4, .9)
    nf <- outer(seq(.8, 1.2, length.out = 12L), sf)
    if (normalization == "size_factors") sizeFactors(dds) <- sf else normalizationFactors(dds) <- nf
    dispersionFunction(dds) <- reloaded
    normalized <- DESeq2::counts(dds, normalized = TRUE)
    grid <- sinh(seq(asinh(0), asinh(max(normalized)), length.out = 1000L))[-1L]
    transformed <- assay(DESeq2::varianceStabilizingTransformation(dds, blind = FALSE, fitType = "local"))
    cases[[length(cases) + 1L]] <- list(
      name = paste(profile, normalization, sep = "_"),
      mi = as.list(fit$mi), nvc = as.list(fit$nvc),
      knots = numbers(fit$eva$xev), residual_values = numbers(fit$eva$coef[, 1L]),
      residual_slopes = numbers(fit$eva$coef[, 2L]), global_polynomial = numbers(fit$eva$pc[1:4]),
      queries = numbers(queries), predictions = numbers(predictions),
      genes = nrow(counts), samples = ncol(counts), counts = as.list(as.integer(t(counts))),
      size_factors = if (normalization == "size_factors") numbers(sf) else NULL,
      normalization_factors = if (normalization == "normalization_factors") numbers(t(nf)) else NULL,
      grid = numbers(grid), dispersions = numbers(reloaded(grid)),
      transformed = numbers(t(transformed)), serialization_identical = TRUE)
  }
}
floor_means <- c(1, 2, 3)
floor_disps <- c(0, 1e-10, 9e-8)
floor_values <- DESeq2:::localDispersionFit(floor_means, floor_disps, 1e-8)
jsonlite::write_json(list(source = list(R = as.character(getRversion()),
  DESeq2 = as.character(packageVersion("DESeq2")), locfit = as.character(packageVersion("locfit"))),
  minimum_floor = list(means = numbers(floor_means), values = numbers(floor_values)),
  cases = cases), args[[1L]], auto_unbox = TRUE, pretty = TRUE, null = "null")
