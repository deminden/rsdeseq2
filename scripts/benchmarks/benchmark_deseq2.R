#!/usr/bin/env Rscript

parse_args <- function(args) {
  out <- list()
  idx <- 1L
  while (idx <= length(args)) {
    key <- args[[idx]]
    if (!startsWith(key, "--")) {
      stop(sprintf("unexpected argument: %s", key), call. = FALSE)
    }
    if (idx == length(args)) {
      stop(sprintf("missing value for %s", key), call. = FALSE)
    }
    out[[substring(key, 3L)]] <- args[[idx + 1L]]
    idx <- idx + 2L
  }
  out
}

read_counts <- function(path) {
  table <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (ncol(table) < 2L) {
    stop("count table must have a gene column and at least one sample column", call. = FALSE)
  }
  gene <- table[[1L]]
  counts <- as.matrix(table[, -1L, drop = FALSE])
  storage.mode(counts) <- "numeric"
  rownames(counts) <- gene
  counts
}

estimate_size_factors <- function(counts, method) {
  method <- match.arg(method, c("ratio", "poscounts"))
  DESeq2::estimateSizeFactorsForMatrix(counts, type = method)
}

write_size_factors <- function(path, size_factors) {
  write.table(
    data.frame(sample = names(size_factors), size_factor = unname(size_factors)),
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

write_base_mean <- function(path, base_mean) {
  write.table(
    data.frame(gene = names(base_mean), base_mean = unname(base_mean)),
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  operation <- args[["operation"]]
  counts_path <- args[["counts"]]
  output_path <- args[["output"]]
  method <- args[["method"]] %||% "ratio"

  if (is.null(operation) || is.null(counts_path) || is.null(output_path)) {
    stop("--operation, --counts, and --output are required", call. = FALSE)
  }
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("DESeq2 is required for the DESeq2 benchmark", call. = FALSE)
  }

  counts <- read_counts(counts_path)
  size_factors <- estimate_size_factors(counts, method)
  names(size_factors) <- colnames(counts)

  if (operation == "size-factors") {
    write_size_factors(output_path, size_factors)
  } else if (operation == "base-mean") {
    normalized <- t(t(counts) / size_factors)
    write_base_mean(output_path, rowMeans(normalized))
  } else {
    stop(sprintf("unsupported operation: %s", operation), call. = FALSE)
  }
}

`%||%` <- function(lhs, rhs) {
  if (is.null(lhs)) rhs else lhs
}

main()
