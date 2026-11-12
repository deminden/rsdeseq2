#!/usr/bin/env Rscript
# Independent original-DESeq2 oracles; no production numerical implementation.

suppressPackageStartupMessages({
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("DESeq2 is required. Install it with BiocManager::install('DESeq2').")
  }
})

file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", gsub("~+~", " ", file_arg[[1L]], fixed = TRUE)), mustWork = FALSE)
} else {
  normalizePath(file.path("scripts", "references", "generate_deseq2_references.R"), mustWork = FALSE)
}
repo_root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = FALSE)
module_dir <- file.path(dirname(script_path), "deseq2")
if (!dir.exists(module_dir)) stop("Missing reference generator modules: ", module_dir)

parse_reference_directories <- function(args, repo_root) {
  directories <- list(
    data_dir = file.path(repo_root, "crates", "rsdeseq2", "tests", "data", "deseq2_reference"),
    out_dir = file.path(repo_root, "results", "parity")
  )
  flags <- c("--output-dir" = "data_dir", "--parity-output-dir" = "out_dir")
  if (identical(args, "--help")) {
    cat("Usage: generate_deseq2_references.R [--output-dir FRESH_DIR] [--parity-output-dir FRESH_DIR]\n")
    return(NULL)
  }
  if (length(args) %% 2L) stop("Each output option requires a directory")
  seen <- character()
  for (pair in seq_len(length(args) %/% 2L)) {
    i <- 2L * pair - 1L
    flag <- args[[i]]
    if (!flag %in% names(flags) || flag %in% seen) stop("Unknown or repeated option: ", flag)
    path <- args[[i + 1L]]
    if (!nzchar(path) || startsWith(path, "--")) stop("Missing directory for ", flag)
    if (file.exists(path)) stop("Explicit output directory must be fresh: ", path)
    directories[[flags[[flag]]]] <- path
    seen <- c(seen, flag)
  }
  directories
}

generate_reference_panels <- function(data_dir, out_dir, module_dir) {
  # A private evaluation environment retains the original panel sequence and
  # fixture dependencies without creating shared variables in the caller.
  context <- environment()
  sys.source(file.path(module_dir, "io.R"), envir = context)
  prepare_reference_directory(data_dir, out_dir)
  sys.source(file.path(module_dir, "helpers.R"), envir = context)
  panels <- c(
    "fixture",
    "normalization",
    "replacement",
    "gene_wise",
    "weighted_mean",
    "mean",
    "local",
    "mean_cox_reid",
    "dispersion_models",
    "fixed",
    "fixed_optim",
    "fixed_weighted"
  )
  for (panel in panels) {
    sys.source(file.path(module_dir, paste0(panel, ".R")), envir = context)
  }
  session_info <- sub("[[:blank:]]+$", "", capture.output(sessionInfo()))
  writeLines(session_info, file.path(data_dir, "sessionInfo.txt"))
  invisible(file.copy(file.path(data_dir, "metadata.tsv"), file.path(out_dir, "metadata.tsv"), overwrite = TRUE))

  message("Wrote DESeq2 references to: ", data_dir)
}

directories <- parse_reference_directories(commandArgs(TRUE), repo_root)
if (!is.null(directories)) {
  generate_reference_panels(directories$data_dir, directories$out_dir, module_dir)
}
