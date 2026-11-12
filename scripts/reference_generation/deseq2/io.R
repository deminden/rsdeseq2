# Output handling for independent original-DESeq2 reference generation.
prepare_reference_directory <- function(data_dir, out_dir) {
  preserve_pattern <- "^hard_se_covariance_.*[.]tsv$"
  preserved_dir <- tempfile("rsdeseq2-preserved-reference-")
  preserved_files <- character()
  if (dir.exists(data_dir)) {
    preserved_files <- list.files(data_dir, pattern = preserve_pattern, full.names = TRUE)
    if (length(preserved_files) > 0) {
      dir.create(preserved_dir, recursive = TRUE, showWarnings = FALSE)
      copied <- file.copy(preserved_files, preserved_dir, overwrite = TRUE)
      if (!all(copied)) {
        stop("failed to preserve separately generated hard SE/covariance fixtures")
      }
    }
    unlink(data_dir, recursive = TRUE)
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
  if (length(preserved_files) > 0) {
    restored <- file.copy(
      file.path(preserved_dir, basename(preserved_files)),
      data_dir,
      overwrite = TRUE
    )
    if (!all(restored)) {
      stop("failed to restore separately generated hard SE/covariance fixtures")
    }
    unlink(preserved_dir, recursive = TRUE)
  }
  invisible(data_dir)
}

write_tsv <- function(x, path) {
  write.table(
    x,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
  )
}

write_matrix_tsv <- function(matrix_value, row_name_col, path) {
  row_ids <- rownames(matrix_value)
  if (is.null(row_ids)) {
    row_ids <- seq_len(nrow(matrix_value))
  }
  write_tsv(
    data.frame(
      setNames(data.frame(row_ids, check.names = FALSE), row_name_col),
      matrix_value,
      check.names = FALSE
    ),
    path
  )
}
