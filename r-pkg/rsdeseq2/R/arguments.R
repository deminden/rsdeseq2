# Keep argument rejection identical across public dispatch and constructor paths.
.rsdeseq2_check_unused_arguments <- function(dots, caller) {
  if (!length(dots)) return(invisible(NULL))
  labels <- names(dots)
  if (is.null(labels)) labels <- rep("...", length(dots))
  labels[is.na(labels) | !nzchar(labels)] <- "..."
  stop(sprintf("unsupported %s arguments: %s", caller, paste(labels, collapse = ", ")),
       call. = FALSE)
}

# DESeq2's estimateDispersionsFit replaces a failed parametric trend with the
# local fit and says so. The Rust core does the replacement and reports the
# trend it fitted (0 parametric, 1 local, 2 mean); return that fitType.
.rsdeseq2_note_fit_type_used <- function(fitType, fitTypeUsed) {
  if (is.null(fitTypeUsed)) return(fitType)
  used <- c("parametric", "local", "mean")[fitTypeUsed + 1L]
  if (identical(fitType, "parametric") && identical(used, "local")) {
    message("-- note: fitType='parametric', but the dispersion trend was not well captured by the
   function: y = a/x + b, and a local regression fit was automatically substituted.
   specify fitType='local' or 'mean' to avoid this message next time.")
  }
  used
}
