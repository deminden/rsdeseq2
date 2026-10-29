# Keep argument rejection identical across public dispatch and constructor paths.
.rsdeseq2_check_unused_arguments <- function(dots, caller) {
  if (!length(dots)) return(invisible(NULL))
  labels <- names(dots)
  if (is.null(labels)) labels <- rep("...", length(dots))
  labels[is.na(labels) | !nzchar(labels)] <- "..."
  stop(sprintf("unsupported %s arguments: %s", caller, paste(labels, collapse = ", ")),
       call. = FALSE)
}
