lfcShrink <- function(dds, coef = NULL, type = "normal", res = NULL,
                      contrast = NULL,
                      lfcThreshold = 0, svalue = FALSE,
                      independentFiltering = TRUE, alpha = 0.1,
                      apeAdapt = TRUE, apeMethod = "nbinomCR",
                      quiet = FALSE, ...) {
  dots <- list(...)
  if (length(dots) && identical(type, "normal")) {
    argumentNames <- names(dots)
    if (is.null(argumentNames)) argumentNames <- rep("...", length(dots))
    argumentNames[!nzchar(argumentNames)] <- "..."
    stop(sprintf("unsupported lfcShrink arguments: %s", paste(argumentNames, collapse = ", ")), call. = FALSE)
  }
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      !type %in% c("normal", "ashr", "apeglm")) {
    stop("type must be 'normal', 'ashr', or 'apeglm'", call. = FALSE)
  }
  if (is.null(coef) == is.null(contrast)) {
    stop("supply exactly one of coef or contrast", call. = FALSE)
  }
  if (!is.logical(svalue) || length(svalue) != 1L || is.na(svalue)) {
    stop("svalue must be TRUE or FALSE", call. = FALSE)
  }
  if (svalue && identical(type, "normal")) {
    stop("svalue is not supported by type='normal'", call. = FALSE)
  }
  if (length(lfcThreshold) != 1L || is.na(lfcThreshold) || !is.finite(lfcThreshold) || lfcThreshold < 0) {
    stop("lfcThreshold must be a finite non-negative number", call. = FALSE)
  }
  if (!is.logical(independentFiltering) || length(independentFiltering) != 1L || is.na(independentFiltering)) {
    stop("independentFiltering must be TRUE or FALSE", call. = FALSE)
  }
  if (length(alpha) != 1L || is.na(alpha) || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must be a finite number between 0 and 1", call. = FALSE)
  }
  if (!is.logical(apeAdapt) || length(apeAdapt) != 1L || is.na(apeAdapt)) {
    stop("apeAdapt must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.character(apeMethod) || length(apeMethod) != 1L || is.na(apeMethod) || !nzchar(apeMethod)) {
    stop("apeMethod must be a non-empty string", call. = FALSE)
  }
  if (!is.logical(quiet) || length(quiet) != 1L || is.na(quiet)) {
    stop("quiet must be TRUE or FALSE", call. = FALSE)
  }
  if (identical(type, "apeglm")) {
    if (!quiet) {
      message("using 'apeglm' for LFC shrinkage")
    }
    return(.rsdeseq2_lfc_shrink_apeglm(
      dds, coef, contrast, res, lfcThreshold, svalue,
      independentFiltering, alpha, apeAdapt, apeMethod, quiet, dots
    ))
  }
  if (identical(type, "ashr")) {
    return(.rsdeseq2_lfc_shrink_ashr(
      dds, coef, contrast, res, lfcThreshold, svalue,
      independentFiltering, alpha, dots
    ))
  }
  .rsdeseq2_lfc_shrink_dataset(
    dds, coef, contrast, res, lfcThreshold, independentFiltering, alpha
  )
}
