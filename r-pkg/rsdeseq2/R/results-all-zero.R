# All-zero contrast handling and observed sample-level aliases.

.rsdeseq2_apply_lrt_results_contrast_all_zero <- function(lrt,
                                                          fit,
                                                          numericContrast,
                                                          allZero) {
  flags <- .rsdeseq2_results_contrast_all_zero_flags(
    fit = fit,
    numericContrast = numericContrast,
    allZero = allZero
  )
  if (is.null(flags)) {
    return(lrt)
  }
  allZeroRows <- if (is.null(fit$counts)) rep(FALSE, length(flags)) else rowSums(fit$counts) == 0
  zeroRows <- which(flags & !allZeroRows)
  if (length(zeroRows) != 0L) {
    lrt$log2FoldChange[zeroRows] <- 0
  }
  lrt
}

.rsdeseq2_apply_results_contrast_all_zero <- function(wald,
                                                      fit,
                                                      numericContrast,
                                                      allZero) {
  flags <- .rsdeseq2_results_contrast_all_zero_flags(
    fit = fit,
    numericContrast = numericContrast,
    allZero = allZero
  )
  if (is.null(flags)) {
    return(wald)
  }
  allZeroRows <- if (is.null(fit$counts)) rep(FALSE, length(flags)) else rowSums(fit$counts) == 0
  zeroRows <- which(flags & !allZeroRows)
  if (length(zeroRows) != 0L) {
    wald$log2FoldChange[zeroRows] <- 0
    wald$stat[zeroRows] <- 0
    wald$pvalue[zeroRows] <- 1
  }
  wald
}

.rsdeseq2_results_contrast_all_zero_flags <- function(fit,
                                                      numericContrast,
                                                      allZero) {
  if (is.null(fit$counts)) {
    return(NULL)
  }
  if (is.null(allZero) || identical(allZero$type, "none")) {
    return(NULL)
  }
  if (identical(allZero$type, "character")) {
    sampleLevels <- .rsdeseq2_character_contrast_sample_levels(fit, allZero)
    numerator <- .rsdeseq2_resolve_sample_level_alias(
      allZero$numerator,
      sampleLevels,
      "numerator"
    )
    denominator <- .rsdeseq2_resolve_sample_level_alias(
      allZero$denominator,
      sampleLevels,
      "denominator"
    )
    if (identical(numerator, denominator)) {
      stop("character contrast all-zero handling requires distinct numerator and denominator levels", call. = FALSE)
    }
    selected <- sampleLevels %in% c(numerator, denominator)
    return(rowSums(fit$counts[, selected, drop = FALSE] == 0) == sum(selected))
  }
  if (!identical(allZero$type, "numeric")) {
    stop("unknown contrast all-zero rule", call. = FALSE)
  }
  if (all(numericContrast >= 0) || all(numericContrast <= 0)) {
    return(rep(FALSE, nrow(fit$counts)))
  }
  if (is.null(fit$modelMatrix)) {
    stop("numeric contrast all-zero handling requires modelMatrix", call. = FALSE)
  }
  contrastBinary <- ifelse(numericContrast == 0, 0, 1)
  selected <- as.numeric(fit$modelMatrix %*% contrastBinary) != 0
  rowSums(fit$counts[, selected, drop = FALSE]) == 0
}

.rsdeseq2_character_contrast_sample_levels <- function(fit, allZero) {
  if (!is.null(fit$sampleLevels)) {
    return(fit$sampleLevels)
  }
  colData <- fit$colData
  factor <- allZero$factor
  field <- .rsdeseq2_resolve_named_metadata_field(
    factor,
    colnames(colData) %||% names(colData),
    "colData"
  )
  if (is.null(colData) || is.null(factor) || is.null(field) || is.null(colData[[field]])) {
    stop("character contrast all-zero handling requires sampleLevels or matching colData factor", call. = FALSE)
  }
  sampleLevels <- colData[[field]]
  if (is.factor(sampleLevels)) {
    sampleLevels <- as.character(sampleLevels)
  }
  if (!is.character(sampleLevels) ||
    length(sampleLevels) != ncol(fit$counts) ||
    anyNA(sampleLevels) ||
    any(!nzchar(sampleLevels))) {
    stop("character contrast all-zero handling requires colData factor values with one non-missing value per sample", call. = FALSE)
  }
  sampleLevels
}

.rsdeseq2_resolve_sample_level_alias <- function(level, sampleLevels, role) {
  uniqueLevels <- unique(sampleLevels)
  exact <- which(uniqueLevels == level)
  if (length(exact) == 1L) {
    return(uniqueLevels[[exact]])
  }
  if (length(exact) > 1L) {
    stop(sprintf("character contrast %s level %s is duplicated", role, level), call. = FALSE)
  }
  alias <- make.names(level)
  aliases <- make.names(uniqueLevels)
  matches <- which(aliases == alias)
  if (length(matches) == 1L) {
    return(uniqueLevels[[matches]])
  }
  if (length(matches) > 1L) {
    stop(sprintf("character contrast %s level alias %s is ambiguous", role, level), call. = FALSE)
  }
  stop(sprintf("character contrast all-zero handling requires sampleLevels or colData to contain %s level %s", role, level), call. = FALSE)
}
