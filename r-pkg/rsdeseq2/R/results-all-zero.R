# Contrast-zero context contains identifiers and borrowed primitive inputs only.
# Every count reduction, design score and result override runs in Rust.

.rsdeseq2_contrast_zero_context <- function(fit, allZero) {
  if (is.null(fit$counts) || is.null(allZero) || identical(allZero$type, "none")) {
    return(NULL)
  }
  if (identical(allZero$type, "character")) {
    levels <- .rsdeseq2_character_contrast_sample_levels(fit, allZero)
    numerator <- .rsdeseq2_resolve_sample_level_alias(allZero$numerator, levels, "numerator")
    denominator <- .rsdeseq2_resolve_sample_level_alias(allZero$denominator, levels, "denominator")
    if (identical(numerator, denominator)) {
      stop("character contrast all-zero handling requires distinct numerator and denominator levels", call. = FALSE)
    }
    return(list(rule = 2L, counts = fit$counts, sampleLevels = levels,
      numerator = numerator, denominator = denominator))
  }
  if (!identical(allZero$type, "numeric")) {
    stop("unknown contrast all-zero rule", call. = FALSE)
  }
  list(rule = 1L, counts = fit$counts, design = fit$modelMatrix)
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
