# Contrast parsing and canonicalization.

.rsdeseq2_parse_results_contrast <- function(contrast,
                                             resultsNames,
                                             listValues,
                                             reference) {
  if (is.numeric(contrast) && !is.list(contrast)) {
    numeric <- .rsdeseq2_validate_numeric_results_contrast(contrast, resultsNames)
    return(list(
      type = "numeric",
      numeric = numeric,
      resultName = "contrast",
      comparison = "primitive numeric contrast",
      allZero = list(type = "numeric")
    ))
  }

  if (is.character(contrast) && !is.list(contrast)) {
    if (length(contrast) != 3L || anyNA(contrast) || any(!nzchar(contrast))) {
      stop("character contrast must be c(factor, numerator, denominator)", call. = FALSE)
    }
    if (!is.null(reference)) {
      if (!is.character(reference) || length(reference) != 1L || is.na(reference) || !nzchar(reference)) {
        stop("reference must be NULL or a single non-empty character value", call. = FALSE)
      }
    }
    factor <- contrast[[1L]]
    numerator <- contrast[[2L]]
    denominator <- contrast[[3L]]
    numeric <- .rsdeseq2_resolve_factor_level_results_contrast(
      resultsNames = resultsNames,
      factor = factor,
      numerator = numerator,
      denominator = denominator,
      reference = reference
    )
    return(list(
      type = "character",
      factor = factor,
      numerator = numerator,
      denominator = denominator,
      reference = reference,
      numeric = numeric,
      resultName = sprintf("%s_%s_vs_%s", make.names(factor), make.names(numerator), make.names(denominator)),
      comparison = sprintf("factor-level contrast: %s %s vs %s", factor, numerator, denominator),
      allZero = list(type = "character", factor = factor, numerator = numerator, denominator = denominator)
    ))
  }

  if (is.list(contrast)) {
    if (length(contrast) < 1L || length(contrast) > 2L) {
      stop("list contrast must contain one or two character vectors", call. = FALSE)
    }
    positive <- contrast[[1L]]
    negative <- if (length(contrast) == 2L) contrast[[2L]] else character()
    if (!is.character(positive) || anyNA(positive) || any(!nzchar(positive))) {
      stop("list contrast numerator must be a character vector", call. = FALSE)
    }
    if (!is.character(negative) || anyNA(negative) || any(!nzchar(negative))) {
      stop("list contrast denominator must be a character vector", call. = FALSE)
    }
    listValues <- .rsdeseq2_validate_list_values(listValues)
    numeric <- .rsdeseq2_resolve_list_results_contrast(
      resultsNames = resultsNames,
      positive = positive,
      negative = negative,
      positiveWeight = listValues[[1L]],
      negativeWeight = listValues[[2L]]
    )
    return(list(
      type = "list",
      positive = positive,
      negative = negative,
      listValues = listValues,
      numeric = numeric,
      resultName = "contrast",
      comparison = .rsdeseq2_list_contrast_comparison(
        positive = positive,
        negative = negative,
        positiveWeight = listValues[[1L]],
        negativeWeight = listValues[[2L]]
      ),
      allZero = list(type = "numeric")
    ))
  }

  stop("contrast must be numeric, character, or list", call. = FALSE)
}

.rsdeseq2_reference_for_character_contrast <- function(object, contrast) {
  if (!is.character(contrast) || length(contrast) != 3L) {
    return(NULL)
  }
  factor <- contrast[[1L]]
  fields <- names(object$factorReferences)
  if (is.null(fields)) {
    fields <- names(object$factorLevels)
  }
  factor <- .rsdeseq2_resolve_named_metadata_field(factor, fields, "factor reference")
  if (is.null(factor)) {
    return(NULL)
  }
  if (!is.null(object$factorReferences) && !is.null(object$factorReferences[[factor]])) {
    return(object$factorReferences[[factor]])
  }
  if (is.null(object$factorLevels)) {
    return(NULL)
  }
  levels <- object$factorLevels[[factor]]
  if (is.null(levels) || length(levels) == 0L) {
    return(NULL)
  }
  levels[[1L]]
}

.rsdeseq2_canonicalize_character_contrast <- function(object, contrast) {
  if (!is.character(contrast) || length(contrast) != 3L) {
    return(contrast)
  }
  factor <- contrast[[1L]]
  fields <- names(object$factorLevels)
  if (is.null(fields) && !is.null(object$colData)) {
    fields <- colnames(object$colData) %||% names(object$colData)
  }
  canonical <- .rsdeseq2_resolve_named_metadata_field(factor, fields, "character contrast factor")
  if (is.null(canonical)) {
    if (is.null(object$sampleLevels)) {
      return(contrast)
    }
    canonical <- factor
  }
  contrast[[1L]] <- canonical
  levelCandidates <- .rsdeseq2_character_contrast_level_candidates(object, canonical)
  if (!is.null(levelCandidates)) {
    contrast[[2L]] <- .rsdeseq2_resolve_sample_level_alias(contrast[[2L]], levelCandidates, "numerator")
    contrast[[3L]] <- .rsdeseq2_resolve_sample_level_alias(contrast[[3L]], levelCandidates, "denominator")
  }
  contrast
}

.rsdeseq2_character_contrast_level_candidates <- function(object, factor) {
  if (!is.null(object$factorLevels) && !is.null(object$factorLevels[[factor]])) {
    return(object$factorLevels[[factor]])
  }
  if (!is.null(object$colData)) {
    fields <- colnames(object$colData) %||% names(object$colData)
    column <- .rsdeseq2_resolve_named_metadata_field(factor, fields, "character contrast factor")
    if (!is.null(column)) {
      values <- object$colData[[column]]
      if (is.factor(values)) {
        return(levels(values))
      }
      return(unique(as.character(values)))
    }
  }
  if (!is.null(object$sampleLevels)) {
    return(unique(as.character(object$sampleLevels)))
  }
  NULL
}

.rsdeseq2_canonicalize_character_contrast_reference <- function(object, contrast, reference) {
  if (is.null(reference) || !is.character(contrast) || length(contrast) != 3L) {
    return(reference)
  }
  levelCandidates <- .rsdeseq2_character_contrast_level_candidates(object, contrast[[1L]])
  if (is.null(levelCandidates)) {
    return(reference)
  }
  .rsdeseq2_resolve_sample_level_alias(reference, levelCandidates, "reference")
}
