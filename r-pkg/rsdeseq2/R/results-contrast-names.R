# Named, numeric, list and factor-level contrast resolution.

.rsdeseq2_resolve_named_metadata_field <- function(name, fields, source) {
  if (is.null(name) || is.null(fields)) {
    return(NULL)
  }
  exact <- which(fields == name)
  if (length(exact) == 1L) {
    return(fields[[exact]])
  }
  if (length(exact) > 1L) {
    stop(sprintf("%s field %s is duplicated", source, name), call. = FALSE)
  }
  alias <- make.names(name)
  aliases <- make.names(fields)
  matches <- which(aliases == alias)
  if (length(matches) == 1L) {
    return(fields[[matches]])
  }
  if (length(matches) > 1L) {
    stop(sprintf("%s field alias %s is ambiguous", source, name), call. = FALSE)
  }
  NULL
}

.rsdeseq2_validate_results_names <- function(resultsNames) {
  if (!is.character(resultsNames) || length(resultsNames) == 0L || anyNA(resultsNames) || any(!nzchar(resultsNames))) {
    stop("resultsNames must be a non-empty character vector", call. = FALSE)
  }
  if (anyDuplicated(resultsNames)) {
    stop("resultsNames must be unique", call. = FALSE)
  }
  resultsNames
}

.rsdeseq2_validate_numeric_results_contrast <- function(contrast, resultsNames) {
  if (length(contrast) != length(resultsNames)) {
    stop("numeric contrast must have one value per results name", call. = FALSE)
  }
  if (any(!is.finite(contrast))) {
    stop("numeric contrast must contain finite values", call. = FALSE)
  }
  if (!any(contrast != 0)) {
    stop("numeric contrast vector cannot be all zero", call. = FALSE)
  }
  out <- as.numeric(contrast)
  names(out) <- resultsNames
  out
}

.rsdeseq2_validate_lrt_stat <- function(lrtStat, nGenes) {
  lrtStat <- .rsdeseq2_validate_result_numeric_vector(lrtStat, "lrtStat", allowNA = TRUE)
  if (length(lrtStat) != nGenes) {
    stop("lrtStat must have one value per beta row", call. = FALSE)
  }
  if (any(!is.na(lrtStat) & lrtStat < 0)) {
    stop("lrtStat must contain non-negative values or NA", call. = FALSE)
  }
  lrtStat
}

.rsdeseq2_validate_lrt_df <- function(lrtDf, nGenes) {
  if (is.null(lrtDf)) {
    return(NULL)
  }
  if (!is.numeric(lrtDf) || !(length(lrtDf) %in% c(1L, nGenes))) {
    stop("lrtDf must be a scalar or one value per gene", call. = FALSE)
  }
  out <- rep(as.numeric(lrtDf), length.out = nGenes)
  if (any(!is.na(out) & (!is.finite(out) | out <= 0))) {
    stop("lrtDf must contain positive finite values or NA", call. = FALSE)
  }
  out
}

.rsdeseq2_validate_list_values <- function(listValues) {
  if (!is.numeric(listValues) || length(listValues) != 2L || any(!is.finite(listValues))) {
    stop("listValues must be a finite numeric vector of length two", call. = FALSE)
  }
  if (listValues[[1L]] <= 0 || listValues[[2L]] >= 0) {
    stop("listValues must have a positive numerator weight and negative denominator weight", call. = FALSE)
  }
  as.numeric(listValues)
}

.rsdeseq2_resolve_list_results_contrast <- function(resultsNames,
                                                    positive,
                                                    negative,
                                                    positiveWeight,
                                                    negativeWeight) {
  if (length(positive) == 0L && length(negative) == 0L) {
    stop("list contrast must contain at least one coefficient name", call. = FALSE)
  }
  positiveIdx <- .rsdeseq2_resolve_results_name_list(resultsNames, positive)
  negativeIdx <- .rsdeseq2_resolve_results_name_list(resultsNames, negative)
  if (length(intersect(positiveIdx, negativeIdx)) != 0L) {
    stop("contrast list entries must not appear in both numerator and denominator", call. = FALSE)
  }
  numeric <- numeric(length(resultsNames))
  numeric[positiveIdx] <- positiveWeight
  numeric[negativeIdx] <- negativeWeight
  names(numeric) <- resultsNames
  .rsdeseq2_validate_numeric_results_contrast(numeric, resultsNames)
}

.rsdeseq2_resolve_results_name_list <- function(resultsNames, names) {
  if (length(names) == 0L) {
    return(integer())
  }
  idx <- vapply(names, .rsdeseq2_resolve_results_name_index, integer(1L), resultsNames = resultsNames)
  unique(idx)
}

.rsdeseq2_resolve_factor_level_results_contrast <- function(resultsNames,
                                                            factor,
                                                            numerator,
                                                            denominator,
                                                            reference) {
  if (identical(numerator, denominator)) {
    stop("contrast numerator and denominator levels must differ", call. = FALSE)
  }
  if (is.null(reference)) {
    direct <- .rsdeseq2_find_first_results_name(
      resultsNames,
      .rsdeseq2_standard_results_names(factor, numerator, denominator)
    )
    if (!is.na(direct)) {
      return(.rsdeseq2_unit_contrast(resultsNames, direct, 1))
    }
    reverse <- .rsdeseq2_find_first_results_name(
      resultsNames,
      .rsdeseq2_standard_results_names(factor, denominator, numerator)
    )
    if (!is.na(reverse)) {
      return(.rsdeseq2_unit_contrast(resultsNames, reverse, -1))
    }
    shared <- .rsdeseq2_find_shared_reference_results_names(resultsNames, factor, numerator, denominator)
    if (!is.null(shared)) {
      out <- numeric(length(resultsNames))
      out[[shared$numerator]] <- 1
      out[[shared$denominator]] <- -1
      names(out) <- resultsNames
      return(out)
    }
  } else {
    if (.rsdeseq2_same_factor_level_name(numerator, reference)) {
      denominatorIdx <- .rsdeseq2_find_first_results_name(
        resultsNames,
        c(
          .rsdeseq2_standard_results_names(factor, denominator, reference),
          .rsdeseq2_expanded_results_names(factor, denominator)
        )
      )
      if (!is.na(denominatorIdx)) {
        return(.rsdeseq2_unit_contrast(resultsNames, denominatorIdx, -1))
      }
    } else if (.rsdeseq2_same_factor_level_name(denominator, reference)) {
      numeratorIdx <- .rsdeseq2_find_first_results_name(
        resultsNames,
        c(
          .rsdeseq2_standard_results_names(factor, numerator, reference),
          .rsdeseq2_expanded_results_names(factor, numerator)
        )
      )
      if (!is.na(numeratorIdx)) {
        return(.rsdeseq2_unit_contrast(resultsNames, numeratorIdx, 1))
      }
    } else {
      numeratorIdx <- .rsdeseq2_find_first_results_name(
        resultsNames,
        .rsdeseq2_standard_results_names(factor, numerator, reference)
      )
      denominatorIdx <- .rsdeseq2_find_first_results_name(
        resultsNames,
        .rsdeseq2_standard_results_names(factor, denominator, reference)
      )
      if (!is.na(numeratorIdx) && !is.na(denominatorIdx)) {
        out <- numeric(length(resultsNames))
        out[[numeratorIdx]] <- 1
        out[[denominatorIdx]] <- -1
        names(out) <- resultsNames
        return(out)
      }
    }
  }

  numeratorIdx <- .rsdeseq2_find_first_results_name(
    resultsNames,
    .rsdeseq2_expanded_results_names(factor, numerator)
  )
  denominatorIdx <- .rsdeseq2_find_first_results_name(
    resultsNames,
    .rsdeseq2_expanded_results_names(factor, denominator)
  )
  if (!is.na(numeratorIdx) && !is.na(denominatorIdx)) {
    out <- numeric(length(resultsNames))
    out[[numeratorIdx]] <- 1
    out[[denominatorIdx]] <- -1
    names(out) <- resultsNames
    return(out)
  }

  stop(
    sprintf(
      "factor-level contrast %s: %s vs %s could not be resolved from resultsNames",
      factor,
      numerator,
      denominator
    ),
    call. = FALSE
  )
}

.rsdeseq2_same_factor_level_name <- function(left, right) {
  identical(left, right) || identical(make.names(left), make.names(right))
}

.rsdeseq2_unit_contrast <- function(resultsNames, index, value) {
  out <- numeric(length(resultsNames))
  out[[index]] <- value
  names(out) <- resultsNames
  out
}

.rsdeseq2_standard_results_names <- function(factor, level, reference) {
  c(
    sprintf("%s_%s_vs_%s", factor, level, reference),
    sprintf("%s%s_vs_%s", factor, level, reference),
    sprintf("%s_%s_vs_%s", make.names(factor), make.names(level), make.names(reference)),
    sprintf("%s%s_vs_%s", make.names(factor), make.names(level), make.names(reference))
  )
}

.rsdeseq2_expanded_results_names <- function(factor, level) {
  c(
    sprintf("%s%s", factor, level),
    sprintf("%s_%s", factor, level),
    sprintf("%s%s", make.names(factor), make.names(level)),
    sprintf("%s_%s", make.names(factor), make.names(level))
  )
}

.rsdeseq2_find_shared_reference_results_names <- function(resultsNames, factor, numerator, denominator) {
  numeratorRows <- .rsdeseq2_standard_results_rows_for_level(resultsNames, factor, numerator)
  denominatorRows <- .rsdeseq2_standard_results_rows_for_level(resultsNames, factor, denominator)
  shared <- list()
  for (numeratorRow in numeratorRows) {
    for (denominatorRow in denominatorRows) {
      if (identical(numeratorRow$reference, denominatorRow$reference)) {
        key <- sprintf("%d:%d", numeratorRow$index, denominatorRow$index)
        shared[[key]] <- list(numerator = numeratorRow$index, denominator = denominatorRow$index)
      }
    }
  }
  if (length(shared) == 0L) {
    return(NULL)
  }
  if (length(shared) > 1L) {
    stop(
      sprintf(
        "factor-level contrast %s: %s vs %s resolves ambiguously through shared-reference result names",
        factor,
        numerator,
        denominator
      ),
      call. = FALSE
    )
  }
  shared[[1L]]
}

.rsdeseq2_standard_results_rows_for_level <- function(resultsNames, factor, level) {
  prefixes <- .rsdeseq2_standard_results_prefixes(factor, level)
  rows <- list()
  for (idx in seq_along(resultsNames)) {
    name <- resultsNames[[idx]]
    for (prefix in prefixes) {
      if (startsWith(name, prefix)) {
        reference <- substring(name, nchar(prefix) + 1L)
        key <- sprintf("%d:%s", idx, reference)
        rows[[key]] <- list(index = idx, reference = reference)
      }
    }
  }
  rows
}

.rsdeseq2_standard_results_prefixes <- function(factor, level) {
  unique(c(
    sprintf("%s_%s_vs_", factor, level),
    sprintf("%s%s_vs_", factor, level),
    sprintf("%s_%s_vs_", make.names(factor), make.names(level)),
    sprintf("%s%s_vs_", make.names(factor), make.names(level))
  ))
}

.rsdeseq2_find_first_results_name <- function(resultsNames, candidates) {
  for (candidate in unique(candidates)) {
    idx <- tryCatch(
      .rsdeseq2_resolve_results_name_index(candidate, resultsNames),
      error = function(error) NA_integer_
    )
    if (!is.na(idx)) {
      return(idx)
    }
  }
  NA_integer_
}

.rsdeseq2_resolve_results_name_index <- function(name, resultsNames) {
  exact <- which(resultsNames == name)
  if (length(exact) == 1L) {
    return(exact)
  }
  if (length(exact) > 1L) {
    stop(sprintf("results name %s is duplicated", name), call. = FALSE)
  }

  alias <- .rsdeseq2_clean_results_name(name)
  aliases <- vapply(resultsNames, .rsdeseq2_clean_results_name, character(1L))
  matches <- which(aliases == alias)
  if (length(matches) == 1L) {
    return(matches)
  }
  if (length(matches) > 1L) {
    stop(sprintf("results name alias %s is ambiguous", name), call. = FALSE)
  }
  stop(sprintf("results name %s not found", name), call. = FALSE)
}

.rsdeseq2_clean_results_name <- function(name) {
  cleaned <- make.names(name)
  cleaned <- gsub("^X.Intercept\\.$", "Intercept", cleaned)
  cleaned <- gsub("^Intercept$", "Intercept", cleaned)
  cleaned
}

.rsdeseq2_list_contrast_comparison <- function(positive,
                                               negative,
                                               positiveWeight,
                                               negativeWeight) {
  if (length(negative) == 0L) {
    return(sprintf(
      "coefficient list contrast: %s effect",
      .rsdeseq2_weighted_name_list(positive, positiveWeight)
    ))
  }
  if (length(positive) == 0L) {
    return(sprintf(
      "coefficient list contrast: %s effect",
      .rsdeseq2_weighted_name_list(negative, negativeWeight)
    ))
  }
  sprintf(
    "coefficient list contrast: %s vs %s",
    .rsdeseq2_weighted_name_list(positive, positiveWeight),
    .rsdeseq2_weighted_name_list(negative, abs(negativeWeight))
  )
}

.rsdeseq2_weighted_name_list <- function(names, weight) {
  prefix <- if (identical(unname(weight), 1)) "" else sprintf("%s ", format(weight, trim = TRUE, scientific = FALSE))
  paste0(prefix, paste(names, collapse = " + "))
}
