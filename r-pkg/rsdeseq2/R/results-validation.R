# Validation and alignment of fitted fields and design metadata.

.rsdeseq2_validate_beta_matrix <- function(beta, resultsNames) {
  if (!is.matrix(beta) || !is.numeric(beta) || nrow(beta) == 0L || ncol(beta) == 0L) {
    stop("beta must be a non-empty numeric matrix", call. = FALSE)
  }
  if (any(!is.na(beta) & !is.finite(beta))) {
    stop("beta must contain finite values or NA", call. = FALSE)
  }
  if (is.null(resultsNames)) {
    resultsNames <- colnames(beta)
  }
  resultsNames <- .rsdeseq2_validate_results_names(resultsNames)
  if (length(resultsNames) != ncol(beta)) {
    stop("resultsNames must have one value per beta column", call. = FALSE)
  }
  betaColumns <- colnames(beta)
  if (!is.null(betaColumns) && .rsdeseq2_has_explicit_row_names(betaColumns, ncol(beta))) {
    if (length(betaColumns) != ncol(beta) || anyNA(betaColumns) || any(!nzchar(betaColumns))) {
      stop("beta column names must be non-empty when supplied", call. = FALSE)
    }
    if (anyDuplicated(betaColumns)) {
      stop("beta column names must be unique when supplied", call. = FALSE)
    }
    order <- vapply(
      resultsNames,
      .rsdeseq2_resolve_results_name_index,
      integer(1L),
      resultsNames = betaColumns
    )
    beta <- beta[, order, drop = FALSE]
  }
  colnames(beta) <- resultsNames
  storage.mode(beta) <- "double"
  beta
}

.rsdeseq2_validate_beta_covariance <- function(betaCovariance, nGenes, resultsNames) {
  if (is.null(betaCovariance)) {
    return(NULL)
  }
  nCoef <- length(resultsNames)
  if (!is.array(betaCovariance) || !is.numeric(betaCovariance) || length(dim(betaCovariance)) != 3L) {
    stop("betaCovariance must be a numeric array with dimensions genes x coefficients x coefficients", call. = FALSE)
  }
  if (!identical(dim(betaCovariance), c(nGenes, nCoef, nCoef))) {
    stop("betaCovariance must have dimensions genes x coefficients x coefficients", call. = FALSE)
  }
  if (any(!is.na(betaCovariance) & !is.finite(betaCovariance))) {
    stop("betaCovariance must contain finite values or NA", call. = FALSE)
  }
  betaCovariance <- .rsdeseq2_align_beta_covariance_axis(betaCovariance, resultsNames, 2L)
  betaCovariance <- .rsdeseq2_align_beta_covariance_axis(betaCovariance, resultsNames, 3L)
  names <- dimnames(betaCovariance)
  if (is.null(names)) {
    names <- vector("list", 3L)
  }
  names[[2L]] <- resultsNames
  names[[3L]] <- resultsNames
  dimnames(betaCovariance) <- names
  betaCovariance
}

.rsdeseq2_align_beta_covariance_axis <- function(betaCovariance, resultsNames, axis) {
  axisNames <- dimnames(betaCovariance)[[axis]]
  if (is.null(axisNames) || !.rsdeseq2_has_explicit_row_names(axisNames, dim(betaCovariance)[[axis]])) {
    return(betaCovariance)
  }
  if (length(axisNames) != dim(betaCovariance)[[axis]] || anyNA(axisNames) || any(!nzchar(axisNames))) {
    stop("betaCovariance coefficient names must be non-empty when supplied", call. = FALSE)
  }
  if (anyDuplicated(axisNames)) {
    stop("betaCovariance coefficient names must be unique when supplied", call. = FALSE)
  }
  order <- vapply(
    resultsNames,
    .rsdeseq2_resolve_results_name_index,
    integer(1L),
    resultsNames = axisNames
  )
  if (axis == 2L) {
    return(betaCovariance[, order, , drop = FALSE])
  }
  betaCovariance[, , order, drop = FALSE]
}

.rsdeseq2_validate_beta_se_matrix <- function(betaSE, nGenes, resultsNames) {
  if (is.null(betaSE)) {
    return(NULL)
  }
  if (!is.matrix(betaSE) || !is.numeric(betaSE) || nrow(betaSE) != nGenes || ncol(betaSE) != length(resultsNames)) {
    stop("betaSE must be a numeric matrix with the same shape as beta", call. = FALSE)
  }
  if (any(!is.na(betaSE) & (!is.finite(betaSE) | betaSE < 0))) {
    stop("betaSE must contain non-negative finite values or NA", call. = FALSE)
  }
  betaSEColumns <- colnames(betaSE)
  if (!is.null(betaSEColumns)) {
    if (length(betaSEColumns) != ncol(betaSE) || anyNA(betaSEColumns) || any(!nzchar(betaSEColumns))) {
      stop("betaSE column names must be non-empty when supplied", call. = FALSE)
    }
    if (anyDuplicated(betaSEColumns)) {
      stop("betaSE column names must be unique when supplied", call. = FALSE)
    }
    order <- vapply(
      resultsNames,
      .rsdeseq2_resolve_results_name_index,
      integer(1L),
      resultsNames = betaSEColumns
    )
    betaSE <- betaSE[, order, drop = FALSE]
  }
  colnames(betaSE) <- resultsNames
  storage.mode(betaSE) <- "double"
  betaSE
}

.rsdeseq2_validate_optional_contrast_counts <- function(counts, nGenes) {
  if (is.null(counts)) {
    return(NULL)
  }
  counts <- .rsdeseq2_validate_counts(counts)
  if (nrow(counts) != nGenes) {
    stop("counts must have one row per beta row", call. = FALSE)
  }
  counts
}

.rsdeseq2_validate_optional_sample_levels <- function(sampleLevels, counts) {
  if (is.null(sampleLevels)) {
    return(NULL)
  }
  if (is.null(counts)) {
    stop("sampleLevels require counts", call. = FALSE)
  }
  if (!is.character(sampleLevels) || length(sampleLevels) != ncol(counts) || anyNA(sampleLevels) || any(!nzchar(sampleLevels))) {
    stop("sampleLevels must be a character vector with one value per sample", call. = FALSE)
  }
  sampleNames <- colnames(counts)
  levelNames <- names(sampleLevels)
  if (!is.null(sampleNames) && !is.null(levelNames)) {
    if (anyNA(levelNames) || any(!nzchar(levelNames)) || anyDuplicated(levelNames)) {
      stop("sampleLevels names must be unique non-empty sample names", call. = FALSE)
    }
    missing <- setdiff(sampleNames, levelNames)
    if (length(missing) > 0L) {
      stop("sampleLevels names must contain all count sample names", call. = FALSE)
    }
    sampleLevels <- sampleLevels[sampleNames]
  }
  sampleLevels
}

.rsdeseq2_align_named_gene_vector <- function(values, geneNames, name, allowScalar = FALSE) {
  valueNames <- names(values)
  if (is.null(geneNames) || is.null(valueNames)) {
    return(values)
  }
  if (isTRUE(allowScalar) && length(values) == 1L && length(geneNames) != 1L) {
    return(values)
  }
  if (anyNA(valueNames) || any(!nzchar(valueNames)) || anyDuplicated(valueNames)) {
    stop(sprintf("%s names must be unique non-empty gene names", name), call. = FALSE)
  }
  missing <- setdiff(geneNames, valueNames)
  if (length(missing) > 0L) {
    stop(sprintf("%s names must contain all beta row names", name), call. = FALSE)
  }
  values[geneNames]
}

.rsdeseq2_align_named_gene_matrix <- function(values, geneNames, name) {
  valueNames <- rownames(values)
  if (is.null(values) || is.null(geneNames) || is.null(valueNames) ||
    !.rsdeseq2_has_explicit_row_names(valueNames, nrow(values))) {
    return(values)
  }
  if (anyNA(valueNames) || any(!nzchar(valueNames)) || anyDuplicated(valueNames)) {
    stop(sprintf("%s row names must be unique non-empty gene names", name), call. = FALSE)
  }
  missing <- setdiff(geneNames, valueNames)
  if (length(missing) > 0L) {
    stop(sprintf("%s row names must contain all beta row names", name), call. = FALSE)
  }
  values[geneNames, , drop = FALSE]
}

.rsdeseq2_align_named_gene_array <- function(values, geneNames, name) {
  valueNames <- dimnames(values)[[1L]]
  if (is.null(values) || is.null(geneNames) || is.null(valueNames) ||
    !.rsdeseq2_has_explicit_row_names(valueNames, dim(values)[[1L]])) {
    return(values)
  }
  if (anyNA(valueNames) || any(!nzchar(valueNames)) || anyDuplicated(valueNames)) {
    stop(sprintf("%s row names must be unique non-empty gene names", name), call. = FALSE)
  }
  missing <- setdiff(geneNames, valueNames)
  if (length(missing) > 0L) {
    stop(sprintf("%s row names must contain all beta row names", name), call. = FALSE)
  }
  values[geneNames, , , drop = FALSE]
}

.rsdeseq2_validate_optional_col_data <- function(colData, counts) {
  if (is.null(colData)) {
    return(NULL)
  }
  colData <- .rsdeseq2_coerce_col_data_frame(colData)
  if (is.null(colData)) {
    stop("colData must be a data.frame or coercible with as.data.frame", call. = FALSE)
  }
  if (!is.null(counts) && nrow(colData) != ncol(counts)) {
    stop("colData must have one row per sample", call. = FALSE)
  }
  if (!is.null(counts)) {
    sampleNames <- colnames(counts)
    colDataNames <- rownames(colData)
    if (.rsdeseq2_has_explicit_row_names(colDataNames, nrow(colData)) &&
      !is.null(sampleNames)) {
      if (anyNA(colDataNames) || any(!nzchar(colDataNames)) || anyDuplicated(colDataNames)) {
        stop("colData row names must be unique non-empty sample names", call. = FALSE)
      }
      missing <- setdiff(sampleNames, colDataNames)
      if (length(missing) > 0L) {
        stop("colData row names must contain all count sample names", call. = FALSE)
      }
      colData <- colData[sampleNames, , drop = FALSE]
    }
  }
  colData
}

.rsdeseq2_coerce_col_data_frame <- function(colData) {
  if (is.null(colData)) {
    return(NULL)
  }
  if (is.data.frame(colData)) {
    return(colData)
  }
  coerced <- tryCatch(
    as.data.frame(colData),
    error = function(e) NULL
  )
  if (is.null(coerced) || !is.data.frame(coerced)) {
    return(NULL)
  }
  coerced
}

.rsdeseq2_has_explicit_row_names <- function(rowNames, nRows) {
  if (is.null(rowNames) || length(rowNames) != nRows) {
    return(FALSE)
  }
  !identical(rowNames, as.character(seq_len(nRows)))
}

.rsdeseq2_validate_optional_factor_levels <- function(factorLevels) {
  if (is.null(factorLevels)) {
    return(NULL)
  }
  if (!is.list(factorLevels) || is.null(names(factorLevels)) || anyNA(names(factorLevels)) || any(!nzchar(names(factorLevels)))) {
    stop("factorLevels must be a named list of character vectors", call. = FALSE)
  }
  out <- lapply(factorLevels, function(levels) {
    if (!is.character(levels) || length(levels) == 0L || anyNA(levels) || any(!nzchar(levels))) {
      stop("factorLevels entries must be non-empty character vectors", call. = FALSE)
    }
    if (anyDuplicated(levels)) {
      stop("factorLevels entries must contain unique levels", call. = FALSE)
    }
    if (anyDuplicated(make.names(levels))) {
      stop("factorLevels entries must resolve to unique R-cleaned level aliases", call. = FALSE)
    }
    levels
  })
  if (anyDuplicated(names(out))) {
    stop("factorLevels names must be unique", call. = FALSE)
  }
  cleanedNames <- make.names(names(out))
  if (anyDuplicated(cleanedNames)) {
    stop("factorLevels names must resolve to unique R-cleaned aliases", call. = FALSE)
  }
  out
}

.rsdeseq2_validate_optional_factor_references <- function(factorReferences, factorLevels = NULL) {
  if (is.null(factorReferences)) {
    return(NULL)
  }
  if (is.list(factorReferences)) {
    if (is.null(names(factorReferences)) || anyNA(names(factorReferences)) || any(!nzchar(names(factorReferences)))) {
      stop("factorReferences must be a named character vector or named list", call. = FALSE)
    }
    out <- vapply(factorReferences, function(reference) {
      if (!is.character(reference) || length(reference) != 1L || is.na(reference) || !nzchar(reference)) {
        stop("factorReferences entries must be single non-empty character values", call. = FALSE)
      }
      reference
    }, character(1L))
  } else {
    if (!is.character(factorReferences) || is.null(names(factorReferences)) || anyNA(names(factorReferences)) || any(!nzchar(names(factorReferences)))) {
      stop("factorReferences must be a named character vector or named list", call. = FALSE)
    }
    if (anyNA(factorReferences) || any(!nzchar(factorReferences))) {
      stop("factorReferences entries must be single non-empty character values", call. = FALSE)
    }
    out <- factorReferences
  }
  if (anyDuplicated(names(out))) {
    stop("factorReferences names must be unique", call. = FALSE)
  }
  if (!is.null(factorLevels)) {
    canonicalNames <- vapply(names(out), function(factor) {
      .rsdeseq2_resolve_named_metadata_field(factor, names(factorLevels), "factor reference")
    }, character(1L))
    matchedCanonicalNames <- canonicalNames[!is.na(canonicalNames)]
    if (anyDuplicated(matchedCanonicalNames)) {
      stop("factorReferences aliases must resolve to unique factor names", call. = FALSE)
    }
    for (idx in seq_along(out)) {
      factor <- canonicalNames[[idx]]
      if (is.na(factor) || is.null(factorLevels[[factor]])) {
        next
      }
      out[[idx]] <- .rsdeseq2_resolve_sample_level_alias(
        out[[idx]],
        factorLevels[[factor]],
        "reference"
      )
      names(out)[[idx]] <- factor
    }
  }
  out
}

.rsdeseq2_validate_optional_model_matrix <- function(modelMatrix, counts, resultsNames) {
  if (is.null(modelMatrix)) {
    return(NULL)
  }
  if (is.null(counts)) {
    stop("modelMatrix requires counts", call. = FALSE)
  }
  if (!is.matrix(modelMatrix) || !is.numeric(modelMatrix) || nrow(modelMatrix) != ncol(counts) || ncol(modelMatrix) != length(resultsNames)) {
    stop("modelMatrix must be a numeric matrix with samples in rows and one column per result name", call. = FALSE)
  }
  if (any(!is.finite(modelMatrix))) {
    stop("modelMatrix must contain finite values", call. = FALSE)
  }
  sampleNames <- colnames(counts)
  modelMatrixNames <- rownames(modelMatrix)
  if (.rsdeseq2_has_explicit_row_names(modelMatrixNames, nrow(modelMatrix)) &&
    !is.null(sampleNames)) {
    if (anyNA(modelMatrixNames) || any(!nzchar(modelMatrixNames)) || anyDuplicated(modelMatrixNames)) {
      stop("modelMatrix row names must be unique non-empty sample names", call. = FALSE)
    }
    missing <- setdiff(sampleNames, modelMatrixNames)
    if (length(missing) > 0L) {
      stop("modelMatrix row names must contain all count sample names", call. = FALSE)
    }
    modelMatrix <- modelMatrix[sampleNames, , drop = FALSE]
  }
  modelMatrixColumns <- colnames(modelMatrix)
  if (!is.null(modelMatrixColumns)) {
    if (length(modelMatrixColumns) != ncol(modelMatrix) || anyNA(modelMatrixColumns) || any(!nzchar(modelMatrixColumns))) {
      stop("modelMatrix column names must be non-empty when supplied", call. = FALSE)
    }
    if (anyDuplicated(modelMatrixColumns)) {
      stop("modelMatrix column names must be unique when supplied", call. = FALSE)
    }
    order <- vapply(
      resultsNames,
      .rsdeseq2_resolve_results_name_index,
      integer(1L),
      resultsNames = modelMatrixColumns
    )
    modelMatrix <- modelMatrix[, order, drop = FALSE]
  }
  colnames(modelMatrix) <- resultsNames
  storage.mode(modelMatrix) <- "double"
  modelMatrix
}
