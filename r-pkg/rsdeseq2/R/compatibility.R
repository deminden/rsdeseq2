.rsdeseq2_not_available <- function(feature) {
  stop(
    sprintf("%s is not implemented in the R wrapper yet; use the Rust crate for current initial stages.", feature),
    call. = FALSE
  )
}

.rsdeseq2_validate_counts <- function(counts) {
  if (!is.matrix(counts)) {
    stop("counts must be a matrix with genes in rows and samples in columns", call. = FALSE)
  }
  if (nrow(counts) == 0L) {
    stop("counts must contain at least one gene", call. = FALSE)
  }
  if (ncol(counts) == 0L) {
    stop("counts must contain at least one sample", call. = FALSE)
  }
  if (!is.numeric(counts) && !is.integer(counts)) {
    stop("counts must be an integer or numeric matrix", call. = FALSE)
  }
  if (any(!is.finite(counts))) {
    stop("counts must be finite", call. = FALSE)
  }
  if (any(counts < 0)) {
    stop("counts must be non-negative", call. = FALSE)
  }
  if (any(counts != floor(counts))) {
    stop("counts must contain integer-valued entries", call. = FALSE)
  }
  storage.mode(counts) <- "double"
  counts
}

.rsdeseq2_validate_size_factors <- function(sizeFactors, nSamples) {
  if (!is.numeric(sizeFactors) || length(sizeFactors) != nSamples) {
    stop("sizeFactors must be a numeric vector with one value per sample", call. = FALSE)
  }
  if (any(!is.finite(sizeFactors)) || any(sizeFactors <= 0)) {
    stop("sizeFactors must be finite and positive", call. = FALSE)
  }
  as.numeric(sizeFactors)
}

.rsdeseq2_validate_normalization_factors <- function(normalizationFactors, counts) {
  if (is.null(normalizationFactors)) {
    return(NULL)
  }
  if (!is.matrix(normalizationFactors) || !is.numeric(normalizationFactors)) {
    stop("normalizationFactors must be a numeric matrix", call. = FALSE)
  }
  if (nrow(normalizationFactors) != nrow(counts) || ncol(normalizationFactors) != ncol(counts)) {
    stop("normalizationFactors must have the same dimensions as counts", call. = FALSE)
  }
  if (any(!is.finite(normalizationFactors)) || any(normalizationFactors <= 0)) {
    stop("normalizationFactors must be finite and positive", call. = FALSE)
  }
  storage.mode(normalizationFactors) <- "double"
  normalizationFactors
}

.rsdeseq2_validate_observation_weights <- function(weights, counts) {
  if (is.null(weights)) {
    return(NULL)
  }
  if (!is.matrix(weights) || !is.numeric(weights)) {
    stop("weights must be a numeric matrix", call. = FALSE)
  }
  if (nrow(weights) != nrow(counts) || ncol(weights) != ncol(counts)) {
    stop("weights must have the same dimensions as counts", call. = FALSE)
  }
  if (any(!is.finite(weights)) || any(weights < 0)) {
    stop("weights must be finite and non-negative", call. = FALSE)
  }
  storage.mode(weights) <- "double"
  weights
}

.rsdeseq2_control_gene_indices <- function(controlGenes, counts) {
  if (is.null(controlGenes)) {
    return(seq_len(nrow(counts)))
  }
  if (is.logical(controlGenes)) {
    if (length(controlGenes) != nrow(counts)) {
      stop("logical controlGenes must have one value per gene", call. = FALSE)
    }
    if (anyNA(controlGenes)) {
      stop("logical controlGenes must not contain NA", call. = FALSE)
    }
    idx <- which(controlGenes)
  } else if (is.numeric(controlGenes)) {
    if (any(!is.finite(controlGenes)) || any(controlGenes != floor(controlGenes))) {
      stop("numeric controlGenes must contain integer row indices", call. = FALSE)
    }
    if (any(controlGenes < 1L) || any(controlGenes > nrow(counts))) {
      stop("numeric controlGenes must be 1-based row indices", call. = FALSE)
    }
    idx <- as.integer(controlGenes)
  } else if (is.character(controlGenes)) {
    if (anyNA(controlGenes)) {
      stop("character controlGenes must not contain NA", call. = FALSE)
    }
    rowNames <- rownames(counts)
    if (is.null(rowNames)) {
      stop("character controlGenes require count row names", call. = FALSE)
    }
    idx <- match(controlGenes, rowNames)
    if (any(is.na(idx))) {
      stop("all character controlGenes must be present in rownames(counts)", call. = FALSE)
    }
  } else {
    stop("controlGenes must be NULL, logical, numeric, or character", call. = FALSE)
  }
  if (length(idx) == 0L) {
    stop("controlGenes selects no genes", call. = FALSE)
  }
  idx
}

.rsdeseq2_validate_geo_means <- function(geoMeans, nGenes) {
  if (is.null(geoMeans)) {
    return(NULL)
  }
  if (!is.numeric(geoMeans) || length(geoMeans) != nGenes) {
    stop("geoMeans must be a numeric vector with one value per gene", call. = FALSE)
  }
  if (anyNA(geoMeans) || any(is.nan(geoMeans)) || any(geoMeans < 0)) {
    stop("geoMeans must be non-negative or infinite", call. = FALSE)
  }
  logGeoMeans <- rep(-Inf, nGenes)
  positive <- geoMeans > 0
  logGeoMeans[positive] <- log(geoMeans[positive])
  logGeoMeans
}

.rsdeseq2_stabilize_size_factors <- function(sizeFactors) {
  scale <- exp(mean(log(sizeFactors)))
  if (!is.finite(scale) || scale <= 0) {
    stop("cannot stabilize size factors to geometric mean one", call. = FALSE)
  }
  sizeFactors / scale
}

.rsdeseq2_validate_n_genes <- function(nGenes) {
  if (!is.numeric(nGenes) || length(nGenes) != 1L || !is.finite(nGenes) || nGenes != floor(nGenes) || nGenes < 0L) {
    stop("nGenes must be a single non-negative integer", call. = FALSE)
  }
  as.integer(nGenes)
}

.rsdeseq2_validate_row_names <- function(rowNames, nGenes) {
  if (is.null(rowNames)) {
    return(NULL)
  }
  if (!is.character(rowNames) || length(rowNames) != nGenes || anyNA(rowNames)) {
    stop("rowNames must be a character vector with one non-missing value per gene", call. = FALSE)
  }
  rowNames
}

.rsdeseq2_optional_vector <- function(value, nGenes, name, type = c("logical", "integer", "numeric")) {
  type <- match.arg(type)
  if (is.null(value)) {
    return(NULL)
  }
  if (length(value) != nGenes) {
    stop(sprintf("%s must have one value per gene", name), call. = FALSE)
  }
  if (type == "logical") {
    if (!is.logical(value)) {
      stop(sprintf("%s must be logical", name), call. = FALSE)
    }
    return(value)
  }
  if (type == "integer") {
    if (!is.numeric(value) || any(!is.na(value) & (value != floor(value) | value < 0))) {
      stop(sprintf("%s must contain non-negative integer values or NA", name), call. = FALSE)
    }
    return(as.integer(value))
  }
  if (!is.numeric(value)) {
    stop(sprintf("%s must be numeric", name), call. = FALSE)
  }
  as.numeric(value)
}

.rsdeseq2_validate_pvalues <- function(pvalue) {
  if (!is.numeric(pvalue) || length(pvalue) == 0L) {
    stop("pvalue must be a non-empty numeric vector", call. = FALSE)
  }
  invalid <- !is.na(pvalue) & (!is.finite(pvalue) | pvalue < 0 | pvalue > 1)
  if (any(invalid)) {
    stop("pvalue must contain values in [0, 1] or NA", call. = FALSE)
  }
  out <- as.numeric(pvalue)
  names(out) <- names(pvalue)
  out
}

.rsdeseq2_validate_result_numeric_vector <- function(value, name, allowNA = TRUE) {
  if (!is.numeric(value) || length(value) == 0L) {
    stop(sprintf("%s must be a non-empty numeric vector", name), call. = FALSE)
  }
  if (any(!is.na(value) & !is.finite(value))) {
    stop(sprintf("%s must contain finite values or NA", name), call. = FALSE)
  }
  if (!allowNA && anyNA(value)) {
    stop(sprintf("%s must not contain NA", name), call. = FALSE)
  }
  out <- as.numeric(value)
  names(out) <- names(value)
  out
}

.rsdeseq2_validate_result_optional_numeric <- function(value, nGenes, name) {
  if (is.null(value)) {
    return(rep(NA_real_, nGenes))
  }
  if (!is.numeric(value) || length(value) != nGenes) {
    stop(sprintf("%s must be a numeric vector with one value per gene", name), call. = FALSE)
  }
  if (any(!is.na(value) & !is.finite(value))) {
    stop(sprintf("%s must contain finite values or NA", name), call. = FALSE)
  }
  as.numeric(value)
}

.rsdeseq2_validate_result_optional_logical <- function(value, nGenes, name) {
  if (!is.logical(value) || length(value) != nGenes) {
    stop(sprintf("%s must be a logical vector with one value per gene", name), call. = FALSE)
  }
  value
}

.rsdeseq2_validate_result_pvalues <- function(value, nGenes, name = "pvalue") {
  if (!is.numeric(value) || length(value) != nGenes) {
    stop(sprintf("%s must be a numeric vector with one value per gene", name), call. = FALSE)
  }
  invalid <- !is.na(value) & (!is.finite(value) | value < 0 | value > 1)
  if (any(invalid)) {
    stop(sprintf("%s must contain values in [0, 1] or NA", name), call. = FALSE)
  }
  as.numeric(value)
}

.rsdeseq2_validate_result_table <- function(results) {
  if (!is.data.frame(results)) {
    stop("results must be a data frame", call. = FALSE)
  }
  required <- c("baseMean", "pvalue", "padj")
  missing <- setdiff(required, colnames(results))
  if (length(missing) != 0L) {
    stop(sprintf("results is missing required columns: %s", paste(missing, collapse = ", ")), call. = FALSE)
  }
  nGenes <- nrow(results)
  if (nGenes == 0L) {
    stop("results must contain at least one row", call. = FALSE)
  }
  results$baseMean <- .rsdeseq2_validate_result_numeric_vector(results$baseMean, "baseMean", allowNA = FALSE)
  results$pvalue <- .rsdeseq2_validate_result_pvalues(results$pvalue, nGenes)
  results$padj <- .rsdeseq2_validate_result_pvalues(results$padj, nGenes, name = "padj")
  results
}

.rsdeseq2_validate_independent_filter_alpha <- function(alpha) {
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must be a single finite numeric value in (0, 1)", call. = FALSE)
  }
  as.numeric(alpha)
}

.rsdeseq2_validate_independent_filter_theta <- function(theta, filter) {
  if (any(!is.finite(filter))) {
    stop("baseMean must contain finite values for independent filtering", call. = FALSE)
  }
  if (is.null(theta)) {
    lower <- sum(filter == 0) / length(filter)
    upper <- if (lower < 0.95) 0.95 else 1.0
    return(seq(lower, upper, length.out = 50L))
  }
  if (!is.numeric(theta) || length(theta) <= 1L || any(!is.finite(theta)) || any(theta < 0 | theta > 1)) {
    stop("theta must be a numeric vector of at least two values in [0, 1]", call. = FALSE)
  }
  as.numeric(theta)
}

.rsdeseq2_select_independent_filter_index <- function(numRej, lowessFit) {
  if (length(numRej) == 0L || max(numRej) <= 10L) {
    return(1L)
  }
  residuals <- numRej[numRej > 0L] - lowessFit[numRej > 0L]
  rmse <- if (length(residuals) == 0L) 0 else sqrt(mean(residuals^2))
  maxFit <- max(lowessFit)
  selected <- which(numRej > (maxFit - rmse))[1]
  if (!is.na(selected)) {
    return(selected)
  }
  selected <- which(numRej > (0.9 * maxFit))[1]
  if (!is.na(selected)) {
    return(selected)
  }
  selected <- which(numRej > (0.8 * maxFit))[1]
  if (!is.na(selected)) {
    return(selected)
  }
  1L
}

.rsdeseq2_validate_optional_numeric_vector <- function(value, n, name) {
  if (!is.numeric(value) || length(value) != n) {
    stop(sprintf("%s must be a numeric vector with one value per p-value", name), call. = FALSE)
  }
  if (any(!is.na(value) & !is.finite(value))) {
    stop(sprintf("%s must contain finite values or NA", name), call. = FALSE)
  }
  as.numeric(value)
}

.rsdeseq2_validate_cooks_cutoff <- function(cooksCutoff) {
  if (is.null(cooksCutoff) || identical(cooksCutoff, FALSE)) {
    return(NULL)
  }
  if (!is.numeric(cooksCutoff) || length(cooksCutoff) != 1L || !is.finite(cooksCutoff)) {
    stop("cooksCutoff must be a single finite numeric value, FALSE, or NULL", call. = FALSE)
  }
  as.numeric(cooksCutoff)
}

.rsdeseq2_validate_cooks_matrix <- function(cooks, nGenes, nSamples) {
  if (!is.matrix(cooks) || !is.numeric(cooks)) {
    stop("cooks must be a numeric matrix", call. = FALSE)
  }
  if (nrow(cooks) != nGenes || ncol(cooks) != nSamples) {
    stop("cooks must have the same dimensions as counts", call. = FALSE)
  }
  if (any(!is.na(cooks) & !is.finite(cooks))) {
    stop("cooks must contain finite values or NA", call. = FALSE)
  }
  storage.mode(cooks) <- "double"
  cooks
}

.rsdeseq2_low_count_cooks_spares_row <- function(countsRow, cooksRow) {
  if (all(is.na(cooksRow))) {
    return(FALSE)
  }
  maxIdx <- which.max(cooksRow)
  if (length(maxIdx) == 0L || is.na(cooksRow[[maxIdx]]) || !is.finite(cooksRow[[maxIdx]])) {
    return(FALSE)
  }
  outCount <- countsRow[[maxIdx]]
  sum(countsRow > outCount) >= 3L
}

.rsdeseq2_cooks_cutoff_output <- function(pvalue, maxCooks, cooksOutlier) {
  out <- data.frame(
    pvalue = pvalue,
    padj = stats::p.adjust(pvalue, method = "BH"),
    maxCooks = maxCooks,
    cooksOutlier = cooksOutlier,
    check.names = FALSE
  )
  rowNames <- names(pvalue)
  if (!is.null(rowNames)) {
    rownames(out) <- rowNames
  }
  out
}

.rsdeseq2_diagnostic_schema_names <- function(native = FALSE) {
  fallback <- c(
    "betaConv",
    "fullBetaConv",
    "reducedBetaConv",
    "betaIter",
    "reducedBetaIter",
    "deviance",
    "maxCooks"
  )
  if (isTRUE(native)) {
    nativeSchema <- tryCatch(
      .Call("rsdeseq2_diagnostic_schema", PACKAGE = "rsdeseq2"),
      error = function(error) NULL
    )
    if (!is.null(nativeSchema)) {
      return(nativeSchema)
    }
    nativeSchema <- tryCatch(
      .Call("rsdeseq2_diagnostic_schema"),
      error = function(error) NULL
    )
    if (!is.null(nativeSchema)) {
      return(nativeSchema)
    }
  }
  fallback
}

rsdeseq2DiagnosticSchemaRust <- function(native = FALSE) {
  .rsdeseq2_diagnostic_schema_names(native = native)
}

deseq2McolsDiagnosticsRust <- function(nGenes,
                                       test = c("none", "Wald", "LRT"),
                                       rowNames = NULL,
                                       betaConv = NULL,
                                       fullBetaConv = NULL,
                                       reducedBetaConv = NULL,
                                       betaIter = NULL,
                                       reducedBetaIter = NULL,
                                       deviance = NULL,
                                       maxCooks = NULL) {
  nGenes <- .rsdeseq2_validate_n_genes(nGenes)
  test <- match.arg(test)
  rowNames <- .rsdeseq2_validate_row_names(rowNames, nGenes)

  columns <- list()
  if (test == "Wald") {
    columns$betaConv <- .rsdeseq2_optional_vector(betaConv, nGenes, "betaConv", "logical")
  } else if (test == "LRT") {
    columns$fullBetaConv <- .rsdeseq2_optional_vector(fullBetaConv, nGenes, "fullBetaConv", "logical")
    columns$reducedBetaConv <- .rsdeseq2_optional_vector(reducedBetaConv, nGenes, "reducedBetaConv", "logical")
  }
  columns$betaIter <- .rsdeseq2_optional_vector(betaIter, nGenes, "betaIter", "integer")
  columns$reducedBetaIter <- .rsdeseq2_optional_vector(reducedBetaIter, nGenes, "reducedBetaIter", "integer")
  columns$deviance <- .rsdeseq2_optional_vector(deviance, nGenes, "deviance", "numeric")
  columns$maxCooks <- .rsdeseq2_optional_vector(maxCooks, nGenes, "maxCooks", "numeric")
  columns <- columns[!vapply(columns, is.null, logical(1))]

  if (length(columns) == 0L) {
    return(data.frame(row.names = rowNames %||% seq_len(nGenes)))
  }
  out <- as.data.frame(columns, optional = TRUE)
  if (!is.null(rowNames)) {
    rownames(out) <- rowNames
  }
  out
}

`%||%` <- function(lhs, rhs) {
  if (is.null(lhs)) {
    rhs
  } else {
    lhs
  }
}
