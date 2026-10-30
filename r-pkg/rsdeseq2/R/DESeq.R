.rsdeseq2_align_fit_samples <- function(matrix, counts, label) {
  if (!is.null(rownames(matrix)) && !is.null(colnames(counts))) {
    automaticNames <- identical(rownames(matrix), as.character(seq_len(nrow(matrix))))
    if (automaticNames) return(matrix)
    if (anyDuplicated(rownames(matrix)) || anyDuplicated(colnames(counts))) {
      stop(sprintf("%s row names and count column names must be unique", label), call. = FALSE)
    }
    order <- match(colnames(counts), rownames(matrix))
    if (anyNA(order)) stop(sprintf("%s row names must match count columns", label), call. = FALSE)
    matrix <- matrix[order, , drop = FALSE]
  }
  matrix
}

.rsdeseq2_align_fit_observations <- function(matrix, counts, label) {
  if (!is.null(rownames(matrix)) && !is.null(rownames(counts))) {
    if (anyDuplicated(rownames(matrix)) || anyDuplicated(rownames(counts))) {
      stop(sprintf("%s row names and count row names must be unique", label), call. = FALSE)
    }
    order <- match(rownames(counts), rownames(matrix))
    if (anyNA(order)) stop(sprintf("%s row names must match count rows", label), call. = FALSE)
    matrix <- matrix[order, , drop = FALSE]
  }
  if (!is.null(colnames(matrix)) && !is.null(colnames(counts))) {
    if (anyDuplicated(colnames(matrix)) || anyDuplicated(colnames(counts))) {
      stop(sprintf("%s column names and count column names must be unique", label), call. = FALSE)
    }
    order <- match(colnames(counts), colnames(matrix))
    if (anyNA(order)) stop(sprintf("%s column names must match count columns", label), call. = FALSE)
    matrix <- matrix[, order, drop = FALSE]
  }
  matrix
}

DESeq <- function(object,
                  fitType = c("parametric", "local", "mean"),
                  test = c("Wald", "LRT"),
                  engine = "rust",
                  parityMode = "strict",
                  nproc = 1L,
                  design = NULL,
                  reduced = NULL,
                  sizeFactors = NULL,
                  normalizationFactors = NULL,
                  weights = NULL,
                  minReplicatesForReplace = 7L,
                  ...) {
  explicitSize <- !missing(sizeFactors)
  explicitFactors <- !missing(normalizationFactors)
  explicitWeights <- !missing(weights)
  normalizationState <- NULL
  isDataSet <- inherits(object, "DESeqDataSet")
  fitType <- match.arg(fitType)
  test <- match.arg(test)
  if (!identical(engine, "rust")) {
    stop("engine must be 'rust'", call. = FALSE)
  }
  parityMode <- match.arg(parityMode, c("strict", "fast", "diagnostic"))
  .rsdeseq2_check_unused_arguments(list(...), "DESeq")
  if (length(nproc) != 1L || is.na(nproc) || !is.finite(nproc) || nproc < 1L || nproc > .Machine$integer.max || nproc != floor(nproc)) {
    stop("nproc must be a positive integer", call. = FALSE)
  }
  disableReplacement <- is.numeric(minReplicatesForReplace) &&
    length(minReplicatesForReplace) == 1L &&
    is.infinite(minReplicatesForReplace) && minReplicatesForReplace > 0
  validReplacementMinimum <- is.numeric(minReplicatesForReplace) &&
    length(minReplicatesForReplace) == 1L &&
    is.finite(minReplicatesForReplace) && minReplicatesForReplace >= 3 &&
    minReplicatesForReplace <= .Machine$integer.max &&
    minReplicatesForReplace == floor(minReplicatesForReplace)
  if (!disableReplacement && !validReplacementMinimum) {
    stop("minReplicatesForReplace must be an integer of at least 3 or Inf", call. = FALSE)
  }
  colData <- NULL
  designOrigin <- list(kind="matrix",variables=character())
  if (inherits(object, "DESeqDataSet")) {
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
      stop("SummarizedExperiment is required for DESeqDataSet input", call. = FALSE)
    }
    # Anonymous genes need stable identities for fitted-state/subset mapping.
    # This is interface metadata; the caller's dataset remains unchanged.
    if (is.null(rownames(object))) rownames(object) <- as.character(seq_len(nrow(object)))
    counts <- SummarizedExperiment::assay(object, "counts")
    colData <- as.data.frame(SummarizedExperiment::colData(object))
    if (requireNamespace("DESeq2", quietly = TRUE)) {
      if (!explicitSize) sizeFactors <- DESeq2::sizeFactors(object)
      if (!explicitFactors) normalizationFactors <- DESeq2::normalizationFactors(object)
    }
    if (!explicitWeights && "weights" %in% SummarizedExperiment::assayNames(object)) {
      weights <- SummarizedExperiment::assay(object, "weights")
    }
    formula <- methods::slot(object, "design")
    if (!inherits(formula, "formula")) {
      stop("DESeqDataSet design must be a formula", call. = FALSE)
    }
    designOrigin <- .rsdeseq2_formula_origin(formula)
    modelFrame <- .rsdeseq2_model_frame(formula, data = colData)
    design <- .rsdeseq2_model_matrix(formula, data = modelFrame)
    formulaOffsets <- .rsdeseq2_model_offsets(modelFrame)
    prepared <- .rsdeseq2_prepare_normalization(
      object, .rsdeseq2_validate_counts(counts), sizeFactors, normalizationFactors,
      formulaOffsets, explicitSize, explicitFactors)
    sizeFactors <- prepared$sizeFactors
    normalizationFactors <- prepared$normalizationFactors
    normalizationState <- prepared$state
    if (test == "LRT") {
      if (is.null(reduced) || !inherits(reduced, "formula")) {
        stop("LRT requires reduced = a formula", call. = FALSE)
      }
      reducedFrame <- .rsdeseq2_model_frame(reduced, data = colData)
      reducedOffsets <- .rsdeseq2_model_offsets(reducedFrame)
      if (!is.null(reducedOffsets) && !identical(as.numeric(reducedOffsets), as.numeric(formulaOffsets))) {
        stop("LRT full and reduced formulas must use the same normalization offsets", call. = FALSE)
      }
      reduced <- .rsdeseq2_model_matrix(reduced, data = reducedFrame)
    }
  } else if (is.list(object) && !is.null(object$counts)) {
    counts <- object$counts
    colData <- object$colData
    if (is.null(design)) {
      design <- object$design
      if(!is.null(object$designOrigin)) designOrigin <- object$designOrigin
    }
    if (is.null(reduced)) reduced <- object$reduced
    if (!explicitSize) sizeFactors <- object$sizeFactors
    if (!explicitFactors) normalizationFactors <- object$normalizationFactors
    if (!explicitWeights) weights <- object$weights
  } else {
    counts <- object
  }
  counts <- .rsdeseq2_validate_counts(counts)
  designOrigin <- .rsdeseq2_result_design_origin(list(schemaVersion=2L,designOrigin=designOrigin))
  if (!is.null(rownames(counts)) && anyDuplicated(rownames(counts))) {
    stop("count row names must be unique", call. = FALSE)
  }
  if (!is.null(colnames(counts)) && anyDuplicated(colnames(counts))) {
    stop("count column names must be unique", call. = FALSE)
  }
  design <- .rsdeseq2_align_fit_samples(design, counts, "design")
  if (test == "LRT") reduced <- .rsdeseq2_align_fit_samples(reduced, counts, "reduced design")
  normalizationFactors <- .rsdeseq2_validate_normalization_factors(normalizationFactors, counts)
  weights <- .rsdeseq2_validate_observation_weights(weights, counts)
  if (!is.null(normalizationFactors)) {
    normalizationFactors <- .rsdeseq2_align_fit_observations(normalizationFactors, counts, "normalizationFactors")
  }
  if (!is.null(weights)) {
    weights <- .rsdeseq2_align_fit_observations(weights, counts, "weights")
  }
  if (is.null(normalizationFactors)) {
    if (!is.null(sizeFactors)) {
      if (!is.numeric(sizeFactors) || length(sizeFactors) != ncol(counts)) {
        stop("sizeFactors must be a numeric vector with one value per sample", call. = FALSE)
      }
      if (!is.null(names(sizeFactors)) && !is.null(colnames(counts))) {
        if (anyDuplicated(names(sizeFactors)) || anyDuplicated(colnames(counts))) {
          stop("named sizeFactors and count columns must have unique names", call. = FALSE)
        }
        sampleOrder <- match(colnames(counts), names(sizeFactors))
        if (anyNA(sampleOrder)) stop("sizeFactors names must match count columns", call. = FALSE)
        sizeFactors <- sizeFactors[sampleOrder]
      }
      sizeFactors <- .rsdeseq2_validate_size_factors(sizeFactors, ncol(counts))
    }
  } else {
    sizeFactors <- NULL
  }
  if (!is.matrix(design) || nrow(design) != ncol(counts) || !is.numeric(design)) {
    stop("design must be a numeric model matrix with one row per sample", call. = FALSE)
  }
  storage.mode(design) <- "double"
  if (any(!is.finite(design))) stop("design must contain finite values", call. = FALSE)
  if (is.null(colnames(design))) colnames(design) <- paste0("coefficient", seq_len(ncol(design)) - 1L)
  .rsdeseq2_validate_results_names(colnames(design))
  if (test == "LRT") {
    if (!is.matrix(reduced) || !is.numeric(reduced) || nrow(reduced) != ncol(counts)) {
      stop("LRT requires a numeric reduced design matrix with one row per sample", call. = FALSE)
    }
    storage.mode(reduced) <- "double"
    if (any(!is.finite(reduced))) stop("reduced design must contain finite values", call. = FALSE)
  } else {
    reduced <- NULL
  }
  results <- .Call(
    "rsdeseq2_fit",
    counts,
    design,
    reduced,
    sizeFactors,
    normalizationFactors,
    weights,
    as.integer(match(fitType, c("parametric", "local", "mean")) - 1L),
    as.integer(match(test, c("Wald", "LRT")) - 1L),
    if (disableReplacement) 0L else as.integer(minReplicatesForReplace),
    .rsdeseq2_default_cooks_cutoff(design) %||% NA_real_,
    as.integer(nproc),
    as.integer(match(parityMode, c("strict", "fast", "diagnostic")) - 1L),
    PACKAGE = "rsdeseq2"
  )
  rawPvalue <- results$results$rawPvalue
  rawStat <- results$results$rawStat
  results$results$rawPvalue <- NULL
  results$results$rawStat <- NULL
  rownames(results$beta) <- rownames(counts)
  rownames(results$betaSE) <- rownames(counts)
  names(results$dispersions) <- rownames(counts)
  names(results$dispFit) <- rownames(counts)
  dimnames(results$mu) <- dimnames(counts)
  dimnames(results$cooks) <- dimnames(counts)
  dimnames(results$betaCovariance)[[1L]] <- rownames(counts)
  fit <- structure(
    list(
      schemaVersion = 2L,
      designOrigin = designOrigin,
      rawPvalue = rawPvalue,
      rawStat = rawStat,
      results = results$results,
      baseMean = results$results$baseMean,
      beta = results$beta,
      betaSE = results$betaSE,
      betaCovariance = results$betaCovariance,
      mu = results$mu,
      cooks = results$cooks,
      replaceCounts = if (results$nRefit > 0L) results$replaceCounts else NULL,
      replace = results$replace,
      replaceable = results$replaceable,
      refitRows = results$refitRows,
      newAllZeroRows = results$newAllZeroRows,
      nRefit = results$nRefit,
      baseVar = results$baseVar,
      allZero = results$allZero,
      minReplicatesForReplace = minReplicatesForReplace,
      sizeFactors = results$sizeFactors,
      dispersions = results$dispersions,
      dispFit = results$dispFit,
      normalizationFactors = normalizationFactors,
      weights = weights,
      counts = counts,
      design = design,
      reduced = reduced,
      colData = colData,
      isDataSet = isDataSet,
      fitType = fitType,
      test = test,
      resultsNames = colnames(design),
      lrtStat = if (test == "LRT") rawStat else NULL,
      lrtPvalue = if (test == "LRT") rawPvalue else NULL,
      lrtDf = if (test == "LRT") ncol(design) - ncol(reduced) else NULL
    ),
    class = "rsdeseq2RustFit"
  )
  if (isDataSet) {
    orderedFactors <- names(modelFrame)[vapply(modelFrame, is.ordered, logical(1))]
    if (length(orderedFactors)) fit$orderedFactors <- orderedFactors
    rowMetadata <- S4Vectors::mcols(object)
    previousFit <- S4Vectors::metadata(object)$rsdeseq2RustFit
    if (!is.null(previousFit$resultsNames)) {
      previousColumns <- c(paste0("MLE_", previousFit$resultsNames), paste0("SE_", previousFit$resultsNames))
      for (column in intersect(previousColumns, names(rowMetadata))) rowMetadata[[column]] <- NULL
    }
    rowMetadata$baseMean <- fit$results$baseMean
    rowMetadata$baseVar <- fit$baseVar
    rowMetadata$allZero <- fit$allZero
    rowMetadata$dispFit <- fit$dispFit
    rowMetadata$dispersion <- fit$dispersions
    rowMetadata$betaConv <- fit$results$converged
    rowMetadata$maxCooks <- fit$results$maxCooks
    if (any(fit$replaceable)) {
      rowMetadata$replace <- fit$replace
      colMetadata <- SummarizedExperiment::colData(object)
      colMetadata$replaceable <- fit$replaceable
      SummarizedExperiment::colData(object) <- colMetadata
    } else {
      rowMetadata$replace <- NULL
      colMetadata <- SummarizedExperiment::colData(object)
      colMetadata$replaceable <- NULL
      SummarizedExperiment::colData(object) <- colMetadata
    }
    for (coefficient in seq_along(fit$resultsNames)) {
      coefficient_name <- fit$resultsNames[[coefficient]]
      rowMetadata[[paste0("MLE_", coefficient_name)]] <- fit$beta[, coefficient]
      rowMetadata[[paste0("SE_", coefficient_name)]] <- fit$betaSE[, coefficient]
    }
    S4Vectors::mcols(object) <- rowMetadata
    SummarizedExperiment::assay(object, "mu") <- fit$mu
    SummarizedExperiment::assay(object, "cooks") <- fit$cooks
    if (!is.null(fit$weights)) {
      SummarizedExperiment::assay(object, "weights", withDimnames = FALSE) <- fit$weights
    } else if ("weights" %in% SummarizedExperiment::assayNames(object)) {
      assays <- SummarizedExperiment::assays(object)
      assays[["weights"]] <- NULL
      SummarizedExperiment::assays(object) <- assays
    }
    if (fit$nRefit > 0L) {
      replaceCounts <- fit$replaceCounts
      storage.mode(replaceCounts) <- storage.mode(counts)
      SummarizedExperiment::assay(object, "replaceCounts") <- replaceCounts
      SummarizedExperiment::assay(object, "replaceCooks") <- fit$cooks
    }
    if (fit$nRefit == 0L) {
      assays <- SummarizedExperiment::assays(object)
      for (name in intersect(c("replaceCounts", "replaceCooks"), names(assays))) assays[[name]] <- NULL
      SummarizedExperiment::assays(object) <- assays
    }
    objectMetadata <- S4Vectors::metadata(object)
    normalizationState$writtenSizeFactors <- fit$sizeFactors
    if (is.null(normalizationState$baseFactors) && is.null(normalizationState$baseSizeFactors)) {
      normalizationState$baseSizeFactors <- fit$sizeFactors
    }
    objectMetadata$rsdeseq2Normalization <- normalizationState
    objectMetadata$rsdeseq2RustFit <- fit
    S4Vectors::metadata(object) <- objectMetadata
    if (requireNamespace("DESeq2", quietly = TRUE)) {
      DESeq2::sizeFactors(object) <- fit$sizeFactors
      if (!is.null(fit$normalizationFactors)) {
        DESeq2::normalizationFactors(object) <- fit$normalizationFactors
      } else if ("normalizationFactors" %in% SummarizedExperiment::assayNames(object)) {
        assays <- SummarizedExperiment::assays(object)
        assays[["normalizationFactors"]] <- NULL
        SummarizedExperiment::assays(object) <- assays
      }
    }
    return(object)
  }
  fit
}

fitRust <- DESeq
