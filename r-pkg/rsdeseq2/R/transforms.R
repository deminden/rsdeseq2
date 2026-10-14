.rsdeseq2_transform_dataset <- function(object,
                                         transformType,
                                         blind,
                                         fitType,
                                         nsub = 0L,
                                         nproc = 1L,
                                         intercept = NULL,
                                         betaPriorVar = NULL) {
  if (!inherits(object, "DESeqDataSet")) {
    stop("object must be a DESeqDataSet", call. = FALSE)
  }
  for (package in c("SummarizedExperiment", "S4Vectors", "DESeq2")) {
    if (!requireNamespace(package, quietly = TRUE)) {
      stop(sprintf("%s is required for DESeqDataSet transformations", package), call. = FALSE)
    }
  }
  if (!is.logical(blind) || length(blind) != 1L || is.na(blind)) {
    stop("blind must be TRUE or FALSE", call. = FALSE)
  }
  blind <- isTRUE(blind)
  counts <- SummarizedExperiment::assay(object, "counts")
  counts <- .rsdeseq2_validate_counts(counts)
  if (length(nproc) != 1L || is.na(nproc) || !is.finite(nproc) || nproc < 1L || nproc != floor(nproc)) {
    stop("nproc must be a positive integer", call. = FALSE)
  }
  fitType <- match.arg(fitType, c("parametric", "local", "mean"))
  requestedFitType <- fitType
  frozenTrendType <- -1L
  frozenTrendParameters <- NULL
  frozenTrendSource <- "estimated"
  frozenDispersionFunction <- NULL
  frozenDispFit <- NULL
  frozenLocalFunction <- NULL
  frozenLocalCurve <- NULL
  if (identical(transformType, "vstFull") && !blind) {
    candidate <- DESeq2::dispersionFunction(object)
    candidateType <- attr(candidate, "fitType", exact = TRUE)
    if (identical(candidateType, "parametric")) {
      coefficients <- attr(candidate, "coefficients", exact = TRUE)
      if (is.numeric(coefficients) && all(c("extraPois", "asymptDisp") %in% names(coefficients))) {
        frozenTrendType <- 0L
        frozenTrendParameters <- unname(coefficients[c("extraPois", "asymptDisp")])
        frozenDispersionFunction <- candidateType
        frozenTrendSource <- "frozen_dispersion_function"
        fitType <- candidateType
      }
    } else if (identical(candidateType, "mean")) {
      meanDispersion <- attr(candidate, "mean", exact = TRUE)
      if (is.numeric(meanDispersion) && length(meanDispersion) == 1L && is.finite(meanDispersion) && meanDispersion > 0) {
        frozenTrendType <- 2L
        frozenTrendParameters <- c(meanDispersion, NA_real_)
        frozenDispersionFunction <- candidateType
        frozenTrendSource <- "frozen_dispersion_function"
        fitType <- candidateType
      }
    } else if (identical(candidateType, "local")) {
      frozenTrendType <- 1L
      frozenLocalFunction <- candidate
      frozenDispersionFunction <- candidateType
      frozenTrendSource <- "frozen_dispersion_function"
      fitType <- candidateType
    }
  }
  if (identical(transformType, "rlog") && !blind) {
    candidateDispFit <- S4Vectors::mcols(object)$dispFit
    if (!is.null(candidateDispFit)) {
      if (!is.numeric(candidateDispFit) || length(candidateDispFit) != nrow(counts)) {
        stop("stored dispFit must be numeric with one value per gene", call. = FALSE)
      }
      if (!is.null(names(candidateDispFit)) && !is.null(rownames(counts))) {
        candidateDispFit <- .rsdeseq2_align_named_gene_vector(
          candidateDispFit, rownames(counts), "dispFit"
        )
      }
      if (any(!is.finite(candidateDispFit)) || any(candidateDispFit < 0)) {
        stop("stored dispFit values must be finite and non-negative", call. = FALSE)
      }
      frozenDispFit <- as.numeric(candidateDispFit)
      frozenTrendSource <- "frozen_dispFit"
      frozenDispersionFunction <- attr(
        DESeq2::dispersionFunction(object), "fitType", exact = TRUE
      )
      if (frozenDispersionFunction %in% c("parametric", "local", "mean")) {
        fitType <- frozenDispersionFunction
      }
    }
  }
  colData <- as.data.frame(SummarizedExperiment::colData(object))
  formula <- if (blind) stats::formula(~ 1) else methods::slot(object, "design")
  if (!inherits(formula, "formula")) {
    stop("DESeqDataSet design must be a formula", call. = FALSE)
  }
  modelFrame <- stats::model.frame(formula, data = colData)
  design <- stats::model.matrix(formula, data = modelFrame)
  formulaOffsets <- stats::model.offset(modelFrame)

  sizeFactors <- DESeq2::sizeFactors(object)
  normalizationFactors <- DESeq2::normalizationFactors(object)
  weights <- if ("weights" %in% SummarizedExperiment::assayNames(object)) {
    SummarizedExperiment::assay(object, "weights")
  } else {
    NULL
  }
  normalizationFactors <- .rsdeseq2_validate_normalization_factors(normalizationFactors, counts)
  weights <- .rsdeseq2_validate_observation_weights(weights, counts)
  if (!is.null(normalizationFactors)) {
    normalizationFactors <- .rsdeseq2_align_fit_observations(
      normalizationFactors, counts, "normalizationFactors"
    )
  }
  if (!is.null(weights)) {
    weights <- .rsdeseq2_align_fit_observations(weights, counts, "weights")
  }
  if (!blind && !is.null(formulaOffsets)) {
    if (any(!is.finite(formulaOffsets))) {
      stop("formula offsets must be finite", call. = FALSE)
    }
    offsetFactors <- matrix(
      rep(exp(formulaOffsets), each = nrow(counts)),
      nrow = nrow(counts),
      ncol = ncol(counts),
      dimnames = dimnames(counts)
    )
    normalizationFactors <- if (is.null(normalizationFactors)) {
      offsetFactors
    } else {
      normalizationFactors * offsetFactors
    }
    normalizationFactors <- .rsdeseq2_validate_normalization_factors(normalizationFactors, counts)
  }
  if (is.null(normalizationFactors) && !is.null(sizeFactors)) {
    if (!is.null(names(sizeFactors)) && !is.null(colnames(counts))) {
      if (anyDuplicated(names(sizeFactors)) || anyDuplicated(colnames(counts))) {
        stop("named sizeFactors and count columns must have unique names", call. = FALSE)
      }
      sampleOrder <- match(colnames(counts), names(sizeFactors))
      if (anyNA(sampleOrder)) stop("sizeFactors names must match count columns", call. = FALSE)
      sizeFactors <- sizeFactors[sampleOrder]
    }
    sizeFactors <- .rsdeseq2_validate_size_factors(sizeFactors, ncol(counts))
  } else if (!is.null(normalizationFactors)) {
    sizeFactors <- NULL
  }
  if (!is.null(frozenLocalFunction)) {
    if (is.null(normalizationFactors) && is.null(sizeFactors)) {
      sizeFactors <- estimateSizeFactorsRust(counts, native = TRUE)
    }
    normalizedForTrend <- if (is.null(normalizationFactors)) {
      sweep(counts, 2L, sizeFactors, "/")
    } else {
      counts / normalizationFactors
    }
    maxNormalizedCount <- max(normalizedForTrend)
    frozenGrid <- sinh(seq(
      asinh(0), asinh(maxNormalizedCount), length.out = 1000L
    ))[-1L]
    frozenDispersions <- as.numeric(frozenLocalFunction(frozenGrid))
    if (length(frozenDispersions) != length(frozenGrid) ||
        any(!is.finite(frozenDispersions)) || any(frozenDispersions <= 0)) {
      stop("frozen local dispersionFunction must return finite positive values on the VST grid", call. = FALSE)
    }
    frozenLocalCurve <- c(frozenGrid, frozenDispersions)
  }
  if (!is.matrix(design) || nrow(design) != ncol(counts) || !is.numeric(design)) {
    stop("transformation design must have one numeric row per sample", call. = FALSE)
  }
  storage.mode(design) <- "double"
  if (any(!is.finite(design))) stop("transformation design must be finite", call. = FALSE)
  if (is.null(colnames(design))) colnames(design) <- paste0("coefficient", seq_len(ncol(design)) - 1L)

  if (transformType == "vst" && (length(nsub) != 1L || is.na(nsub) || !is.finite(nsub) || nsub < 1L || nsub != floor(nsub))) {
    stop("nsub must be a positive integer", call. = FALSE)
  }
  if (transformType == "rlog") {
    if (xor(is.null(intercept), is.null(betaPriorVar))) {
      stop("intercept and betaPriorVar must be supplied together", call. = FALSE)
    }
    if (!is.null(intercept) && (!is.numeric(intercept) || length(intercept) != nrow(counts) || any(!is.finite(intercept)))) {
      stop("intercept must be a finite numeric vector with one value per gene", call. = FALSE)
    }
    if (!is.null(intercept) && !is.null(names(intercept)) && !is.null(rownames(counts))) {
      intercept <- .rsdeseq2_align_named_gene_vector(intercept, rownames(counts), "intercept")
    }
    if (!is.null(betaPriorVar) && (!is.numeric(betaPriorVar) || length(betaPriorVar) != 1L || !is.finite(betaPriorVar) || betaPriorVar <= 0)) {
      stop("betaPriorVar must be a positive finite number", call. = FALSE)
    }
  }

  typeCode <- switch(transformType, vst = 0L, rlog = 1L, vstFull = 2L)
  runTransform <- function(selectedFitType) {
    backendFitType <- if (frozenTrendType >= 0L || !is.null(frozenDispFit)) "mean" else selectedFitType
    .Call(
      "rsdeseq2_transform",
      counts,
      design,
      sizeFactors,
      normalizationFactors,
      weights,
      as.integer(match(backendFitType, c("parametric", "local", "mean")) - 1L),
      as.integer(typeCode),
      blind,
      as.integer(nsub),
      as.integer(nproc),
      intercept,
      betaPriorVar,
      as.integer(frozenTrendType),
      frozenTrendParameters,
      frozenDispFit,
      frozenLocalCurve,
      PACKAGE = "rsdeseq2"
    )
  }
  output <- tryCatch(
    runTransform(fitType),
    error = function(error) {
      parametricFailure <- grepl(
        "dispersion fit did not converge|parametric dispersion fit failed",
        conditionMessage(error)
      )
      if (fitType != "parametric" || !parametricFailure) {
        stop(error)
      }
      fitType <<- "local"
      runTransform(fitType)
    }
  )
  if (is.null(normalizationFactors)) {
    DESeq2::sizeFactors(object) <- output$sizeFactors
  } else {
    DESeq2::normalizationFactors(object) <- normalizationFactors
  }
  SummarizedExperiment::assays(object) <- S4Vectors::SimpleList(counts = output$transformed)
  metadata <- S4Vectors::metadata(object)
  metadata$rsdeseq2Transform <- list(
    type = transformType,
    blind = blind,
    fitType = fitType,
    fitTypeRequested = requestedFitType,
    trendSource = frozenTrendSource,
    frozenFitType = frozenDispersionFunction,
    nsub = if (transformType == "vst") as.integer(nsub) else NULL,
    sizeFactors = output$sizeFactors,
    normalizationFactors = normalizationFactors,
    formulaOffsets = if (is.null(formulaOffsets) || blind) NULL else as.numeric(formulaOffsets),
    weightsUsed = !is.null(weights),
    samplePriorVariance = output$samplePriorVariance
  )
  S4Vectors::metadata(object) <- metadata
  transformed <- DESeq2::DESeqTransform(object)
  if (transformType == "rlog") {
    S4Vectors::mcols(transformed)$rlogIntercept <- output$intercept
    attr(transformed, "betaPriorVar") <- output$samplePriorVariance
  }
  transformed
}

vst <- function(object, blind = TRUE, nsub = 1000L,
                fitType = c("parametric", "local", "mean"), nproc = 1L) {
  .rsdeseq2_transform_dataset(object, "vst", blind, fitType, nsub, nproc)
}

varianceStabilizingTransformation <- function(object, blind = TRUE,
                                               fitType = c("parametric", "local", "mean"),
                                               nproc = 1L) {
  .rsdeseq2_transform_dataset(object, "vstFull", blind, fitType, nproc = nproc)
}

rlog <- function(object, blind = TRUE, intercept, betaPriorVar,
                 fitType = c("parametric", "local", "mean"), nproc = 1L) {
  interceptMissing <- missing(intercept)
  priorMissing <- missing(betaPriorVar)
  if (xor(interceptMissing, priorMissing)) {
    stop("intercept and betaPriorVar must be supplied together", call. = FALSE)
  }
  .rsdeseq2_transform_dataset(
    object,
    "rlog",
    blind,
    fitType,
    nproc = nproc,
    intercept = if (interceptMissing) NULL else intercept,
    betaPriorVar = if (priorMissing) NULL else betaPriorVar
  )
}
