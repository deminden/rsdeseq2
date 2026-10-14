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
                  ...) {
  isDataSet <- inherits(object, "DESeqDataSet")
  fitType <- match.arg(fitType)
  test <- match.arg(test)
  if (!identical(engine, "rust")) {
    stop("engine must be 'rust'", call. = FALSE)
  }
  parityMode <- match.arg(parityMode, c("strict", "fast", "diagnostic"))
  dots <- list(...)
  if (length(dots)) {
    argumentNames <- names(dots)
    if (is.null(argumentNames)) argumentNames <- rep("...", length(dots))
    argumentNames[!nzchar(argumentNames)] <- "..."
    stop(sprintf("unsupported DESeq arguments: %s", paste(argumentNames, collapse = ", ")), call. = FALSE)
  }
  if (length(nproc) != 1L || is.na(nproc) || !is.finite(nproc) || nproc < 1L || nproc != floor(nproc)) {
    stop("nproc must be a positive integer", call. = FALSE)
  }
  colData <- NULL
  if (inherits(object, "DESeqDataSet")) {
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
      stop("SummarizedExperiment is required for DESeqDataSet input", call. = FALSE)
    }
    counts <- SummarizedExperiment::assay(object, "counts")
    colData <- as.data.frame(SummarizedExperiment::colData(object))
    if (requireNamespace("DESeq2", quietly = TRUE)) {
      if (is.null(sizeFactors)) sizeFactors <- DESeq2::sizeFactors(object)
      if (is.null(normalizationFactors)) normalizationFactors <- DESeq2::normalizationFactors(object)
    }
    if (is.null(weights) && "weights" %in% SummarizedExperiment::assayNames(object)) {
      weights <- SummarizedExperiment::assay(object, "weights")
    }
    formula <- methods::slot(object, "design")
    if (!inherits(formula, "formula")) {
      stop("DESeqDataSet design must be a formula", call. = FALSE)
    }
    modelFrame <- stats::model.frame(formula, data = colData)
    design <- stats::model.matrix(formula, data = modelFrame)
    formulaOffsets <- stats::model.offset(modelFrame)
    if (!is.null(formulaOffsets)) {
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
    }
    if (test == "LRT") {
      if (is.null(reduced) || !inherits(reduced, "formula")) {
        stop("LRT requires reduced = a formula", call. = FALSE)
      }
      reduced <- stats::model.matrix(reduced, data = colData)
    }
  } else if (is.list(object) && !is.null(object$counts)) {
    counts <- object$counts
    colData <- object$colData
    if (is.null(design)) design <- object$design
    if (is.null(reduced)) reduced <- object$reduced
    if (is.null(sizeFactors)) sizeFactors <- object$sizeFactors
    if (is.null(normalizationFactors)) normalizationFactors <- object$normalizationFactors
    if (is.null(weights)) weights <- object$weights
  } else {
    counts <- object
  }
  counts <- .rsdeseq2_validate_counts(counts)
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
    as.integer(nproc),
    as.integer(match(parityMode, c("strict", "fast", "diagnostic")) - 1L),
    PACKAGE = "rsdeseq2"
  )
  rownames(results$beta) <- rownames(counts)
  rownames(results$betaSE) <- rownames(counts)
  names(results$dispersions) <- rownames(counts)
  names(results$dispFit) <- rownames(counts)
  dimnames(results$betaCovariance)[[1L]] <- rownames(counts)
  fit <- structure(
    list(
      results = results$results,
      baseMean = results$results$baseMean,
      beta = results$beta,
      betaSE = results$betaSE,
      betaCovariance = results$betaCovariance,
      sizeFactors = results$sizeFactors,
      dispersions = results$dispersions,
      dispFit = results$dispFit,
      normalizationFactors = normalizationFactors,
      weights = weights,
      counts = counts,
      design = design,
      reduced = reduced,
      colData = colData,
      fitType = fitType,
      test = test,
      resultsNames = colnames(design),
      lrtStat = if (test == "LRT") results$results$stat else NULL,
      lrtPvalue = if (test == "LRT") results$results$pvalue else NULL,
      lrtDf = if (test == "LRT") ncol(design) - ncol(reduced) else NULL
    ),
    class = "rsdeseq2RustFit"
  )
  if (isDataSet) {
    objectMetadata <- S4Vectors::metadata(object)
    objectMetadata$rsdeseq2RustFit <- fit
    S4Vectors::metadata(object) <- objectMetadata
    if (requireNamespace("DESeq2", quietly = TRUE)) {
      DESeq2::sizeFactors(object) <- fit$sizeFactors
      if (!is.null(fit$normalizationFactors)) {
        DESeq2::normalizationFactors(object) <- fit$normalizationFactors
      }
    }
    return(object)
  }
  fit
}

fitRust <- DESeq
