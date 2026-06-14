nbinomWaldTestRust <- function(baseMean = NULL,
                               beta = NULL,
                               betaCovariance = NULL,
                               betaSE = NULL,
                               resultsNames = NULL,
                               counts = NULL,
                               sampleLevels = NULL,
                               factorLevels = NULL,
                               factorReferences = NULL,
                               modelMatrix = NULL,
                               colData = NULL,
                               rowNames = NULL,
                               object = NULL,
                               ...) {
  if (!is.null(object)) {
    if (inherits(object, "rsdeseq2PrimitiveWaldFit")) {
      return(object)
    }
    fields <- .rsdeseq2_extract_wald_fit_fields(object)
    baseMean <- baseMean %||% fields$baseMean
    beta <- beta %||% fields$beta
    betaCovariance <- betaCovariance %||% fields$betaCovariance
    betaSE <- betaSE %||% fields$betaSE
    resultsNames <- resultsNames %||% fields$resultsNames
    counts <- counts %||% fields$counts
    sampleLevels <- sampleLevels %||% fields$sampleLevels
    factorLevels <- factorLevels %||% fields$factorLevels
    factorReferences <- factorReferences %||% fields$factorReferences
    modelMatrix <- modelMatrix %||% fields$modelMatrix
    colData <- colData %||% fields$colData
    rowNames <- rowNames %||% fields$rowNames
  }
  if (is.null(baseMean) || is.null(beta)) {
    stop("baseMean and beta are required for primitive nbinomWaldTestRust", call. = FALSE)
  }
  waldFitRust(
    baseMean = baseMean,
    beta = beta,
    betaCovariance = betaCovariance,
    betaSE = betaSE,
    resultsNames = resultsNames,
    counts = counts,
    sampleLevels = sampleLevels,
    factorLevels = factorLevels,
    factorReferences = factorReferences,
    modelMatrix = modelMatrix,
    colData = colData,
    rowNames = rowNames
  )
}
