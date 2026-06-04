nbinomLRTRust <- function(baseMean = NULL,
                          beta = NULL,
                          betaCovariance = NULL,
                          betaSE = NULL,
                          lrtStat = NULL,
                          lrtPvalue = NULL,
                          lrtDf = NULL,
                          resultsNames = NULL,
                          counts = NULL,
                          sampleLevels = NULL,
                          factorLevels = NULL,
                          modelMatrix = NULL,
                          rowNames = NULL,
                          object = NULL,
                          ...) {
  if (!is.null(object)) {
    if (inherits(object, "rsdeseq2PrimitiveLrtFit")) {
      return(object)
    }
    fields <- .rsdeseq2_extract_lrt_fit_fields(object)
    baseMean <- baseMean %||% fields$baseMean
    beta <- beta %||% fields$beta
    betaCovariance <- betaCovariance %||% fields$betaCovariance
    betaSE <- betaSE %||% fields$betaSE
    lrtStat <- lrtStat %||% fields$lrtStat
    lrtPvalue <- lrtPvalue %||% fields$lrtPvalue
    lrtDf <- lrtDf %||% fields$lrtDf
    resultsNames <- resultsNames %||% fields$resultsNames
    counts <- counts %||% fields$counts
    sampleLevels <- sampleLevels %||% fields$sampleLevels
    factorLevels <- factorLevels %||% fields$factorLevels
    modelMatrix <- modelMatrix %||% fields$modelMatrix
    rowNames <- rowNames %||% fields$rowNames
  }
  if (is.null(baseMean) || is.null(beta) || is.null(lrtStat)) {
    stop("baseMean, beta, and lrtStat are required for primitive nbinomLRTRust", call. = FALSE)
  }
  lrtFitRust(
    baseMean = baseMean,
    beta = beta,
    betaCovariance = betaCovariance,
    betaSE = betaSE,
    lrtStat = lrtStat,
    lrtPvalue = lrtPvalue,
    lrtDf = lrtDf,
    resultsNames = resultsNames,
    counts = counts,
    sampleLevels = sampleLevels,
    factorLevels = factorLevels,
    modelMatrix = modelMatrix,
    rowNames = rowNames
  )
}
