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
                          factorReferences = NULL,
                          modelMatrix = NULL,
                          colData = NULL,
                          rowNames = NULL,
                          object = NULL,
                          ...) {
  .rsdeseq2_check_unused_arguments(list(...), "nbinomLRTRust")
  supplied <- names(match.call(expand.dots = FALSE))
  if (!is.null(object)) {
    fieldsSupplied <- setdiff(supplied, c("", "object", "..."))
    if (inherits(object, "rsdeseq2PrimitiveLrtFit") && !length(fieldsSupplied)) {
      return(object)
    }
    fields <- .rsdeseq2_extract_lrt_fit_fields(object)
    if (!"baseMean" %in% supplied) baseMean <- fields$baseMean
    if (!"beta" %in% supplied) beta <- fields$beta
    if (!"betaCovariance" %in% supplied) betaCovariance <- fields$betaCovariance
    if (!"betaSE" %in% supplied) betaSE <- fields$betaSE
    if (!"lrtStat" %in% supplied) lrtStat <- fields$lrtStat
    if (!"lrtPvalue" %in% supplied) lrtPvalue <- fields$lrtPvalue
    if (!"lrtDf" %in% supplied) lrtDf <- fields$lrtDf
    if (!"resultsNames" %in% supplied) resultsNames <- fields$resultsNames
    if (!"counts" %in% supplied) counts <- fields$counts
    if (!"sampleLevels" %in% supplied) sampleLevels <- fields$sampleLevels
    if (!"factorLevels" %in% supplied) factorLevels <- fields$factorLevels
    if (!"factorReferences" %in% supplied) factorReferences <- fields$factorReferences
    if (!"modelMatrix" %in% supplied) modelMatrix <- fields$modelMatrix
    if (!"colData" %in% supplied) colData <- fields$colData
    if (!"rowNames" %in% supplied) rowNames <- fields$rowNames
  }
  if (any(c("lrtStat", "lrtDf") %in% supplied) && !"lrtPvalue" %in% supplied) lrtPvalue <- NULL
  if (is.null(baseMean) || is.null(beta) || is.null(lrtStat)) {
    stop("baseMean, beta, and lrtStat are required for primitive nbinomLRTRust", call. = FALSE)
  }
  fit <- lrtFitRust(
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
    factorReferences = factorReferences,
    modelMatrix = modelMatrix,
    colData = colData,
    rowNames = rowNames
  )
  orderedFactors <- .rsdeseq2_object_field(object, "orderedFactors")
  if (!is.null(orderedFactors)) fit$orderedFactors <- .rsdeseq2_validate_results_names(orderedFactors)
  fit
}
