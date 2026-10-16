# Result dispatch and metadata assembly.

results <- function(object,
                    contrast = NULL,
                    name = NULL,
                    listValues = c(1, -1),
                    reference = NULL,
                    lfcThreshold = 0,
                    altHypothesis = c("greaterAbs", "lessAbs", "greater", "less", "greaterAbs2014", "greaterAbsUPSHOT"),
                    useT = FALSE,
                    degreesOfFreedom = NULL,
                    alpha = 0.1,
                    independentFiltering = TRUE,
                    ...) {
  altHypothesis <- match.arg(
    altHypothesis,
    c("greaterAbs", "lessAbs", "greater", "less", "greaterAbs2014", "greaterAbsUPSHOT")
  )
  if (inherits(object, "DESeqDataSet")) {
    fit <- S4Vectors::metadata(object)$rsdeseq2RustFit
    if (is.null(fit)) {
      stop("DESeqDataSet must be fitted with rsdeseq2::DESeq() before calling results()", call. = FALSE)
    }
    return(results(
      object = fit,
      contrast = contrast,
      name = name,
      listValues = listValues,
      reference = reference,
      lfcThreshold = lfcThreshold,
      altHypothesis = altHypothesis,
      useT = useT,
      degreesOfFreedom = degreesOfFreedom,
      alpha = alpha,
      independentFiltering = independentFiltering,
      ...
    ))
  }
  if (inherits(object, "rsdeseq2RustFit")) {
    if (is.null(contrast) && is.null(name) &&
        identical(as.numeric(lfcThreshold), 0) && identical(altHypothesis, "greaterAbs") &&
        !isTRUE(useT) && is.null(degreesOfFreedom)) {
      out <- applyIndependentFilteringRust(
        object$results,
        alpha = alpha,
        enabled = independentFiltering
      )
      if (!isTRUE(independentFiltering)) out$filtered <- NULL
      if (identical(object$test, "Wald")) {
        resultName <- utils::tail(object$resultsNames, 1L)
        comparison <- sprintf("coefficient %s", resultName)
      } else {
        resultName <- utils::tail(object$resultsNames, 1L)
        comparison <- sprintf("coefficient %s; LRT full vs reduced", resultName)
      }
      metadata <- list(
        resultName = resultName,
        comparison = comparison,
        test = object$test,
        contrastAllZero = "none",
        lfcThreshold = 0,
        altHypothesis = "greaterAbs",
        contrast = NULL,
        alpha = alpha,
        lrtDf = object$lrtDf,
        independentFiltering = attr(out, "independentFiltering"),
        cooksReplacement = .rsdeseq2_cooks_replacement_metadata(object)
      )
      return(.rsdeseq2_finalize_result_metadata(out, metadata, object))
    }
    primitive <- if (identical(object$test, "LRT")) {
      nbinomLRTRust(object = object)
    } else {
      nbinomWaldTestRust(object = object)
    }
    out <- results(
      object = primitive,
      contrast = contrast,
      name = name,
      listValues = listValues,
      reference = reference,
      lfcThreshold = lfcThreshold,
      altHypothesis = altHypothesis,
      useT = useT,
      degreesOfFreedom = degreesOfFreedom,
      alpha = alpha,
      independentFiltering = independentFiltering
    )
    metadata <- attr(out, "metadata")
    if (is.null(metadata)) {
      metadata <- list(
        resultName = attr(out, "resultName"),
        comparison = attr(out, "comparison"),
        test = object$test,
        contrastAllZero = attr(out, "contrastAllZero"),
        contrast = attr(out, "contrast"),
        lfcThreshold = attr(out, "lfcThreshold"),
        altHypothesis = attr(out, "altHypothesis"),
        degreesOfFreedom = attr(out, "degreesOfFreedom"),
        lrtDf = attr(out, "lrtDf"),
        alpha = alpha,
        independentFiltering = attr(out, "independentFiltering")
      )
    }
    metadata$cooksReplacement <- .rsdeseq2_cooks_replacement_metadata(object)
    return(.rsdeseq2_finalize_result_metadata(out, metadata, object))
  }
  if (inherits(object, "rsdeseq2PrimitiveWaldFit")) {
    return(.rsdeseq2_results_primitive_wald_fit(
      object = object,
      contrast = contrast,
      name = name,
      listValues = listValues,
      reference = reference,
      lfcThreshold = lfcThreshold,
      altHypothesis = altHypothesis,
      useT = useT,
      degreesOfFreedom = degreesOfFreedom,
      alpha = alpha,
      independentFiltering = independentFiltering
    ))
  }
  if (inherits(object, "rsdeseq2PrimitiveLrtFit")) {
    return(.rsdeseq2_results_primitive_lrt_fit(
      object = object,
      contrast = contrast,
      name = name,
      listValues = listValues,
      reference = reference,
      lfcThreshold = lfcThreshold,
      altHypothesis = altHypothesis,
      useT = useT,
      degreesOfFreedom = degreesOfFreedom,
      alpha = alpha,
      independentFiltering = independentFiltering
    ))
  }
  if (.rsdeseq2_is_list_like_lrt_input(object)) {
    return(results(
      object = nbinomLRTRust(object = object),
      contrast = contrast,
      name = name,
      listValues = listValues,
      reference = reference,
      lfcThreshold = lfcThreshold,
      altHypothesis = altHypothesis,
      useT = useT,
      degreesOfFreedom = degreesOfFreedom,
      alpha = alpha,
      independentFiltering = independentFiltering
    ))
  }
  if (.rsdeseq2_is_list_like_wald_input(object)) {
    return(results(
      object = nbinomWaldTestRust(object = object),
      contrast = contrast,
      name = name,
      listValues = listValues,
      reference = reference,
      lfcThreshold = lfcThreshold,
      altHypothesis = altHypothesis,
      useT = useT,
      degreesOfFreedom = degreesOfFreedom,
      alpha = alpha,
      independentFiltering = independentFiltering
    ))
  }
  .rsdeseq2_not_available("results")
}

.rsdeseq2_cooks_replacement_metadata <- function(object) {
  eligible <- object$replaceable
  enabled <- !is.null(eligible) && any(eligible)
  refitRows <- object$refitRows
  newAllZeroRows <- object$newAllZeroRows
  refitRows <- if (enabled && !is.null(refitRows)) {
    indices <- which(refitRows)
    if (is.null(names(refitRows))) indices else names(refitRows)[indices]
  } else {
    character()
  }
  newAllZeroRows <- if (enabled && !is.null(newAllZeroRows)) {
    indices <- which(newAllZeroRows)
    if (is.null(names(newAllZeroRows))) indices else names(newAllZeroRows)[indices]
  } else {
    character()
  }
  replaceableSamples <- if (enabled) {
    indices <- which(eligible)
    if (is.null(names(eligible))) indices else names(eligible)[indices]
  } else {
    character()
  }
  list(
    enabled = enabled,
    minReplicatesForReplace = object$minReplicatesForReplace,
    nRefit = if (enabled) object$nRefit else 0L,
    refitRows = refitRows,
    newAllZeroRows = newAllZeroRows,
    replaceableSamples = replaceableSamples
  )
}

.rsdeseq2_finalize_result_metadata <- function(out, resultMetadata, fit) {
  attr(out, "resultName") <- resultMetadata$resultName
  attr(out, "comparison") <- resultMetadata$comparison
  attr(out, "test") <- resultMetadata$test
  attr(out, "contrastAllZero") <- resultMetadata$contrastAllZero
  attr(out, "metadata") <- resultMetadata
  if (!isTRUE(fit$isDataSet) || !requireNamespace("S4Vectors", quietly = TRUE)) {
    return(out)
  }

  biocMetadata <- resultMetadata
  filter <- resultMetadata$independentFiltering
  if (isTRUE(filter$enabled)) {
    biocMetadata$filterTheta <- filter$filterTheta
    biocMetadata$filterThreshold <- filter$filterThreshold
    biocMetadata$filterNumRej <- data.frame(theta = filter$theta, numRej = filter$numRej)
    biocMetadata$lo.fit <- filter$lo.fit
  }
  result <- S4Vectors::DataFrame(out)
  descriptions <- c(
    baseMean = "mean of normalized counts for all samples",
    log2FoldChange = sprintf("log2 fold change (MLE): %s", resultMetadata$resultName %||% "contrast"),
    lfcSE = "standard error: log2 fold change",
    stat = "Wald statistic" ,
    pvalue = if (identical(resultMetadata$test, "LRT")) "LRT p-value" else "Wald test p-value",
    padj = "BH adjusted p-values",
    dispersion = "gene-wise dispersion estimate",
    converged = "Wald/LRT fit convergence",
    maxCooks = "maximum Cook's distance",
    cooksOutlier = "Cook's distance outlier",
    filtered = "independent filtering status"
  )
  columnMetadata <- S4Vectors::DataFrame(
    type = rep("results", ncol(result)),
    description = unname(descriptions[colnames(result)]),
    row.names = colnames(result)
  )
  S4Vectors::mcols(result) <- columnMetadata
  S4Vectors::metadata(result) <- biocMetadata
  if (requireNamespace("DESeq2", quietly = TRUE) && methods::isClass("DESeqResults")) {
    priorInfo <- list(type = "none", package = "rsdeseq2", version = utils::packageVersion("rsdeseq2"))
    result <- methods::new("DESeqResults", result, priorInfo = priorInfo)
  }
  # Reapply after constructing DESeqResults, whose S4 constructor may normalize
  # metadata, so downstream DESeq2-compatible consumers see both our structured
  # metadata and the conventional independent-filtering fields.
  S4Vectors::metadata(result) <- biocMetadata
  for (attribute in c(
    "resultName", "comparison", "test", "contrastAllZero",
    "independentFiltering", "contrast", "lfcThreshold", "altHypothesis",
    "degreesOfFreedom", "lrtDf"
  )) {
    attr(result, attribute) <- attr(out, attribute)
  }
  result
}
