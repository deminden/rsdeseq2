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
                    cooksCutoff = TRUE,
                    ...) {
  .rsdeseq2_check_unused_arguments(list(...), "results")
  if (!identical(cooksCutoff, TRUE)) {
    cutoff <- .rsdeseq2_validate_cooks_cutoff(cooksCutoff)
    if (!is.null(cutoff) && cutoff < 0) stop("cooksCutoff must be non-negative", call. = FALSE)
  }
  if (!inherits(object, c("DESeqDataSet", "rsdeseq2RustFit")) &&
      !missing(cooksCutoff) && !identical(cooksCutoff, FALSE) && !is.null(cooksCutoff)) {
    stop("Cook's filtering requires a fitted dataset with retained Cook's context", call. = FALSE)
  }
  altHypothesis <- match.arg(
    altHypothesis,
    c("greaterAbs", "lessAbs", "greater", "less", "greaterAbs2014", "greaterAbsUPSHOT")
  )
  if (inherits(object, "DESeqDataSet")) {
    fit <- .rsdeseq2_fitted_view(object)
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
      cooksCutoff = cooksCutoff,
      ...
    ))
  }
  if (inherits(object, "rsdeseq2RustFit")) {
    args <- list(contrast=contrast,name=name,listValues=listValues,reference=reference,
      lfcThreshold=lfcThreshold,altHypothesis=altHypothesis,useT=useT,degreesOfFreedom=degreesOfFreedom,
      alpha=alpha,independentFiltering=independentFiltering,cooksCutoff=cooksCutoff)
    if(.rsdeseq2_cached_result_request(args)) {
      return(.rsdeseq2_finish_fitted_result(object,args,object$results,cached=TRUE))
    }
    lrt <- identical(object$test,"LRT")
    primitive <- if(lrt)nbinomLRTRust(object=object) else nbinomWaldTestRust(object=object)
    prepared <- .rsdeseq2_prepare_result_request(primitive,args,lrt)
    raw <- if(lrt) .rsdeseq2_lrt_from_primitive_fit(primitive,prepared$selector$numeric,prepared$selector$allZero) else
      .rsdeseq2_wald_from_primitive_fit(primitive,prepared$selector$numeric,prepared$wald,prepared$selector$allZero)
    out <- .rsdeseq2_format_primitive_result(primitive,prepared,raw,adjust=FALSE)
    return(.rsdeseq2_finish_fitted_result(object,args,out))
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

# Mask presentation p-values before the single final filtering/BH step.
.rsdeseq2_fit_cooks_policy <- function(out, fit, cutoff, cached = FALSE) {
  if (cached) {
    if (!is.null(fit$rawPvalue)) {
      out$pvalue <- fit$rawPvalue
    } else if (!identical(cutoff, TRUE)) {
      stop("This saved fit lacks raw p-values; refit before changing cooksCutoff", call. = FALSE)
    }
  } else if (identical(fit$test, "LRT") && is.null(fit$schemaVersion)) {
    stop("This saved LRT fit lacks verified raw statistics; refit before selecting results", call. = FALSE)
  }
  zero <- which(fit$newAllZeroRows %||% logical())
  if (length(zero)) {
    out$log2FoldChange[zero] <- 0
    out$lfcSE[zero] <- 0
    out$stat[zero] <- 0
    out$pvalue[zero] <- 1
  }
  if (identical(cutoff, TRUE)) {
    cutoff <- .rsdeseq2_default_cooks_cutoff(fit$design)
  }
  if (is.null(cutoff) || identical(cutoff, FALSE)) {
    if ("cooksOutlier" %in% names(out)) out$cooksOutlier <- rep(NA, nrow(out))
    return(out)
  }
  maxCooks <- fit$results$maxCooks
  if (is.null(maxCooks)) stop("The fit lacks Cook's context; refit", call. = FALSE)
  # Successfully refitted rows must not be masked using the original outlier.
  if (!is.null(fit$refitRows)) maxCooks[which(fit$refitRows)] <- 0
  if (length(zero)) maxCooks[zero] <- 0
  lowCount <- FALSE
  if (ncol(fit$design) == 2L && !is.null(fit$colData)) {
    groups <- fit$design[, 2L]
    lowCount <- any(vapply(fit$colData, function(x) {
      is.factor(x) && nlevels(x) == 2L &&
        identical(outer(unname(as.character(x)), unname(as.character(x)), "=="), outer(unname(groups), unname(groups), "=="))
    }, logical(1)))
  }
  masked <- applyCooksCutoffRust(out$pvalue, maxCooks, cutoff,
                                counts = fit$counts, cooks = fit$cooks,
                                lowCountHeuristic = lowCount)
  out$pvalue <- masked$pvalue
  if ("cooksOutlier" %in% names(out)) out$cooksOutlier <- masked$cooksOutlier
  out
}
