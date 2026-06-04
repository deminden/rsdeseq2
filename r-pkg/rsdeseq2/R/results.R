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

waldFitRust <- function(baseMean,
                        beta,
                        betaCovariance = NULL,
                        betaSE = NULL,
                        resultsNames = NULL,
                        counts = NULL,
                        sampleLevels = NULL,
                        factorLevels = NULL,
                        modelMatrix = NULL,
                        rowNames = NULL) {
  beta <- .rsdeseq2_validate_beta_matrix(beta, resultsNames)
  resultsNames <- colnames(beta)
  nGenes <- nrow(beta)
  baseMean <- .rsdeseq2_validate_result_numeric_vector(baseMean, "baseMean", allowNA = FALSE)
  if (length(baseMean) != nGenes) {
    stop("baseMean must have one value per beta row", call. = FALSE)
  }
  rowNames <- .rsdeseq2_validate_row_names(rowNames %||% rownames(beta) %||% names(baseMean), nGenes)
  betaCovariance <- .rsdeseq2_validate_beta_covariance(betaCovariance, nGenes, length(resultsNames))
  betaSE <- .rsdeseq2_validate_beta_se_matrix(betaSE, nGenes, resultsNames)
  counts <- .rsdeseq2_validate_optional_contrast_counts(counts, nGenes)
  sampleLevels <- .rsdeseq2_validate_optional_sample_levels(sampleLevels, counts)
  factorLevels <- .rsdeseq2_validate_optional_factor_levels(factorLevels)
  modelMatrix <- .rsdeseq2_validate_optional_model_matrix(modelMatrix, counts, resultsNames)

  out <- list(
    baseMean = baseMean,
    beta = beta,
    betaCovariance = betaCovariance,
    betaSE = betaSE,
    resultsNames = resultsNames,
    counts = counts,
    sampleLevels = sampleLevels,
    factorLevels = factorLevels,
    modelMatrix = modelMatrix,
    rowNames = rowNames
  )
  class(out) <- c("rsdeseq2PrimitiveWaldFit", "list")
  out
}

lrtFitRust <- function(baseMean,
                       beta,
                       betaCovariance = NULL,
                       betaSE = NULL,
                       lrtStat,
                       lrtPvalue = NULL,
                       lrtDf = NULL,
                       resultsNames = NULL,
                       counts = NULL,
                       sampleLevels = NULL,
                       factorLevels = NULL,
                       modelMatrix = NULL,
                       rowNames = NULL) {
  beta <- .rsdeseq2_validate_beta_matrix(beta, resultsNames)
  resultsNames <- colnames(beta)
  nGenes <- nrow(beta)
  baseMean <- .rsdeseq2_validate_result_numeric_vector(baseMean, "baseMean", allowNA = FALSE)
  if (length(baseMean) != nGenes) {
    stop("baseMean must have one value per beta row", call. = FALSE)
  }
  rowNames <- .rsdeseq2_validate_row_names(rowNames %||% rownames(beta) %||% names(baseMean), nGenes)
  betaCovariance <- .rsdeseq2_validate_beta_covariance(betaCovariance, nGenes, length(resultsNames))
  betaSE <- .rsdeseq2_validate_beta_se_matrix(betaSE, nGenes, resultsNames)
  lrtStat <- .rsdeseq2_validate_lrt_stat(lrtStat, nGenes)
  lrtDf <- .rsdeseq2_validate_lrt_df(lrtDf, nGenes)
  if (is.null(lrtPvalue)) {
    if (is.null(lrtDf)) {
      stop("lrtDf is required when lrtPvalue is NULL", call. = FALSE)
    }
    lrtPvalue <- stats::pchisq(lrtStat, df = lrtDf, lower.tail = FALSE)
    lrtPvalue[is.na(lrtStat) | is.na(lrtDf)] <- NA_real_
  } else {
    lrtPvalue <- .rsdeseq2_validate_result_pvalues(lrtPvalue, nGenes)
  }
  counts <- .rsdeseq2_validate_optional_contrast_counts(counts, nGenes)
  sampleLevels <- .rsdeseq2_validate_optional_sample_levels(sampleLevels, counts)
  factorLevels <- .rsdeseq2_validate_optional_factor_levels(factorLevels)
  modelMatrix <- .rsdeseq2_validate_optional_model_matrix(modelMatrix, counts, resultsNames)

  out <- list(
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
  class(out) <- c("rsdeseq2PrimitiveLrtFit", "list")
  out
}

resolveResultsContrastRust <- function(contrast,
                                       resultsNames,
                                       listValues = c(1, -1),
                                       reference = NULL) {
  resultsNames <- .rsdeseq2_validate_results_names(resultsNames)
  parsed <- .rsdeseq2_parse_results_contrast(
    contrast = contrast,
    resultsNames = resultsNames,
    listValues = listValues,
    reference = reference
  )
  class(parsed) <- c("rsdeseq2ResultsContrast", "list")
  parsed
}

resultsNamesRust <- function(object) {
  if (inherits(object, c("rsdeseq2PrimitiveWaldFit", "rsdeseq2PrimitiveLrtFit"))) {
    return(.rsdeseq2_validate_results_names(object$resultsNames))
  }
  if (is.character(object)) {
    return(.rsdeseq2_validate_results_names(object))
  }
  if (is.matrix(object)) {
    beta <- .rsdeseq2_validate_beta_matrix(object, colnames(object))
    return(colnames(beta))
  }
  if (is.list(object)) {
    if (!is.null(object$resultsNames)) {
      return(.rsdeseq2_validate_results_names(object$resultsNames))
    }
    if (!is.null(object$beta)) {
      return(resultsNamesRust(object$beta))
    }
  }
  fields <- .rsdeseq2_extract_wald_fit_fields(object)
  if (!is.null(fields$resultsNames)) {
    return(.rsdeseq2_validate_results_names(fields$resultsNames))
  }
  .rsdeseq2_not_available("resultsNamesRust object integration")
}

resultsNames <- function(object, ...) {
  resultsNamesRust(object)
}

.rsdeseq2_is_list_like_lrt_input <- function(object) {
  !inherits(object, c("rsdeseq2PrimitiveWaldFit", "rsdeseq2PrimitiveLrtFit")) &&
    {
      fields <- .rsdeseq2_extract_lrt_fit_fields(object)
      !is.null(fields$baseMean) &&
        !is.null(fields$beta) &&
        !is.null(fields$lrtStat)
    }
}

.rsdeseq2_is_list_like_wald_input <- function(object) {
  !inherits(object, c("rsdeseq2PrimitiveWaldFit", "rsdeseq2PrimitiveLrtFit")) &&
    {
      fields <- .rsdeseq2_extract_wald_fit_fields(object)
      !is.null(fields$baseMean) &&
        !is.null(fields$beta)
    }
}

.rsdeseq2_extract_wald_fit_fields <- function(object) {
  mcols <- .rsdeseq2_object_field(object, "mcols")
  assays <- .rsdeseq2_object_field(object, "assays")
  colData <- .rsdeseq2_object_field(object, "colData")
  beta <- .rsdeseq2_object_field(object, "beta") %||%
    .rsdeseq2_object_field(object, "betaMatrix") %||%
    .rsdeseq2_mcols_beta_matrix(mcols)
  resultsNames <- .rsdeseq2_object_field(object, "resultsNames") %||% colnames(beta)
  list(
    baseMean = .rsdeseq2_object_field(object, "baseMean") %||%
      .rsdeseq2_mcols_column(mcols, "baseMean"),
    beta = beta,
    betaCovariance = .rsdeseq2_object_field(object, "betaCovariance"),
    betaSE = .rsdeseq2_object_field(object, "betaSE") %||%
      .rsdeseq2_mcols_prefixed_matrix(mcols, "SE_", resultsNames),
    resultsNames = resultsNames,
    counts = .rsdeseq2_object_field(object, "counts") %||%
      .rsdeseq2_object_field(assays, "counts"),
    sampleLevels = .rsdeseq2_object_field(object, "sampleLevels"),
    factorLevels = .rsdeseq2_object_field(object, "factorLevels") %||%
      .rsdeseq2_factor_levels_from_col_data(colData),
    modelMatrix = .rsdeseq2_object_field(object, "modelMatrix"),
    rowNames = .rsdeseq2_object_field(object, "rowNames") %||%
      rownames(beta) %||%
      rownames(mcols)
  )
}

.rsdeseq2_extract_lrt_fit_fields <- function(object) {
  fields <- .rsdeseq2_extract_wald_fit_fields(object)
  mcols <- .rsdeseq2_object_field(object, "mcols")
  fields$lrtStat <- .rsdeseq2_object_field(object, "lrtStat") %||%
    .rsdeseq2_object_field(object, "stat") %||%
    .rsdeseq2_mcols_column(mcols, "stat")
  fields$lrtPvalue <- .rsdeseq2_object_field(object, "lrtPvalue") %||%
    .rsdeseq2_object_field(object, "pvalue") %||%
    .rsdeseq2_mcols_column(mcols, "pvalue")
  fields$lrtDf <- .rsdeseq2_object_field(object, "lrtDf")
  fields
}

.rsdeseq2_object_field <- function(object, name) {
  if (is.null(object)) {
    return(NULL)
  }
  if (is.list(object) && !is.null(object[[name]])) {
    return(object[[name]])
  }
  accessed <- .rsdeseq2_optional_object_accessor(object, name)
  if (!is.null(accessed)) {
    return(accessed)
  }
  if (isS4(object)) {
    slots <- methods::slotNames(object)
    if (name %in% slots) {
      return(methods::slot(object, name))
    }
    aliases <- switch(
      name,
      mcols = c("elementMetadata", "rowData"),
      colData = c("colData"),
      rowNames = c("NAMES"),
      character()
    )
    for (alias in aliases) {
      if (alias %in% slots) {
        return(methods::slot(object, alias))
      }
    }
  }
  NULL
}

.rsdeseq2_optional_object_accessor <- function(object, name) {
  accessors <- switch(
    name,
    mcols = list(c("BiocGenerics", "mcols"), c("S4Vectors", "mcols")),
    assays = list(c("SummarizedExperiment", "assays")),
    colData = list(c("SummarizedExperiment", "colData")),
    counts = list(c("DESeq2", "counts")),
    sizeFactors = list(c("DESeq2", "sizeFactors")),
    rowNames = list(c("BiocGenerics", "rownames")),
    list()
  )
  for (accessor in accessors) {
    namespace <- accessor[[1L]]
    symbol <- accessor[[2L]]
    if (requireNamespace(namespace, quietly = TRUE) && exists(symbol, envir = asNamespace(namespace), mode = "function")) {
      value <- tryCatch(
        get(symbol, envir = asNamespace(namespace), mode = "function")(object),
        error = function(e) NULL
      )
      if (!is.null(value)) {
        return(value)
      }
    }
  }
  NULL
}

.rsdeseq2_mcols_column <- function(mcols, name) {
  if (is.null(mcols) || is.null(mcols[[name]])) {
    return(NULL)
  }
  mcols[[name]]
}

.rsdeseq2_factor_levels_from_col_data <- function(colData) {
  if (is.null(colData)) {
    return(NULL)
  }
  names <- colnames(colData) %||% names(colData) %||% character()
  out <- list()
  for (name in names) {
    column <- colData[[name]]
    if (is.factor(column)) {
      out[[name]] <- levels(column)
    }
  }
  if (length(out) == 0L) {
    return(NULL)
  }
  out
}

.rsdeseq2_mcols_beta_matrix <- function(mcols) {
  if (is.null(mcols)) {
    return(NULL)
  }
  betaMatrix <- .rsdeseq2_mcols_column(mcols, "betaMatrix")
  if (!is.null(betaMatrix)) {
    return(betaMatrix)
  }
  if (is.null(colnames(mcols))) {
    return(NULL)
  }
  seColumns <- grep("^SE_", colnames(mcols), value = TRUE)
  betaColumns <- sub("^SE_", "", seColumns)
  betaColumns <- betaColumns[betaColumns %in% colnames(mcols)]
  if (length(betaColumns) == 0L) {
    return(NULL)
  }
  as.matrix(mcols[, betaColumns, drop = FALSE])
}

.rsdeseq2_mcols_prefixed_matrix <- function(mcols, prefix, resultsNames) {
  if (is.null(mcols) || is.null(resultsNames)) {
    return(NULL)
  }
  columns <- paste0(prefix, resultsNames)
  if (is.null(colnames(mcols)) || !all(columns %in% colnames(mcols))) {
    return(NULL)
  }
  out <- as.matrix(mcols[, columns, drop = FALSE])
  colnames(out) <- resultsNames
  out
}

.rsdeseq2_results_primitive_wald_fit <- function(object,
                                                 contrast,
                                                 name,
                                                 listValues,
                                                 reference,
                                                 lfcThreshold,
                                                 altHypothesis,
                                                 useT,
                                                 degreesOfFreedom,
                                                 alpha,
                                                 independentFiltering) {
  resultsNames <- object$resultsNames
  reference <- reference %||% .rsdeseq2_reference_for_character_contrast(object, contrast)
  if (is.null(contrast)) {
    if (is.null(name)) {
      if (length(resultsNames) < 2L) {
        stop("default results require at least two result names", call. = FALSE)
      }
      name <- resultsNames[[length(resultsNames)]]
    }
    coefficient <- .rsdeseq2_resolve_results_name_index(name, resultsNames)
    numericContrast <- numeric(length(resultsNames))
    numericContrast[[coefficient]] <- 1
    names(numericContrast) <- resultsNames
    resultName <- resultsNames[[coefficient]]
    comparison <- sprintf("coefficient %s", resultName)
    allZero <- list(type = "none")
  } else {
    resolved <- resolveResultsContrastRust(
      contrast = contrast,
      resultsNames = resultsNames,
      listValues = listValues,
      reference = reference
    )
    numericContrast <- resolved$numeric
    resultName <- resolved$resultName
    comparison <- resolved$comparison
    allZero <- resolved$allZero
  }

  waldOptions <- .rsdeseq2_validate_wald_result_options(
    lfcThreshold = lfcThreshold,
    altHypothesis = altHypothesis,
    useT = useT,
    degreesOfFreedom = degreesOfFreedom,
    nGenes = nrow(object$beta)
  )
  wald <- .rsdeseq2_wald_from_primitive_fit(object, numericContrast, waldOptions)
  wald <- .rsdeseq2_apply_results_contrast_all_zero(
    wald = wald,
    fit = object,
    numericContrast = numericContrast,
    allZero = allZero
  )
  out <- resultsTableRust(
    baseMean = object$baseMean,
    log2FoldChange = wald$log2FoldChange,
    lfcSE = wald$lfcSE,
    stat = wald$stat,
    pvalue = wald$pvalue,
    rowNames = object$rowNames
  )
  attr(out, "resultName") <- resultName
  attr(out, "comparison") <- comparison
  attr(out, "contrast") <- numericContrast
  attr(out, "contrastAllZero") <- allZero$type
  attr(out, "lfcThreshold") <- waldOptions$lfcThreshold
  attr(out, "altHypothesis") <- waldOptions$altHypothesis
  if (!is.null(waldOptions$degreesOfFreedom)) {
    attr(out, "degreesOfFreedom") <- waldOptions$degreesOfFreedom
  }
  applyIndependentFilteringRust(
    out,
    alpha = alpha,
    enabled = independentFiltering
  )
}

.rsdeseq2_results_primitive_lrt_fit <- function(object,
                                                contrast,
                                                name,
                                                listValues,
                                                reference,
                                                lfcThreshold,
                                                altHypothesis,
                                                useT,
                                                degreesOfFreedom,
                                                alpha,
                                                independentFiltering) {
  if (!is.numeric(lfcThreshold) || length(lfcThreshold) != 1L || !is.finite(lfcThreshold) || lfcThreshold != 0 ||
    !identical(altHypothesis, "greaterAbs") || isTRUE(useT) || !is.null(degreesOfFreedom)) {
    stop("lfcThreshold, altHypothesis, useT, and degreesOfFreedom are only supported for Wald results", call. = FALSE)
  }
  resultsNames <- object$resultsNames
  reference <- reference %||% .rsdeseq2_reference_for_character_contrast(object, contrast)
  if (is.null(contrast)) {
    if (is.null(name)) {
      if (length(resultsNames) < 2L) {
        stop("default results require at least two result names", call. = FALSE)
      }
      name <- resultsNames[[length(resultsNames)]]
    }
    coefficient <- .rsdeseq2_resolve_results_name_index(name, resultsNames)
    numericContrast <- numeric(length(resultsNames))
    numericContrast[[coefficient]] <- 1
    names(numericContrast) <- resultsNames
    resultName <- resultsNames[[coefficient]]
    comparison <- sprintf("coefficient %s; LRT full vs reduced", resultName)
    allZero <- list(type = "none")
  } else {
    resolved <- resolveResultsContrastRust(
      contrast = contrast,
      resultsNames = resultsNames,
      listValues = listValues,
      reference = reference
    )
    numericContrast <- resolved$numeric
    resultName <- resolved$resultName
    comparison <- sprintf("%s; LRT full vs reduced", resolved$comparison)
    allZero <- resolved$allZero
  }

  lrt <- .rsdeseq2_lrt_from_primitive_fit(object, numericContrast)
  lrt <- .rsdeseq2_apply_lrt_results_contrast_all_zero(
    lrt = lrt,
    fit = object,
    numericContrast = numericContrast,
    allZero = allZero
  )
  out <- resultsTableRust(
    baseMean = object$baseMean,
    log2FoldChange = lrt$log2FoldChange,
    lfcSE = lrt$lfcSE,
    stat = object$lrtStat,
    pvalue = object$lrtPvalue,
    rowNames = object$rowNames
  )
  attr(out, "resultName") <- resultName
  attr(out, "comparison") <- comparison
  attr(out, "contrast") <- numericContrast
  attr(out, "contrastAllZero") <- allZero$type
  if (!is.null(object$lrtDf)) {
    attr(out, "lrtDf") <- object$lrtDf
  }
  applyIndependentFilteringRust(
    out,
    alpha = alpha,
    enabled = independentFiltering
  )
}

resultsTableRust <- function(baseMean,
                             log2FoldChange = NULL,
                             lfcSE = NULL,
                             stat = NULL,
                             pvalue = NULL,
                             padj = NULL,
                             dispersion = NULL,
                             converged = NULL,
                             rowNames = NULL) {
  baseMean <- .rsdeseq2_validate_result_numeric_vector(baseMean, "baseMean", allowNA = FALSE)
  nGenes <- length(baseMean)
  rowNames <- .rsdeseq2_validate_row_names(rowNames %||% names(baseMean), nGenes)

  if (is.null(pvalue)) {
    if (!is.null(padj)) {
      stop("padj cannot be supplied when pvalue is NULL", call. = FALSE)
    }
    padj <- NULL
  } else {
    pvalue <- .rsdeseq2_validate_result_pvalues(pvalue, nGenes)
    if (is.null(padj)) {
      padj <- stats::p.adjust(pvalue, method = "BH")
    } else {
      padj <- .rsdeseq2_validate_result_pvalues(padj, nGenes, name = "padj")
    }
  }

  out <- data.frame(
    baseMean = baseMean,
    log2FoldChange = .rsdeseq2_validate_result_optional_numeric(log2FoldChange, nGenes, "log2FoldChange"),
    lfcSE = .rsdeseq2_validate_result_optional_numeric(lfcSE, nGenes, "lfcSE"),
    stat = .rsdeseq2_validate_result_optional_numeric(stat, nGenes, "stat"),
    pvalue = pvalue %||% rep(NA_real_, nGenes),
    padj = padj %||% rep(NA_real_, nGenes),
    check.names = FALSE
  )

  if (!is.null(dispersion)) {
    out$dispersion <- .rsdeseq2_validate_result_optional_numeric(dispersion, nGenes, "dispersion")
  }
  if (!is.null(converged)) {
    out$converged <- .rsdeseq2_validate_result_optional_logical(converged, nGenes, "converged")
  }
  if (!is.null(rowNames)) {
    rownames(out) <- rowNames
  }
  out
}

applyIndependentFilteringRust <- function(results,
                                          alpha = 0.1,
                                          theta = NULL,
                                          enabled = TRUE) {
  results <- .rsdeseq2_validate_result_table(results)
  alpha <- .rsdeseq2_validate_independent_filter_alpha(alpha)

  if (!isTRUE(enabled)) {
    results$padj <- stats::p.adjust(results$pvalue, method = "BH")
    results$filtered <- rep(NA, nrow(results))
    attr(results, "independentFiltering") <- list(
      enabled = FALSE,
      theta = numeric(),
      numRej = integer(),
      selected = NA_integer_,
      filterTheta = NA_real_,
      filterThreshold = NA_real_,
      lo.fit = NULL,
      alpha = alpha
    )
    return(results)
  }

  filter <- results$baseMean
  theta <- .rsdeseq2_validate_independent_filter_theta(theta, filter)
  cutoffs <- as.numeric(stats::quantile(filter, probs = theta, names = FALSE, type = 7))
  filteredPadj <- vector("list", length(cutoffs))
  numRej <- integer(length(cutoffs))

  for (idx in seq_along(cutoffs)) {
    selected <- filter >= cutoffs[[idx]]
    adjusted <- rep(NA_real_, length(filter))
    adjusted[selected] <- stats::p.adjust(results$pvalue[selected], method = "BH")
    filteredPadj[[idx]] <- adjusted
    numRej[[idx]] <- sum(adjusted < alpha, na.rm = TRUE)
  }

  loFit <- stats::lowess(theta, numRej, f = 1 / 5)
  selectedIndex <- .rsdeseq2_select_independent_filter_index(numRej, loFit$y)
  filterThreshold <- cutoffs[[selectedIndex]]

  results$padj <- filteredPadj[[selectedIndex]]
  results$filtered <- ifelse(is.na(results$pvalue), NA, filter < filterThreshold)
  attr(results, "independentFiltering") <- list(
    enabled = TRUE,
    theta = theta,
    numRej = numRej,
    selected = selectedIndex,
    filterTheta = theta[[selectedIndex]],
    filterThreshold = filterThreshold,
    lo.fit = loFit,
    alpha = alpha
  )
  results
}

applyCooksCutoffRust <- function(pvalue,
                                 maxCooks,
                                 cooksCutoff,
                                 counts = NULL,
                                 cooks = NULL,
                                 lowCountHeuristic = FALSE,
                                 native = FALSE) {
  pvalue <- .rsdeseq2_validate_pvalues(pvalue)
  maxCooks <- .rsdeseq2_validate_optional_numeric_vector(maxCooks, length(pvalue), "maxCooks")
  cooksCutoff <- .rsdeseq2_validate_cooks_cutoff(cooksCutoff)

  if (is.null(cooksCutoff)) {
    return(.rsdeseq2_cooks_cutoff_output(pvalue, maxCooks, rep(NA, length(pvalue))))
  }

  if (isTRUE(lowCountHeuristic)) {
    counts <- .rsdeseq2_validate_counts(counts)
    cooks <- .rsdeseq2_validate_cooks_matrix(cooks, nrow(counts), ncol(counts))
    if (nrow(counts) != length(pvalue)) {
      stop("counts must have one row per p-value", call. = FALSE)
    }
  } else {
    counts <- NULL
    cooks <- NULL
  }

  if (isTRUE(native)) {
    nativeMasked <- .rsdeseq2_try_native_cooks_cutoff(
      pvalue = pvalue,
      maxCooks = maxCooks,
      cooksCutoff = cooksCutoff,
      counts = counts,
      cooks = cooks,
      lowCountHeuristic = lowCountHeuristic
    )
    if (!is.null(nativeMasked)) {
      return(.rsdeseq2_cooks_cutoff_output(
        nativeMasked$pvalue,
        maxCooks,
        nativeMasked$cooksOutlier
      ))
    }
  }

  .rsdeseq2_apply_cooks_cutoff_fallback(
    pvalue = pvalue,
    maxCooks = maxCooks,
    cooksCutoff = cooksCutoff,
    counts = counts,
    cooks = cooks,
    lowCountHeuristic = lowCountHeuristic
  )
}

.rsdeseq2_parse_results_contrast <- function(contrast,
                                             resultsNames,
                                             listValues,
                                             reference) {
  if (is.numeric(contrast) && !is.list(contrast)) {
    numeric <- .rsdeseq2_validate_numeric_results_contrast(contrast, resultsNames)
    return(list(
      type = "numeric",
      numeric = numeric,
      resultName = "contrast",
      comparison = "primitive numeric contrast",
      allZero = list(type = "numeric")
    ))
  }

  if (is.character(contrast) && !is.list(contrast)) {
    if (length(contrast) != 3L || anyNA(contrast) || any(!nzchar(contrast))) {
      stop("character contrast must be c(factor, numerator, denominator)", call. = FALSE)
    }
    if (!is.null(reference)) {
      if (!is.character(reference) || length(reference) != 1L || is.na(reference) || !nzchar(reference)) {
        stop("reference must be NULL or a single non-empty character value", call. = FALSE)
      }
    }
    factor <- contrast[[1L]]
    numerator <- contrast[[2L]]
    denominator <- contrast[[3L]]
    numeric <- .rsdeseq2_resolve_factor_level_results_contrast(
      resultsNames = resultsNames,
      factor = factor,
      numerator = numerator,
      denominator = denominator,
      reference = reference
    )
    return(list(
      type = "character",
      factor = factor,
      numerator = numerator,
      denominator = denominator,
      reference = reference,
      numeric = numeric,
      resultName = sprintf("%s_%s_vs_%s", factor, numerator, denominator),
      comparison = sprintf("factor-level contrast: %s %s vs %s", factor, numerator, denominator),
      allZero = list(type = "character", numerator = numerator, denominator = denominator)
    ))
  }

  if (is.list(contrast)) {
    if (length(contrast) < 1L || length(contrast) > 2L) {
      stop("list contrast must contain one or two character vectors", call. = FALSE)
    }
    positive <- contrast[[1L]]
    negative <- if (length(contrast) == 2L) contrast[[2L]] else character()
    if (!is.character(positive) || anyNA(positive) || any(!nzchar(positive))) {
      stop("list contrast numerator must be a character vector", call. = FALSE)
    }
    if (!is.character(negative) || anyNA(negative) || any(!nzchar(negative))) {
      stop("list contrast denominator must be a character vector", call. = FALSE)
    }
    listValues <- .rsdeseq2_validate_list_values(listValues)
    numeric <- .rsdeseq2_resolve_list_results_contrast(
      resultsNames = resultsNames,
      positive = positive,
      negative = negative,
      positiveWeight = listValues[[1L]],
      negativeWeight = listValues[[2L]]
    )
    return(list(
      type = "list",
      positive = positive,
      negative = negative,
      listValues = listValues,
      numeric = numeric,
      resultName = "contrast",
      comparison = .rsdeseq2_list_contrast_comparison(
        positive = positive,
        negative = negative,
        positiveWeight = listValues[[1L]],
        negativeWeight = listValues[[2L]]
      ),
      allZero = list(type = "numeric")
    ))
  }

  stop("contrast must be numeric, character, or list", call. = FALSE)
}

.rsdeseq2_reference_for_character_contrast <- function(object, contrast) {
  if (!is.character(contrast) || length(contrast) != 3L || is.null(object$factorLevels)) {
    return(NULL)
  }
  factor <- contrast[[1L]]
  levels <- object$factorLevels[[factor]]
  if (is.null(levels) || length(levels) == 0L) {
    return(NULL)
  }
  levels[[1L]]
}

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
  colnames(beta) <- resultsNames
  storage.mode(beta) <- "double"
  beta
}

.rsdeseq2_validate_beta_covariance <- function(betaCovariance, nGenes, nCoef) {
  if (is.null(betaCovariance)) {
    return(NULL)
  }
  if (!is.array(betaCovariance) || !is.numeric(betaCovariance) || length(dim(betaCovariance)) != 3L) {
    stop("betaCovariance must be a numeric array with dimensions genes x coefficients x coefficients", call. = FALSE)
  }
  if (!identical(dim(betaCovariance), c(nGenes, nCoef, nCoef))) {
    stop("betaCovariance must have dimensions genes x coefficients x coefficients", call. = FALSE)
  }
  if (any(!is.na(betaCovariance) & !is.finite(betaCovariance))) {
    stop("betaCovariance must contain finite values or NA", call. = FALSE)
  }
  betaCovariance
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
  sampleLevels
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
    levels
  })
  if (anyDuplicated(names(out))) {
    stop("factorLevels names must be unique", call. = FALSE)
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
  colnames(modelMatrix) <- resultsNames
  storage.mode(modelMatrix) <- "double"
  modelMatrix
}

.rsdeseq2_wald_from_primitive_fit <- function(fit, numericContrast, options) {
  beta <- fit$beta
  lfc <- as.numeric(beta %*% numericContrast)
  lfcSE <- .rsdeseq2_contrast_se_from_primitive_fit(fit, numericContrast)
  thresholded <- .rsdeseq2_wald_stat_pvalue_with_options(lfc, lfcSE, options)
  list(
    log2FoldChange = lfc,
    lfcSE = lfcSE,
    stat = thresholded$stat,
    pvalue = thresholded$pvalue
  )
}

.rsdeseq2_lrt_from_primitive_fit <- function(fit, numericContrast) {
  beta <- fit$beta
  lfc <- as.numeric(beta %*% numericContrast)
  lfcSE <- .rsdeseq2_contrast_se_from_primitive_fit(fit, numericContrast)
  list(
    log2FoldChange = lfc,
    lfcSE = lfcSE
  )
}

.rsdeseq2_validate_wald_result_options <- function(lfcThreshold,
                                                   altHypothesis,
                                                   useT,
                                                   degreesOfFreedom,
                                                   nGenes) {
  if (!is.numeric(lfcThreshold) || length(lfcThreshold) != 1L || !is.finite(lfcThreshold) || lfcThreshold < 0) {
    stop("lfcThreshold must be a single finite non-negative number", call. = FALSE)
  }
  if (!is.character(altHypothesis) || length(altHypothesis) != 1L || is.na(altHypothesis)) {
    stop("altHypothesis must be a single character value", call. = FALSE)
  }
  if (identical(altHypothesis, "lessAbs") && lfcThreshold == 0) {
    stop("altHypothesis='lessAbs' requires a positive lfcThreshold", call. = FALSE)
  }
  if (!is.logical(useT) || length(useT) != 1L || is.na(useT)) {
    stop("useT must be TRUE or FALSE", call. = FALSE)
  }
  useT <- isTRUE(useT)
  if (useT && identical(altHypothesis, "greaterAbsUPSHOT")) {
    stop("greaterAbsUPSHOT with useT=TRUE is not supported", call. = FALSE)
  }
  if (useT) {
    if (is.null(degreesOfFreedom)) {
      stop("degreesOfFreedom is required when useT is TRUE", call. = FALSE)
    }
    if (!is.numeric(degreesOfFreedom) || !(length(degreesOfFreedom) %in% c(1L, nGenes))) {
      stop("degreesOfFreedom must be a scalar or one value per gene", call. = FALSE)
    }
    df <- rep(as.numeric(degreesOfFreedom), length.out = nGenes)
    df[!is.finite(df) | df <= 0] <- NA_real_
  } else {
    df <- NULL
  }
  list(
    lfcThreshold = as.numeric(lfcThreshold),
    altHypothesis = altHypothesis,
    useT = useT,
    degreesOfFreedom = df
  )
}

.rsdeseq2_wald_stat_pvalue_with_options <- function(beta, betaSE, options) {
  stat <- rep(NA_real_, length(beta))
  pvalue <- rep(NA_real_, length(beta))
  valid <- !is.na(beta) & !is.na(betaSE) & is.finite(beta) & is.finite(betaSE) & betaSE > 0
  if (!any(valid)) {
    return(list(stat = stat, pvalue = pvalue))
  }
  for (idx in which(valid)) {
    one <- .rsdeseq2_wald_one_with_options(
      beta = beta[[idx]],
      betaSE = betaSE[[idx]],
      gene = idx,
      options = options
    )
    stat[[idx]] <- one$stat
    pvalue[[idx]] <- one$pvalue
  }
  list(stat = stat, pvalue = pvalue)
}

.rsdeseq2_wald_one_with_options <- function(beta, betaSE, gene, options) {
  defaultStat <- beta / betaSE
  threshold <- options$lfcThreshold
  absBeta <- abs(beta)
  tail <- .rsdeseq2_wald_tail(options, gene)
  alt <- options$altHypothesis

  if (identical(alt, "greaterAbs") && threshold == 0) {
    return(list(stat = defaultStat, pvalue = tail$twoSided(defaultStat)))
  }
  if (identical(alt, "greaterAbs")) {
    q1 <- (-absBeta + threshold) / betaSE
    q2 <- (-absBeta - threshold) / betaSE
    return(list(stat = defaultStat, pvalue = .rsdeseq2_clamp_probability(tail$lower(q1) + tail$lower(q2))))
  }
  if (identical(alt, "greaterAbsUPSHOT") && threshold == 0) {
    return(list(stat = defaultStat, pvalue = tail$twoSided(defaultStat)))
  }
  if (identical(alt, "greaterAbsUPSHOT")) {
    return(list(
      stat = defaultStat,
      pvalue = .rsdeseq2_greater_abs_upshot_normal_pvalue(absBeta, betaSE, threshold)
    ))
  }
  if (identical(alt, "greaterAbs2014")) {
    shifted <- (absBeta - threshold) / betaSE
    return(list(
      stat = sign(beta) * max(shifted, 0),
      pvalue = .rsdeseq2_clamp_probability(2 * tail$upper(shifted))
    ))
  }
  if (identical(alt, "lessAbs")) {
    aboveShift <- (threshold - beta) / betaSE
    belowShift <- (beta + threshold) / betaSE
    return(list(
      stat = min(max(aboveShift, 0), max(belowShift, 0)),
      pvalue = max(tail$upper(aboveShift), tail$upper(belowShift))
    ))
  }
  if (identical(alt, "greater")) {
    shifted <- (beta - threshold) / betaSE
    return(list(stat = max(shifted, 0), pvalue = tail$upper(shifted)))
  }
  shiftedStat <- (beta + threshold) / betaSE
  shiftedPvalue <- (-threshold - beta) / betaSE
  list(stat = min(shiftedStat, 0), pvalue = tail$upper(shiftedPvalue))
}

.rsdeseq2_wald_tail <- function(options, gene) {
  if (isTRUE(options$useT)) {
    df <- options$degreesOfFreedom[[gene]]
    if (is.na(df)) {
      return(list(
        twoSided = function(q) NA_real_,
        upper = function(q) NA_real_,
        lower = function(q) NA_real_
      ))
    }
    return(list(
      twoSided = function(q) .rsdeseq2_clamp_probability(2 * stats::pt(abs(q), df = df, lower.tail = FALSE)),
      upper = function(q) .rsdeseq2_clamp_probability(stats::pt(q, df = df, lower.tail = FALSE)),
      lower = function(q) .rsdeseq2_clamp_probability(stats::pt(q, df = df, lower.tail = TRUE))
    ))
  }
  list(
    twoSided = function(q) .rsdeseq2_clamp_probability(2 * stats::pnorm(abs(q), lower.tail = FALSE)),
    upper = function(q) .rsdeseq2_clamp_probability(stats::pnorm(q, lower.tail = FALSE)),
    lower = function(q) .rsdeseq2_clamp_probability(stats::pnorm(q, lower.tail = TRUE))
  )
}

.rsdeseq2_greater_abs_upshot_normal_pvalue <- function(absBeta, betaSE, threshold) {
  a <- (absBeta - threshold) / betaSE
  b <- (absBeta + threshold) / betaSE
  denominator <- b - a
  if (!is.finite(denominator) || denominator == 0) {
    return(NA_real_)
  }
  value <- (2 / denominator) * (-a * stats::pnorm(a, lower.tail = FALSE) +
    stats::dnorm(a) +
    b * stats::pnorm(b, lower.tail = FALSE) -
    stats::dnorm(b))
  .rsdeseq2_clamp_probability(value)
}

.rsdeseq2_clamp_probability <- function(value) {
  if (!is.finite(value)) {
    return(NA_real_)
  }
  max(0, min(1, value))
}

.rsdeseq2_contrast_se_from_primitive_fit <- function(fit, numericContrast) {
  nonzero <- which(numericContrast != 0)
  if (length(nonzero) == 1L && !is.null(fit$betaSE)) {
    out <- abs(numericContrast[[nonzero]]) * fit$betaSE[, nonzero]
    return(as.numeric(out))
  }
  if (is.null(fit$betaCovariance)) {
    stop("betaCovariance is required for multi-coefficient or numeric contrast SEs", call. = FALSE)
  }
  cov <- fit$betaCovariance
  out <- vapply(seq_len(dim(cov)[[1L]]), function(gene) {
    geneCov <- cov[gene, , ]
    if (anyNA(geneCov)) {
      return(NA_real_)
    }
    variance <- as.numeric(t(numericContrast) %*% geneCov %*% numericContrast)
    if (!is.finite(variance) || variance < 0) {
      return(NA_real_)
    }
    sqrt(variance)
  }, numeric(1L))
  out
}

.rsdeseq2_apply_lrt_results_contrast_all_zero <- function(lrt,
                                                          fit,
                                                          numericContrast,
                                                          allZero) {
  flags <- .rsdeseq2_results_contrast_all_zero_flags(
    fit = fit,
    numericContrast = numericContrast,
    allZero = allZero
  )
  if (is.null(flags)) {
    return(lrt)
  }
  allZeroRows <- if (is.null(fit$counts)) rep(FALSE, length(flags)) else rowSums(fit$counts) == 0
  zeroRows <- which(flags & !allZeroRows)
  if (length(zeroRows) != 0L) {
    lrt$log2FoldChange[zeroRows] <- 0
  }
  lrt
}

.rsdeseq2_apply_results_contrast_all_zero <- function(wald,
                                                      fit,
                                                      numericContrast,
                                                      allZero) {
  flags <- .rsdeseq2_results_contrast_all_zero_flags(
    fit = fit,
    numericContrast = numericContrast,
    allZero = allZero
  )
  if (is.null(flags)) {
    return(wald)
  }
  allZeroRows <- if (is.null(fit$counts)) rep(FALSE, length(flags)) else rowSums(fit$counts) == 0
  zeroRows <- which(flags & !allZeroRows)
  if (length(zeroRows) != 0L) {
    wald$log2FoldChange[zeroRows] <- 0
    wald$stat[zeroRows] <- 0
    wald$pvalue[zeroRows] <- 1
  }
  wald
}

.rsdeseq2_results_contrast_all_zero_flags <- function(fit,
                                                      numericContrast,
                                                      allZero) {
  if (is.null(fit$counts)) {
    return(NULL)
  }
  if (is.null(allZero) || identical(allZero$type, "none")) {
    return(NULL)
  }
  if (identical(allZero$type, "character")) {
    if (is.null(fit$sampleLevels)) {
      stop("character contrast all-zero handling requires sampleLevels", call. = FALSE)
    }
    if (!all(c(allZero$numerator, allZero$denominator) %in% fit$sampleLevels)) {
      stop("character contrast all-zero handling requires sampleLevels to contain numerator and denominator levels", call. = FALSE)
    }
    selected <- fit$sampleLevels %in% c(allZero$numerator, allZero$denominator)
    return(rowSums(fit$counts[, selected, drop = FALSE] == 0) == sum(selected))
  }
  if (!identical(allZero$type, "numeric")) {
    stop("unknown contrast all-zero rule", call. = FALSE)
  }
  if (all(numericContrast >= 0) || all(numericContrast <= 0)) {
    return(rep(FALSE, nrow(fit$counts)))
  }
  if (is.null(fit$modelMatrix)) {
    stop("numeric contrast all-zero handling requires modelMatrix", call. = FALSE)
  }
  contrastBinary <- ifelse(numericContrast == 0, 0, 1)
  selected <- as.numeric(fit$modelMatrix %*% contrastBinary) != 0
  rowSums(fit$counts[, selected, drop = FALSE]) == 0
}

.rsdeseq2_validate_results_names <- function(resultsNames) {
  if (!is.character(resultsNames) || length(resultsNames) == 0L || anyNA(resultsNames) || any(!nzchar(resultsNames))) {
    stop("resultsNames must be a non-empty character vector", call. = FALSE)
  }
  resultsNames
}

.rsdeseq2_validate_numeric_results_contrast <- function(contrast, resultsNames) {
  if (length(contrast) != length(resultsNames)) {
    stop("numeric contrast must have one value per results name", call. = FALSE)
  }
  if (any(!is.finite(contrast))) {
    stop("numeric contrast must contain finite values", call. = FALSE)
  }
  if (!any(contrast != 0)) {
    stop("numeric contrast vector cannot be all zero", call. = FALSE)
  }
  out <- as.numeric(contrast)
  names(out) <- resultsNames
  out
}

.rsdeseq2_validate_lrt_stat <- function(lrtStat, nGenes) {
  lrtStat <- .rsdeseq2_validate_result_numeric_vector(lrtStat, "lrtStat", allowNA = TRUE)
  if (length(lrtStat) != nGenes) {
    stop("lrtStat must have one value per beta row", call. = FALSE)
  }
  if (any(!is.na(lrtStat) & lrtStat < 0)) {
    stop("lrtStat must contain non-negative values or NA", call. = FALSE)
  }
  lrtStat
}

.rsdeseq2_validate_lrt_df <- function(lrtDf, nGenes) {
  if (is.null(lrtDf)) {
    return(NULL)
  }
  if (!is.numeric(lrtDf) || !(length(lrtDf) %in% c(1L, nGenes))) {
    stop("lrtDf must be a scalar or one value per gene", call. = FALSE)
  }
  out <- rep(as.numeric(lrtDf), length.out = nGenes)
  if (any(!is.na(out) & (!is.finite(out) | out <= 0))) {
    stop("lrtDf must contain positive finite values or NA", call. = FALSE)
  }
  out
}

.rsdeseq2_validate_list_values <- function(listValues) {
  if (!is.numeric(listValues) || length(listValues) != 2L || any(!is.finite(listValues))) {
    stop("listValues must be a finite numeric vector of length two", call. = FALSE)
  }
  if (listValues[[1L]] <= 0 || listValues[[2L]] >= 0) {
    stop("listValues must have a positive numerator weight and negative denominator weight", call. = FALSE)
  }
  as.numeric(listValues)
}

.rsdeseq2_resolve_list_results_contrast <- function(resultsNames,
                                                    positive,
                                                    negative,
                                                    positiveWeight,
                                                    negativeWeight) {
  if (length(positive) == 0L && length(negative) == 0L) {
    stop("list contrast must contain at least one coefficient name", call. = FALSE)
  }
  positiveIdx <- .rsdeseq2_resolve_results_name_list(resultsNames, positive)
  negativeIdx <- .rsdeseq2_resolve_results_name_list(resultsNames, negative)
  if (length(intersect(positiveIdx, negativeIdx)) != 0L) {
    stop("contrast list entries must not appear in both numerator and denominator", call. = FALSE)
  }
  numeric <- numeric(length(resultsNames))
  numeric[positiveIdx] <- positiveWeight
  numeric[negativeIdx] <- negativeWeight
  names(numeric) <- resultsNames
  .rsdeseq2_validate_numeric_results_contrast(numeric, resultsNames)
}

.rsdeseq2_resolve_results_name_list <- function(resultsNames, names) {
  if (length(names) == 0L) {
    return(integer())
  }
  idx <- vapply(names, .rsdeseq2_resolve_results_name_index, integer(1L), resultsNames = resultsNames)
  unique(idx)
}

.rsdeseq2_resolve_factor_level_results_contrast <- function(resultsNames,
                                                            factor,
                                                            numerator,
                                                            denominator,
                                                            reference) {
  if (identical(numerator, denominator)) {
    stop("contrast numerator and denominator levels must differ", call. = FALSE)
  }
  if (is.null(reference)) {
    direct <- .rsdeseq2_find_first_results_name(
      resultsNames,
      .rsdeseq2_standard_results_names(factor, numerator, denominator)
    )
    if (!is.na(direct)) {
      return(.rsdeseq2_unit_contrast(resultsNames, direct, 1))
    }
    reverse <- .rsdeseq2_find_first_results_name(
      resultsNames,
      .rsdeseq2_standard_results_names(factor, denominator, numerator)
    )
    if (!is.na(reverse)) {
      return(.rsdeseq2_unit_contrast(resultsNames, reverse, -1))
    }
    shared <- .rsdeseq2_find_shared_reference_results_names(resultsNames, factor, numerator, denominator)
    if (!is.null(shared)) {
      out <- numeric(length(resultsNames))
      out[[shared$numerator]] <- 1
      out[[shared$denominator]] <- -1
      names(out) <- resultsNames
      return(out)
    }
  } else {
    if (identical(numerator, reference)) {
      denominatorIdx <- .rsdeseq2_find_first_results_name(
        resultsNames,
        c(
          .rsdeseq2_standard_results_names(factor, denominator, reference),
          .rsdeseq2_expanded_results_names(factor, denominator)
        )
      )
      if (!is.na(denominatorIdx)) {
        return(.rsdeseq2_unit_contrast(resultsNames, denominatorIdx, -1))
      }
    } else if (identical(denominator, reference)) {
      numeratorIdx <- .rsdeseq2_find_first_results_name(
        resultsNames,
        c(
          .rsdeseq2_standard_results_names(factor, numerator, reference),
          .rsdeseq2_expanded_results_names(factor, numerator)
        )
      )
      if (!is.na(numeratorIdx)) {
        return(.rsdeseq2_unit_contrast(resultsNames, numeratorIdx, 1))
      }
    } else {
      numeratorIdx <- .rsdeseq2_find_first_results_name(
        resultsNames,
        .rsdeseq2_standard_results_names(factor, numerator, reference)
      )
      denominatorIdx <- .rsdeseq2_find_first_results_name(
        resultsNames,
        .rsdeseq2_standard_results_names(factor, denominator, reference)
      )
      if (!is.na(numeratorIdx) && !is.na(denominatorIdx)) {
        out <- numeric(length(resultsNames))
        out[[numeratorIdx]] <- 1
        out[[denominatorIdx]] <- -1
        names(out) <- resultsNames
        return(out)
      }
    }
  }

  numeratorIdx <- .rsdeseq2_find_first_results_name(
    resultsNames,
    .rsdeseq2_expanded_results_names(factor, numerator)
  )
  denominatorIdx <- .rsdeseq2_find_first_results_name(
    resultsNames,
    .rsdeseq2_expanded_results_names(factor, denominator)
  )
  if (!is.na(numeratorIdx) && !is.na(denominatorIdx)) {
    out <- numeric(length(resultsNames))
    out[[numeratorIdx]] <- 1
    out[[denominatorIdx]] <- -1
    names(out) <- resultsNames
    return(out)
  }

  stop(
    sprintf(
      "factor-level contrast %s: %s vs %s could not be resolved from resultsNames",
      factor,
      numerator,
      denominator
    ),
    call. = FALSE
  )
}

.rsdeseq2_unit_contrast <- function(resultsNames, index, value) {
  out <- numeric(length(resultsNames))
  out[[index]] <- value
  names(out) <- resultsNames
  out
}

.rsdeseq2_standard_results_names <- function(factor, level, reference) {
  c(
    sprintf("%s_%s_vs_%s", factor, level, reference),
    sprintf("%s%s_vs_%s", factor, level, reference),
    sprintf("%s_%s_vs_%s", make.names(factor), make.names(level), make.names(reference)),
    sprintf("%s%s_vs_%s", make.names(factor), make.names(level), make.names(reference))
  )
}

.rsdeseq2_expanded_results_names <- function(factor, level) {
  c(
    sprintf("%s%s", factor, level),
    sprintf("%s_%s", factor, level),
    sprintf("%s%s", make.names(factor), make.names(level)),
    sprintf("%s_%s", make.names(factor), make.names(level))
  )
}

.rsdeseq2_find_shared_reference_results_names <- function(resultsNames, factor, numerator, denominator) {
  numeratorRows <- .rsdeseq2_standard_results_rows_for_level(resultsNames, factor, numerator)
  denominatorRows <- .rsdeseq2_standard_results_rows_for_level(resultsNames, factor, denominator)
  shared <- list()
  for (numeratorRow in numeratorRows) {
    for (denominatorRow in denominatorRows) {
      if (identical(numeratorRow$reference, denominatorRow$reference)) {
        key <- sprintf("%d:%d", numeratorRow$index, denominatorRow$index)
        shared[[key]] <- list(numerator = numeratorRow$index, denominator = denominatorRow$index)
      }
    }
  }
  if (length(shared) == 0L) {
    return(NULL)
  }
  if (length(shared) > 1L) {
    stop(
      sprintf(
        "factor-level contrast %s: %s vs %s resolves ambiguously through shared-reference result names",
        factor,
        numerator,
        denominator
      ),
      call. = FALSE
    )
  }
  shared[[1L]]
}

.rsdeseq2_standard_results_rows_for_level <- function(resultsNames, factor, level) {
  prefixes <- .rsdeseq2_standard_results_prefixes(factor, level)
  rows <- list()
  for (idx in seq_along(resultsNames)) {
    name <- resultsNames[[idx]]
    for (prefix in prefixes) {
      if (startsWith(name, prefix)) {
        reference <- substring(name, nchar(prefix) + 1L)
        key <- sprintf("%d:%s", idx, reference)
        rows[[key]] <- list(index = idx, reference = reference)
      }
    }
  }
  rows
}

.rsdeseq2_standard_results_prefixes <- function(factor, level) {
  unique(c(
    sprintf("%s_%s_vs_", factor, level),
    sprintf("%s%s_vs_", factor, level),
    sprintf("%s_%s_vs_", make.names(factor), make.names(level)),
    sprintf("%s%s_vs_", make.names(factor), make.names(level))
  ))
}

.rsdeseq2_find_first_results_name <- function(resultsNames, candidates) {
  for (candidate in unique(candidates)) {
    idx <- tryCatch(
      .rsdeseq2_resolve_results_name_index(candidate, resultsNames),
      error = function(error) NA_integer_
    )
    if (!is.na(idx)) {
      return(idx)
    }
  }
  NA_integer_
}

.rsdeseq2_resolve_results_name_index <- function(name, resultsNames) {
  exact <- which(resultsNames == name)
  if (length(exact) == 1L) {
    return(exact)
  }
  if (length(exact) > 1L) {
    stop(sprintf("results name %s is duplicated", name), call. = FALSE)
  }

  alias <- .rsdeseq2_clean_results_name(name)
  aliases <- vapply(resultsNames, .rsdeseq2_clean_results_name, character(1L))
  matches <- which(aliases == alias)
  if (length(matches) == 1L) {
    return(matches)
  }
  if (length(matches) > 1L) {
    stop(sprintf("results name alias %s is ambiguous", name), call. = FALSE)
  }
  stop(sprintf("results name %s not found", name), call. = FALSE)
}

.rsdeseq2_clean_results_name <- function(name) {
  cleaned <- make.names(name)
  cleaned <- gsub("^X.Intercept\\.$", "Intercept", cleaned)
  cleaned <- gsub("^Intercept$", "Intercept", cleaned)
  cleaned
}

.rsdeseq2_list_contrast_comparison <- function(positive,
                                               negative,
                                               positiveWeight,
                                               negativeWeight) {
  if (length(negative) == 0L) {
    return(sprintf(
      "coefficient list contrast: %s effect",
      .rsdeseq2_weighted_name_list(positive, positiveWeight)
    ))
  }
  if (length(positive) == 0L) {
    return(sprintf(
      "coefficient list contrast: %s effect",
      .rsdeseq2_weighted_name_list(negative, negativeWeight)
    ))
  }
  sprintf(
    "coefficient list contrast: %s vs %s",
    .rsdeseq2_weighted_name_list(positive, positiveWeight),
    .rsdeseq2_weighted_name_list(negative, abs(negativeWeight))
  )
}

.rsdeseq2_weighted_name_list <- function(names, weight) {
  prefix <- if (identical(unname(weight), 1)) "" else sprintf("%s ", format(weight, trim = TRUE, scientific = FALSE))
  paste0(prefix, paste(names, collapse = " + "))
}

.rsdeseq2_apply_cooks_cutoff_fallback <- function(pvalue,
                                                  maxCooks,
                                                  cooksCutoff,
                                                  counts,
                                                  cooks,
                                                  lowCountHeuristic) {
  cooksOutlier <- ifelse(is.na(maxCooks), NA, maxCooks > cooksCutoff)
  if (isTRUE(lowCountHeuristic)) {
    for (gene in which(cooksOutlier %in% TRUE)) {
      if (.rsdeseq2_low_count_cooks_spares_row(counts[gene, ], cooks[gene, ])) {
        cooksOutlier[gene] <- FALSE
      }
    }
  }

  masked <- pvalue
  masked[which(cooksOutlier %in% TRUE)] <- NA_real_
  .rsdeseq2_cooks_cutoff_output(masked, maxCooks, cooksOutlier)
}

.rsdeseq2_try_native_cooks_cutoff <- function(pvalue,
                                              maxCooks,
                                              cooksCutoff,
                                              counts,
                                              cooks,
                                              lowCountHeuristic) {
  args <- list(
    pvalue,
    maxCooks,
    cooksCutoff,
    counts %||% NULL,
    cooks %||% NULL,
    as.logical(lowCountHeuristic)
  )
  out <- tryCatch(
    do.call(.Call, c(list("rsdeseq2_apply_cooks_cutoff", PACKAGE = "rsdeseq2"), args)),
    error = function(error) NULL
  )
  if (!is.null(out)) {
    return(out)
  }
  tryCatch(
    do.call(.Call, c(list("rsdeseq2_apply_cooks_cutoff"), args)),
    error = function(error) NULL
  )
}
