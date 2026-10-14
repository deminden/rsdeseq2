.rsdeseq2_expanded_contrast_design <- function(dds, standardDesign, contrast) {
  formula <- methods::slot(dds, "design")
  if (!inherits(formula, "formula")) {
    stop("normal contrast shrinkage requires a formula design", call. = FALSE)
  }
  if (any(attr(stats::terms(formula), "order") > 1L)) {
    stop("normal contrast shrinkage does not support formula interactions", call. = FALSE)
  }
  colData <- as.data.frame(SummarizedExperiment::colData(dds))
  modelFrame <- stats::model.frame(formula, data = colData)
  factorNames <- names(modelFrame)[vapply(modelFrame, is.factor, logical(1))]
  if (!all(factorNames %in% names(colData))) {
    stop("normal contrast shrinkage currently requires factors used directly in the formula", call. = FALSE)
  }
  if (is.character(contrast)) {
    if (length(contrast) != 3L || anyNA(contrast) || any(!nzchar(contrast))) {
      stop("character contrast must be c(factor, numerator, denominator)", call. = FALSE)
    }
    factorName <- contrast[[1L]]
    if (!factorName %in% factorNames) {
      stop(sprintf("contrast factor '%s' must name a factor in the formula", factorName), call. = FALSE)
    }
    selectedLevels <- base::levels(modelFrame[[factorName]])
    if (!contrast[[2L]] %in% selectedLevels || !contrast[[3L]] %in% selectedLevels) {
      stop("contrast numerator and denominator must be declared factor levels", call. = FALSE)
    }
    if (identical(contrast[[2L]], contrast[[3L]])) {
      stop("contrast numerator and denominator must differ", call. = FALSE)
    }
  } else if (is.numeric(contrast) && is.null(dim(contrast))) {
    factorName <- NULL
  } else if (is.list(contrast)) {
    if (length(contrast) == 1L) contrast <- list(contrast[[1L]], character())
    if (length(contrast) != 2L || !all(vapply(contrast, is.character, logical(1)))) {
      stop("list contrast must contain one or two character vectors", call. = FALSE)
    }
    if (anyNA(unlist(contrast, use.names = FALSE)) || any(!nzchar(unlist(contrast, use.names = FALSE)))) {
      stop("list contrast coefficient names must be non-empty", call. = FALSE)
    }
    contrastNames <- unlist(contrast, use.names = FALSE)
    if (anyDuplicated(contrastNames) || length(contrastNames) == 0L) {
      stop("list contrast must select distinct expanded-design coefficients", call. = FALSE)
    }
    factorName <- NULL
  } else {
    stop("contrast must be a factor-level triplet, list, or numeric expanded-design vector", call. = FALSE)
  }
  nullLevel <- "_null_level_"
  if (any(vapply(factorNames, function(name) nullLevel %in% levels(modelFrame[[name]]), logical(1)))) {
    stop("a formula factor already uses the reserved expanded-model level '_null_level_'", call. = FALSE)
  }
  expandedData <- rbind(colData, colData[nrow(colData), , drop = FALSE])
  for (name in factorNames) {
    value <- expandedData[[name]]
    levels(value) <- c(levels(value), nullLevel)
    value <- stats::relevel(value, ref = nullLevel)
    value[nrow(value)] <- nullLevel
    expandedData[[name]] <- value
  }
  expanded <- stats::model.matrix(formula, data = expandedData)
  nullColumns <- grepl(nullLevel, colnames(expanded), fixed = TRUE)
  expanded <- expanded[-nrow(expanded), !nullColumns, drop = FALSE]
  colnames(expanded)[colnames(expanded) == "(Intercept)"] <- "Intercept"
  colnames(expanded) <- make.names(colnames(expanded))
  rowOrder <- match(rownames(SummarizedExperiment::colData(dds)), rownames(expanded))
  if (anyNA(rowOrder)) stop("expanded design rows do not match dataset samples", call. = FALSE)
  expanded <- expanded[rowOrder, , drop = FALSE]

  expandedGroups <- integer(ncol(expanded))
  standardGroups <- integer(ncol(standardDesign))
  expandedPriorMap <- integer(ncol(expanded))
  standardNames <- make.names(ifelse(colnames(standardDesign) == "(Intercept)", "Intercept", colnames(standardDesign)))
  expandedFactorColumns <- vector("list", length(factorNames))
  names(expandedFactorColumns) <- factorNames
  for (group in seq_along(factorNames)) {
    name <- factorNames[[group]]
    factorLevels <- levels(modelFrame[[name]])
    expandedNames <- make.names(paste0(name, factorLevels))
    expandedColumns <- match(expandedNames, colnames(expanded))
    if (anyNA(expandedColumns)) {
      stop(sprintf("could not map expanded coefficients for factor '%s'", name), call. = FALSE)
    }
    expandedGroups[expandedColumns] <- group
    expandedFactorColumns[[name]] <- expandedColumns
    standardFactorLevels <- if (attr(stats::terms(formula), "intercept") == 1L) {
      factorLevels[-1L]
    } else {
      factorLevels
    }
    standardNamesForLevels <- make.names(paste0(name, standardFactorLevels))
    standardColumns <- match(standardNamesForLevels, standardNames)
    standardColumns <- standardColumns[!is.na(standardColumns)]
    if (length(standardColumns)) standardGroups[standardColumns] <- group
  }
  for (column in seq_len(ncol(expanded))) {
    if (expandedGroups[[column]] > 0L) next
    standardColumn <- match(colnames(expanded)[[column]], standardNames)
    if (is.na(standardColumn)) {
      stop(sprintf("could not map expanded coefficient '%s' to the fitted design", colnames(expanded)[[column]]), call. = FALSE)
    }
    expandedPriorMap[[column]] <- standardColumn
  }
  if (is.numeric(contrast)) {
    if (length(contrast) != ncol(expanded)) {
      stop("numeric contrast vector must have one value per expanded-design coefficient", call. = FALSE)
    }
    if (any(!is.finite(contrast)) || all(contrast == 0)) {
      stop("numeric contrast vector must be finite and contain a non-zero value", call. = FALSE)
    }
    contrastVector <- as.numeric(contrast)
  } else if (is.list(contrast)) {
    contrastVector <- numeric(ncol(expanded))
    positive <- match(contrast[[1L]], colnames(expanded))
    negative <- match(contrast[[2L]], colnames(expanded))
    if (anyNA(positive) || anyNA(negative)) {
      stop("all list-contrast coefficients must name expanded-design columns", call. = FALSE)
    }
    if (length(intersect(positive, negative))) {
      stop("list contrast coefficients must not appear in both numerator and denominator", call. = FALSE)
    }
    if (length(positive)) contrastVector[positive] <- 1 / length(positive)
    if (length(negative)) contrastVector[negative] <- -1 / length(negative)
  } else {
    contrastVector <- numeric(ncol(expanded))
    selectedColumns <- expandedFactorColumns[[factorName]]
    contrastVector[[selectedColumns[match(contrast[[2L]], selectedLevels)]]] <- 1
    contrastVector[[selectedColumns[match(contrast[[3L]], selectedLevels)]]] <- -1
  }
  list(
    design = expanded,
    expandedGroups = as.integer(expandedGroups),
    standardGroups = as.integer(standardGroups),
    expandedPriorMap = as.integer(expandedPriorMap),
    contrast = contrastVector,
    interceptIndex = match("Intercept", colnames(expanded)),
    standardInterceptIndex = match("Intercept", standardNames)
  )
}

.rsdeseq2_lfc_shrink_dataset <- function(dds, coef, contrast, res, lfcThreshold,
                                         independentFiltering, alpha) {
  if (!inherits(dds, "DESeqDataSet")) {
    stop("object must be a DESeqDataSet", call. = FALSE)
  }
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE) ||
      !requireNamespace("S4Vectors", quietly = TRUE) ||
      !requireNamespace("DESeq2", quietly = TRUE)) {
    stop("SummarizedExperiment, S4Vectors, and DESeq2 are required for DESeqDataSet shrinkage", call. = FALSE)
  }
  fit <- S4Vectors::metadata(dds)$rsdeseq2RustFit
  if (is.null(fit) || !inherits(fit, "rsdeseq2RustFit")) {
    stop("DESeqDataSet must be fitted with rsdeseq2::DESeq() before calling lfcShrink()", call. = FALSE)
  }
  if (!identical(fit$test, "Wald")) {
    stop("lfcShrink() requires a Wald fit", call. = FALSE)
  }
  formula <- methods::slot(dds, "design")
  if (inherits(formula, "formula") && any(attr(stats::terms(formula), "order") > 1L)) {
    stop("normal lfcShrink does not support formula interactions", call. = FALSE)
  }
  resultsNames <- fit$resultsNames
  if (is.null(contrast)) {
    if (length(coef) != 1L || is.na(coef)) {
      stop("supply either coef or contrast", call. = FALSE)
    }
    if (is.numeric(coef)) {
      if (!is.finite(coef) || coef != floor(coef) || coef < 1L || coef > length(resultsNames)) {
        stop("numeric coef must be a one-based coefficient index", call. = FALSE)
      }
      coefficient <- as.integer(coef)
    } else if (is.character(coef)) {
      coefficient <- tryCatch(
        .rsdeseq2_resolve_results_name_index(coef, resultsNames),
        error = function(error) NA_integer_
      )
      if (is.na(coefficient)) {
        stop(sprintf("coef '%s' is not in resultsNames(object)", coef), call. = FALSE)
      }
    } else {
      stop("coef must be a coefficient name or one-based integer index", call. = FALSE)
    }
    resultName <- resultsNames[[coefficient]]
  } else {
    if (!is.null(coef)) stop("supply either coef or contrast, not both", call. = FALSE)
    if (is.character(contrast)) {
      if (length(contrast) != 3L || anyNA(contrast) || any(!nzchar(contrast))) {
        stop("character contrast must be c(factor, numerator, denominator)", call. = FALSE)
      }
      resultName <- sprintf("%s_%s_vs_%s", contrast[[1L]], contrast[[2L]], contrast[[3L]])
    } else if (is.numeric(contrast) && is.null(dim(contrast))) {
      if (is.null(res)) {
        stop("normal shrinkage with a numeric contrast requires res, matching DESeq2", call. = FALSE)
      }
      resultName <- "numeric contrast"
    } else if (is.list(contrast)) {
      resultName <- "list contrast"
    } else {
      stop("normal contrast shrinkage supports factor-level, list, or numeric expanded-design contrasts", call. = FALSE)
    }
    coefficient <- NA_integer_
  }
  counts <- SummarizedExperiment::assay(dds, "counts")
  counts <- .rsdeseq2_validate_counts(counts)
  design <- fit$design
  if (!is.matrix(design) || nrow(design) != ncol(counts)) {
    stop("the stored Rust fit has no compatible design matrix", call. = FALSE)
  }
  if (is.null(fit$dispersions) || is.null(fit$dispFit) ||
      length(fit$dispersions) != nrow(counts) || length(fit$dispFit) != nrow(counts)) {
    stop("the stored Rust fit has no compatible final dispersions and dispersion trend", call. = FALSE)
  }
  interceptMatches <- which(colnames(design) %in% c("(Intercept)", "Intercept"))
  interceptIndex <- if (length(interceptMatches)) interceptMatches[[1L]] else NA_integer_
  if (!is.null(res)) {
    if (!is.data.frame(res) && !methods::is(res, "DataFrame")) {
      stop("res must be a DESeqResults or data-frame-like result table", call. = FALSE)
    }
    res <- as.data.frame(res)
  } else if (is.null(contrast)) {
    res <- as.data.frame(results(
      dds, name = resultName, lfcThreshold = lfcThreshold,
      independentFiltering = independentFiltering, alpha = alpha
    ))
  } else if (is.numeric(contrast)) {
    stop("normal shrinkage with a numeric contrast requires a supplied result table", call. = FALSE)
  } else {
    res <- as.data.frame(results(
      dds, contrast = contrast, lfcThreshold = lfcThreshold,
      independentFiltering = independentFiltering, alpha = alpha
    ))
  }
  if (!all(c("log2FoldChange", "lfcSE") %in% names(res))) {
    stop("res must contain log2FoldChange and lfcSE columns", call. = FALSE)
  }
  geneOrder <- match(rownames(res), rownames(counts))
  if (is.null(rownames(res)) || anyNA(geneOrder) || anyDuplicated(rownames(res))) {
    stop("res row names must uniquely match count rows", call. = FALSE)
  }
  if (!identical(geneOrder, seq_len(nrow(counts)))) {
    counts <- counts[geneOrder, , drop = FALSE]
    fit$beta <- fit$beta[geneOrder, , drop = FALSE]
    fit$dispersions <- fit$dispersions[geneOrder]
    fit$dispFit <- fit$dispFit[geneOrder]
    fit$baseMean <- fit$baseMean[geneOrder]
    if (!is.null(fit$normalizationFactors)) fit$normalizationFactors <- fit$normalizationFactors[geneOrder, , drop = FALSE]
    if (!is.null(fit$weights)) fit$weights <- fit$weights[geneOrder, , drop = FALSE]
  }
  allZero <- rowSums(counts) == 0
  invalidDispersion <- !is.finite(fit$dispersions) | fit$dispersions <= 0 |
    !is.finite(fit$dispFit) | fit$dispFit <= 0
  if (any(invalidDispersion & !allZero)) {
    stop("normal lfcShrink requires finite positive dispersions and fitted trends for non-zero genes", call. = FALSE)
  }
  useGene <- which(!allZero)
  if (!length(useGene)) {
    stop("normal lfcShrink cannot estimate a beta prior when every gene is all zero", call. = FALSE)
  }
  shrinkCounts <- counts[useGene, , drop = FALSE]
  shrinkNormalizationFactors <- if (is.null(fit$normalizationFactors)) NULL else {
    fit$normalizationFactors[useGene, , drop = FALSE]
  }
  shrinkWeights <- if (is.null(fit$weights)) NULL else fit$weights[useGene, , drop = FALSE]
  if (is.null(contrast)) {
    shrunk <- .Call(
      "rsdeseq2_lfc_shrink_normal",
      shrinkCounts,
      design,
      if (is.null(fit$normalizationFactors)) fit$sizeFactors else NULL,
      shrinkNormalizationFactors,
      shrinkWeights,
      as.numeric(fit$dispersions[useGene]),
      as.numeric(fit$baseMean[useGene]),
      as.numeric(fit$dispFit[useGene]),
      as.integer(coefficient),
      as.integer(if (is.na(interceptIndex)) -1L else interceptIndex - 1L),
      as.numeric(lfcThreshold),
      PACKAGE = "rsdeseq2"
    )
    priorNames <- resultsNames
  } else {
    expanded <- .rsdeseq2_expanded_contrast_design(dds, design, contrast)
    shrunk <- .Call(
      "rsdeseq2_lfc_shrink_normal_contrast",
      shrinkCounts,
      expanded$design,
      if (is.null(fit$normalizationFactors)) fit$sizeFactors else NULL,
      shrinkNormalizationFactors,
      shrinkWeights,
      as.numeric(fit$dispersions[useGene]),
      as.numeric(fit$baseMean[useGene]),
      as.numeric(fit$dispFit[useGene]),
      fit$beta[useGene, , drop = FALSE],
      expanded$standardGroups,
      expanded$expandedPriorMap,
      expanded$expandedGroups,
      expanded$contrast,
      as.integer(if (is.na(expanded$interceptIndex)) 0L else expanded$interceptIndex),
      as.integer(if (is.na(expanded$standardInterceptIndex)) 0L else expanded$standardInterceptIndex),
      as.numeric(lfcThreshold),
      PACKAGE = "rsdeseq2"
    )
    priorNames <- colnames(expanded$design)
  }
  posteriorBeta <- stats::setNames(rep(NA_real_, nrow(res)), rownames(res))
  posteriorSE <- posteriorBeta
  names(shrunk$beta) <- rownames(shrinkCounts)
  names(shrunk$betaSE) <- rownames(shrinkCounts)
  posteriorBeta[names(shrunk$beta)] <- shrunk$beta
  posteriorSE[names(shrunk$betaSE)] <- shrunk$betaSE
  res$log2FoldChange <- posteriorBeta
  res$lfcSE <- posteriorSE
  independentFilteringInfo <- attr(res, "independentFiltering")
  independentFilteringInfo <- attr(res, "independentFiltering")
  if (lfcThreshold > 0) {
    posteriorStat <- stats::setNames(rep(NA_real_, nrow(res)), rownames(res))
    posteriorPvalue <- posteriorStat
    names(shrunk$stat) <- rownames(shrinkCounts)
    names(shrunk$pvalue) <- rownames(shrinkCounts)
    posteriorStat[names(shrunk$stat)] <- shrunk$stat
    posteriorPvalue[names(shrunk$pvalue)] <- shrunk$pvalue
    res$stat <- posteriorStat
    res$pvalue <- posteriorPvalue
    res <- applyIndependentFilteringRust(
      res, alpha = alpha, enabled = independentFiltering
    )
    independentFilteringInfo <- attr(res, "independentFiltering")
  }
  resultData <- S4Vectors::DataFrame(res)
  resultMetadata <- S4Vectors::metadata(resultData)
  resultMetadata$lfcThreshold <- lfcThreshold
  resultMetadata$resultName <- resultName
  if (!is.null(independentFilteringInfo)) {
    resultMetadata$independentFiltering <- independentFilteringInfo
  }
  S4Vectors::metadata(resultData) <- resultMetadata
  descriptions <- stats::setNames(rep("", ncol(resultData)), names(resultData))
  types <- stats::setNames(rep("results", ncol(resultData)), names(resultData))
  if ("baseMean" %in% names(descriptions)) {
    descriptions["baseMean"] <- "mean of normalized counts for all samples"
    types["baseMean"] <- "intermediate"
  }
  if ("log2FoldChange" %in% names(descriptions)) descriptions["log2FoldChange"] <- paste("log2 fold change (MAP):", resultName)
  if ("lfcSE" %in% names(descriptions)) descriptions["lfcSE"] <- paste("posterior standard deviation:", resultName)
  if ("stat" %in% names(descriptions)) descriptions["stat"] <- paste("Wald statistic:", resultName)
  if ("pvalue" %in% names(descriptions)) descriptions["pvalue"] <- paste("Wald test p-value:", resultName)
  if ("padj" %in% names(descriptions)) descriptions["padj"] <- "BH adjusted p-values"
  S4Vectors::mcols(resultData) <- S4Vectors::DataFrame(
    type = unname(types),
    description = unname(descriptions),
    row.names = names(resultData)
  )
  DESeq2::DESeqResults(
    resultData,
    priorInfo = list(
      type = "normal",
      package = "rsdeseq2",
      version = as.character(utils::packageVersion("rsdeseq2")),
      betaPriorVar = stats::setNames(shrunk$betaPriorVar, priorNames)
    )
  )
}

.rsdeseq2_lfc_shrink_apeglm <- function(dds, coef, contrast, res,
                                        lfcThreshold, svalue,
                                        independentFiltering, alpha,
                                        apeAdapt, apeMethod, quiet, dots) {
  if (!is.null(contrast)) {
    stop("type='apeglm' shrinkage supports coefficients only", call. = FALSE)
  }
  if (!inherits(dds, "DESeqDataSet")) {
    stop("object must be a DESeqDataSet", call. = FALSE)
  }
  for (package in c("SummarizedExperiment", "S4Vectors", "DESeq2", "apeglm")) {
    if (!requireNamespace(package, quietly = TRUE)) {
      stop(sprintf("type='apeglm' requires package '%s'", package), call. = FALSE)
    }
  }
  fit <- S4Vectors::metadata(dds)$rsdeseq2RustFit
  if (is.null(fit) || !inherits(fit, "rsdeseq2RustFit")) {
    stop("DESeqDataSet must be fitted with rsdeseq2::DESeq() before calling lfcShrink()", call. = FALSE)
  }
  if (!identical(fit$test, "Wald")) {
    stop("lfcShrink() requires a Wald fit", call. = FALSE)
  }
  resultsNames <- fit$resultsNames
  if (length(coef) != 1L || is.na(coef)) {
    stop("type='apeglm' requires one coefficient", call. = FALSE)
  }
  if (is.numeric(coef)) {
    if (!is.finite(coef) || coef != floor(coef) || coef < 1L || coef > length(resultsNames)) {
      stop("numeric coef must be a one-based coefficient index", call. = FALSE)
    }
    coefficient <- as.integer(coef)
  } else if (is.character(coef)) {
    coefficient <- match(coef, resultsNames)
    if (is.na(coefficient)) {
      stop(sprintf("coef '%s' is not in resultsNames(object)", coef), call. = FALSE)
    }
  } else {
    stop("coef must be a coefficient name or one-based integer index", call. = FALSE)
  }
  resultName <- resultsNames[[coefficient]]
  if (is.null(res)) {
    res <- results(
      dds, name = resultName, lfcThreshold = lfcThreshold,
      independentFiltering = independentFiltering, alpha = alpha
    )
  }
  if (!is.data.frame(res) && !methods::is(res, "DataFrame")) {
    stop("res must be a DESeqResults or data-frame-like result table", call. = FALSE)
  }
  resObject <- res
  res <- as.data.frame(res)
  if (!all(c("log2FoldChange", "lfcSE") %in% names(res))) {
    stop("res must contain log2FoldChange and lfcSE columns", call. = FALSE)
  }
  counts <- .rsdeseq2_validate_counts(SummarizedExperiment::assay(dds, "counts"))
  if (is.null(rownames(res)) || anyDuplicated(rownames(res)) ||
      !setequal(rownames(res), rownames(counts))) {
    stop("res row names must uniquely match count rows", call. = FALSE)
  }
  res <- res[match(rownames(counts), rownames(res)), , drop = FALSE]
  if (all(is.na(res$lfcSE))) {
    stop("lfcShrink requires standard errors, use a supported dispersion fit", call. = FALSE)
  }
  if (apeAdapt && methods::is(resObject, "DESeqResults")) {
    descriptions <- S4Vectors::mcols(resObject)$description
    if (length(descriptions) >= 2L) {
      incomingCoef <- gsub(" ", "_", sub("log2 fold change \\(MLE\\): ", "", descriptions[[2L]]))
      if (!identical(resultName, incomingCoef)) {
        stop("'coef' should specify the same coefficient as in results 'res'", call. = FALSE)
      }
    }
  }
  design <- fit$design
  dispersions <- fit$dispersions
  if (!is.matrix(design) || nrow(design) != ncol(counts) ||
      length(dispersions) != nrow(counts)) {
    stop("the stored Rust fit has no compatible design and dispersions", call. = FALSE)
  }
  if (is.null(fit$normalizationFactors)) {
    sizeFactors <- fit$sizeFactors
    if (length(sizeFactors) != ncol(counts) || any(!is.finite(sizeFactors)) || any(sizeFactors <= 0)) {
      stop("the stored Rust fit has no compatible size factors", call. = FALSE)
    }
    offset <- matrix(rep(log(sizeFactors), each = nrow(counts)), nrow = nrow(counts))
  } else {
    normalizationFactors <- fit$normalizationFactors
    if (!all(dim(normalizationFactors) == dim(counts)) ||
        any(!is.finite(normalizationFactors)) || any(normalizationFactors <= 0)) {
      stop("the stored Rust fit has no compatible normalization factors", call. = FALSE)
    }
    offset <- log(normalizationFactors)
  }
  weights <- fit$weights
  if (is.null(weights)) {
    weights <- matrix(1, nrow = nrow(counts), ncol = ncol(counts))
  } else if (!all(dim(weights) == dim(counts)) || any(!is.finite(weights)) || any(weights <= 0)) {
    stop("the stored Rust fit has no compatible observation weights", call. = FALSE)
  }
  apeArguments <- utils::modifyList(
    list(
      Y = counts,
      x = design,
      log.lik = if (identical(apeMethod, "general")) apeglm::logLikNB else NULL,
      param = as.numeric(dispersions),
      coef = coefficient,
      mle = if (apeAdapt) log(2) * cbind(res$log2FoldChange, res$lfcSE) else NULL,
      threshold = if (lfcThreshold > 0) log(2) * lfcThreshold else NULL,
      weights = weights,
      offset = offset,
      method = apeMethod
    ),
    dots
  )
  apeFit <- do.call(apeglm::apeglm, apeArguments)
  if (nrow(apeFit$map) != nrow(counts)) {
    stop("apeglm returned a coefficient matrix with an incompatible row count", call. = FALSE)
  }
  if (!quiet && !is.null(apeFit$diag) && "conv" %in% colnames(apeFit$diag) &&
      any(apeFit$diag[, "conv"][!is.na(apeFit$diag[, "conv"])] != 0)) {
    message("some rows did not converge in finding the MAP")
  }
  res$log2FoldChange <- log2(exp(1)) * apeFit$map[, coefficient]
  names(res$log2FoldChange) <- rownames(res)
  res$lfcSE <- if (is.null(apeFit$sd)) {
    rep(NA_real_, nrow(res))
  } else {
    log2(exp(1)) * apeFit$sd[, coefficient]
  }
  names(res$lfcSE) <- rownames(res)
  returnSvalue <- isTRUE(svalue) || lfcThreshold > 0
  resultColumns <- if (returnSvalue) {
    c("baseMean", "log2FoldChange", "lfcSE")
  } else {
    c("baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj")
  }
  resultColumns <- intersect(resultColumns, names(res))
  res <- res[, resultColumns, drop = FALSE]
  if (returnSvalue) {
    res$svalue <- if (lfcThreshold > 0) {
      as.numeric(apeglm::svalue(apeFit$thresh))
    } else {
      as.numeric(apeFit$svalue)
    }
    names(res$svalue) <- rownames(res)
  }
  resultData <- S4Vectors::DataFrame(res)
  resultMetadata <- S4Vectors::metadata(resultData)
  resultMetadata$lfcThreshold <- lfcThreshold
  resultMetadata$resultName <- resultName
  S4Vectors::metadata(resultData) <- resultMetadata
  descriptions <- stats::setNames(rep("", ncol(resultData)), names(resultData))
  types <- stats::setNames(rep("results", ncol(resultData)), names(resultData))
  if ("baseMean" %in% names(descriptions)) {
    descriptions["baseMean"] <- "mean of normalized counts for all samples"
    types["baseMean"] <- "intermediate"
  }
  if ("log2FoldChange" %in% names(descriptions)) descriptions["log2FoldChange"] <- paste("log2 fold change (MAP):", resultName)
  if ("lfcSE" %in% names(descriptions)) descriptions["lfcSE"] <- paste("posterior standard deviation:", resultName)
  if ("svalue" %in% names(descriptions)) descriptions["svalue"] <- if (lfcThreshold > 0) paste0("FSOS s-value (T=", round(lfcThreshold, 3), "): ", gsub("_", " ", resultName)) else paste("s-value:", gsub("_", " ", resultName))
  S4Vectors::mcols(resultData) <- S4Vectors::DataFrame(
    type = unname(types),
    description = unname(descriptions),
    row.names = names(resultData)
  )
  DESeq2::DESeqResults(
    resultData,
    priorInfo = list(
      type = "apeglm",
      package = "apeglm",
      version = as.character(utils::packageVersion("apeglm")),
      prior.control = apeFit$prior.control
    )
  )
}

.rsdeseq2_lfc_shrink_ashr <- function(dds, coef, contrast, res, lfcThreshold,
                                     svalue, independentFiltering, alpha, dots) {
  if (!inherits(dds, "DESeqDataSet")) {
    stop("object must be a DESeqDataSet", call. = FALSE)
  }
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE) ||
      !requireNamespace("S4Vectors", quietly = TRUE) ||
      !requireNamespace("DESeq2", quietly = TRUE)) {
    stop("SummarizedExperiment, S4Vectors, and DESeq2 are required for DESeqDataSet shrinkage", call. = FALSE)
  }
  if (!requireNamespace("ashr", quietly = TRUE)) {
    stop("type='ashr' requires installing the optional package 'ashr'", call. = FALSE)
  }
  fit <- S4Vectors::metadata(dds)$rsdeseq2RustFit
  if (is.null(fit) || !inherits(fit, "rsdeseq2RustFit")) {
    stop("DESeqDataSet must be fitted with rsdeseq2::DESeq() before calling lfcShrink()", call. = FALSE)
  }
  if (!identical(fit$test, "Wald")) {
    stop("lfcShrink() requires a Wald fit", call. = FALSE)
  }
  resultsNames <- fit$resultsNames
  if (is.null(contrast)) {
    if (length(coef) != 1L || is.na(coef)) {
      stop("supply either coef or contrast", call. = FALSE)
    }
    if (is.numeric(coef)) {
      if (!is.finite(coef) || coef != floor(coef) || coef < 1L || coef > length(resultsNames)) {
        stop("numeric coef must be a one-based coefficient index", call. = FALSE)
      }
      resultName <- resultsNames[[as.integer(coef)]]
    } else if (is.character(coef)) {
      resultName <- coef
      if (!resultName %in% resultsNames) {
        stop(sprintf("coef '%s' is not in resultsNames(object)", resultName), call. = FALSE)
      }
    } else {
      stop("coef must be a coefficient name or one-based integer index", call. = FALSE)
    }
  } else {
    if (!is.null(coef)) stop("supply either coef or contrast, not both", call. = FALSE)
    if (is.character(contrast)) {
      if (length(contrast) != 3L || anyNA(contrast) || any(!nzchar(contrast))) {
        stop("character contrast must be c(factor, numerator, denominator)", call. = FALSE)
      }
      resultName <- sprintf("%s_%s_vs_%s", contrast[[1L]], contrast[[2L]], contrast[[3L]])
    } else if (is.numeric(contrast) && is.null(dim(contrast))) {
      if (is.null(res)) {
        stop("ashr shrinkage with a numeric contrast requires res, matching DESeq2", call. = FALSE)
      }
      resultName <- "numeric contrast"
    } else if (is.list(contrast)) {
      resultName <- "list contrast"
    } else {
      stop("contrast must be a factor-level, list, or numeric contrast", call. = FALSE)
    }
  }
  if (is.null(res)) {
    res <- if (is.null(contrast)) {
      results(
        dds, name = resultName, lfcThreshold = lfcThreshold,
        independentFiltering = independentFiltering, alpha = alpha
      )
    } else {
      results(
        dds, contrast = contrast, lfcThreshold = lfcThreshold,
        independentFiltering = independentFiltering, alpha = alpha
      )
    }
  }
  if (!is.data.frame(res) && !methods::is(res, "DataFrame")) {
    stop("res must be a DESeqResults or data-frame-like result table", call. = FALSE)
  }
  res <- as.data.frame(res)
  if (!all(c("log2FoldChange", "lfcSE") %in% names(res))) {
    stop("res must contain log2FoldChange and lfcSE columns", call. = FALSE)
  }
  counts <- SummarizedExperiment::assay(dds, "counts")
  if (is.null(rownames(res)) || anyDuplicated(rownames(res)) ||
      !setequal(rownames(res), rownames(counts))) {
    stop("res row names must uniquely match count rows", call. = FALSE)
  }
  res <- res[match(rownames(counts), rownames(res)), , drop = FALSE]
  if (all(is.na(res$lfcSE))) {
    stop("lfcShrink requires standard errors, use a supported dispersion fit", call. = FALSE)
  }
  betaHat <- res$log2FoldChange
  seBetaHat <- res$lfcSE
  ashArguments <- utils::modifyList(
    list(
      betahat = betaHat,
      sebetahat = seBetaHat,
      mixcompdist = "normal",
      method = "shrink"
    ),
    dots
  )
  ashFit <- do.call(ashr::ash, ashArguments)
  res$log2FoldChange <- ashFit$result$PosteriorMean
  res$lfcSE <- ashFit$result$PosteriorSD
  independentFilteringInfo <- attr(res, "independentFiltering")
  if (lfcThreshold > 0) {
    cdfPositive <- ashr::cdf_post(
      ashFit$fitted_g, lfcThreshold,
      ashr::set_data(betaHat, seBetaHat)
    )
    cdfNegative <- ashr::cdf_post(
      ashFit$fitted_g, -lfcThreshold,
      ashr::set_data(betaHat, seBetaHat)
    )
    lfsr <- ifelse(res$log2FoldChange > 0, cdfPositive, 1 - cdfNegative)
    sortedLfsr <- sort(lfsr, na.last = TRUE)
    res$svalue <- cumsum(sortedLfsr) / seq_along(sortedLfsr)
    res$svalue <- res$svalue[rank(lfsr, ties.method = "first", na.last = TRUE)]
  } else if (svalue) {
    res$svalue <- ashFit$result$svalue
    res <- res[, c("baseMean", "log2FoldChange", "lfcSE", "svalue"), drop = FALSE]
  } else {
    outputColumns <- intersect(
      c("baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj"),
      names(res)
    )
    res <- res[, outputColumns, drop = FALSE]
  }
  resultData <- S4Vectors::DataFrame(res)
  resultMetadata <- S4Vectors::metadata(resultData)
  resultMetadata$lfcThreshold <- lfcThreshold
  resultMetadata$resultName <- resultName
  if (!is.null(independentFilteringInfo)) {
    resultMetadata$independentFiltering <- independentFilteringInfo
  }
  S4Vectors::metadata(resultData) <- resultMetadata
  descriptions <- stats::setNames(rep("", ncol(resultData)), names(resultData))
  types <- stats::setNames(rep("results", ncol(resultData)), names(resultData))
  if ("baseMean" %in% names(descriptions)) {
    descriptions["baseMean"] <- "mean of normalized counts for all samples"
    types["baseMean"] <- "intermediate"
  }
  if ("log2FoldChange" %in% names(descriptions)) descriptions["log2FoldChange"] <- paste("log2 fold change (MMSE):", resultName)
  if ("lfcSE" %in% names(descriptions)) descriptions["lfcSE"] <- paste("posterior standard deviation:", resultName)
  if ("stat" %in% names(descriptions)) descriptions["stat"] <- paste("Wald statistic:", resultName)
  if ("pvalue" %in% names(descriptions)) descriptions["pvalue"] <- paste("Wald test p-value:", resultName)
  if ("padj" %in% names(descriptions)) descriptions["padj"] <- "BH adjusted p-values"
  if ("svalue" %in% names(descriptions)) descriptions["svalue"] <- "s-value"
  S4Vectors::mcols(resultData) <- S4Vectors::DataFrame(
    type = unname(types),
    description = unname(descriptions),
    row.names = names(resultData)
  )
  DESeq2::DESeqResults(
    resultData,
    priorInfo = list(
      type = "ashr",
      package = "ashr",
      version = as.character(utils::packageVersion("ashr")),
      fitted_g = ashFit$fitted_g
    )
  )
}

lfcShrink <- function(dds, coef = NULL, type = "normal", res = NULL,
                      contrast = NULL,
                      lfcThreshold = 0, svalue = FALSE,
                      independentFiltering = TRUE, alpha = 0.1,
                      apeAdapt = TRUE, apeMethod = "nbinomCR",
                      quiet = FALSE, ...) {
  dots <- list(...)
  if (length(dots) && identical(type, "normal")) {
    argumentNames <- names(dots)
    if (is.null(argumentNames)) argumentNames <- rep("...", length(dots))
    argumentNames[!nzchar(argumentNames)] <- "..."
    stop(sprintf("unsupported lfcShrink arguments: %s", paste(argumentNames, collapse = ", ")), call. = FALSE)
  }
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      !type %in% c("normal", "ashr", "apeglm")) {
    stop("type must be 'normal', 'ashr', or 'apeglm'", call. = FALSE)
  }
  if (is.null(coef) == is.null(contrast)) {
    stop("supply exactly one of coef or contrast", call. = FALSE)
  }
  if (!is.logical(svalue) || length(svalue) != 1L || is.na(svalue)) {
    stop("svalue must be TRUE or FALSE", call. = FALSE)
  }
  if (svalue && identical(type, "normal")) {
    stop("svalue is not supported by type='normal'", call. = FALSE)
  }
  if (length(lfcThreshold) != 1L || is.na(lfcThreshold) || !is.finite(lfcThreshold) || lfcThreshold < 0) {
    stop("lfcThreshold must be a finite non-negative number", call. = FALSE)
  }
  if (!is.logical(independentFiltering) || length(independentFiltering) != 1L || is.na(independentFiltering)) {
    stop("independentFiltering must be TRUE or FALSE", call. = FALSE)
  }
  if (length(alpha) != 1L || is.na(alpha) || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must be a finite number between 0 and 1", call. = FALSE)
  }
  if (!is.logical(apeAdapt) || length(apeAdapt) != 1L || is.na(apeAdapt)) {
    stop("apeAdapt must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.character(apeMethod) || length(apeMethod) != 1L || is.na(apeMethod) || !nzchar(apeMethod)) {
    stop("apeMethod must be a non-empty string", call. = FALSE)
  }
  if (!is.logical(quiet) || length(quiet) != 1L || is.na(quiet)) {
    stop("quiet must be TRUE or FALSE", call. = FALSE)
  }
  if (identical(type, "apeglm")) {
    if (!quiet) {
      message("using 'apeglm' for LFC shrinkage")
    }
    return(.rsdeseq2_lfc_shrink_apeglm(
      dds, coef, contrast, res, lfcThreshold, svalue,
      independentFiltering, alpha, apeAdapt, apeMethod, quiet, dots
    ))
  }
  if (identical(type, "ashr")) {
    return(.rsdeseq2_lfc_shrink_ashr(
      dds, coef, contrast, res, lfcThreshold, svalue,
      independentFiltering, alpha, dots
    ))
  }
  .rsdeseq2_lfc_shrink_dataset(
    dds, coef, contrast, res, lfcThreshold, independentFiltering, alpha
  )
}
