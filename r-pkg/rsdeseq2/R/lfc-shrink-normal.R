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
  fit <- .rsdeseq2_fitted_view(dds)
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
    if (!is.null(fit$allZero)) fit$allZero <- fit$allZero[geneOrder]
    if (!is.null(fit$normalizationFactors)) fit$normalizationFactors <- fit$normalizationFactors[geneOrder, , drop = FALSE]
    if (!is.null(fit$weights)) fit$weights <- fit$weights[geneOrder, , drop = FALSE]
  }
  useGene <- .Call("rsdeseq2_normal_shrinkage_rows", counts,
    as.numeric(fit$dispersions), as.numeric(fit$dispFit), fit$allZero,
    PACKAGE="rsdeseq2")
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

