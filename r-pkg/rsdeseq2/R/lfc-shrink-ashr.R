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
  fit <- .rsdeseq2_fitted_view(dds)
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

