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
  for (package in c("SummarizedExperiment", "S4Vectors", "DESeq2")) {
    if (!requireNamespace(package, quietly = TRUE)) {
      stop(sprintf("type='apeglm' requires package '%s'", package), call. = FALSE)
    }
  }
  fit <- .rsdeseq2_fitted_view(dds)
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
    offset <- list(size_factors = sizeFactors)
  } else {
    normalizationFactors <- fit$normalizationFactors
    if (!all(dim(normalizationFactors) == dim(counts)) ||
        any(!is.finite(normalizationFactors)) || any(normalizationFactors <= 0)) {
      stop("the stored Rust fit has no compatible normalization factors", call. = FALSE)
    }
    offset <- list(factors = normalizationFactors)
  }
  weights <- fit$weights
  if (!is.null(weights) && (!all(dim(weights) == dim(counts)) || any(!is.finite(weights)) || any(weights <= 0))) {
    stop("the stored Rust fit has no compatible observation weights", call. = FALSE)
  }
  apeArguments <- utils::modifyList(
    list(
      Y = counts,
      x = design,
      log.lik = if (identical(apeMethod, "general")) list(type = "negative_binomial") else NULL,
      param = as.numeric(dispersions),
      coef = coefficient,
      mle = if (apeAdapt) list(estimates = res$log2FoldChange, standard_errors = res$lfcSE, scale = "log2") else NULL,
      threshold = if (lfcThreshold > 0) list(value = lfcThreshold, scale = "log2") else NULL,
      weights = weights,
      offset = offset,
      method = apeMethod
    ),
    dots
  )
  nativeMethod <- match.arg(apeArguments$method, c("general", "nbinomR", "nbinomCR", "nbinomC", "nbinomC*", "betabinR", "betabinCR", "betabinC", "betabinC*"))
  if (identical(nativeMethod, "general")) {
    if (is.function(apeArguments$log.lik)) apeArguments$log.lik <- .rsdeseq2_apeglm_likelihood(apeArguments$log.lik)
    if (is.null(apeArguments$log.lik)) stop("general apeglm requires a likelihood function", call. = FALSE)
  } else apeArguments$log.lik <- NULL
  apeFit <- .rsdeseq2_formula_request(list(operation = "apeglm_shrinkage", arguments = apeArguments))
  if (length(apeFit$posterior_mode) != nrow(counts)) {
    stop("apeglm returned a coefficient matrix with an incompatible row count", call. = FALSE)
  }
  if (!quiet && any(apeFit$convergence[!is.na(apeFit$convergence)] != 0)) {
    message("some rows did not converge in finding the MAP")
  }
  res$log2FoldChange <- apeFit$posterior_mode
  names(res$log2FoldChange) <- rownames(res)
  res$lfcSE <- apeFit$posterior_sd
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
      if (is.null(apeFit$threshold_svalue)) stop("thresholded apeglm output requires a non-NULL threshold", call. = FALSE)
      apeFit$threshold_svalue
    } else {
      apeFit$svalue
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
      package = "rsdeseq2",
      version = as.character(utils::packageVersion("rsdeseq2")),
      engine = "Rust",
      prior.control = apeFit$prior_control
    )
  )
}
