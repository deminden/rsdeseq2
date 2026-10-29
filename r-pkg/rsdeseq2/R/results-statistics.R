# Wald/LRT statistics, thresholds and contrast standard errors.

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
  .Call("rsdeseq2_wald_statistics", as.numeric(beta), as.numeric(betaSE),
    options$lfcThreshold,
    as.integer(match(options$altHypothesis,
      c("greaterAbs", "greaterAbsUPSHOT", "greaterAbs2014", "lessAbs", "greater", "less")) - 1L),
    options$useT, options$degreesOfFreedom, PACKAGE = "rsdeseq2")
}

.rsdeseq2_default_cooks_cutoff <- function(design) {
  .Call("rsdeseq2_default_cooks_cutoff", as.integer(nrow(design)),
    as.integer(ncol(design)), PACKAGE = "rsdeseq2")
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
