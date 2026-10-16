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
