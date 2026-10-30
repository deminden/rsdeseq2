# All numerical inference and projection runs in the standalone Rust engine.

.rsdeseq2_wald_from_primitive_fit <- function(fit, numericContrast, options, allZero = NULL) {
  .Call("rsdeseq2_wald_contrast_zero", fit$beta, fit$betaSE, fit$betaCovariance,
    as.numeric(numericContrast), options$lfcThreshold,
    .rsdeseq2_wald_alternative_code(options$altHypothesis),
    options$useT, options$degreesOfFreedom,
    .rsdeseq2_contrast_zero_context(fit, allZero), PACKAGE = "rsdeseq2")
}

.rsdeseq2_lrt_from_primitive_fit <- function(fit, numericContrast, allZero = NULL) {
  .Call("rsdeseq2_project_contrast_zero", fit$beta, fit$betaSE, fit$betaCovariance,
    as.numeric(numericContrast), .rsdeseq2_contrast_zero_context(fit, allZero), PACKAGE = "rsdeseq2")
}

.rsdeseq2_wald_alternative_code <- function(alternative) {
  as.integer(match(alternative,
    c("greaterAbs", "greaterAbsUPSHOT", "greaterAbs2014", "lessAbs", "greater", "less")) - 1L)
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
    .rsdeseq2_wald_alternative_code(options$altHypothesis),
    options$useT, options$degreesOfFreedom, PACKAGE = "rsdeseq2")
}

.rsdeseq2_default_cooks_cutoff <- function(design) {
  .Call("rsdeseq2_default_cooks_cutoff", as.integer(nrow(design)),
    as.integer(ncol(design)), PACKAGE = "rsdeseq2")
}

.rsdeseq2_evaluate_result_batch <- function(object,requests) {
  .Call("rsdeseq2_contrast_batch",object$beta,object$betaSE,object$betaCovariance,
    requests,PACKAGE="rsdeseq2")
}
