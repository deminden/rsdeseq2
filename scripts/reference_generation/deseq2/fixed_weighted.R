weighted_fixed_reference <- tryCatch({
  dds_weighted_fixed <- DESeq2:::getBaseMeansAndVariances(dds_weighted)
  dds_weighted_fixed <- DESeq2::`dispersions<-`(dds_weighted_fixed, value = fixed_dispersions)
  mcols_weighted_fixed <- S4Vectors::mcols(dds_weighted_fixed)
  mcols_weighted_fixed$dispersion <- fixed_dispersions
  dds_weighted_fixed <- S4Vectors::`mcols<-`(dds_weighted_fixed, value = mcols_weighted_fixed)
  checked_weighted_fixed <- suppressWarnings(DESeq2:::getAndCheckWeights(
    dds_weighted_fixed,
    full_design
  )$object)
  fixed_all_zero <- S4Vectors::mcols(checked_weighted_fixed)$allZero
  object_nz <- checked_weighted_fixed[!fixed_all_zero, , drop = FALSE]

  full_fit <- suppressWarnings(DESeq2:::fitNbinomGLMs(
    object_nz,
    modelMatrix = full_design,
    renameCols = FALSE,
    betaTol = 1e-8,
    maxit = 100,
    useOptim = FALSE,
    useQR = FALSE,
    warnNonposVar = FALSE,
    minmu = 0.5
  ))
  reduced_fit <- suppressWarnings(DESeq2:::fitNbinomGLMs(
    object_nz,
    modelMatrix = reduced_design,
    renameCols = FALSE,
    betaTol = 1e-8,
    maxit = 100,
    useOptim = FALSE,
    useQR = FALSE,
    warnNonposVar = FALSE,
    minmu = 0.5
  ))

  coefficient <- 2L
  wald_stat <- full_fit$betaMatrix[, coefficient] / full_fit$betaSE[, coefficient]
  wald_pvalue <- 2 * stats::pnorm(abs(wald_stat), lower.tail = FALSE)
  lrt_stat <- 2 * (full_fit$logLike - reduced_fit$logLike)
  lrt_df <- ncol(full_design) - ncol(reduced_design)
  lrt_pvalue <- stats::pchisq(lrt_stat, df = lrt_df, lower.tail = FALSE)

  expand_numeric <- function(value) {
    expanded <- rep(NA_real_, length(fixed_all_zero))
    expanded[!fixed_all_zero] <- value
    expanded
  }
  expand_integer <- function(value) {
    expanded <- integer(length(fixed_all_zero))
    expanded[!fixed_all_zero] <- value
    expanded
  }
  expand_logical <- function(value) {
    expanded <- rep(FALSE, length(fixed_all_zero))
    expanded[!fixed_all_zero] <- value
    expanded
  }
  expand_matrix <- function(value) {
    expanded <- DESeq2:::buildMatrixWithNARows(value, fixed_all_zero)
    dimnames(expanded) <- dimnames(counts)
    expanded
  }

  weighted_meta <- S4Vectors::mcols(checked_weighted_fixed)
  weights_fail <- if ("weightsFail" %in% names(weighted_meta)) {
    weighted_meta$weightsFail
  } else {
    rep(FALSE, nrow(counts))
  }

  list(
    full_fit = full_fit,
    reduced_fit = reduced_fit,
    wald_stat = expand_numeric(wald_stat),
    wald_pvalue = expand_numeric(wald_pvalue),
    lrt_stat = expand_numeric(lrt_stat),
    lrt_df = lrt_df,
    lrt_pvalue = expand_numeric(lrt_pvalue),
    row_meta = weighted_meta,
    weights_fail = weights_fail,
    beta_intercept = expand_numeric(full_fit$betaMatrix[, 1]),
    beta_conditionB = expand_numeric(full_fit$betaMatrix[, 2]),
    beta_se_intercept = expand_numeric(full_fit$betaSE[, 1]),
    beta_se_conditionB = expand_numeric(full_fit$betaSE[, 2]),
    log_like = expand_numeric(full_fit$logLike),
    converged = expand_logical(full_fit$betaConv),
    iterations = expand_integer(full_fit$betaIter),
    reduced_log_like = expand_numeric(reduced_fit$logLike),
    reduced_converged = expand_logical(reduced_fit$betaConv),
    reduced_iterations = expand_integer(reduced_fit$betaIter),
    reduced_mu = expand_matrix(reduced_fit$mu),
    reduced_hat = expand_matrix(reduced_fit$hat_diagonals)
  )
}, error = function(error) {
  writeLines(conditionMessage(error), file.path(data_dir, "fixed_weighted_glm_reference_error.txt"))
  NULL
})

if (!is.null(weighted_fixed_reference)) {
  write_tsv(
    data.frame(
      gene = rownames(counts),
      dispersion = fixed_dispersions,
      baseMean = weighted_fixed_reference$row_meta$baseMean,
      allZero = weighted_fixed_reference$row_meta$allZero,
      weightsFail = weighted_fixed_reference$weights_fail,
      beta_intercept = weighted_fixed_reference$beta_intercept,
      beta_conditionB = weighted_fixed_reference$beta_conditionB,
      beta_se_intercept = weighted_fixed_reference$beta_se_intercept,
      beta_se_conditionB = weighted_fixed_reference$beta_se_conditionB,
      stat_conditionB = weighted_fixed_reference$wald_stat,
      pvalue_conditionB = weighted_fixed_reference$wald_pvalue,
      log_like = weighted_fixed_reference$log_like,
      converged = weighted_fixed_reference$converged,
      iterations = weighted_fixed_reference$iterations,
      check.names = FALSE
    ),
    file.path(data_dir, "fixed_wald_weighted_reference.tsv")
  )

  write_tsv(
    data.frame(
      gene = rownames(counts),
      dispersion = fixed_dispersions,
      baseMean = weighted_fixed_reference$row_meta$baseMean,
      allZero = weighted_fixed_reference$row_meta$allZero,
      weightsFail = weighted_fixed_reference$weights_fail,
      beta_intercept = weighted_fixed_reference$beta_intercept,
      beta_conditionB = weighted_fixed_reference$beta_conditionB,
      beta_se_intercept = weighted_fixed_reference$beta_se_intercept,
      beta_se_conditionB = weighted_fixed_reference$beta_se_conditionB,
      log_like_full = weighted_fixed_reference$log_like,
      log_like_reduced = weighted_fixed_reference$reduced_log_like,
      lrt_stat = weighted_fixed_reference$lrt_stat,
      pvalue = weighted_fixed_reference$lrt_pvalue,
      df = weighted_fixed_reference$lrt_df,
      full_converged = weighted_fixed_reference$converged,
      reduced_converged = weighted_fixed_reference$reduced_converged,
      full_iterations = weighted_fixed_reference$iterations,
      reduced_iterations = weighted_fixed_reference$reduced_iterations,
      check.names = FALSE
    ),
    file.path(data_dir, "fixed_lrt_weighted_reference.tsv")
  )
  write_matrix_tsv(
    weighted_fixed_reference$reduced_mu,
    "gene",
    file.path(data_dir, "fixed_weighted_mu_reduced.tsv")
  )
  write_matrix_tsv(
    weighted_fixed_reference$reduced_hat,
    "gene",
    file.path(data_dir, "fixed_weighted_hat_reduced.tsv")
  )
}
