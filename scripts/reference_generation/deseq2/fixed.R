full_results <- tryCatch({
  dds_full <- DESeq2::DESeq(dds_ratio, quiet = TRUE)
  DESeq2::results(dds_full)
}, error = function(error) {
  writeLines(conditionMessage(error), file.path(data_dir, "results_wald_ratio_error.txt"))
  NULL
})

if (!is.null(full_results)) {
  write_tsv(
    data.frame(gene = rownames(full_results), as.data.frame(full_results), check.names = FALSE),
    file.path(data_dir, "results_wald_ratio.tsv")
  )
}

fixed_reference <- tryCatch({
  dds_fixed <- dds_ratio
  dds_fixed <- DESeq2::`dispersions<-`(dds_fixed, value = fixed_dispersions)
  mcols_fixed <- S4Vectors::mcols(dds_fixed)
  mcols_fixed$dispersion <- fixed_dispersions
  mcols_fixed$allZero <- rowSums(counts) == 0
  dds_fixed <- S4Vectors::`mcols<-`(dds_fixed, value = mcols_fixed)

  full_fit <- DESeq2:::fitNbinomGLMs(
    dds_fixed,
    modelMatrix = full_design,
    renameCols = FALSE,
    betaTol = 1e-8,
    maxit = 100,
    useOptim = FALSE,
    useQR = FALSE,
    warnNonposVar = FALSE,
    minmu = 0.5
  )
  reduced_fit <- DESeq2:::fitNbinomGLMs(
    dds_fixed,
    modelMatrix = reduced_design,
    renameCols = FALSE,
    betaTol = 1e-8,
    maxit = 100,
    useOptim = FALSE,
    useQR = FALSE,
    warnNonposVar = FALSE,
    minmu = 0.5
  )

  coefficient <- 2L
  wald_stat <- full_fit$betaMatrix[, coefficient] / full_fit$betaSE[, coefficient]
  wald_pvalue <- 2 * stats::pnorm(abs(wald_stat), lower.tail = FALSE)
  wald_t_df <- ncol(counts) - ncol(full_design)
  wald_t_pvalue <- 2 * stats::pt(abs(wald_stat), df = wald_t_df, lower.tail = FALSE)
  lrt_stat <- 2 * (full_fit$logLike - reduced_fit$logLike)
  lrt_df <- ncol(full_design) - ncol(reduced_design)
  lrt_pvalue <- stats::pchisq(lrt_stat, df = lrt_df, lower.tail = FALSE)

  assays_fixed <- SummarizedExperiment::assays(dds_fixed, withDimnames = FALSE)
  assays_fixed[["mu"]] <- full_fit$mu
  dds_fixed <- SummarizedExperiment::`assays<-`(
    dds_fixed,
    withDimnames = FALSE,
    value = assays_fixed
  )
  cooks <- DESeq2:::calculateCooksDistance(dds_fixed, full_fit$hat_diagonals, full_design)

  list(
    full_fit = full_fit,
    reduced_fit = reduced_fit,
    wald_stat = wald_stat,
    wald_pvalue = wald_pvalue,
    wald_t_df = wald_t_df,
    wald_t_pvalue = wald_t_pvalue,
    lrt_stat = lrt_stat,
    lrt_df = lrt_df,
    lrt_pvalue = lrt_pvalue,
    cooks = cooks
  )
}, error = function(error) {
  writeLines(conditionMessage(error), file.path(data_dir, "fixed_glm_reference_error.txt"))
  NULL
})

if (!is.null(fixed_reference)) {
  full_fit <- fixed_reference$full_fit
  reduced_fit <- fixed_reference$reduced_fit
  dds_beta_prior <- dds_ratio
  dds_beta_prior <- DESeq2::`dispersions<-`(dds_beta_prior, value = fixed_dispersions)
  attr(dds_beta_prior, "modelMatrixType") <- "standard"
  beta_prior_mcols <- S4Vectors::mcols(dds_beta_prior)
  beta_prior_mcols$dispersion <- fixed_dispersions
  beta_prior_mcols$allZero <- rowSums(counts) == 0
  beta_prior_mcols$baseMean <- base_mean
  beta_prior_mcols$dispFit <- fixed_dispersions
  for (coefficient in seq_len(ncol(full_fit$betaMatrix))) {
    beta_prior_mcols[[paste0("MLE_", colnames(full_fit$betaMatrix)[coefficient])]] <-
      full_fit$betaMatrix[, coefficient]
  }
  dds_beta_prior <- S4Vectors::`mcols<-`(dds_beta_prior, value = beta_prior_mcols)
  beta_prior_weighted <- DESeq2:::estimateBetaPriorVar(
    dds_beta_prior,
    betaPriorMethod = "weighted",
    upperQuantile = 0.05,
    modelMatrix = full_design
  )
  beta_prior_quantile <- DESeq2:::estimateBetaPriorVar(
    dds_beta_prior,
    betaPriorMethod = "quantile",
    upperQuantile = 0.05,
    modelMatrix = full_design
  )
  beta_prior_refit <- DESeq2:::fitNbinomGLMs(
    dds_beta_prior,
    modelMatrix = full_design,
    lambda = 1 / beta_prior_weighted,
    renameCols = FALSE,
    betaTol = 1e-8,
    maxit = 100,
    useOptim = FALSE,
    useQR = FALSE,
    warnNonposVar = FALSE,
    minmu = 0.5
  )

  write_tsv(
    data.frame(
      coefficient = names(beta_prior_weighted),
      weighted = as.numeric(beta_prior_weighted),
      quantile = as.numeric(beta_prior_quantile),
      check.names = FALSE
    ),
    file.path(data_dir, "beta_prior_variance_reference.tsv")
  )
  write_tsv(
    data.frame(
      gene = rownames(counts),
      beta_intercept = beta_prior_refit$betaMatrix[, 1],
      beta_conditionB = beta_prior_refit$betaMatrix[, 2],
      beta_se_intercept = beta_prior_refit$betaSE[, 1],
      beta_se_conditionB = beta_prior_refit$betaSE[, 2],
      log_like = beta_prior_refit$logLike,
      converged = beta_prior_refit$betaConv,
      iterations = beta_prior_refit$betaIter,
      check.names = FALSE
    ),
    file.path(data_dir, "beta_prior_refit_reference.tsv")
  )
  dimnames(beta_prior_refit$mu) <- dimnames(counts)
  dimnames(beta_prior_refit$hat_diagonals) <- dimnames(counts)
  write_matrix_tsv(
    beta_prior_refit$mu,
    "gene",
    file.path(data_dir, "beta_prior_refit_mu.tsv")
  )
  write_matrix_tsv(
    beta_prior_refit$hat_diagonals,
    "gene",
    file.path(data_dir, "beta_prior_refit_hat.tsv")
  )

  write_tsv(
    data.frame(
      gene = rownames(counts),
      dispersion = fixed_dispersions,
      beta_intercept = full_fit$betaMatrix[, 1],
      beta_conditionB = full_fit$betaMatrix[, 2],
      beta_se_intercept = full_fit$betaSE[, 1],
      beta_se_conditionB = full_fit$betaSE[, 2],
      stat_conditionB = fixed_reference$wald_stat,
      pvalue_conditionB = fixed_reference$wald_pvalue,
      log_like = full_fit$logLike,
      converged = full_fit$betaConv,
      iterations = full_fit$betaIter,
      check.names = FALSE
    ),
    file.path(data_dir, "fixed_wald_reference.tsv")
  )

  write_tsv(
    data.frame(
      gene = rownames(counts),
      dispersion = fixed_dispersions,
      stat_conditionB = fixed_reference$wald_stat,
      df_conditionB = fixed_reference$wald_t_df,
      pvalue_conditionB = fixed_reference$wald_t_pvalue,
      check.names = FALSE
    ),
    file.path(data_dir, "fixed_wald_t_reference.tsv")
  )

  write_tsv(
    data.frame(
      gene = rownames(counts),
      dispersion = fixed_dispersions,
      beta_intercept = full_fit$betaMatrix[, 1],
      beta_conditionB = full_fit$betaMatrix[, 2],
      beta_se_intercept = full_fit$betaSE[, 1],
      beta_se_conditionB = full_fit$betaSE[, 2],
      log_like_full = full_fit$logLike,
      log_like_reduced = reduced_fit$logLike,
      lrt_stat = fixed_reference$lrt_stat,
      pvalue = fixed_reference$lrt_pvalue,
      df = fixed_reference$lrt_df,
      full_converged = full_fit$betaConv,
      reduced_converged = reduced_fit$betaConv,
      check.names = FALSE
    ),
    file.path(data_dir, "fixed_lrt_reference.tsv")
  )

  write_matrix_tsv(full_fit$mu, "gene", file.path(data_dir, "fixed_mu_full.tsv"))
  write_matrix_tsv(full_fit$hat_diagonals, "gene", file.path(data_dir, "fixed_hat_full.tsv"))
  dimnames(reduced_fit$mu) <- dimnames(counts)
  dimnames(reduced_fit$hat_diagonals) <- dimnames(counts)
  write_matrix_tsv(
    reduced_fit$mu,
    "gene",
    file.path(data_dir, "fixed_mu_reduced.tsv")
  )
  write_matrix_tsv(
    reduced_fit$hat_diagonals,
    "gene",
    file.path(data_dir, "fixed_hat_reduced.tsv")
  )
  write_matrix_tsv(fixed_reference$cooks, "gene", file.path(data_dir, "fixed_cooks_full.tsv"))
}
