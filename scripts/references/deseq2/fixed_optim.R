fixed_force_optim_reference <- tryCatch({
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
    useOptim = TRUE,
    useQR = FALSE,
    forceOptim = TRUE,
    warnNonposVar = FALSE,
    minmu = 0.5
  )

  coefficient <- 2L
  wald_stat <- full_fit$betaMatrix[, coefficient] / full_fit$betaSE[, coefficient]
  wald_pvalue <- 2 * stats::pnorm(abs(wald_stat), lower.tail = FALSE)

  list(
    full_fit = full_fit,
    wald_stat = wald_stat,
    wald_pvalue = wald_pvalue
  )
}, error = function(error) {
  writeLines(conditionMessage(error), file.path(data_dir, "fixed_force_optim_reference_error.txt"))
  NULL
})

if (!is.null(fixed_force_optim_reference)) {
  full_fit <- fixed_force_optim_reference$full_fit
  write_tsv(
    data.frame(
      gene = rownames(counts),
      dispersion = fixed_dispersions,
      beta_intercept = full_fit$betaMatrix[, 1],
      beta_conditionB = full_fit$betaMatrix[, 2],
      beta_se_intercept = full_fit$betaSE[, 1],
      beta_se_conditionB = full_fit$betaSE[, 2],
      stat_conditionB = fixed_force_optim_reference$wald_stat,
      pvalue_conditionB = fixed_force_optim_reference$wald_pvalue,
      log_like = full_fit$logLike,
      converged = full_fit$betaConv,
      iterations = full_fit$betaIter,
      check.names = FALSE
    ),
    file.path(data_dir, "fixed_force_optim_wald_reference.tsv")
  )
  write_matrix_tsv(full_fit$mu, "gene", file.path(data_dir, "fixed_force_optim_mu_full.tsv"))
  dimnames(full_fit$hat_diagonals) <- dimnames(counts)
  write_matrix_tsv(
    full_fit$hat_diagonals,
    "gene",
    file.path(data_dir, "fixed_force_optim_hat_full.tsv")
  )
}
