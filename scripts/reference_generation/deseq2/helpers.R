write_result_reference <- function(reference, filename, test) {
  if (identical(test, "Wald")) {
    table <- data.frame(
      gene = rownames(counts),
      baseMean = reference$row_meta$baseMean,
      log2FoldChange = reference$beta_conditionB,
      lfcSE = reference$beta_se_conditionB,
      stat = reference$stat_conditionB,
      pvalue = reference$pvalue_conditionB,
      padj = reference$padj_conditionB,
      dispersion = reference$row_meta$dispersion,
      betaConv = reference$beta_converged,
      check.names = FALSE
    )
  } else if (identical(test, "LRT")) {
    table <- data.frame(
      gene = rownames(counts),
      baseMean = reference$row_meta$baseMean,
      log2FoldChange = reference$beta_conditionB,
      lfcSE = reference$beta_se_conditionB,
      stat = reference$lrt_stat,
      pvalue = reference$lrt_pvalue,
      padj = reference$lrt_padj,
      dispersion = reference$row_meta$dispersion,
      betaConv = reference$beta_converged,
      check.names = FALSE
    )
  } else {
    stop("unsupported result reference test: ", test)
  }
  write_tsv(table, file.path(data_dir, filename))
}

build_glm_mu_reference <- function(dds_input, fit_type, use_cr) {
  dds_gene <- DESeq2:::estimateDispersionsGeneEst(
    dds_input,
    modelMatrix = full_design,
    linearMu = FALSE,
    niter = 2,
    useCR = use_cr,
    quiet = TRUE
  )
  dds_fit <- DESeq2:::estimateDispersionsFit(
    dds_gene,
    fitType = fit_type,
    quiet = TRUE
  )
  disp_prior_var <- DESeq2:::estimateDispersionsPriorVar(
    dds_fit,
    modelMatrix = full_design
  )
  dds_map <- DESeq2:::estimateDispersionsMAP(
    dds_fit,
    dispPriorVar = disp_prior_var,
    modelMatrix = full_design,
    useCR = use_cr,
    quiet = TRUE
  )

  row_meta <- S4Vectors::mcols(dds_map)
  weights_fail <- if ("weightsFail" %in% names(row_meta)) {
    row_meta$weightsFail
  } else {
    rep(FALSE, nrow(counts))
  }
  all_zero <- row_meta$allZero
  object_nz <- dds_map[!all_zero, , drop = FALSE]
  full_fit_nz <- suppressWarnings(DESeq2:::fitNbinomGLMs(
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
  reduced_fit_nz <- suppressWarnings(DESeq2:::fitNbinomGLMs(
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

  expand_vector <- function(value, missing_rows) {
    expanded <- rep(NA_real_, length(missing_rows))
    expanded[!missing_rows] <- value
    expanded
  }
  expand_logical <- function(value, missing_rows) {
    expanded <- rep(NA, length(missing_rows))
    expanded[!missing_rows] <- value
    expanded
  }
  expand_matrix <- function(value, missing_rows, dim_names) {
    expanded <- DESeq2:::buildMatrixWithNARows(value, missing_rows)
    dimnames(expanded) <- dim_names
    expanded
  }

  coefficient <- 2L
  wald_stat_nz <- full_fit_nz$betaMatrix[, coefficient] / full_fit_nz$betaSE[, coefficient]
  wald_pvalue_nz <- 2 * stats::pnorm(abs(wald_stat_nz), lower.tail = FALSE)
  wald_padj_nz <- stats::p.adjust(wald_pvalue_nz, method = "BH")
  lrt_stat_nz <- 2 * (full_fit_nz$logLike - reduced_fit_nz$logLike)
  lrt_df <- ncol(full_design) - ncol(reduced_design)
  lrt_pvalue_nz <- stats::pchisq(lrt_stat_nz, df = lrt_df, lower.tail = FALSE)
  lrt_padj_nz <- stats::p.adjust(lrt_pvalue_nz, method = "BH")
  list(
    row_meta = row_meta,
    weights_fail = weights_fail,
    disp_prior_var = disp_prior_var,
    var_log_disp_estimates = attr(DESeq2::dispersionFunction(dds_map), "varLogDispEsts"),
    dispersion_mu = SummarizedExperiment::assays(dds_map)[["mu"]],
    beta_intercept = expand_vector(full_fit_nz$betaMatrix[, 1], all_zero),
    beta_conditionB = expand_vector(full_fit_nz$betaMatrix[, 2], all_zero),
    beta_se_intercept = expand_vector(full_fit_nz$betaSE[, 1], all_zero),
    beta_se_conditionB = expand_vector(full_fit_nz$betaSE[, 2], all_zero),
    stat_conditionB = expand_vector(wald_stat_nz, all_zero),
    pvalue_conditionB = expand_vector(wald_pvalue_nz, all_zero),
    padj_conditionB = expand_vector(wald_padj_nz, all_zero),
    lrt_stat = expand_vector(lrt_stat_nz, all_zero),
    lrt_df = lrt_df,
    lrt_pvalue = expand_vector(lrt_pvalue_nz, all_zero),
    lrt_padj = expand_vector(lrt_padj_nz, all_zero),
    log_like = expand_vector(full_fit_nz$logLike, all_zero),
    reduced_log_like = expand_vector(reduced_fit_nz$logLike, all_zero),
    beta_converged = expand_logical(full_fit_nz$betaConv, all_zero),
    beta_iterations = expand_vector(full_fit_nz$betaIter, all_zero),
    reduced_beta_converged = expand_logical(reduced_fit_nz$betaConv, all_zero),
    reduced_beta_iterations = expand_vector(reduced_fit_nz$betaIter, all_zero),
    wald_mu = expand_matrix(full_fit_nz$mu, all_zero, dimnames(counts)),
    wald_hat = expand_matrix(full_fit_nz$hat_diagonals, all_zero, dimnames(counts))
  )
}
