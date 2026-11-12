native_weighted_glm_mu_reference <- tryCatch({
  dds_weighted_gene <- DESeq2:::estimateDispersionsGeneEst(
    dds_weighted,
    modelMatrix = full_design,
    linearMu = FALSE,
    niter = 2,
    useCR = FALSE,
    quiet = TRUE
  )
  dds_weighted_fit <- DESeq2:::estimateDispersionsFit(
    dds_weighted_gene,
    fitType = "mean",
    quiet = TRUE
  )
  disp_prior_var <- DESeq2:::estimateDispersionsPriorVar(
    dds_weighted_fit,
    modelMatrix = full_design
  )
  dds_weighted_map <- DESeq2:::estimateDispersionsMAP(
    dds_weighted_fit,
    dispPriorVar = disp_prior_var,
    modelMatrix = full_design,
    useCR = FALSE,
    quiet = TRUE
  )

  row_meta <- S4Vectors::mcols(dds_weighted_map)
  weights_fail <- if ("weightsFail" %in% names(row_meta)) {
    row_meta$weightsFail
  } else {
    rep(FALSE, nrow(counts))
  }
  all_zero <- row_meta$allZero
  object_nz <- dds_weighted_map[!all_zero, , drop = FALSE]
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
    var_log_disp_estimates = attr(DESeq2::dispersionFunction(dds_weighted_map), "varLogDispEsts"),
    dispersion_mu = SummarizedExperiment::assays(dds_weighted_map)[["mu"]],
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
}, error = function(error) {
  writeLines(conditionMessage(error), file.path(data_dir, "native_weighted_glm_mu_reference_error.txt"))
  NULL
})

if (!is.null(native_weighted_glm_mu_reference)) {
  native_weighted_rows <- data.frame(
    gene = rownames(counts),
    baseMean = native_weighted_glm_mu_reference$row_meta$baseMean,
    baseVar = native_weighted_glm_mu_reference$row_meta$baseVar,
    allZero = native_weighted_glm_mu_reference$row_meta$allZero,
    weightsFail = native_weighted_glm_mu_reference$weights_fail,
    dispGeneEst = native_weighted_glm_mu_reference$row_meta$dispGeneEst,
    dispGeneIter = native_weighted_glm_mu_reference$row_meta$dispGeneIter,
    dispFit = native_weighted_glm_mu_reference$row_meta$dispFit,
    dispPriorVar = native_weighted_glm_mu_reference$disp_prior_var,
    varLogDispEsts = native_weighted_glm_mu_reference$var_log_disp_estimates,
    dispMAP = native_weighted_glm_mu_reference$row_meta$dispMAP,
    dispersion = native_weighted_glm_mu_reference$row_meta$dispersion,
    dispIter = native_weighted_glm_mu_reference$row_meta$dispIter,
    dispOutlier = native_weighted_glm_mu_reference$row_meta$dispOutlier,
    beta_intercept = native_weighted_glm_mu_reference$beta_intercept,
    beta_conditionB = native_weighted_glm_mu_reference$beta_conditionB,
    beta_se_intercept = native_weighted_glm_mu_reference$beta_se_intercept,
    beta_se_conditionB = native_weighted_glm_mu_reference$beta_se_conditionB,
    stat_conditionB = native_weighted_glm_mu_reference$stat_conditionB,
    pvalue_conditionB = native_weighted_glm_mu_reference$pvalue_conditionB,
    padj_conditionB = native_weighted_glm_mu_reference$padj_conditionB,
    log_like = native_weighted_glm_mu_reference$log_like,
    beta_converged = native_weighted_glm_mu_reference$beta_converged,
    beta_iterations = native_weighted_glm_mu_reference$beta_iterations,
    check.names = FALSE
  )
  write_tsv(
    native_weighted_rows,
    file.path(data_dir, "native_weighted_glm_mu_reference.tsv")
  )
  write_result_reference(
    native_weighted_glm_mu_reference,
    "native_weighted_glm_mu_results_wald.tsv",
    "Wald"
  )
  native_weighted_lrt_rows <- data.frame(
    gene = rownames(counts),
    baseMean = native_weighted_glm_mu_reference$row_meta$baseMean,
    baseVar = native_weighted_glm_mu_reference$row_meta$baseVar,
    allZero = native_weighted_glm_mu_reference$row_meta$allZero,
    weightsFail = native_weighted_glm_mu_reference$weights_fail,
    dispGeneEst = native_weighted_glm_mu_reference$row_meta$dispGeneEst,
    dispGeneIter = native_weighted_glm_mu_reference$row_meta$dispGeneIter,
    dispFit = native_weighted_glm_mu_reference$row_meta$dispFit,
    dispMAP = native_weighted_glm_mu_reference$row_meta$dispMAP,
    dispersion = native_weighted_glm_mu_reference$row_meta$dispersion,
    beta_intercept = native_weighted_glm_mu_reference$beta_intercept,
    beta_conditionB = native_weighted_glm_mu_reference$beta_conditionB,
    beta_se_intercept = native_weighted_glm_mu_reference$beta_se_intercept,
    beta_se_conditionB = native_weighted_glm_mu_reference$beta_se_conditionB,
    log_like_full = native_weighted_glm_mu_reference$log_like,
    log_like_reduced = native_weighted_glm_mu_reference$reduced_log_like,
    lrt_stat = native_weighted_glm_mu_reference$lrt_stat,
    pvalue = native_weighted_glm_mu_reference$lrt_pvalue,
    padj = native_weighted_glm_mu_reference$lrt_padj,
    df = native_weighted_glm_mu_reference$lrt_df,
    full_converged = native_weighted_glm_mu_reference$beta_converged,
    reduced_converged = native_weighted_glm_mu_reference$reduced_beta_converged,
    full_iterations = native_weighted_glm_mu_reference$beta_iterations,
    reduced_iterations = native_weighted_glm_mu_reference$reduced_beta_iterations,
    check.names = FALSE
  )
  write_tsv(
    native_weighted_lrt_rows,
    file.path(data_dir, "native_weighted_glm_mu_lrt_reference.tsv")
  )
  write_result_reference(
    native_weighted_glm_mu_reference,
    "native_weighted_glm_mu_results_lrt.tsv",
    "LRT"
  )
  write_matrix_tsv(
    native_weighted_glm_mu_reference$dispersion_mu,
    "gene",
    file.path(data_dir, "native_weighted_glm_mu_dispersion_mu.tsv")
  )
  write_matrix_tsv(
    native_weighted_glm_mu_reference$wald_mu,
    "gene",
    file.path(data_dir, "native_weighted_glm_mu_wald_mu.tsv")
  )
  write_matrix_tsv(
    native_weighted_glm_mu_reference$wald_hat,
    "gene",
    file.path(data_dir, "native_weighted_glm_mu_wald_hat.tsv")
  )
}

native_weighted_glm_mu_mean_cr_map_reference <- tryCatch({
  dds_weighted_cr_gene <- DESeq2:::estimateDispersionsGeneEst(
    dds_weighted,
    modelMatrix = full_design,
    linearMu = FALSE,
    niter = 2,
    useCR = TRUE,
    quiet = TRUE
  )
  dds_weighted_cr_fit <- DESeq2:::estimateDispersionsFit(
    dds_weighted_cr_gene,
    fitType = "mean",
    quiet = TRUE
  )
  disp_prior_var <- DESeq2:::estimateDispersionsPriorVar(
    dds_weighted_cr_fit,
    modelMatrix = full_design
  )
  dds_weighted_cr_map <- DESeq2:::estimateDispersionsMAP(
    dds_weighted_cr_fit,
    dispPriorVar = disp_prior_var,
    modelMatrix = full_design,
    useCR = TRUE,
    quiet = TRUE
  )

  row_meta <- S4Vectors::mcols(dds_weighted_cr_map)
  weights_fail <- if ("weightsFail" %in% names(row_meta)) {
    row_meta$weightsFail
  } else {
    rep(FALSE, nrow(counts))
  }
  all_zero <- row_meta$allZero
  object_nz <- dds_weighted_cr_map[!all_zero, , drop = FALSE]
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
    var_log_disp_estimates = attr(DESeq2::dispersionFunction(dds_weighted_cr_map), "varLogDispEsts"),
    dispersion_mu = SummarizedExperiment::assays(dds_weighted_cr_map)[["mu"]],
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
}, error = function(error) {
  writeLines(conditionMessage(error), file.path(data_dir, "native_weighted_glm_mu_mean_cr_map_reference_error.txt"))
  NULL
})

if (!is.null(native_weighted_glm_mu_mean_cr_map_reference)) {
  native_weighted_glm_mu_mean_cr_map_rows <- data.frame(
    gene = rownames(counts),
    baseMean = native_weighted_glm_mu_mean_cr_map_reference$row_meta$baseMean,
    baseVar = native_weighted_glm_mu_mean_cr_map_reference$row_meta$baseVar,
    allZero = native_weighted_glm_mu_mean_cr_map_reference$row_meta$allZero,
    weightsFail = native_weighted_glm_mu_mean_cr_map_reference$weights_fail,
    dispGeneEst = native_weighted_glm_mu_mean_cr_map_reference$row_meta$dispGeneEst,
    dispGeneIter = native_weighted_glm_mu_mean_cr_map_reference$row_meta$dispGeneIter,
    dispFit = native_weighted_glm_mu_mean_cr_map_reference$row_meta$dispFit,
    dispPriorVar = native_weighted_glm_mu_mean_cr_map_reference$disp_prior_var,
    varLogDispEsts = native_weighted_glm_mu_mean_cr_map_reference$var_log_disp_estimates,
    dispMAP = native_weighted_glm_mu_mean_cr_map_reference$row_meta$dispMAP,
    dispersion = native_weighted_glm_mu_mean_cr_map_reference$row_meta$dispersion,
    dispIter = native_weighted_glm_mu_mean_cr_map_reference$row_meta$dispIter,
    dispOutlier = native_weighted_glm_mu_mean_cr_map_reference$row_meta$dispOutlier,
    beta_intercept = native_weighted_glm_mu_mean_cr_map_reference$beta_intercept,
    beta_conditionB = native_weighted_glm_mu_mean_cr_map_reference$beta_conditionB,
    beta_se_intercept = native_weighted_glm_mu_mean_cr_map_reference$beta_se_intercept,
    beta_se_conditionB = native_weighted_glm_mu_mean_cr_map_reference$beta_se_conditionB,
    stat_conditionB = native_weighted_glm_mu_mean_cr_map_reference$stat_conditionB,
    pvalue_conditionB = native_weighted_glm_mu_mean_cr_map_reference$pvalue_conditionB,
    padj_conditionB = native_weighted_glm_mu_mean_cr_map_reference$padj_conditionB,
    log_like = native_weighted_glm_mu_mean_cr_map_reference$log_like,
    beta_converged = native_weighted_glm_mu_mean_cr_map_reference$beta_converged,
    beta_iterations = native_weighted_glm_mu_mean_cr_map_reference$beta_iterations,
    check.names = FALSE
  )
  write_tsv(
    native_weighted_glm_mu_mean_cr_map_rows,
    file.path(data_dir, "native_weighted_glm_mu_mean_cr_map_reference.tsv")
  )
  write_result_reference(
    native_weighted_glm_mu_mean_cr_map_reference,
    "native_weighted_glm_mu_mean_cr_results_wald.tsv",
    "Wald"
  )
  native_weighted_glm_mu_mean_cr_lrt_rows <- data.frame(
    gene = rownames(counts),
    baseMean = native_weighted_glm_mu_mean_cr_map_reference$row_meta$baseMean,
    baseVar = native_weighted_glm_mu_mean_cr_map_reference$row_meta$baseVar,
    allZero = native_weighted_glm_mu_mean_cr_map_reference$row_meta$allZero,
    weightsFail = native_weighted_glm_mu_mean_cr_map_reference$weights_fail,
    dispGeneEst = native_weighted_glm_mu_mean_cr_map_reference$row_meta$dispGeneEst,
    dispGeneIter = native_weighted_glm_mu_mean_cr_map_reference$row_meta$dispGeneIter,
    dispFit = native_weighted_glm_mu_mean_cr_map_reference$row_meta$dispFit,
    dispMAP = native_weighted_glm_mu_mean_cr_map_reference$row_meta$dispMAP,
    dispersion = native_weighted_glm_mu_mean_cr_map_reference$row_meta$dispersion,
    beta_intercept = native_weighted_glm_mu_mean_cr_map_reference$beta_intercept,
    beta_conditionB = native_weighted_glm_mu_mean_cr_map_reference$beta_conditionB,
    beta_se_intercept = native_weighted_glm_mu_mean_cr_map_reference$beta_se_intercept,
    beta_se_conditionB = native_weighted_glm_mu_mean_cr_map_reference$beta_se_conditionB,
    log_like_full = native_weighted_glm_mu_mean_cr_map_reference$log_like,
    log_like_reduced = native_weighted_glm_mu_mean_cr_map_reference$reduced_log_like,
    lrt_stat = native_weighted_glm_mu_mean_cr_map_reference$lrt_stat,
    pvalue = native_weighted_glm_mu_mean_cr_map_reference$lrt_pvalue,
    padj = native_weighted_glm_mu_mean_cr_map_reference$lrt_padj,
    df = native_weighted_glm_mu_mean_cr_map_reference$lrt_df,
    full_converged = native_weighted_glm_mu_mean_cr_map_reference$beta_converged,
    reduced_converged = native_weighted_glm_mu_mean_cr_map_reference$reduced_beta_converged,
    full_iterations = native_weighted_glm_mu_mean_cr_map_reference$beta_iterations,
    reduced_iterations = native_weighted_glm_mu_mean_cr_map_reference$reduced_beta_iterations,
    check.names = FALSE
  )
  write_tsv(
    native_weighted_glm_mu_mean_cr_lrt_rows,
    file.path(data_dir, "native_weighted_glm_mu_mean_cr_lrt_reference.tsv")
  )
  write_result_reference(
    native_weighted_glm_mu_mean_cr_map_reference,
    "native_weighted_glm_mu_mean_cr_results_lrt.tsv",
    "LRT"
  )
  write_matrix_tsv(
    native_weighted_glm_mu_mean_cr_map_reference$dispersion_mu,
    "gene",
    file.path(data_dir, "native_weighted_glm_mu_mean_cr_map_dispersion_mu.tsv")
  )
  write_matrix_tsv(
    native_weighted_glm_mu_mean_cr_map_reference$wald_mu,
    "gene",
    file.path(data_dir, "native_weighted_glm_mu_mean_cr_wald_mu.tsv")
  )
  write_matrix_tsv(
    native_weighted_glm_mu_mean_cr_map_reference$wald_hat,
    "gene",
    file.path(data_dir, "native_weighted_glm_mu_mean_cr_wald_hat.tsv")
  )
}
