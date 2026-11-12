native_glm_mu_local_reference <- tryCatch({
  build_glm_mu_reference(dds_ratio, "local", FALSE)
}, error = function(error) {
  writeLines(conditionMessage(error), file.path(data_dir, "native_glm_mu_local_reference_error.txt"))
  NULL
})

if (!is.null(native_glm_mu_local_reference)) {
  native_glm_mu_local_rows <- data.frame(
    gene = rownames(counts),
    baseMean = native_glm_mu_local_reference$row_meta$baseMean,
    baseVar = native_glm_mu_local_reference$row_meta$baseVar,
    allZero = native_glm_mu_local_reference$row_meta$allZero,
    dispGeneEst = native_glm_mu_local_reference$row_meta$dispGeneEst,
    dispGeneIter = native_glm_mu_local_reference$row_meta$dispGeneIter,
    dispFit = native_glm_mu_local_reference$row_meta$dispFit,
    dispPriorVar = native_glm_mu_local_reference$disp_prior_var,
    varLogDispEsts = native_glm_mu_local_reference$var_log_disp_estimates,
    dispMAP = native_glm_mu_local_reference$row_meta$dispMAP,
    dispersion = native_glm_mu_local_reference$row_meta$dispersion,
    dispIter = native_glm_mu_local_reference$row_meta$dispIter,
    dispOutlier = native_glm_mu_local_reference$row_meta$dispOutlier,
    beta_intercept = native_glm_mu_local_reference$beta_intercept,
    beta_conditionB = native_glm_mu_local_reference$beta_conditionB,
    beta_se_intercept = native_glm_mu_local_reference$beta_se_intercept,
    beta_se_conditionB = native_glm_mu_local_reference$beta_se_conditionB,
    stat_conditionB = native_glm_mu_local_reference$stat_conditionB,
    pvalue_conditionB = native_glm_mu_local_reference$pvalue_conditionB,
    padj_conditionB = native_glm_mu_local_reference$padj_conditionB,
    log_like = native_glm_mu_local_reference$log_like,
    beta_converged = native_glm_mu_local_reference$beta_converged,
    beta_iterations = native_glm_mu_local_reference$beta_iterations,
    check.names = FALSE
  )
  write_tsv(
    native_glm_mu_local_rows,
    file.path(data_dir, "native_glm_mu_local_reference.tsv")
  )
  write_result_reference(
    native_glm_mu_local_reference,
    "native_glm_mu_local_results_wald.tsv",
    "Wald"
  )
  native_glm_mu_local_lrt_rows <- data.frame(
    gene = rownames(counts),
    baseMean = native_glm_mu_local_reference$row_meta$baseMean,
    baseVar = native_glm_mu_local_reference$row_meta$baseVar,
    allZero = native_glm_mu_local_reference$row_meta$allZero,
    dispGeneEst = native_glm_mu_local_reference$row_meta$dispGeneEst,
    dispGeneIter = native_glm_mu_local_reference$row_meta$dispGeneIter,
    dispFit = native_glm_mu_local_reference$row_meta$dispFit,
    dispMAP = native_glm_mu_local_reference$row_meta$dispMAP,
    dispersion = native_glm_mu_local_reference$row_meta$dispersion,
    beta_intercept = native_glm_mu_local_reference$beta_intercept,
    beta_conditionB = native_glm_mu_local_reference$beta_conditionB,
    beta_se_intercept = native_glm_mu_local_reference$beta_se_intercept,
    beta_se_conditionB = native_glm_mu_local_reference$beta_se_conditionB,
    log_like_full = native_glm_mu_local_reference$log_like,
    log_like_reduced = native_glm_mu_local_reference$reduced_log_like,
    lrt_stat = native_glm_mu_local_reference$lrt_stat,
    pvalue = native_glm_mu_local_reference$lrt_pvalue,
    padj = native_glm_mu_local_reference$lrt_padj,
    df = native_glm_mu_local_reference$lrt_df,
    full_converged = native_glm_mu_local_reference$beta_converged,
    reduced_converged = native_glm_mu_local_reference$reduced_beta_converged,
    full_iterations = native_glm_mu_local_reference$beta_iterations,
    reduced_iterations = native_glm_mu_local_reference$reduced_beta_iterations,
    check.names = FALSE
  )
  write_tsv(
    native_glm_mu_local_lrt_rows,
    file.path(data_dir, "native_glm_mu_local_lrt_reference.tsv")
  )
  write_result_reference(
    native_glm_mu_local_reference,
    "native_glm_mu_local_results_lrt.tsv",
    "LRT"
  )
  write_matrix_tsv(
    native_glm_mu_local_reference$dispersion_mu,
    "gene",
    file.path(data_dir, "native_glm_mu_local_dispersion_mu.tsv")
  )
  write_matrix_tsv(
    native_glm_mu_local_reference$wald_mu,
    "gene",
    file.path(data_dir, "native_glm_mu_local_wald_mu.tsv")
  )
  write_matrix_tsv(
    native_glm_mu_local_reference$wald_hat,
    "gene",
    file.path(data_dir, "native_glm_mu_local_wald_hat.tsv")
  )
}

native_glm_mu_local_cr_map_reference <- tryCatch({
  build_glm_mu_reference(dds_ratio, "local", TRUE)
}, error = function(error) {
  writeLines(conditionMessage(error), file.path(data_dir, "native_glm_mu_local_cr_map_reference_error.txt"))
  NULL
})

if (!is.null(native_glm_mu_local_cr_map_reference)) {
  native_glm_mu_local_cr_map_rows <- data.frame(
    gene = rownames(counts),
    baseMean = native_glm_mu_local_cr_map_reference$row_meta$baseMean,
    baseVar = native_glm_mu_local_cr_map_reference$row_meta$baseVar,
    allZero = native_glm_mu_local_cr_map_reference$row_meta$allZero,
    dispGeneEst = native_glm_mu_local_cr_map_reference$row_meta$dispGeneEst,
    dispGeneIter = native_glm_mu_local_cr_map_reference$row_meta$dispGeneIter,
    dispFit = native_glm_mu_local_cr_map_reference$row_meta$dispFit,
    dispPriorVar = native_glm_mu_local_cr_map_reference$disp_prior_var,
    varLogDispEsts = native_glm_mu_local_cr_map_reference$var_log_disp_estimates,
    dispMAP = native_glm_mu_local_cr_map_reference$row_meta$dispMAP,
    dispersion = native_glm_mu_local_cr_map_reference$row_meta$dispersion,
    dispIter = native_glm_mu_local_cr_map_reference$row_meta$dispIter,
    dispOutlier = native_glm_mu_local_cr_map_reference$row_meta$dispOutlier,
    beta_intercept = native_glm_mu_local_cr_map_reference$beta_intercept,
    beta_conditionB = native_glm_mu_local_cr_map_reference$beta_conditionB,
    beta_se_intercept = native_glm_mu_local_cr_map_reference$beta_se_intercept,
    beta_se_conditionB = native_glm_mu_local_cr_map_reference$beta_se_conditionB,
    stat_conditionB = native_glm_mu_local_cr_map_reference$stat_conditionB,
    pvalue_conditionB = native_glm_mu_local_cr_map_reference$pvalue_conditionB,
    padj_conditionB = native_glm_mu_local_cr_map_reference$padj_conditionB,
    log_like = native_glm_mu_local_cr_map_reference$log_like,
    beta_converged = native_glm_mu_local_cr_map_reference$beta_converged,
    beta_iterations = native_glm_mu_local_cr_map_reference$beta_iterations,
    check.names = FALSE
  )
  write_tsv(
    native_glm_mu_local_cr_map_rows,
    file.path(data_dir, "native_glm_mu_local_cr_map_reference.tsv")
  )
  write_result_reference(
    native_glm_mu_local_cr_map_reference,
    "native_glm_mu_local_cr_results_wald.tsv",
    "Wald"
  )
  native_glm_mu_local_cr_lrt_rows <- data.frame(
    gene = rownames(counts),
    baseMean = native_glm_mu_local_cr_map_reference$row_meta$baseMean,
    baseVar = native_glm_mu_local_cr_map_reference$row_meta$baseVar,
    allZero = native_glm_mu_local_cr_map_reference$row_meta$allZero,
    dispGeneEst = native_glm_mu_local_cr_map_reference$row_meta$dispGeneEst,
    dispGeneIter = native_glm_mu_local_cr_map_reference$row_meta$dispGeneIter,
    dispFit = native_glm_mu_local_cr_map_reference$row_meta$dispFit,
    dispMAP = native_glm_mu_local_cr_map_reference$row_meta$dispMAP,
    dispersion = native_glm_mu_local_cr_map_reference$row_meta$dispersion,
    beta_intercept = native_glm_mu_local_cr_map_reference$beta_intercept,
    beta_conditionB = native_glm_mu_local_cr_map_reference$beta_conditionB,
    beta_se_intercept = native_glm_mu_local_cr_map_reference$beta_se_intercept,
    beta_se_conditionB = native_glm_mu_local_cr_map_reference$beta_se_conditionB,
    log_like_full = native_glm_mu_local_cr_map_reference$log_like,
    log_like_reduced = native_glm_mu_local_cr_map_reference$reduced_log_like,
    lrt_stat = native_glm_mu_local_cr_map_reference$lrt_stat,
    pvalue = native_glm_mu_local_cr_map_reference$lrt_pvalue,
    padj = native_glm_mu_local_cr_map_reference$lrt_padj,
    df = native_glm_mu_local_cr_map_reference$lrt_df,
    full_converged = native_glm_mu_local_cr_map_reference$beta_converged,
    reduced_converged = native_glm_mu_local_cr_map_reference$reduced_beta_converged,
    full_iterations = native_glm_mu_local_cr_map_reference$beta_iterations,
    reduced_iterations = native_glm_mu_local_cr_map_reference$reduced_beta_iterations,
    check.names = FALSE
  )
  write_tsv(
    native_glm_mu_local_cr_lrt_rows,
    file.path(data_dir, "native_glm_mu_local_cr_lrt_reference.tsv")
  )
  write_result_reference(
    native_glm_mu_local_cr_map_reference,
    "native_glm_mu_local_cr_results_lrt.tsv",
    "LRT"
  )
  write_matrix_tsv(
    native_glm_mu_local_cr_map_reference$dispersion_mu,
    "gene",
    file.path(data_dir, "native_glm_mu_local_cr_map_dispersion_mu.tsv")
  )
  write_matrix_tsv(
    native_glm_mu_local_cr_map_reference$wald_mu,
    "gene",
    file.path(data_dir, "native_glm_mu_local_cr_wald_mu.tsv")
  )
  write_matrix_tsv(
    native_glm_mu_local_cr_map_reference$wald_hat,
    "gene",
    file.path(data_dir, "native_glm_mu_local_cr_wald_hat.tsv")
  )
}

native_weighted_glm_mu_local_reference <- tryCatch({
  build_glm_mu_reference(dds_weighted, "local", FALSE)
}, error = function(error) {
  writeLines(conditionMessage(error), file.path(data_dir, "native_weighted_glm_mu_local_reference_error.txt"))
  NULL
})

if (!is.null(native_weighted_glm_mu_local_reference)) {
  native_weighted_glm_mu_local_rows <- data.frame(
    gene = rownames(counts),
    baseMean = native_weighted_glm_mu_local_reference$row_meta$baseMean,
    baseVar = native_weighted_glm_mu_local_reference$row_meta$baseVar,
    allZero = native_weighted_glm_mu_local_reference$row_meta$allZero,
    weightsFail = native_weighted_glm_mu_local_reference$weights_fail,
    dispGeneEst = native_weighted_glm_mu_local_reference$row_meta$dispGeneEst,
    dispGeneIter = native_weighted_glm_mu_local_reference$row_meta$dispGeneIter,
    dispFit = native_weighted_glm_mu_local_reference$row_meta$dispFit,
    dispPriorVar = native_weighted_glm_mu_local_reference$disp_prior_var,
    varLogDispEsts = native_weighted_glm_mu_local_reference$var_log_disp_estimates,
    dispMAP = native_weighted_glm_mu_local_reference$row_meta$dispMAP,
    dispersion = native_weighted_glm_mu_local_reference$row_meta$dispersion,
    dispIter = native_weighted_glm_mu_local_reference$row_meta$dispIter,
    dispOutlier = native_weighted_glm_mu_local_reference$row_meta$dispOutlier,
    beta_intercept = native_weighted_glm_mu_local_reference$beta_intercept,
    beta_conditionB = native_weighted_glm_mu_local_reference$beta_conditionB,
    beta_se_intercept = native_weighted_glm_mu_local_reference$beta_se_intercept,
    beta_se_conditionB = native_weighted_glm_mu_local_reference$beta_se_conditionB,
    stat_conditionB = native_weighted_glm_mu_local_reference$stat_conditionB,
    pvalue_conditionB = native_weighted_glm_mu_local_reference$pvalue_conditionB,
    padj_conditionB = native_weighted_glm_mu_local_reference$padj_conditionB,
    log_like = native_weighted_glm_mu_local_reference$log_like,
    beta_converged = native_weighted_glm_mu_local_reference$beta_converged,
    beta_iterations = native_weighted_glm_mu_local_reference$beta_iterations,
    check.names = FALSE
  )
  write_tsv(
    native_weighted_glm_mu_local_rows,
    file.path(data_dir, "native_weighted_glm_mu_local_reference.tsv")
  )
  write_result_reference(
    native_weighted_glm_mu_local_reference,
    "native_weighted_glm_mu_local_results_wald.tsv",
    "Wald"
  )
  native_weighted_glm_mu_local_lrt_rows <- data.frame(
    gene = rownames(counts),
    baseMean = native_weighted_glm_mu_local_reference$row_meta$baseMean,
    baseVar = native_weighted_glm_mu_local_reference$row_meta$baseVar,
    allZero = native_weighted_glm_mu_local_reference$row_meta$allZero,
    weightsFail = native_weighted_glm_mu_local_reference$weights_fail,
    dispGeneEst = native_weighted_glm_mu_local_reference$row_meta$dispGeneEst,
    dispGeneIter = native_weighted_glm_mu_local_reference$row_meta$dispGeneIter,
    dispFit = native_weighted_glm_mu_local_reference$row_meta$dispFit,
    dispMAP = native_weighted_glm_mu_local_reference$row_meta$dispMAP,
    dispersion = native_weighted_glm_mu_local_reference$row_meta$dispersion,
    beta_intercept = native_weighted_glm_mu_local_reference$beta_intercept,
    beta_conditionB = native_weighted_glm_mu_local_reference$beta_conditionB,
    beta_se_intercept = native_weighted_glm_mu_local_reference$beta_se_intercept,
    beta_se_conditionB = native_weighted_glm_mu_local_reference$beta_se_conditionB,
    log_like_full = native_weighted_glm_mu_local_reference$log_like,
    log_like_reduced = native_weighted_glm_mu_local_reference$reduced_log_like,
    lrt_stat = native_weighted_glm_mu_local_reference$lrt_stat,
    pvalue = native_weighted_glm_mu_local_reference$lrt_pvalue,
    padj = native_weighted_glm_mu_local_reference$lrt_padj,
    df = native_weighted_glm_mu_local_reference$lrt_df,
    full_converged = native_weighted_glm_mu_local_reference$beta_converged,
    reduced_converged = native_weighted_glm_mu_local_reference$reduced_beta_converged,
    full_iterations = native_weighted_glm_mu_local_reference$beta_iterations,
    reduced_iterations = native_weighted_glm_mu_local_reference$reduced_beta_iterations,
    check.names = FALSE
  )
  write_tsv(
    native_weighted_glm_mu_local_lrt_rows,
    file.path(data_dir, "native_weighted_glm_mu_local_lrt_reference.tsv")
  )
  write_result_reference(
    native_weighted_glm_mu_local_reference,
    "native_weighted_glm_mu_local_results_lrt.tsv",
    "LRT"
  )
  write_matrix_tsv(
    native_weighted_glm_mu_local_reference$dispersion_mu,
    "gene",
    file.path(data_dir, "native_weighted_glm_mu_local_dispersion_mu.tsv")
  )
  write_matrix_tsv(
    native_weighted_glm_mu_local_reference$wald_mu,
    "gene",
    file.path(data_dir, "native_weighted_glm_mu_local_wald_mu.tsv")
  )
  write_matrix_tsv(
    native_weighted_glm_mu_local_reference$wald_hat,
    "gene",
    file.path(data_dir, "native_weighted_glm_mu_local_wald_hat.tsv")
  )
}

native_weighted_glm_mu_local_cr_map_reference <- tryCatch({
  build_glm_mu_reference(dds_weighted, "local", TRUE)
}, error = function(error) {
  writeLines(conditionMessage(error), file.path(data_dir, "native_weighted_glm_mu_local_cr_map_reference_error.txt"))
  NULL
})

if (!is.null(native_weighted_glm_mu_local_cr_map_reference)) {
  native_weighted_glm_mu_local_cr_map_rows <- data.frame(
    gene = rownames(counts),
    baseMean = native_weighted_glm_mu_local_cr_map_reference$row_meta$baseMean,
    baseVar = native_weighted_glm_mu_local_cr_map_reference$row_meta$baseVar,
    allZero = native_weighted_glm_mu_local_cr_map_reference$row_meta$allZero,
    weightsFail = native_weighted_glm_mu_local_cr_map_reference$weights_fail,
    dispGeneEst = native_weighted_glm_mu_local_cr_map_reference$row_meta$dispGeneEst,
    dispGeneIter = native_weighted_glm_mu_local_cr_map_reference$row_meta$dispGeneIter,
    dispFit = native_weighted_glm_mu_local_cr_map_reference$row_meta$dispFit,
    dispPriorVar = native_weighted_glm_mu_local_cr_map_reference$disp_prior_var,
    varLogDispEsts = native_weighted_glm_mu_local_cr_map_reference$var_log_disp_estimates,
    dispMAP = native_weighted_glm_mu_local_cr_map_reference$row_meta$dispMAP,
    dispersion = native_weighted_glm_mu_local_cr_map_reference$row_meta$dispersion,
    dispIter = native_weighted_glm_mu_local_cr_map_reference$row_meta$dispIter,
    dispOutlier = native_weighted_glm_mu_local_cr_map_reference$row_meta$dispOutlier,
    beta_intercept = native_weighted_glm_mu_local_cr_map_reference$beta_intercept,
    beta_conditionB = native_weighted_glm_mu_local_cr_map_reference$beta_conditionB,
    beta_se_intercept = native_weighted_glm_mu_local_cr_map_reference$beta_se_intercept,
    beta_se_conditionB = native_weighted_glm_mu_local_cr_map_reference$beta_se_conditionB,
    stat_conditionB = native_weighted_glm_mu_local_cr_map_reference$stat_conditionB,
    pvalue_conditionB = native_weighted_glm_mu_local_cr_map_reference$pvalue_conditionB,
    padj_conditionB = native_weighted_glm_mu_local_cr_map_reference$padj_conditionB,
    log_like = native_weighted_glm_mu_local_cr_map_reference$log_like,
    beta_converged = native_weighted_glm_mu_local_cr_map_reference$beta_converged,
    beta_iterations = native_weighted_glm_mu_local_cr_map_reference$beta_iterations,
    check.names = FALSE
  )
  write_tsv(
    native_weighted_glm_mu_local_cr_map_rows,
    file.path(data_dir, "native_weighted_glm_mu_local_cr_map_reference.tsv")
  )
  write_result_reference(
    native_weighted_glm_mu_local_cr_map_reference,
    "native_weighted_glm_mu_local_cr_results_wald.tsv",
    "Wald"
  )
  native_weighted_glm_mu_local_cr_lrt_rows <- data.frame(
    gene = rownames(counts),
    baseMean = native_weighted_glm_mu_local_cr_map_reference$row_meta$baseMean,
    baseVar = native_weighted_glm_mu_local_cr_map_reference$row_meta$baseVar,
    allZero = native_weighted_glm_mu_local_cr_map_reference$row_meta$allZero,
    weightsFail = native_weighted_glm_mu_local_cr_map_reference$weights_fail,
    dispGeneEst = native_weighted_glm_mu_local_cr_map_reference$row_meta$dispGeneEst,
    dispGeneIter = native_weighted_glm_mu_local_cr_map_reference$row_meta$dispGeneIter,
    dispFit = native_weighted_glm_mu_local_cr_map_reference$row_meta$dispFit,
    dispMAP = native_weighted_glm_mu_local_cr_map_reference$row_meta$dispMAP,
    dispersion = native_weighted_glm_mu_local_cr_map_reference$row_meta$dispersion,
    beta_intercept = native_weighted_glm_mu_local_cr_map_reference$beta_intercept,
    beta_conditionB = native_weighted_glm_mu_local_cr_map_reference$beta_conditionB,
    beta_se_intercept = native_weighted_glm_mu_local_cr_map_reference$beta_se_intercept,
    beta_se_conditionB = native_weighted_glm_mu_local_cr_map_reference$beta_se_conditionB,
    log_like_full = native_weighted_glm_mu_local_cr_map_reference$log_like,
    log_like_reduced = native_weighted_glm_mu_local_cr_map_reference$reduced_log_like,
    lrt_stat = native_weighted_glm_mu_local_cr_map_reference$lrt_stat,
    pvalue = native_weighted_glm_mu_local_cr_map_reference$lrt_pvalue,
    padj = native_weighted_glm_mu_local_cr_map_reference$lrt_padj,
    df = native_weighted_glm_mu_local_cr_map_reference$lrt_df,
    full_converged = native_weighted_glm_mu_local_cr_map_reference$beta_converged,
    reduced_converged = native_weighted_glm_mu_local_cr_map_reference$reduced_beta_converged,
    full_iterations = native_weighted_glm_mu_local_cr_map_reference$beta_iterations,
    reduced_iterations = native_weighted_glm_mu_local_cr_map_reference$reduced_beta_iterations,
    check.names = FALSE
  )
  write_tsv(
    native_weighted_glm_mu_local_cr_lrt_rows,
    file.path(data_dir, "native_weighted_glm_mu_local_cr_lrt_reference.tsv")
  )
  write_result_reference(
    native_weighted_glm_mu_local_cr_map_reference,
    "native_weighted_glm_mu_local_cr_results_lrt.tsv",
    "LRT"
  )
  write_matrix_tsv(
    native_weighted_glm_mu_local_cr_map_reference$dispersion_mu,
    "gene",
    file.path(data_dir, "native_weighted_glm_mu_local_cr_map_dispersion_mu.tsv")
  )
  write_matrix_tsv(
    native_weighted_glm_mu_local_cr_map_reference$wald_mu,
    "gene",
    file.path(data_dir, "native_weighted_glm_mu_local_cr_wald_mu.tsv")
  )
  write_matrix_tsv(
    native_weighted_glm_mu_local_cr_map_reference$wald_hat,
    "gene",
    file.path(data_dir, "native_weighted_glm_mu_local_cr_wald_hat.tsv")
  )
}
