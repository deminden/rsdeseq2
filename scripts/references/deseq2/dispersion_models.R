trend_means <- c(10, 20, 40, 80, 160, 320)
trend_disps <- 0.05 + 2 / trend_means
trend_fit <- DESeq2:::parametricDispersionFit(trend_means, trend_disps)
trend_coefs <- attr(trend_fit, "coefficients")
write_tsv(
  data.frame(
    gene = paste0("trend_gene", seq_along(trend_means)),
    baseMean = trend_means,
    dispGeneEst = trend_disps,
    dispFit = trend_fit(trend_means),
    useForFit = trend_disps > 100 * 1e-8,
    asymptDisp = unname(trend_coefs[["asymptDisp"]]),
    extraPois = unname(trend_coefs[["extraPois"]]),
    check.names = FALSE
  ),
  file.path(data_dir, "parametric_trend_reference.tsv")
)

trend_prediction_means <- c(12, 30, 120, 500)
write_tsv(
  data.frame(
    mean = trend_prediction_means,
    dispFit = trend_fit(trend_prediction_means),
    check.names = FALSE
  ),
  file.path(data_dir, "parametric_trend_prediction_reference.tsv")
)

mean_trend_means <- c(10, 20, 30, 40, 50, 60)
mean_trend_disps <- c(5e-8, 2e-7, 0.1, 0.2, 0.3, 1000)
dds_mean_trend <- DESeq2::DESeqDataSetFromMatrix(
  countData = matrix(1L, nrow = length(mean_trend_means), ncol = ncol(counts)),
  colData = col_data,
  design = ~ condition
)
S4Vectors::mcols(dds_mean_trend)$baseMean <- mean_trend_means
S4Vectors::mcols(dds_mean_trend)$allZero <- rep(FALSE, length(mean_trend_means))
S4Vectors::mcols(dds_mean_trend)$dispGeneEst <- mean_trend_disps
dds_mean_trend_fit <- DESeq2:::estimateDispersionsFit(
  dds_mean_trend,
  fitType = "mean",
  minDisp = 1e-8,
  quiet = TRUE
)
mean_trend_function <- DESeq2::dispersionFunction(dds_mean_trend_fit)
write_tsv(
  data.frame(
    gene = paste0("mean_trend_gene", seq_along(mean_trend_means)),
    baseMean = mean_trend_means,
    dispGeneEst = mean_trend_disps,
    dispFit = mean_trend_function(mean_trend_means),
    meanDisp = unname(attr(mean_trend_function, "mean")),
    useForFit = mean_trend_disps > 100 * 1e-8,
    useForMean = mean_trend_disps > 10 * 1e-8,
    check.names = FALSE
  ),
  file.path(data_dir, "mean_trend_reference.tsv")
)

local_trend_fit <- DESeq2:::localDispersionFit(trend_means, trend_disps, minDisp = 1e-8)
write_tsv(
  data.frame(
    gene = paste0("trend_gene", seq_along(trend_means)),
    baseMean = trend_means,
    dispGeneEst = trend_disps,
    dispFit = local_trend_fit(trend_means),
    useForFit = trend_disps >= 10 * 1e-8,
    usedMinDispFloor = all(trend_disps < 10 * 1e-8),
    check.names = FALSE
  ),
  file.path(data_dir, "local_trend_reference.tsv")
)

local_prediction_means <- c(12, 30, 120, 500)
write_tsv(
  data.frame(
    mean = local_prediction_means,
    dispFit = local_trend_fit(local_prediction_means),
    check.names = FALSE
  ),
  file.path(data_dir, "local_trend_prediction_reference.tsv")
)

local_floor_means <- c(5, 10, 20, 40)
local_floor_disps <- c(1e-8, 2e-8, 5e-8, 9e-8)
local_floor_fit <- DESeq2:::localDispersionFit(local_floor_means, local_floor_disps, minDisp = 1e-8)
local_floor_disp_fit <- if (is.function(local_floor_fit)) {
  local_floor_fit(local_floor_means)
} else {
  local_floor_fit
}
write_tsv(
  data.frame(
    gene = paste0("local_floor_gene", seq_along(local_floor_means)),
    baseMean = local_floor_means,
    dispGeneEst = local_floor_disps,
    dispFit = local_floor_disp_fit,
    useForFit = local_floor_disps >= 10 * 1e-8,
    usedMinDispFloor = all(local_floor_disps < 10 * 1e-8),
    check.names = FALSE
  ),
  file.path(data_dir, "local_trend_floor_reference.tsv")
)

local_mixed_means <- c(6, 9, 13, 20, 32, 50, 80, 125)
local_mixed_disps <- c(1e-8, 5e-8, 1.1e-7, 0.045, 0.09, 0.055, 0.08, 0.04)
local_mixed_fit <- DESeq2:::localDispersionFit(local_mixed_means, local_mixed_disps, minDisp = 1e-8)
write_tsv(
  data.frame(
    gene = paste0("local_mixed_gene", seq_along(local_mixed_means)),
    baseMean = local_mixed_means,
    dispGeneEst = local_mixed_disps,
    dispFit = local_mixed_fit(local_mixed_means),
    useForFit = local_mixed_disps >= 10 * 1e-8,
    usedMinDispFloor = all(local_mixed_disps < 10 * 1e-8),
    check.names = FALSE
  ),
  file.path(data_dir, "local_trend_mixed_threshold_reference.tsv")
)

prior_residuals <- c(-2, -1, 0, 1, 2)
prior_disp_fit <- rep(1, length(prior_residuals))
prior_disp_gene_est <- exp(prior_residuals) * prior_disp_fit
prior_min_disp <- 1e-8
prior_n_samples <- 10L
prior_n_coefficients <- 2L
prior_residual_df <- prior_n_samples - prior_n_coefficients
prior_above_min <- prior_disp_gene_est >= 100 * prior_min_disp
prior_var_log <- mad(log(prior_disp_gene_est) - log(prior_disp_fit), na.rm = TRUE)^2
prior_expected_sampling <- trigamma(prior_residual_df / 2)
prior_var <- pmax(prior_var_log - prior_expected_sampling, 0.25)
write_tsv(
  data.frame(
    gene = paste0("prior_gene", seq_along(prior_residuals)),
    dispGeneEst = prior_disp_gene_est,
    dispFit = prior_disp_fit,
    aboveMinDisp = prior_above_min,
    nSamples = prior_n_samples,
    nCoefficients = prior_n_coefficients,
    residualDf = prior_residual_df,
    varLogDispEsts = prior_var_log,
    expectedLogDispVariance = prior_expected_sampling,
    dispPriorVar = prior_var,
    check.names = FALSE
  ),
  file.path(data_dir, "dispersion_prior_variance_reference.tsv")
)

low_df_prior_residuals <- c(-1, 0, 1, 2)
low_df_prior_disp_fit <- rep(1, length(low_df_prior_residuals))
low_df_prior_disp_gene_est <- exp(low_df_prior_residuals) * low_df_prior_disp_fit
low_df_prior_counts <- matrix(1L, nrow = length(low_df_prior_residuals), ncol = ncol(counts))
rownames(low_df_prior_counts) <- paste0("low_df_prior_gene", seq_along(low_df_prior_residuals))
colnames(low_df_prior_counts) <- colnames(counts)
dds_low_df_prior <- DESeq2::DESeqDataSetFromMatrix(
  countData = low_df_prior_counts,
  colData = col_data,
  design = ~ condition
)
DESeq2::sizeFactors(dds_low_df_prior) <- rep(1, ncol(counts))
S4Vectors::mcols(dds_low_df_prior)$baseMean <- rep(1, length(low_df_prior_residuals))
S4Vectors::mcols(dds_low_df_prior)$allZero <- rep(FALSE, length(low_df_prior_residuals))
S4Vectors::mcols(dds_low_df_prior)$dispGeneEst <- low_df_prior_disp_gene_est
S4Vectors::mcols(dds_low_df_prior)$dispFit <- low_df_prior_disp_fit
low_df_prior_var_log <- mad(
  log(low_df_prior_disp_gene_est) - log(low_df_prior_disp_fit),
  na.rm = TRUE
)^2
low_df_prior_function <- function(means) rep(1, length(means))
attr(low_df_prior_function, "varLogDispEsts") <- low_df_prior_var_log
DESeq2::dispersionFunction(dds_low_df_prior) <- low_df_prior_function
low_df_prior_design <- matrix(
  c(1, 0, 1, 0, 1, 1, 1, 1),
  nrow = ncol(counts),
  byrow = TRUE
)
low_df_prior_var <- DESeq2:::estimateDispersionsPriorVar(
  dds_low_df_prior,
  minDisp = prior_min_disp,
  modelMatrix = low_df_prior_design
)
write_tsv(
  data.frame(
    gene = rownames(low_df_prior_counts),
    dispGeneEst = low_df_prior_disp_gene_est,
    dispFit = low_df_prior_disp_fit,
    aboveMinDisp = low_df_prior_disp_gene_est >= 100 * prior_min_disp,
    nSamples = nrow(low_df_prior_design),
    nCoefficients = ncol(low_df_prior_design),
    residualDf = nrow(low_df_prior_design) - ncol(low_df_prior_design),
    varLogDispEsts = low_df_prior_var_log,
    expectedLogDispVariance = trigamma((nrow(low_df_prior_design) - ncol(low_df_prior_design)) / 2),
    dispPriorVar = low_df_prior_var,
    check.names = FALSE
  ),
  file.path(data_dir, "dispersion_prior_variance_low_df_reference.tsv")
)

map_counts <- matrix(c(10L, 30L, 10L, 30L), nrow = 1)
rownames(map_counts) <- "map_gene1"
colnames(map_counts) <- colnames(counts)
map_mu <- matrix(rep(20, 4), nrow = 1)
map_disp_gene_est <- 0.5
map_disp_fit <- 0.02
map_prior_var <- 0.05
map_var_log_disp_est <- 100
map_disp_init <- ifelse(map_disp_gene_est > 0.1 * map_disp_fit, map_disp_gene_est, map_disp_fit)
map_fit <- DESeq2:::fitDisp(
  ySEXP = map_counts,
  xSEXP = full_design,
  mu_hatSEXP = map_mu,
  log_alphaSEXP = log(map_disp_init),
  log_alpha_prior_meanSEXP = log(map_disp_fit),
  log_alpha_prior_sigmasqSEXP = map_prior_var,
  min_log_alphaSEXP = log(1e-8 / 10),
  kappa_0SEXP = 1,
  tolSEXP = 1e-6,
  maxitSEXP = 100,
  usePriorSEXP = TRUE,
  weightsSEXP = matrix(1, nrow = 1, ncol = 4),
  useWeightsSEXP = FALSE,
  weightThresholdSEXP = 1e-2,
  useCRSEXP = FALSE
)
map_disp_map <- exp(map_fit$log_alpha)
map_disp_outlier <- log(map_disp_gene_est) > log(map_disp_fit) + 2 * sqrt(map_var_log_disp_est)
map_dispersion <- ifelse(map_disp_outlier, map_disp_gene_est, map_disp_map)
write_tsv(
  data.frame(
    gene = "map_gene1",
    dispGeneEst = map_disp_gene_est,
    dispFit = map_disp_fit,
    dispInit = map_disp_init,
    dispPriorVar = map_prior_var,
    varLogDispEsts = map_var_log_disp_est,
    dispMAP = map_disp_map,
    dispersion = map_dispersion,
    dispIter = map_fit$iter,
    dispOutlier = map_disp_outlier,
    converged = map_fit$iter < 100,
    useCR = FALSE,
    check.names = FALSE
  ),
  file.path(data_dir, "map_dispersion_reference.tsv")
)
