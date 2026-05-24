//! Common public imports for users of the `rsdeseq2` crate.

pub use crate::contrasts::{
    contrast_all_zero_factor_levels, contrast_all_zero_numeric, resolve_contrast, ContrastSpec,
    FactorLevelContrast,
};
pub use crate::cooks::{
    calculate_cooks_distance, max_cooks_after_replacement_refit, prepare_cooks_replacement_refit,
    record_max_cooks, replace_outlier_counts, robust_method_of_moments_dispersion,
    samples_for_cooks, CooksOutput, CooksRefitPlan, CooksReplacementOptions,
    CooksReplacementOutput,
};
pub use crate::core::{
    CooksReplacementLrtOutput, CooksReplacementWaldOutput, CountMatrix, CountsSummary,
    DeseqBuilder, DeseqFit, FastVstGlmMuMetadata, FastVstGlmMuOutput, VstFullDataReason,
    VstGlmMuMetadata, VstGlmMuOutput, VstTrendSource, DEFAULT_FAST_VST_NSUB,
};
pub use crate::design::DesignMatrix;
pub use crate::diagnostics::{
    Deseq2McolsDiagnosticColumn, Deseq2McolsDiagnosticValues, Deseq2McolsDiagnostics,
    Deseq2McolsDiagnosticsDataFrame, DiagnosticSummary,
};
pub use crate::dispersion::{
    cox_reid_adjustment, cox_reid_adjustment_derivative, cox_reid_adjustment_derivative_weighted,
    cox_reid_adjustment_second_derivative, cox_reid_adjustment_second_derivative_weighted,
    cox_reid_adjustment_weighted, dispersion_log_posterior, dispersion_log_posterior_derivative,
    dispersion_log_posterior_derivative_with_prior,
    dispersion_log_posterior_derivative_with_prior_and_weights,
    dispersion_log_posterior_second_derivative,
    dispersion_log_posterior_second_derivative_with_prior,
    dispersion_log_posterior_second_derivative_with_prior_and_weights,
    dispersion_log_posterior_with_prior, dispersion_log_posterior_with_prior_and_weights,
    dispersion_nb_log_likelihood_kernel, dispersion_nb_log_likelihood_kernel_derivative,
    dispersion_nb_log_likelihood_kernel_derivative_weighted,
    dispersion_nb_log_likelihood_kernel_second_derivative,
    dispersion_nb_log_likelihood_kernel_second_derivative_weighted,
    dispersion_nb_log_likelihood_kernel_weighted, dispersion_prior_derivative,
    dispersion_prior_log_density, dispersion_prior_second_derivative, estimate_dispersion_prior,
    estimate_dispersion_prior_variance, estimate_gene_wise_dispersions_glm_mu,
    estimate_gene_wise_dispersions_linear_mu, estimate_low_df_prior_variance,
    estimate_map_dispersions, fit_dispersion_grid, fit_dispersion_grid_no_cr,
    fit_dispersion_grid_no_cr_with_prior, fit_dispersion_grid_no_cr_with_prior_and_weights,
    fit_dispersion_grid_with_prior, fit_dispersion_grid_with_prior_and_weights,
    fit_dispersion_line_search, fit_dispersion_line_search_no_cr,
    fit_dispersion_line_search_no_cr_with_prior,
    fit_dispersion_line_search_no_cr_with_prior_and_weights, fit_dispersion_line_search_with_prior,
    fit_dispersion_line_search_with_prior_and_weights, fit_dispersion_trend,
    fit_local_dispersion_trend, fit_mean_dispersion_trend, fit_parametric_dispersion_trend,
    initial_dispersion_estimates, linear_model_mu, local_trend_use_for_fit,
    log_dispersion_residuals_above_min, mad_squared, map_dispersion_initial_value,
    map_dispersion_outlier, mean_trend_use_for_mean, moments_dispersion_estimates,
    moments_dispersion_estimates_with_normalization_factors, parametric_trend_use_for_fit,
    rough_dispersion_estimates, DispersionLineSearchOutput, DispersionPrior,
    DispersionPriorVarianceOutput, DispersionTrendFit, GeneWiseDispersionFitMethod,
    GeneWiseDispersionInput, GeneWiseDispersionOptions, GeneWiseDispersionOutput,
    LocalDispersionTrend, LocalDispersionTrendFit, LocalDispersionTrendOptions, MapDispersionInput,
    MapDispersionOptions, MapDispersionOutput, MeanDispersionTrend, MeanDispersionTrendFit,
    MeanDispersionTrendOptions, ParametricDispersionTrend, ParametricDispersionTrendFit,
    ParametricDispersionTrendOptions, WeightedDispersionFitInput,
};
pub use crate::errors::DeseqError;
pub use crate::glm::{
    beta_prior_variance_to_ridge_lambda, estimate_beta, estimate_beta_prior_variance,
    fit_fixed_dispersion_irls, fit_fixed_dispersion_irls_with_normalization_factors,
    fit_fixed_dispersion_irls_with_normalization_factors_and_weights,
    fit_fixed_dispersion_irls_with_weights, fit_glms_with_beta_prior_variance,
    fit_glms_with_beta_prior_variance_and_normalization_factors,
    fit_glms_with_beta_prior_variance_and_normalization_factors_and_weights,
    fit_glms_with_estimated_beta_prior_variance, fit_intercept_only_fixed_dispersion,
    fit_intercept_only_fixed_dispersion_with_normalization_factors,
    fit_intercept_only_fixed_dispersion_with_weights, fit_irls, fit_with_dispersion, lrt_test,
    match_upper_quantile_for_variance, match_weighted_upper_quantile_for_variance,
    nbinom_log_likelihood, nbinom_log_likelihood_matrix, nbinom_log_likelihood_weighted,
    nbinom_log_pmf, nbinom_negative_twice_log_likelihood, optim_fallback_rows,
    preprocess_observation_weights, preprocess_observation_weights_with_options,
    two_sided_normal_pvalue, two_sided_t_pvalue, wald_stat_and_pvalue,
    wald_stat_and_pvalue_with_options, wald_test, wald_test_coefficient,
    wald_test_coefficient_with_options, wald_test_contrast, wald_test_contrast_with_options,
    BetaPriorGlmFit, BetaPriorRefitOptions, BetaPriorVarianceMethod, BetaPriorVarianceOptions,
    IrlsOptions, IrlsSolver, LrtOutput, NbinomGlmFit, ObservationWeightOptions, ObservationWeights,
    OptimFallbackRows, WaldAlternative, WaldContrastOutput, WaldDegreesOfFreedom, WaldOutput,
    WaldPvalueType, WaldTestOptions,
};
pub use crate::independent_filtering::{
    apply_independent_filtering, default_theta, filtered_p_adjustments, lowess_fitted_values,
    select_filter_index, select_filter_index_with_lowess, IndependentFilterLowessRow,
    IndependentFilterMetadataEntry, IndependentFilterNumRejRow, IndependentFilteringOptions,
    IndependentFilteringOutput,
};
pub use crate::io::{
    read_count_matrix_tsv, write_base_mean_tsv, write_base_metadata_tsv, write_count_matrix_tsv,
    write_deseq_mcols_diagnostics_tsv, write_deseq_result_column_metadata_tsv,
    write_deseq_result_table_metadata_tsv, write_deseq_results_tidy_tsv, write_deseq_results_tsv,
    write_independent_filter_lowess_tsv, write_independent_filter_metadata_tsv,
    write_independent_filter_num_rej_tsv, write_normalization_factors_tsv,
    write_normalized_counts_tsv, write_size_factors_tsv,
};
pub use crate::math::{
    negative_binomial_helpers, negative_binomial_log_likelihood,
    negative_binomial_log_likelihood_matrix, negative_binomial_log_likelihood_weighted,
    negative_binomial_log_pmf, negative_binomial_negative_twice_log_likelihood, trigamma,
    NegativeBinomialHelpers,
};
pub use crate::matrix::RowMajorMatrix;
pub use crate::multiple_testing::{bh_adjust, bh_adjust_f64};
pub use crate::normalization::{
    base_mean, base_mean_with_weights, base_variance, base_variance_with_weights,
    estimate_size_factors, estimate_size_factors_poscounts,
    estimate_size_factors_poscounts_with_options, estimate_size_factors_ratio,
    estimate_size_factors_ratio_with_options, estimate_size_factors_with_options,
    normalization_factors_from_size_factors, normalized_counts, normalized_counts_with_factors,
    validate_normalization_factors,
};
pub use crate::options::{
    ControlGenes, CooksCutoff, ExecutionMode, FitType, SizeFactorMethod, SizeFactorOptions,
    TestType,
};
pub use crate::results::{
    apply_cooks_cutoff, apply_cooks_cutoff_with_low_count_heuristic, build_lrt_results,
    build_wald_contrast_results, build_wald_results, build_wald_results_from_wald,
    default_cooks_cutoff, deseq2_result_core_column_names, recompute_padj, resolve_cooks_cutoff,
    rsdeseq2_result_diagnostic_column_names, DeseqResultColumn, DeseqResultColumnMetadata,
    DeseqResultColumnValues, DeseqResultRow, DeseqResults, DeseqResultsDataFrame,
    DeseqResultsMetadata, DeseqResultsTableMetadata, DeseqResultsTableMetadataEntry,
    DESEQ2_RESULT_CORE_COLUMNS, RSDESEQ2_RESULT_DIAGNOSTIC_COLUMNS,
};
pub use crate::transform::{
    fast_vst_eligible_count, fast_vst_subset, fast_vst_subset_indices, fast_vst_subset_matrix_rows,
    fast_vst_subset_normalized_counts, local_vst_inverse_size_factor_mean,
    local_vst_inverse_size_factor_mean_from_normalization_factors, norm_transform,
    norm_transform_value, rlog, vst, vst_local, vst_mean, vst_mean_value, vst_parametric,
    vst_parametric_value, vst_with_dispersion_trend,
    vst_with_dispersion_trend_and_normalization_factors,
    vst_with_dispersion_trend_and_size_factors, FastVstSubset, FastVstSubsetMetadata,
};
