use crate::all_zero::{
    expand_gene_values_with_fill_rows, expand_gene_values_with_nan_rows,
    expand_matrix_with_nan_rows, mask_all_zero_values_with_nan_rows, nonzero_gene_indices,
};
use crate::contrasts::{
    contrast_all_zero_factor_levels, contrast_all_zero_numeric, resolve_coefficient_index,
    resolve_contrast, resolve_results_contrast, factor_level_contrast_from_model_frame,
    ContrastSpec, FactorLevelContrast, ResultsContrast,
};
use crate::cooks::{
    calculate_cooks_distance, prepare_cooks_replacement_refit_with_weights, CooksOutput, CooksRefitPlan,
    CooksReplacementOptions,
};
use crate::design::{r_like_name_candidates as design_name_candidates, DesignMatrix, FormulaModelFrame};
use crate::design::{
    expanded_formula_design_from_model_frame, expanded_formula_design_with_offsets_from_model_frame,
    ExpandedAdditiveFactorDesign, ExpandedFormulaDesignWithOffsets,
};
use crate::dispersion::{
    estimate_dispersion_prior_variance, estimate_gene_wise_dispersions_auto,
    estimate_gene_wise_dispersions_glm_mu,
    estimate_gene_wise_dispersions_linear_mu, estimate_map_dispersions, fit_local_dispersion_trend,
    fit_mean_dispersion_trend, fit_parametric_dispersion_trend, DispersionTrendFit,
    GeneWiseDispersionInput, GeneWiseDispersionOptions, LocalDispersionTrendOptions,
    MapDispersionInput, MapDispersionOptions, MeanDispersionTrendOptions,
    ParametricDispersionTrendOptions,
};
use crate::errors::{invalid_dimensions, DeseqError};
use crate::glm::{
    fit_fixed_dispersion_irls_with_normalization_factors_and_weights,
    fit_fixed_dispersion_irls_with_weights,
    fit_intercept_only_fixed_dispersion_with_normalization_factors,
    fit_intercept_only_fixed_dispersion_with_weights, lrt_test,
    preprocess_observation_weights_with_options, wald_test_coefficient_with_options,
    wald_test_contrast_with_options, IrlsOptions, LrtOutput, NbinomGlmFit,
    ObservationWeightOptions, WaldAlternative, WaldContrastOutput, WaldOutput, WaldTestOptions,
};
use crate::independent_filtering::{apply_independent_filtering, IndependentFilteringOptions};
use crate::matrix::RowMajorMatrix;
use crate::normalization::{
    base_mean, base_mean_with_weights, base_variance, base_variance_with_weights,
    estimate_size_factors_with_options, normalized_counts, normalized_counts_with_factors,
    validate_normalization_factors,
};
use crate::options::{
    ControlGenes, CooksCutoff, FitType, SizeFactorMethod, SizeFactorOptions,
    TestType,
};
use crate::results::{
    build_lrt_projection_results_unadjusted_view, build_lrt_results_unadjusted, build_wald_contrast_results_unadjusted,
    build_wald_results_from_wald_unadjusted, resolve_cooks_cutoff, DeseqResultRow, DeseqResults,
};
use crate::transform::{
    fast_vst_eligible_count as count_fast_vst_eligible_rows,
    fast_vst_subset as build_fast_vst_subset, norm_transform,
    rlog_frozen_with_normalization_factors, rlog_frozen_with_size_factors,
    rlog_with_estimated_prior_and_normalization_factors,
    rlog_with_estimated_prior_and_size_factors,
    vst_with_dispersion_trend_and_normalization_factors,
    vst_with_dispersion_trend_and_size_factors, FastVstSubset, RlogMetadata, RlogOffsetMode,
    RlogOutput,
};

/// DESeq2's default number of genes used to fit the fast `vst()` trend.
pub const DEFAULT_FAST_VST_NSUB: usize = 1000;

/// Count storage lives in [`crate::matrix`]; these paths stay for compatibility.
pub use crate::matrix::{CountMatrix, CountsSummary};
