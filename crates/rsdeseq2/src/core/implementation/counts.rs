pub(super) use crate::all_zero::{
    expand_gene_values_with_fill_rows, expand_gene_values_with_nan_rows,
    expand_matrix_with_nan_rows, mask_all_zero_values_with_nan_rows, nonzero_gene_indices,
};
pub(super) use crate::contrasts::{
    ContrastSpec, FactorLevelContrast, ResultsContrast, contrast_all_zero_factor_levels,
    contrast_all_zero_numeric, factor_level_contrast_from_model_frame, resolve_coefficient_index,
    resolve_contrast,
};
pub(super) use crate::cooks::{
    CooksOutput, CooksRefitPlan, CooksReplacementOptions, calculate_cooks_distance,
    prepare_cooks_replacement_refit_with_weights,
};
pub(super) use crate::design::{
    DesignMatrix, FormulaModelFrame, r_like_name_candidates as design_name_candidates,
};
pub(super) use crate::design::{
    ExpandedAdditiveFactorDesign, ExpandedFormulaDesignWithOffsets,
    expanded_formula_design_from_model_frame,
    expanded_formula_design_with_offsets_from_model_frame,
};
pub(super) use crate::dispersion::{
    DispersionTrendFit, GeneWiseDispersionInput, GeneWiseDispersionOptions, MapDispersionInput,
    MapDispersionOptions, estimate_dispersion_prior_variance, estimate_dispersion_trend,
    estimate_gene_wise_dispersions_auto, estimate_gene_wise_dispersions_glm_mu,
    estimate_gene_wise_dispersions_linear_mu, estimate_map_dispersions,
};
pub(super) use crate::errors::{DeseqError, invalid_dimensions};
pub(super) use crate::glm::{
    GlmInput, IrlsOptions, LrtOutput, NbinomGlmFit, ObservationWeightOptions, WaldAlternative,
    WaldContrastOutput, WaldOutput, WaldTestOptions, fit_fixed_dispersion_glm, lrt_test,
    preprocess_observation_weights_with_options, wald_test_coefficient_with_options,
    wald_test_contrast_with_options,
};
pub(super) use crate::independent_filtering::{
    IndependentFilteringOptions, apply_independent_filtering,
};
pub(super) use crate::matrix::RowMajorMatrix;
pub(super) use crate::normalization::{
    Offsets, base_mean, base_mean_with_weights, base_variance, base_variance_with_weights,
    estimate_size_factors_for_dataset, normalized_counts, normalized_counts_with_factors,
    validate_normalization_factors,
};
pub(super) use crate::options::{
    ControlGenes, CooksCutoff, FitType, SizeFactorMethod, SizeFactorOptions, TestType,
};
pub(super) use crate::results::{
    DeseqResultRow, DeseqResults, build_lrt_projection_results_unadjusted_view,
    build_lrt_results_unadjusted, build_wald_contrast_results_unadjusted,
    build_wald_results_from_wald_unadjusted, resolve_cooks_cutoff,
};
pub(super) use crate::transform::{
    FastVstSubset, RlogMetadata, RlogOffsetMode, RlogOutput,
    fast_vst_eligible_count as count_fast_vst_eligible_rows,
    fast_vst_subset as build_fast_vst_subset, norm_transform, rlog_frozen,
    rlog_with_estimated_prior, vst_with_dispersion_trend_and_offsets,
};

/// DESeq2's default number of genes used to fit the fast `vst()` trend.
pub const DEFAULT_FAST_VST_NSUB: usize = 1000;

/// Count storage lives in [`crate::matrix`]; these paths stay for compatibility.
pub use crate::matrix::{CountMatrix, CountsSummary};
