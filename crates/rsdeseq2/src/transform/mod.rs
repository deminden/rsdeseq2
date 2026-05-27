//! Implemented and planned DESeq2-compatible transformations.

pub mod norm;
pub mod rlog;
pub mod vst;

pub use norm::{norm_transform, norm_transform_value};
pub use rlog::{
    estimate_rlog_sample_prior_variance, estimate_rlog_sample_prior_variance_with_quantile,
    rlog_beta_prior_variance, rlog_fit_with_normalization_factors, rlog_fit_with_size_factors,
    rlog_frozen_with_normalization_factors, rlog_frozen_with_size_factors, rlog_sample_design,
    rlog_sample_effect_design, rlog_sample_effect_prior_variance,
    rlog_with_estimated_prior_and_normalization_factors,
    rlog_with_estimated_prior_and_size_factors, rlog_with_normalization_factors,
    rlog_with_size_factors, RlogFitOutput, RlogMetadata, RlogOffsetMode, RlogOutput,
    RLOG_INTERCEPT_PRIOR_VARIANCE, RLOG_PRIOR_UPPER_QUANTILE,
};
pub use vst::{
    fast_vst_eligible_count, fast_vst_subset, fast_vst_subset_indices, fast_vst_subset_matrix_rows,
    fast_vst_subset_normalized_counts, local_vst_inverse_size_factor_mean,
    local_vst_inverse_size_factor_mean_from_normalization_factors, vst, vst_local, vst_mean,
    vst_mean_value, vst_parametric, vst_parametric_value, vst_with_dispersion_trend,
    vst_with_dispersion_trend_and_normalization_factors,
    vst_with_dispersion_trend_and_size_factors, FastVstSubset, FastVstSubsetMetadata,
};
