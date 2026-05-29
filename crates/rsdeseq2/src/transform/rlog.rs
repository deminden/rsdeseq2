use crate::core::CountMatrix;
use crate::design::DesignMatrix;
use crate::errors::{invalid_dimensions, DeseqError};
use crate::glm::{
    fit_glms_with_beta_prior_variance, fit_glms_with_beta_prior_variance_and_normalization_factors,
    match_weighted_upper_quantile_for_variance, IrlsOptions, NbinomGlmFit,
};
use crate::matrix::RowMajorMatrix;
use crate::normalization::{normalized_counts, normalized_counts_with_factors};

/// Wide log2-scale prior variance used for the rlog intercept.
///
/// DESeq2 keeps the intercept effectively unregularized and applies the sample
/// prior to sample effects. This finite value maps to a tiny ridge while still
/// keeping the GLM system numerically regularized.
pub const RLOG_INTERCEPT_PRIOR_VARIANCE: f64 = 1.0e6;

/// Upper-tail quantile used by DESeq2 when estimating rlog sample prior variance.
pub const RLOG_PRIOR_UPPER_QUANTILE: f64 = 0.05;

/// Result of an rlog fit that estimates the sample-effect prior variance.
#[derive(Clone, Debug, PartialEq)]
pub struct RlogOutput {
    /// Genes x samples rlog values on the log2 normalized-expression scale.
    pub transformed: RowMajorMatrix<f64>,
    /// Per-gene fitted rlog intercepts on the log2 scale.
    pub intercept: Vec<f64>,
    /// Estimated shared prior variance for sample effects.
    pub sample_prior_variance: f64,
    /// Offset source used by the rlog fit.
    pub offset_mode: RlogOffsetMode,
}

/// Low-level rlog output retaining the fitted sample-effect GLM.
#[derive(Clone, Debug, PartialEq)]
pub struct RlogFitOutput {
    /// Genes x samples rlog values on the log2 normalized-expression scale.
    pub transformed: RowMajorMatrix<f64>,
    /// Per-gene fitted rlog intercepts on the log2 scale.
    pub intercept: Vec<f64>,
    /// Fitted GLM for the intercept-plus-sample-effect rlog design.
    pub fit: NbinomGlmFit,
}

/// Offset source used by an rlog fit.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum RlogOffsetMode {
    /// Ordinary sample-specific size factors.
    SizeFactors,
    /// Gene/sample normalization factors.
    NormalizationFactors,
}

impl RlogOffsetMode {
    /// Stable label for wrappers and benchmark logs.
    pub fn as_str(self) -> &'static str {
        match self {
            Self::SizeFactors => "sizeFactors",
            Self::NormalizationFactors => "normalizationFactors",
        }
    }
}

/// Metadata summary for an rlog output.
#[derive(Clone, Debug, PartialEq)]
pub struct RlogMetadata {
    /// Number of rows in the transformed matrix.
    pub transformed_rows: usize,
    /// Number of columns in the transformed matrix.
    pub transformed_cols: usize,
    /// Number of fitted rlog intercepts.
    pub intercept_len: usize,
    /// Estimated shared prior variance for sample effects.
    pub sample_prior_variance: f64,
    /// Stable offset source label.
    pub offset_mode: &'static str,
}

impl RlogOutput {
    /// Metadata view for wrappers, diagnostics, and benchmark logs.
    pub fn metadata(&self) -> RlogMetadata {
        RlogMetadata {
            transformed_rows: self.transformed.n_rows(),
            transformed_cols: self.transformed.n_cols(),
            intercept_len: self.intercept.len(),
            sample_prior_variance: self.sample_prior_variance,
            offset_mode: self.offset_mode.as_str(),
        }
    }
}

/// Fit the low-level regularized-log sample-effect model with size factors.
///
/// This is the core numeric shape behind DESeq2's rlog: one intercept plus one
/// indicator column per sample, with a wide intercept prior and a shared prior
/// variance for sample effects. The returned matrix is genes x samples on the
/// log2 normalized-expression scale, assembled as intercept plus the matching
/// sample effect for each sample.
///
/// This helper expects already-estimated dispersions and the sample-effect
/// prior variance. Full high-level DESeq2 rlog object semantics, including
/// automatic dispersion/prior estimation and blind/design-aware wrapper
/// dispatch, remain separate workflow work.
pub fn rlog_with_size_factors(
    counts: &CountMatrix,
    size_factors: &[f64],
    dispersions: &[f64],
    sample_prior_variance: f64,
    options: IrlsOptions,
) -> Result<RowMajorMatrix<f64>, DeseqError> {
    rlog_fit_with_size_factors(
        counts,
        size_factors,
        dispersions,
        sample_prior_variance,
        options,
    )
    .map(|output| output.transformed)
}

/// Fit the low-level rlog sample-effect model with size factors and retain the GLM fit.
pub fn rlog_fit_with_size_factors(
    counts: &CountMatrix,
    size_factors: &[f64],
    dispersions: &[f64],
    sample_prior_variance: f64,
    options: IrlsOptions,
) -> Result<RlogFitOutput, DeseqError> {
    let design = rlog_sample_design(counts.n_samples())?;
    let beta_prior_variance = rlog_beta_prior_variance(counts.n_samples(), sample_prior_variance)?;
    let fit = fit_glms_with_beta_prior_variance(
        counts,
        &design,
        size_factors,
        dispersions,
        &beta_prior_variance,
        options,
    )?;
    rlog_output_from_fit(fit, counts.n_samples())
}

/// Fit an rlog transform using caller-supplied frozen intercepts and size factors.
///
/// The supplied intercepts are treated as gene-specific log2 offsets. Only the
/// sample-effect coefficients are fit, using the shared sample prior variance,
/// and the output is assembled as `frozen_intercept_i + sample_effect_ij`.
pub fn rlog_frozen_with_size_factors(
    counts: &CountMatrix,
    size_factors: &[f64],
    dispersions: &[f64],
    sample_prior_variance: f64,
    frozen_intercept: &[f64],
    options: IrlsOptions,
) -> Result<RlogOutput, DeseqError> {
    if size_factors.len() != counts.n_samples() {
        return Err(invalid_dimensions(
            "rlog size factors",
            counts.n_samples(),
            size_factors.len(),
        ));
    }
    let normalization_factors = rlog_normalization_factors_from_size_factors(counts, size_factors)?;
    rlog_frozen_with_normalization_factors(
        counts,
        &normalization_factors,
        dispersions,
        sample_prior_variance,
        frozen_intercept,
        options,
    )
    .map(|mut output| {
        output.offset_mode = RlogOffsetMode::SizeFactors;
        output
    })
}

/// Estimate the rlog sample prior from normalized counts and fit with size factors.
///
/// This keeps the caller in control of dispersion inputs while avoiding a
/// duplicated manual prior-estimation step for workflows that already have
/// `baseMean` and fitted dispersion-trend values.
pub fn rlog_with_estimated_prior_and_size_factors(
    counts: &CountMatrix,
    size_factors: &[f64],
    base_mean: &[f64],
    disp_fit: &[f64],
    dispersions: &[f64],
    options: IrlsOptions,
) -> Result<RlogOutput, DeseqError> {
    let normalized = normalized_counts(counts, size_factors)?;
    let sample_prior_variance =
        estimate_rlog_sample_prior_variance(&normalized, base_mean, disp_fit)?;
    let fit_output = rlog_fit_with_size_factors(
        counts,
        size_factors,
        dispersions,
        sample_prior_variance,
        options,
    )?;
    Ok(RlogOutput {
        transformed: fit_output.transformed,
        intercept: fit_output.intercept,
        sample_prior_variance,
        offset_mode: RlogOffsetMode::SizeFactors,
    })
}

/// Estimate DESeq2's rlog sample-effect prior variance from normalized counts.
///
/// This mirrors the rlog-specific prior estimate:
///
/// ```text
/// log2(normalized_count + 0.5) - log2(baseMean + 0.5)
/// weight = 1 / (1 / baseMean + dispFit)
/// ```
///
/// The flattened log-fold-change values are matched to a zero-centered Normal
/// using the same weighted upper-quantile helper as the beta-prior code.
pub fn estimate_rlog_sample_prior_variance(
    normalized_counts: &RowMajorMatrix<f64>,
    base_mean: &[f64],
    disp_fit: &[f64],
) -> Result<f64, DeseqError> {
    estimate_rlog_sample_prior_variance_with_quantile(
        normalized_counts,
        base_mean,
        disp_fit,
        RLOG_PRIOR_UPPER_QUANTILE,
    )
}

/// Estimate rlog sample prior variance with an explicit upper-tail quantile.
pub fn estimate_rlog_sample_prior_variance_with_quantile(
    normalized_counts: &RowMajorMatrix<f64>,
    base_mean: &[f64],
    disp_fit: &[f64],
    upper_quantile: f64,
) -> Result<f64, DeseqError> {
    if base_mean.len() != normalized_counts.n_rows() {
        return Err(invalid_dimensions(
            "rlog baseMean rows",
            normalized_counts.n_rows(),
            base_mean.len(),
        ));
    }
    if disp_fit.len() != normalized_counts.n_rows() {
        return Err(invalid_dimensions(
            "rlog dispFit rows",
            normalized_counts.n_rows(),
            disp_fit.len(),
        ));
    }
    let mut log_fold_changes = Vec::with_capacity(normalized_counts.len());
    let mut weights = Vec::with_capacity(normalized_counts.len());
    for gene in 0..normalized_counts.n_rows() {
        let mean = base_mean[gene];
        let dispersion = disp_fit[gene];
        if !mean.is_finite() || mean <= 0.0 {
            return Err(DeseqError::NonFiniteValue {
                context: "rlog baseMean".to_string(),
                index: Some(gene),
                value: mean,
            });
        }
        if !dispersion.is_finite() || dispersion < 0.0 {
            return Err(DeseqError::NonFiniteValue {
                context: "rlog dispFit".to_string(),
                index: Some(gene),
                value: dispersion,
            });
        }
        let variance = 1.0 / mean + dispersion;
        if !variance.is_finite() || variance <= 0.0 {
            return Err(DeseqError::InvalidOptions {
                reason: format!("rlog prior weight variance is invalid for gene {gene}"),
            });
        }
        let weight = 1.0 / variance;
        let mean_log = (mean + 0.5).log2();
        for sample in 0..normalized_counts.n_cols() {
            let value = normalized_counts
                .get(gene, sample)
                .copied()
                .unwrap_or(f64::NAN);
            if !value.is_finite() || value < 0.0 {
                return Err(DeseqError::NonFiniteValue {
                    context: "rlog normalized count".to_string(),
                    index: Some(gene * normalized_counts.n_cols() + sample),
                    value,
                });
            }
            log_fold_changes.push((value + 0.5).log2() - mean_log);
            weights.push(weight);
        }
    }
    match_weighted_upper_quantile_for_variance(&log_fold_changes, &weights, upper_quantile)
}

/// Fit the low-level regularized-log sample-effect model with normalization factors.
pub fn rlog_with_normalization_factors(
    counts: &CountMatrix,
    normalization_factors: &RowMajorMatrix<f64>,
    dispersions: &[f64],
    sample_prior_variance: f64,
    options: IrlsOptions,
) -> Result<RowMajorMatrix<f64>, DeseqError> {
    rlog_fit_with_normalization_factors(
        counts,
        normalization_factors,
        dispersions,
        sample_prior_variance,
        options,
    )
    .map(|output| output.transformed)
}

/// Fit the low-level rlog sample-effect model with normalization factors and retain the GLM fit.
pub fn rlog_fit_with_normalization_factors(
    counts: &CountMatrix,
    normalization_factors: &RowMajorMatrix<f64>,
    dispersions: &[f64],
    sample_prior_variance: f64,
    options: IrlsOptions,
) -> Result<RlogFitOutput, DeseqError> {
    let design = rlog_sample_design(counts.n_samples())?;
    let beta_prior_variance = rlog_beta_prior_variance(counts.n_samples(), sample_prior_variance)?;
    let fit = fit_glms_with_beta_prior_variance_and_normalization_factors(
        counts,
        &design,
        normalization_factors,
        dispersions,
        &beta_prior_variance,
        options,
    )?;
    rlog_output_from_fit(fit, counts.n_samples())
}

/// Fit an rlog transform using caller-supplied frozen intercepts and normalization factors.
pub fn rlog_frozen_with_normalization_factors(
    counts: &CountMatrix,
    normalization_factors: &RowMajorMatrix<f64>,
    dispersions: &[f64],
    sample_prior_variance: f64,
    frozen_intercept: &[f64],
    options: IrlsOptions,
) -> Result<RlogOutput, DeseqError> {
    validate_rlog_frozen_inputs(
        counts,
        normalization_factors,
        dispersions,
        sample_prior_variance,
        frozen_intercept,
    )?;
    let design = rlog_sample_effect_design(counts.n_samples())?;
    let beta_prior_variance =
        rlog_sample_effect_prior_variance(counts.n_samples(), sample_prior_variance)?;
    let frozen_factors =
        rlog_frozen_normalization_factors(normalization_factors, frozen_intercept)?;
    let fit = fit_glms_with_beta_prior_variance_and_normalization_factors(
        counts,
        &design,
        &frozen_factors,
        dispersions,
        &beta_prior_variance,
        options,
    )?;
    let transformed = rlog_transform_from_sample_effect_fit(&fit.beta, frozen_intercept)?;
    Ok(RlogOutput {
        transformed,
        intercept: frozen_intercept.to_vec(),
        sample_prior_variance,
        offset_mode: RlogOffsetMode::NormalizationFactors,
    })
}

/// Estimate the rlog sample prior and fit with gene/sample normalization factors.
pub fn rlog_with_estimated_prior_and_normalization_factors(
    counts: &CountMatrix,
    normalization_factors: &RowMajorMatrix<f64>,
    base_mean: &[f64],
    disp_fit: &[f64],
    dispersions: &[f64],
    options: IrlsOptions,
) -> Result<RlogOutput, DeseqError> {
    let normalized = normalized_counts_with_factors(counts, normalization_factors)?;
    let sample_prior_variance =
        estimate_rlog_sample_prior_variance(&normalized, base_mean, disp_fit)?;
    let fit_output = rlog_fit_with_normalization_factors(
        counts,
        normalization_factors,
        dispersions,
        sample_prior_variance,
        options,
    )?;
    Ok(RlogOutput {
        transformed: fit_output.transformed,
        intercept: fit_output.intercept,
        sample_prior_variance,
        offset_mode: RlogOffsetMode::NormalizationFactors,
    })
}

/// Build the rlog sample design: intercept plus one indicator per sample.
pub fn rlog_sample_design(n_samples: usize) -> Result<DesignMatrix, DeseqError> {
    if n_samples == 0 {
        return Err(DeseqError::InvalidOptions {
            reason: "rlog requires at least one sample".to_string(),
        });
    }
    let n_coefficients = n_samples + 1;
    let mut values = Vec::with_capacity(n_samples * n_coefficients);
    for sample in 0..n_samples {
        values.push(1.0);
        for candidate in 0..n_samples {
            values.push((sample == candidate) as u8 as f64);
        }
    }
    let mut names = Vec::with_capacity(n_coefficients);
    names.push("Intercept".to_string());
    names.extend((0..n_samples).map(|sample| format!("sample_{sample}")));
    DesignMatrix::from_row_major(n_samples, n_coefficients, values, Some(names))
}

/// Build the frozen-rlog sample-effect design: one indicator per sample.
pub fn rlog_sample_effect_design(n_samples: usize) -> Result<DesignMatrix, DeseqError> {
    if n_samples == 0 {
        return Err(DeseqError::InvalidOptions {
            reason: "rlog requires at least one sample".to_string(),
        });
    }
    let mut values = Vec::with_capacity(n_samples * n_samples);
    for sample in 0..n_samples {
        for candidate in 0..n_samples {
            values.push((sample == candidate) as u8 as f64);
        }
    }
    let names = (0..n_samples)
        .map(|sample| format!("sample_{sample}"))
        .collect();
    DesignMatrix::from_row_major(n_samples, n_samples, values, Some(names))
}

/// Build the log2-scale rlog prior vector for one intercept plus sample effects.
pub fn rlog_beta_prior_variance(
    n_samples: usize,
    sample_prior_variance: f64,
) -> Result<Vec<f64>, DeseqError> {
    if n_samples == 0 {
        return Err(DeseqError::InvalidOptions {
            reason: "rlog requires at least one sample".to_string(),
        });
    }
    validate_rlog_sample_prior_variance(sample_prior_variance)?;
    let mut beta_prior_variance = Vec::with_capacity(n_samples + 1);
    beta_prior_variance.push(RLOG_INTERCEPT_PRIOR_VARIANCE);
    beta_prior_variance.extend(std::iter::repeat_n(sample_prior_variance, n_samples));
    Ok(beta_prior_variance)
}

/// Build the log2-scale prior vector for frozen-rlog sample effects.
pub fn rlog_sample_effect_prior_variance(
    n_samples: usize,
    sample_prior_variance: f64,
) -> Result<Vec<f64>, DeseqError> {
    if n_samples == 0 {
        return Err(DeseqError::InvalidOptions {
            reason: "rlog requires at least one sample".to_string(),
        });
    }
    validate_rlog_sample_prior_variance(sample_prior_variance)?;
    Ok(std::iter::repeat_n(sample_prior_variance, n_samples).collect())
}

fn rlog_transform_from_fit(
    beta: &RowMajorMatrix<f64>,
    n_samples: usize,
) -> Result<RowMajorMatrix<f64>, DeseqError> {
    if beta.n_cols() != n_samples + 1 {
        return Err(invalid_dimensions(
            "rlog beta columns",
            n_samples + 1,
            beta.n_cols(),
        ));
    }
    let mut values = Vec::with_capacity(beta.n_rows() * n_samples);
    for gene in 0..beta.n_rows() {
        let intercept = beta.get(gene, 0).copied().unwrap_or(f64::NAN);
        if !intercept.is_finite() {
            return Err(DeseqError::NonFiniteValue {
                context: "rlog intercept beta".to_string(),
                index: Some(gene),
                value: intercept,
            });
        }
        for sample in 0..n_samples {
            let effect = beta.get(gene, sample + 1).copied().unwrap_or(f64::NAN);
            let transformed = intercept + effect;
            if !transformed.is_finite() {
                return Err(DeseqError::NonFiniteValue {
                    context: "rlog transformed value".to_string(),
                    index: Some(gene * n_samples + sample),
                    value: transformed,
                });
            }
            values.push(transformed);
        }
    }
    RowMajorMatrix::from_row_major(beta.n_rows(), n_samples, values)
}

fn rlog_transform_from_sample_effect_fit(
    beta: &RowMajorMatrix<f64>,
    frozen_intercept: &[f64],
) -> Result<RowMajorMatrix<f64>, DeseqError> {
    if beta.n_rows() != frozen_intercept.len() {
        return Err(invalid_dimensions(
            "rlog frozen intercept rows",
            beta.n_rows(),
            frozen_intercept.len(),
        ));
    }
    let n_samples = beta.n_cols();
    let mut values = Vec::with_capacity(beta.n_rows() * n_samples);
    for (gene, intercept) in frozen_intercept.iter().enumerate() {
        if !intercept.is_finite() {
            return Err(DeseqError::NonFiniteValue {
                context: "rlog frozen intercept".to_string(),
                index: Some(gene),
                value: *intercept,
            });
        }
        for sample in 0..n_samples {
            let effect = beta.get(gene, sample).copied().unwrap_or(f64::NAN);
            let transformed = intercept + effect;
            if !transformed.is_finite() {
                return Err(DeseqError::NonFiniteValue {
                    context: "rlog frozen transformed value".to_string(),
                    index: Some(gene * n_samples + sample),
                    value: transformed,
                });
            }
            values.push(transformed);
        }
    }
    RowMajorMatrix::from_row_major(beta.n_rows(), n_samples, values)
}

fn rlog_output_from_fit(fit: NbinomGlmFit, n_samples: usize) -> Result<RlogFitOutput, DeseqError> {
    let transformed = rlog_transform_from_fit(&fit.beta, n_samples)?;
    let intercept = rlog_intercept_from_beta(&fit.beta, n_samples)?;
    Ok(RlogFitOutput {
        transformed,
        intercept,
        fit,
    })
}

fn rlog_intercept_from_beta(
    beta: &RowMajorMatrix<f64>,
    n_samples: usize,
) -> Result<Vec<f64>, DeseqError> {
    if beta.n_cols() != n_samples + 1 {
        return Err(invalid_dimensions(
            "rlog beta columns",
            n_samples + 1,
            beta.n_cols(),
        ));
    }
    let mut intercepts = Vec::with_capacity(beta.n_rows());
    for gene in 0..beta.n_rows() {
        let intercept = beta.get(gene, 0).copied().unwrap_or(f64::NAN);
        if !intercept.is_finite() {
            return Err(DeseqError::NonFiniteValue {
                context: "rlog intercept beta".to_string(),
                index: Some(gene),
                value: intercept,
            });
        }
        intercepts.push(intercept);
    }
    Ok(intercepts)
}

fn validate_rlog_sample_prior_variance(sample_prior_variance: f64) -> Result<(), DeseqError> {
    if !sample_prior_variance.is_finite() || sample_prior_variance <= 0.0 {
        return Err(DeseqError::InvalidOptions {
            reason: "rlog sample prior variance must be positive and finite".to_string(),
        });
    }
    Ok(())
}

fn validate_rlog_frozen_inputs(
    counts: &CountMatrix,
    normalization_factors: &RowMajorMatrix<f64>,
    dispersions: &[f64],
    sample_prior_variance: f64,
    frozen_intercept: &[f64],
) -> Result<(), DeseqError> {
    if normalization_factors.n_rows() != counts.n_genes() {
        return Err(invalid_dimensions(
            "rlog normalization factor rows",
            counts.n_genes(),
            normalization_factors.n_rows(),
        ));
    }
    if normalization_factors.n_cols() != counts.n_samples() {
        return Err(invalid_dimensions(
            "rlog normalization factor columns",
            counts.n_samples(),
            normalization_factors.n_cols(),
        ));
    }
    if dispersions.len() != counts.n_genes() {
        return Err(invalid_dimensions(
            "rlog dispersions",
            counts.n_genes(),
            dispersions.len(),
        ));
    }
    if frozen_intercept.len() != counts.n_genes() {
        return Err(invalid_dimensions(
            "rlog frozen intercept rows",
            counts.n_genes(),
            frozen_intercept.len(),
        ));
    }
    validate_rlog_sample_prior_variance(sample_prior_variance)?;
    for (gene, intercept) in frozen_intercept.iter().enumerate() {
        if !intercept.is_finite() {
            return Err(DeseqError::NonFiniteValue {
                context: "rlog frozen intercept".to_string(),
                index: Some(gene),
                value: *intercept,
            });
        }
    }
    Ok(())
}

fn rlog_normalization_factors_from_size_factors(
    counts: &CountMatrix,
    size_factors: &[f64],
) -> Result<RowMajorMatrix<f64>, DeseqError> {
    let mut values = Vec::with_capacity(counts.n_genes() * counts.n_samples());
    for _gene in 0..counts.n_genes() {
        for (sample, size_factor) in size_factors.iter().enumerate() {
            if !size_factor.is_finite() || *size_factor <= 0.0 {
                return Err(DeseqError::NonFiniteValue {
                    context: "rlog size factor".to_string(),
                    index: Some(sample),
                    value: *size_factor,
                });
            }
            values.push(*size_factor);
        }
    }
    RowMajorMatrix::from_row_major(counts.n_genes(), counts.n_samples(), values)
}

fn rlog_frozen_normalization_factors(
    normalization_factors: &RowMajorMatrix<f64>,
    frozen_intercept: &[f64],
) -> Result<RowMajorMatrix<f64>, DeseqError> {
    let mut values = Vec::with_capacity(normalization_factors.len());
    for (gene, intercept) in frozen_intercept.iter().enumerate() {
        let intercept_factor = checked_intercept_factor(*intercept, gene)?;
        for sample in 0..normalization_factors.n_cols() {
            let base_factor = normalization_factors
                .get(gene, sample)
                .copied()
                .unwrap_or(f64::NAN);
            if !base_factor.is_finite() || base_factor <= 0.0 {
                return Err(DeseqError::NonFiniteValue {
                    context: "rlog normalization factor".to_string(),
                    index: Some(gene * normalization_factors.n_cols() + sample),
                    value: base_factor,
                });
            }
            let factor = base_factor * intercept_factor;
            if !factor.is_finite() || factor <= 0.0 {
                return Err(DeseqError::NonFiniteValue {
                    context: "rlog frozen normalization factor".to_string(),
                    index: Some(gene * normalization_factors.n_cols() + sample),
                    value: factor,
                });
            }
            values.push(factor);
        }
    }
    RowMajorMatrix::from_row_major(
        normalization_factors.n_rows(),
        normalization_factors.n_cols(),
        values,
    )
}

fn checked_intercept_factor(intercept: f64, gene: usize) -> Result<f64, DeseqError> {
    if !intercept.is_finite() {
        return Err(DeseqError::NonFiniteValue {
            context: "rlog frozen intercept".to_string(),
            index: Some(gene),
            value: intercept,
        });
    }
    let factor = 2.0_f64.powf(intercept);
    if !factor.is_finite() || factor <= 0.0 {
        return Err(DeseqError::NonFiniteValue {
            context: "rlog frozen intercept factor".to_string(),
            index: Some(gene),
            value: factor,
        });
    }
    Ok(factor)
}
