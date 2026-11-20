use crate::design::DesignMatrix;
use crate::errors::{DeseqError, invalid_dimensions};
use crate::errors::{ErrorKind, missing_prerequisite};
use crate::glm::beta::{fit_prior_for_output, match_row_weighted_upper_quantile_for_variance};
use crate::glm::irls::{IrlsOutput, OutputRequirements};
use crate::glm::{GlmInput, IrlsOptions, NbinomGlmFit};
use crate::matrix::CountMatrix;
use crate::matrix::RowMajorMatrix;
use crate::matrix::try_alloc_matrix;
use crate::normalization::{Offsets, normalized_counts_for_rows};
use crate::numeric::math::ProfileMath;
use rayon::prelude::*;
use std::borrow::Cow;

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
#[must_use]
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

fn rlog_fit_for_output(
    input: GlmInput<'_>,
    sample_prior_variance: f64,
    options: IrlsOptions,
    selected: Option<&[usize]>,
    requirements: OutputRequirements,
) -> Result<IrlsOutput, DeseqError> {
    let samples = input.counts.n_samples();
    let design = rlog_sample_design(samples)?;
    let prior = rlog_beta_prior_variance(samples, sample_prior_variance)?;
    fit_prior_for_output(
        input,
        &design,
        &prior,
        options,
        selected,
        None,
        requirements,
    )
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
    estimate_rlog_prior(
        Cow::Borrowed(normalized_counts),
        base_mean,
        disp_fit,
        None,
        upper_quantile,
    )
}

fn estimate_rlog_prior(
    normalized: Cow<'_, RowMajorMatrix<f64>>,
    base_mean: &[f64],
    disp_fit: &[f64],
    selected: Option<&[usize]>,
    upper_quantile: f64,
) -> Result<f64, DeseqError> {
    let genes = normalized.n_rows();
    let samples = normalized.n_cols();
    let length = normalized.len();
    if selected.is_none() {
        if base_mean.len() != genes {
            return Err(invalid_dimensions(
                "rlog baseMean rows",
                genes,
                base_mean.len(),
            ));
        }
        if disp_fit.len() != genes {
            return Err(invalid_dimensions(
                "rlog dispFit rows",
                genes,
                disp_fit.len(),
            ));
        }
    }
    let (normalized_rows, mut log_fold_changes) = match normalized {
        Cow::Borrowed(matrix) => {
            let mut changes = try_alloc_matrix(genes, samples, "rlog prior log fold changes")?;
            changes.resize(length, 0.0);
            (Some(matrix), changes)
        }
        Cow::Owned(matrix) => (None, matrix.into_values()),
    };
    let mut weights = try_alloc_matrix(genes, 1, "rlog prior row weights")?;
    weights.resize(genes, 0.0);
    let fill_row = |gene: usize, changes: &mut [f64], row_weight: &mut f64| {
        let source = selected.map_or(gene, |selected| selected[gene]);
        let stat = |values: &[f64]| {
            values.get(source).copied().ok_or_else(|| {
                invalid_dimensions("rlog prior source rows", source + 1, values.len())
            })
        };
        let mean = stat(base_mean)?;
        let dispersion = stat(disp_fit)?;
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
        let mean_log = (mean + 0.5).log2_p();
        let normalized_row = normalized_rows.map(|matrix| matrix.row(gene)).transpose()?;
        for (sample, change) in changes.iter_mut().enumerate() {
            let value = normalized_row.map_or(*change, |row| row[sample]);
            if !value.is_finite() || value < 0.0 {
                return Err(DeseqError::NonFiniteValue {
                    context: "rlog normalized count".to_string(),
                    index: Some(gene * samples + sample),
                    value,
                });
            }
            *change = (value + 0.5).log2_p() - mean_log;
        }
        *row_weight = weight;
        Ok(())
    };
    let rows = log_fold_changes
        .par_chunks_mut(samples)
        .zip(weights.par_iter_mut())
        .enumerate()
        .map(|(gene, (changes, weights))| fill_row(gene, changes, weights))
        .collect::<Vec<_>>();
    // Report the first row error in the same order as the serial loop.
    rows.into_iter().collect::<Result<(), DeseqError>>()?;
    match_row_weighted_upper_quantile_for_variance(
        &log_fold_changes,
        &weights,
        samples,
        upper_quantile,
    )
}

fn rlog_coefficient_count(n_samples: usize) -> Result<usize, DeseqError> {
    n_samples.checked_add(1).ok_or_else(|| {
        invalid_dimensions("rlog coefficient count overflow", usize::MAX - 1, n_samples)
    })
}

/// Build DESeq2's rlog sample design: an unpenalized intercept plus one
/// regularized indicator coefficient per sample.
pub fn rlog_sample_design(n_samples: usize) -> Result<DesignMatrix, DeseqError> {
    if n_samples == 0 {
        return Err(missing_prerequisite(
            ErrorKind::Options,
            "rlog",
            "at least one sample",
        ));
    }
    let n_coefficients = rlog_coefficient_count(n_samples)?;
    let mut values = try_alloc_matrix(n_samples, n_coefficients, "rlog sample design")?;
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
        return Err(missing_prerequisite(
            ErrorKind::Options,
            "rlog",
            "at least one sample",
        ));
    }
    let mut values = try_alloc_matrix(n_samples, n_samples, "rlog sample effect design")?;
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
        return Err(missing_prerequisite(
            ErrorKind::Options,
            "rlog",
            "at least one sample",
        ));
    }
    validate_rlog_sample_prior_variance(sample_prior_variance)?;
    let mut beta_prior_variance =
        try_alloc_matrix(1, rlog_coefficient_count(n_samples)?, "rlog prior variance")?;
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
        return Err(missing_prerequisite(
            ErrorKind::Options,
            "rlog",
            "at least one sample",
        ));
    }
    validate_rlog_sample_prior_variance(sample_prior_variance)?;
    let mut variances = try_alloc_matrix(1, n_samples, "rlog sample effect prior variance")?;
    variances.resize(n_samples, sample_prior_variance);
    Ok(variances)
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
    let mut values = try_alloc_matrix(beta.n_rows(), n_samples, "rlog transformed values")?;
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
    let mut values = try_alloc_matrix(beta.n_rows(), n_samples, "rlog transformed values")?;
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
    offsets: Offsets<'_>,
    dispersions: &[f64],
    sample_prior_variance: f64,
    frozen_intercept: &[f64],
    selected: Option<&[usize]>,
) -> Result<(), DeseqError> {
    match offsets {
        Offsets::PerSample(factors) => {
            if factors.len() != counts.n_samples() {
                return Err(invalid_dimensions(
                    "rlog size factors",
                    counts.n_samples(),
                    factors.len(),
                ));
            }
            for (sample, &value) in factors.iter().enumerate() {
                if !value.is_finite() || value <= 0.0 {
                    return Err(DeseqError::NonFiniteValue {
                        context: "rlog size factor".to_string(),
                        index: Some(sample),
                        value,
                    });
                }
            }
        }
        Offsets::PerCell(factors) => {
            if factors.n_rows() != counts.n_genes() {
                return Err(invalid_dimensions(
                    "rlog normalization factor rows",
                    counts.n_genes(),
                    factors.n_rows(),
                ));
            }
            if factors.n_cols() != counts.n_samples() {
                return Err(invalid_dimensions(
                    "rlog normalization factor columns",
                    counts.n_samples(),
                    factors.n_cols(),
                ));
            }
        }
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
    let genes = selected.map_or(counts.n_genes(), <[usize]>::len);
    for gene in 0..genes {
        let source = selected.map_or(gene, |selected| selected[gene]);
        let intercept = frozen_intercept[source];
        if !intercept.is_finite() {
            return Err(DeseqError::NonFiniteValue {
                context: "rlog frozen intercept".to_string(),
                index: Some(gene),
                value: intercept,
            });
        }
    }
    Ok(())
}

fn rlog_frozen_offset_scales(
    counts: &CountMatrix,
    offsets: Offsets<'_>,
    intercepts: &[f64],
    selected: Option<&[usize]>,
) -> Result<Vec<f64>, DeseqError> {
    let samples = counts.n_samples();
    let mut scales = try_alloc_matrix(intercepts.len(), 1, "rlog frozen offset scales")?;
    for (gene, &intercept) in intercepts.iter().enumerate() {
        let source = selected.map_or(gene, |selected| selected[gene]);
        let intercept_factor = checked_intercept_factor(intercept, gene)?;
        for (sample, &base_factor) in offsets.row(source)?.iter().enumerate() {
            if !base_factor.is_finite() || base_factor <= 0.0 {
                return Err(DeseqError::NonFiniteValue {
                    context: "rlog normalization factor".to_string(),
                    index: Some(gene * samples + sample),
                    value: base_factor,
                });
            }
            let factor = base_factor * intercept_factor;
            if !factor.is_finite() || factor <= 0.0 {
                return Err(DeseqError::NonFiniteValue {
                    context: "rlog frozen normalization factor".to_string(),
                    index: Some(gene * samples + sample),
                    value: factor,
                });
            }
        }
        scales.push(intercept_factor);
    }
    Ok(scales)
}

fn offset_mode(offsets: Offsets<'_>) -> RlogOffsetMode {
    match offsets {
        Offsets::PerSample(_) => RlogOffsetMode::SizeFactors,
        Offsets::PerCell(_) => RlogOffsetMode::NormalizationFactors,
    }
}

fn checked_intercept_factor(intercept: f64, gene: usize) -> Result<f64, DeseqError> {
    if !intercept.is_finite() {
        return Err(DeseqError::NonFiniteValue {
            context: "rlog frozen intercept".to_string(),
            index: Some(gene),
            value: intercept,
        });
    }
    let factor = 2.0_f64.powf_p(intercept);
    if !factor.is_finite() || factor <= 0.0 {
        return Err(DeseqError::NonFiniteValue {
            context: "rlog frozen intercept factor".to_string(),
            index: Some(gene),
            value: factor,
        });
    }
    Ok(factor)
}

/// Fit DESeq2's low-level rlog sample-effect model (an unpenalized intercept
/// plus one regularized effect per sample) and return the log2 values.
pub fn rlog_transform(
    counts: &CountMatrix,
    offsets: Offsets<'_>,
    dispersions: &[f64],
    sample_prior_variance: f64,
    options: IrlsOptions,
) -> Result<RowMajorMatrix<f64>, DeseqError> {
    let beta = rlog_fit_for_output(
        GlmInput::new(counts, offsets, dispersions),
        sample_prior_variance,
        options,
        None,
        OutputRequirements::Transform,
    )?
    .into_matrix()?;
    rlog_transform_from_fit(&beta, counts.n_samples())
}

/// Deprecated form of [`rlog_transform`] with `Offsets::PerSample`.
#[deprecated(
    since = "0.6.0",
    note = "use rlog_transform(counts, Offsets::PerSample(size_factors), ..)"
)]
pub fn rlog_with_size_factors(
    counts: &CountMatrix,
    size_factors: &[f64],
    dispersions: &[f64],
    sample_prior_variance: f64,
    options: IrlsOptions,
) -> Result<RowMajorMatrix<f64>, DeseqError> {
    rlog_transform(
        counts,
        Offsets::PerSample(size_factors),
        dispersions,
        sample_prior_variance,
        options,
    )
}

/// Deprecated form of [`rlog_transform`] with `Offsets::PerCell`.
#[deprecated(
    since = "0.6.0",
    note = "use rlog_transform(counts, Offsets::PerCell(normalization_factors), ..)"
)]
pub fn rlog_with_normalization_factors(
    counts: &CountMatrix,
    normalization_factors: &RowMajorMatrix<f64>,
    dispersions: &[f64],
    sample_prior_variance: f64,
    options: IrlsOptions,
) -> Result<RowMajorMatrix<f64>, DeseqError> {
    rlog_transform(
        counts,
        Offsets::PerCell(normalization_factors),
        dispersions,
        sample_prior_variance,
        options,
    )
}

/// Fit the low-level rlog sample-effect model and keep the GLM fit.
pub fn rlog_fit(
    counts: &CountMatrix,
    offsets: Offsets<'_>,
    dispersions: &[f64],
    sample_prior_variance: f64,
    options: IrlsOptions,
) -> Result<RlogFitOutput, DeseqError> {
    let fit = rlog_fit_for_output(
        GlmInput::new(counts, offsets, dispersions),
        sample_prior_variance,
        options,
        None,
        OutputRequirements::Inference,
    )?
    .into_inference()?;
    rlog_output_from_fit(fit, counts.n_samples())
}

/// Deprecated form of [`rlog_fit`] with `Offsets::PerSample`.
#[deprecated(
    since = "0.6.0",
    note = "use rlog_fit(counts, Offsets::PerSample(size_factors), ..)"
)]
pub fn rlog_fit_with_size_factors(
    counts: &CountMatrix,
    size_factors: &[f64],
    dispersions: &[f64],
    sample_prior_variance: f64,
    options: IrlsOptions,
) -> Result<RlogFitOutput, DeseqError> {
    rlog_fit(
        counts,
        Offsets::PerSample(size_factors),
        dispersions,
        sample_prior_variance,
        options,
    )
}

/// Deprecated form of [`rlog_fit`] with `Offsets::PerCell`.
#[deprecated(
    since = "0.6.0",
    note = "use rlog_fit(counts, Offsets::PerCell(normalization_factors), ..)"
)]
pub fn rlog_fit_with_normalization_factors(
    counts: &CountMatrix,
    normalization_factors: &RowMajorMatrix<f64>,
    dispersions: &[f64],
    sample_prior_variance: f64,
    options: IrlsOptions,
) -> Result<RlogFitOutput, DeseqError> {
    rlog_fit(
        counts,
        Offsets::PerCell(normalization_factors),
        dispersions,
        sample_prior_variance,
        options,
    )
}

/// Fit rlog with caller-supplied frozen intercepts: only the sample effects
/// are fit, and the output is `frozen_intercept_i + sample_effect_ij`.
pub fn rlog_frozen(
    counts: &CountMatrix,
    offsets: Offsets<'_>,
    dispersions: &[f64],
    sample_prior_variance: f64,
    frozen_intercept: &[f64],
    options: IrlsOptions,
) -> Result<RlogOutput, DeseqError> {
    rlog_frozen_for_rows(
        GlmInput::new(counts, offsets, dispersions),
        sample_prior_variance,
        frozen_intercept,
        options,
        None,
    )
}

pub(crate) fn rlog_frozen_for_rows(
    input: GlmInput<'_>,
    sample_prior_variance: f64,
    frozen_intercept: &[f64],
    options: IrlsOptions,
    selected: Option<&[usize]>,
) -> Result<RlogOutput, DeseqError> {
    validate_rlog_frozen_inputs(
        input.counts,
        input.offsets,
        input.dispersions,
        sample_prior_variance,
        frozen_intercept,
        selected,
    )?;
    let samples = input.counts.n_samples();
    let design = rlog_sample_effect_design(samples)?;
    let prior = rlog_sample_effect_prior_variance(samples, sample_prior_variance)?;
    let intercepts = match selected {
        Some(genes) => {
            let mut values = try_alloc_matrix(genes.len(), 1, "rlog frozen intercepts")?;
            values.extend(genes.iter().map(|&gene| frozen_intercept[gene]));
            Cow::Owned(values)
        }
        None => Cow::Borrowed(frozen_intercept),
    };
    let scales = rlog_frozen_offset_scales(input.counts, input.offsets, &intercepts, selected)?;
    let beta = fit_prior_for_output(
        input,
        &design,
        &prior,
        options,
        selected,
        Some(&scales),
        OutputRequirements::Transform,
    )?
    .into_matrix()?;
    Ok(RlogOutput {
        transformed: rlog_transform_from_sample_effect_fit(&beta, &intercepts)?,
        intercept: intercepts.into_owned(),
        sample_prior_variance,
        offset_mode: offset_mode(input.offsets),
    })
}

/// Deprecated form of [`rlog_frozen`] with `Offsets::PerSample`.
#[deprecated(
    since = "0.6.0",
    note = "use rlog_frozen(counts, Offsets::PerSample(size_factors), ..)"
)]
pub fn rlog_frozen_with_size_factors(
    counts: &CountMatrix,
    size_factors: &[f64],
    dispersions: &[f64],
    sample_prior_variance: f64,
    frozen_intercept: &[f64],
    options: IrlsOptions,
) -> Result<RlogOutput, DeseqError> {
    rlog_frozen(
        counts,
        Offsets::PerSample(size_factors),
        dispersions,
        sample_prior_variance,
        frozen_intercept,
        options,
    )
}

/// Deprecated form of [`rlog_frozen`] with `Offsets::PerCell`.
#[deprecated(
    since = "0.6.0",
    note = "use rlog_frozen(counts, Offsets::PerCell(normalization_factors), ..)"
)]
pub fn rlog_frozen_with_normalization_factors(
    counts: &CountMatrix,
    normalization_factors: &RowMajorMatrix<f64>,
    dispersions: &[f64],
    sample_prior_variance: f64,
    frozen_intercept: &[f64],
    options: IrlsOptions,
) -> Result<RlogOutput, DeseqError> {
    rlog_frozen(
        counts,
        Offsets::PerCell(normalization_factors),
        dispersions,
        sample_prior_variance,
        frozen_intercept,
        options,
    )
}

/// Estimate the rlog sample prior variance from normalized counts, then fit.
pub fn rlog_with_estimated_prior(
    counts: &CountMatrix,
    offsets: Offsets<'_>,
    base_mean: &[f64],
    disp_fit: &[f64],
    dispersions: &[f64],
    options: IrlsOptions,
) -> Result<RlogOutput, DeseqError> {
    rlog_with_estimated_prior_for_rows(
        GlmInput::new(counts, offsets, dispersions),
        base_mean,
        disp_fit,
        options,
        None,
    )
}

pub(crate) fn rlog_with_estimated_prior_for_rows(
    input: GlmInput<'_>,
    base_mean: &[f64],
    disp_fit: &[f64],
    options: IrlsOptions,
    selected: Option<&[usize]>,
) -> Result<RlogOutput, DeseqError> {
    let normalized = normalized_counts_for_rows(input.counts, input.offsets, selected)?;
    let sample_prior_variance = estimate_rlog_prior(
        Cow::Owned(normalized),
        base_mean,
        disp_fit,
        selected,
        RLOG_PRIOR_UPPER_QUANTILE,
    )?;
    let beta = rlog_fit_for_output(
        input,
        sample_prior_variance,
        options,
        selected,
        OutputRequirements::Transform,
    )?
    .into_matrix()?;
    let samples = input.counts.n_samples();
    Ok(RlogOutput {
        transformed: rlog_transform_from_fit(&beta, samples)?,
        intercept: rlog_intercept_from_beta(&beta, samples)?,
        sample_prior_variance,
        offset_mode: offset_mode(input.offsets),
    })
}

/// Deprecated form of [`rlog_with_estimated_prior`] with `Offsets::PerSample`.
#[deprecated(
    since = "0.6.0",
    note = "use rlog_with_estimated_prior(counts, Offsets::PerSample(size_factors), ..)"
)]
pub fn rlog_with_estimated_prior_and_size_factors(
    counts: &CountMatrix,
    size_factors: &[f64],
    base_mean: &[f64],
    disp_fit: &[f64],
    dispersions: &[f64],
    options: IrlsOptions,
) -> Result<RlogOutput, DeseqError> {
    rlog_with_estimated_prior(
        counts,
        Offsets::PerSample(size_factors),
        base_mean,
        disp_fit,
        dispersions,
        options,
    )
}

/// Deprecated form of [`rlog_with_estimated_prior`] with `Offsets::PerCell`.
#[deprecated(
    since = "0.6.0",
    note = "use rlog_with_estimated_prior(counts, Offsets::PerCell(normalization_factors), ..)"
)]
pub fn rlog_with_estimated_prior_and_normalization_factors(
    counts: &CountMatrix,
    normalization_factors: &RowMajorMatrix<f64>,
    base_mean: &[f64],
    disp_fit: &[f64],
    dispersions: &[f64],
    options: IrlsOptions,
) -> Result<RlogOutput, DeseqError> {
    rlog_with_estimated_prior(
        counts,
        Offsets::PerCell(normalization_factors),
        base_mean,
        disp_fit,
        dispersions,
        options,
    )
}
