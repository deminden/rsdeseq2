//! Dispersion trend dispatch and row-selection rules.

use super::local::{LocalDispersionTrendOptions, fit_local_dispersion_trend};
use super::mean::{MeanDispersionTrendOptions, fit_mean_dispersion_trend};
use super::model::DispersionTrendFit;
use super::parametric::{ParametricDispersionTrendOptions, fit_parametric_dispersion_trend};
use crate::errors::{DeseqError, invalid_dimensions};
use crate::options::FitType;

/// DESeq2's parametric trend row-selection rule.
pub fn parametric_trend_use_for_fit(
    base_mean: &[f64],
    disp_gene_est: &[f64],
    min_disp: f64,
) -> Result<Vec<bool>, DeseqError> {
    validate_trend_selection_inputs(
        base_mean,
        disp_gene_est,
        min_disp,
        "parametric dispersion trend rows",
    )?;
    Ok(base_mean
        .iter()
        .copied()
        .zip(disp_gene_est.iter().copied())
        .map(|(mean, disp)| {
            mean.is_finite() && mean > 0.0 && disp.is_finite() && disp > 100.0 * min_disp
        })
        .collect())
}

/// DESeq2's `fitType="local"` row-selection rule for local regression.
pub fn local_trend_use_for_fit(
    base_mean: &[f64],
    disp_gene_est: &[f64],
    min_disp: f64,
) -> Result<Vec<bool>, DeseqError> {
    validate_trend_selection_inputs(
        base_mean,
        disp_gene_est,
        min_disp,
        "local dispersion trend rows",
    )?;
    Ok(base_mean
        .iter()
        .copied()
        .zip(disp_gene_est.iter().copied())
        .map(|(mean, disp)| {
            mean.is_finite() && mean > 0.0 && disp.is_finite() && disp >= 10.0 * min_disp
        })
        .collect())
}

/// DESeq2's `fitType="mean"` preliminary viability rule.
pub fn mean_trend_use_for_fit(
    base_mean: &[f64],
    disp_gene_est: &[f64],
    min_disp: f64,
) -> Result<Vec<bool>, DeseqError> {
    validate_trend_selection_inputs(
        base_mean,
        disp_gene_est,
        min_disp,
        "mean dispersion trend fit rows",
    )?;
    Ok(base_mean
        .iter()
        .copied()
        .zip(disp_gene_est.iter().copied())
        .map(|(mean, disp)| {
            mean.is_finite() && mean > 0.0 && disp.is_finite() && disp > 100.0 * min_disp
        })
        .collect())
}

/// DESeq2's `fitType="mean"` row-selection rule for the trimmed mean itself.
pub fn mean_trend_use_for_mean(
    base_mean: &[f64],
    disp_gene_est: &[f64],
    min_disp: f64,
) -> Result<Vec<bool>, DeseqError> {
    validate_trend_selection_inputs(
        base_mean,
        disp_gene_est,
        min_disp,
        "mean dispersion trend rows",
    )?;
    Ok(base_mean
        .iter()
        .copied()
        .zip(disp_gene_est.iter().copied())
        .map(|(mean, disp)| {
            mean.is_finite() && mean > 0.0 && disp.is_finite() && disp > 10.0 * min_disp
        })
        .collect())
}

fn validate_trend_selection_inputs(
    base_mean: &[f64],
    disp_gene_est: &[f64],
    min_disp: f64,
    context: &str,
) -> Result<(), DeseqError> {
    if base_mean.len() != disp_gene_est.len() {
        return Err(invalid_dimensions(
            context,
            base_mean.len(),
            disp_gene_est.len(),
        ));
    }
    if !min_disp.is_finite() || min_disp <= 0.0 {
        return Err(DeseqError::InvalidDispersion {
            reason: "min_disp must be finite and positive".to_string(),
        });
    }
    Ok(())
}

/// Fit a dispersion trend using DESeq2-style `fitType` defaults.
///
/// Use the type-specific fitters when custom trend options are needed.
pub fn fit_dispersion_trend(
    base_mean: &[f64],
    disp_gene_est: &[f64],
    fit_type: FitType,
) -> Result<DispersionTrendFit, DeseqError> {
    match fit_type {
        FitType::Parametric => Ok(DispersionTrendFit::Parametric(
            fit_parametric_dispersion_trend(
                base_mean,
                disp_gene_est,
                ParametricDispersionTrendOptions::default(),
            )?,
        )),
        FitType::Local => Ok(DispersionTrendFit::Local(fit_local_dispersion_trend(
            base_mean,
            disp_gene_est,
            LocalDispersionTrendOptions::default(),
        )?)),
        FitType::Mean => Ok(DispersionTrendFit::Mean(fit_mean_dispersion_trend(
            base_mean,
            disp_gene_est,
            MeanDispersionTrendOptions::default(),
        )?)),
        FitType::GlmGamPoi => Err(DeseqError::UnsupportedFeature {
            feature: "glmGamPoi dispersion trend fitting".to_string(),
        }),
    }
}

pub(super) fn validate_trend_fit_rows(
    base_mean: &[f64],
    disp_gene_est: &[f64],
    min_disp: f64,
    context: &str,
) -> Result<(), DeseqError> {
    validate_trend_selection_inputs(base_mean, disp_gene_est, min_disp, context)?;
    if base_mean.is_empty() {
        return Err(DeseqError::InvalidDimensions {
            context: context.to_string(),
            expected: 1,
            actual: 0,
        });
    }
    Ok(())
}
