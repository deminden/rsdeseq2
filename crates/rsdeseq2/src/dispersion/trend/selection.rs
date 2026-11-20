//! Dispersion trend dispatch and row-selection rules.

use super::local::{LocalDispersionTrendOptions, fit_local_dispersion_trend};
use super::mean::{MeanDispersionTrendOptions, fit_mean_dispersion_trend};
use super::model::DispersionTrendFit;
use super::parametric::{ParametricDispersionTrendOptions, fit_parametric_dispersion_trend};
use crate::errors::{DeseqError, ErrorKind, invalid_dimensions};
use crate::errors::{Requirement, invalid_parameter};
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
        return Err(invalid_parameter(
            ErrorKind::Dispersion,
            "min_disp",
            Requirement::FiniteAndPositive,
        ));
    }
    Ok(())
}

/// Fit a dispersion trend using DESeq2-style `fitType` defaults.
///
/// This follows DESeq2's `estimateDispersionsFit` (see
/// [`estimate_dispersion_trend`]) with `minDisp = 1e-8`. Use the type-specific
/// fitters when custom trend options are needed.
pub fn fit_dispersion_trend(
    base_mean: &[f64],
    disp_gene_est: &[f64],
    fit_type: FitType,
) -> Result<DispersionTrendFit, DeseqError> {
    estimate_dispersion_trend(base_mean, disp_gene_est, fit_type, DEFAULT_MIN_DISP)
}

const DEFAULT_MIN_DISP: f64 = 1e-8;

/// Fit the trend the way DESeq2 1.52 `estimateDispersionsFit` does.
///
/// - Only rows with `dispGeneEst > 100 * minDisp` feed any trend; when none
///   remain the call fails, before any trend is attempted.
/// - `fitType="local"` sees only those rows, so `localDispersionFit`'s own
///   `>= 10 * minDisp` rule never widens the set.
/// - A parametric fit that fails is replaced by the local fit, as DESeq2's
///   `try(parametricDispersionFit(...))` does. The returned variant records
///   which trend was fitted. Only dispersion and numeric failures fall back;
///   shape and option errors indicate a caller bug and are returned.
/// - `fitType="mean"` keeps its own `> 10 * minDisp` rows for the trimmed mean.
pub fn estimate_dispersion_trend(
    base_mean: &[f64],
    disp_gene_est: &[f64],
    fit_type: FitType,
    min_disp: f64,
) -> Result<DispersionTrendFit, DeseqError> {
    let use_for_fit = parametric_trend_use_for_fit(base_mean, disp_gene_est, min_disp)?;
    if !use_for_fit.iter().any(|used| *used) {
        return Err(DeseqError::NoUsableDispersionTrendRows {
            min_disp,
            genes_checked: base_mean.len(),
        });
    }
    match fit_type {
        FitType::Parametric => {
            let options = ParametricDispersionTrendOptions {
                min_disp,
                ..ParametricDispersionTrendOptions::default()
            };
            match fit_parametric_dispersion_trend(base_mean, disp_gene_est, options) {
                Ok(fit) => Ok(DispersionTrendFit::Parametric(fit)),
                Err(error)
                    if matches!(error.kind(), ErrorKind::Dispersion | ErrorKind::Numeric) =>
                {
                    fit_masked_local_trend(base_mean, disp_gene_est, &use_for_fit, min_disp)
                }
                Err(error) => Err(error),
            }
        }
        FitType::Local => fit_masked_local_trend(base_mean, disp_gene_est, &use_for_fit, min_disp),
        FitType::Mean => Ok(DispersionTrendFit::Mean(fit_mean_dispersion_trend(
            base_mean,
            disp_gene_est,
            MeanDispersionTrendOptions {
                min_disp,
                ..MeanDispersionTrendOptions::default()
            },
        )?)),
    }
}

/// Local trend on DESeq2's `useForFit` rows only; other rows are passed as
/// missing so they are evaluated but not fitted.
fn fit_masked_local_trend(
    base_mean: &[f64],
    disp_gene_est: &[f64],
    use_for_fit: &[bool],
    min_disp: f64,
) -> Result<DispersionTrendFit, DeseqError> {
    let masked = disp_gene_est
        .iter()
        .zip(use_for_fit)
        .map(|(disp, used)| if *used { *disp } else { f64::NAN })
        .collect::<Vec<_>>();
    Ok(DispersionTrendFit::Local(fit_local_dispersion_trend(
        base_mean,
        &masked,
        LocalDispersionTrendOptions {
            min_disp,
            ..LocalDispersionTrendOptions::default()
        },
    )?))
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
