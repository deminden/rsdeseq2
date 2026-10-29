//! Constant mean trend model and trimmed-mean fitting.

use super::numerics::{stable_midpoint, validate_positive_mean};
use super::selection::{mean_trend_use_for_fit, mean_trend_use_for_mean, validate_trend_fit_rows};
use crate::errors::DeseqError;

/// Options for DESeq2-style mean dispersion trend fitting.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MeanDispersionTrendOptions {
    /// Minimum dispersion used to select rows for trend fitting.
    pub min_disp: f64,
    /// Fraction trimmed from each tail before averaging dispersion estimates.
    pub trim: f64,
}

impl Default for MeanDispersionTrendOptions {
    fn default() -> Self {
        Self {
            min_disp: 1e-8,
            trim: 0.001,
        }
    }
}

/// Constant dispersion trend used by DESeq2's `fitType="mean"` path.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MeanDispersionTrend {
    /// Trimmed mean of usable gene-wise dispersion estimates.
    pub mean_disp: f64,
}

impl MeanDispersionTrend {
    /// Evaluate the constant mean dispersion trend at one positive mean.
    pub fn evaluate(&self, mean: f64) -> Result<f64, DeseqError> {
        validate_positive_mean(mean, None)?;
        self.validate_fit()?;
        Ok(self.evaluate_validated())
    }

    /// Evaluate the trend for all finite positive means, returning `NaN` for missing rows.
    pub fn evaluate_many_allow_missing(&self, means: &[f64]) -> Result<Vec<f64>, DeseqError> {
        self.validate_fit()?;
        means
            .iter()
            .copied()
            .map(|mean| {
                if mean.is_finite() && mean > 0.0 {
                    Ok(self.evaluate_validated())
                } else {
                    Ok(f64::NAN)
                }
            })
            .collect()
    }

    fn evaluate_validated(&self) -> f64 {
        self.mean_disp
    }

    fn validate_fit(&self) -> Result<(), DeseqError> {
        if !self.mean_disp.is_finite() || self.mean_disp <= 0.0 {
            return Err(DeseqError::InvalidDispersion {
                reason: "mean dispersion trend value must be finite and positive".to_string(),
            });
        }
        Ok(())
    }
}

/// Output from a mean dispersion trend fit.
#[derive(Clone, Debug, PartialEq)]
pub struct MeanDispersionTrendFit {
    /// Fitted constant trend.
    pub trend: MeanDispersionTrend,
    /// Fitted dispersion values for every input row; missing rows are `NaN`.
    pub disp_fit: Vec<f64>,
    /// Rows passing DESeq2's preliminary `dispGeneEst > 100 * minDisp` fit rule.
    pub use_for_fit: Vec<bool>,
    /// Rows used by DESeq2's mean rule, `dispGeneEst > 10 * minDisp`.
    pub use_for_mean: Vec<bool>,
    /// Number of rows retained by `use_for_fit`.
    pub genes_used: usize,
    /// Number of rows used for the trimmed mean.
    pub genes_used_for_mean: usize,
}

/// Fit DESeq2's `fitType="mean"` dispersion trend.
///
/// The compatibility behavior follows `estimateDispersionsFit`: first require
/// at least one non-all-zero row with `dispGeneEst > 100 * minDisp`, then compute
/// the constant fitted dispersion as `mean(dispGeneEst[dispGeneEst > 10 *
/// minDisp], trim = 0.001, na.rm = TRUE)`.
pub fn fit_mean_dispersion_trend(
    base_mean: &[f64],
    disp_gene_est: &[f64],
    options: MeanDispersionTrendOptions,
) -> Result<MeanDispersionTrendFit, DeseqError> {
    validate_mean_trend_inputs(base_mean, disp_gene_est, options)?;
    let use_for_fit = mean_trend_use_for_fit(base_mean, disp_gene_est, options.min_disp)?;
    let genes_used = use_for_fit.iter().filter(|value| **value).count();
    if genes_used == 0 {
        return Err(DeseqError::InvalidDispersion {
            reason: "all gene-wise dispersion estimates are within 2 orders of magnitude from the minimum value".to_string(),
        });
    }

    let use_for_mean = mean_trend_use_for_mean(base_mean, disp_gene_est, options.min_disp)?;
    let mean_disps = disp_gene_est
        .iter()
        .copied()
        .zip(use_for_mean.iter().copied())
        .filter_map(|(disp, used)| used.then_some(disp))
        .collect::<Vec<_>>();
    let genes_used_for_mean = mean_disps.len();
    let mean_disp = trimmed_mean(&mean_disps, options.trim)?;
    let trend = MeanDispersionTrend { mean_disp };

    Ok(MeanDispersionTrendFit {
        trend,
        disp_fit: trend.evaluate_many_allow_missing(base_mean)?,
        use_for_fit,
        use_for_mean,
        genes_used,
        genes_used_for_mean,
    })
}

fn validate_mean_trend_inputs(
    base_mean: &[f64],
    disp_gene_est: &[f64],
    options: MeanDispersionTrendOptions,
) -> Result<(), DeseqError> {
    validate_trend_fit_rows(
        base_mean,
        disp_gene_est,
        options.min_disp,
        "mean dispersion trend rows",
    )?;
    if !options.trim.is_finite() || !(0.0..=0.5).contains(&options.trim) {
        return Err(DeseqError::InvalidDispersion {
            reason: "mean dispersion trim must be finite and between 0 and 0.5".to_string(),
        });
    }
    Ok(())
}

fn trimmed_mean(values: &[f64], trim: f64) -> Result<f64, DeseqError> {
    if values.is_empty() {
        return Err(DeseqError::InvalidDispersion {
            reason: "no gene-wise dispersion estimates are usable for mean trend fitting"
                .to_string(),
        });
    }
    let mut sorted = values.to_vec();
    sorted.sort_by(|left, right| left.total_cmp(right));
    if trim == 0.5 {
        let mid = sorted.len() / 2;
        return if sorted.len().is_multiple_of(2) {
            Ok(stable_midpoint(sorted[mid - 1], sorted[mid]))
        } else {
            Ok(sorted[mid])
        };
    }
    let drop_each_tail = (sorted.len() as f64 * trim).floor() as usize;
    let start = drop_each_tail;
    let end = sorted.len().saturating_sub(drop_each_tail);
    if start >= end {
        return Err(DeseqError::InvalidDispersion {
            reason: "mean dispersion trim removed all usable estimates".to_string(),
        });
    }
    Ok(stable_mean(&sorted[start..end]))
}

fn stable_mean(values: &[f64]) -> f64 {
    let mut mean = 0.0;
    for (idx, value) in values.iter().copied().enumerate() {
        mean += (value - mean) / (idx + 1) as f64;
    }
    mean
}
