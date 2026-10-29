//! Local regression trend model, fitting, and validation.

use super::numerics::validate_positive_mean;
use super::selection::{local_trend_use_for_fit, validate_trend_fit_rows};
use crate::errors::{DeseqError, invalid_dimensions};
use rcompat_locfit::{Kernel, LocalFit, LocalRegressionConfig, PredictionMethod};

/// Options for the pure-Rust local log-dispersion trend.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LocalDispersionTrendOptions {
    /// Minimum dispersion used by DESeq2's local trend floor and fit rule.
    pub min_disp: f64,
    /// Fraction of fitted rows used in the adaptive local neighborhood.
    pub span: f64,
    /// Local polynomial degree on log(mean), supported through quadratic.
    pub degree: usize,
}

impl Default for LocalDispersionTrendOptions {
    fn default() -> Self {
        Self {
            min_disp: 1e-8,
            span: 0.7,
            degree: 2,
        }
    }
}

/// Local dispersion trend fit on log mean/log dispersion scale.
///
/// DESeq2 uses the R `locfit` package for `fitType="local"`. This pure-Rust
/// representation keeps the same fitted rows and base-mean weights, then
/// evaluates a locfit-compatible local polynomial smoother.
#[derive(Clone, Debug, PartialEq)]
pub struct LocalDispersionTrend {
    /// DESeq2 minimum dispersion used for the all-near-minimum fallback.
    pub min_disp: f64,
    /// Adaptive neighborhood fraction.
    pub span: f64,
    /// Local polynomial degree.
    pub degree: usize,
    /// Sorted log base means used for fitting.
    pub log_means: Vec<f64>,
    /// Log gene-wise dispersions used for fitting, aligned to `log_means`.
    pub log_disps: Vec<f64>,
    /// Base-mean weights used for local regression, aligned to `log_means`.
    pub weights: Vec<f64>,
}

impl LocalDispersionTrend {
    /// Evaluate the local dispersion trend at one positive mean.
    pub fn evaluate(&self, mean: f64) -> Result<f64, DeseqError> {
        validate_positive_mean(mean, None)?;
        self.validate_fit()?;
        self.evaluate_validated(mean)
    }

    fn evaluate_validated(&self, mean: f64) -> Result<f64, DeseqError> {
        if self.log_means.is_empty() {
            return Ok(self.min_disp);
        }
        if self.log_means.len() == 1 {
            return self.single_row_dispersion();
        }
        let predicted = self
            .rcompat_fit()?
            .predict_one(mean.ln())
            .map_err(rcompat_local_dispersion_error)?;
        let dispersion = predicted.exp();
        if !dispersion.is_finite() || dispersion <= 0.0 {
            return Err(DeseqError::InvalidDispersion {
                reason: "local dispersion trend produced a non-finite or non-positive fitted value"
                    .to_string(),
            });
        }
        Ok(dispersion)
    }

    /// Evaluate the trend for all finite positive means, returning `NaN` for missing rows.
    pub fn evaluate_many_allow_missing(&self, means: &[f64]) -> Result<Vec<f64>, DeseqError> {
        self.validate_fit()?;
        if self.log_means.is_empty() {
            return Ok(means
                .iter()
                .copied()
                .map(|mean| {
                    if mean.is_finite() && mean > 0.0 {
                        self.min_disp
                    } else {
                        f64::NAN
                    }
                })
                .collect());
        }
        if self.log_means.len() == 1 {
            let dispersion = self.single_row_dispersion()?;
            return Ok(means
                .iter()
                .copied()
                .map(|mean| {
                    if mean.is_finite() && mean > 0.0 {
                        dispersion
                    } else {
                        f64::NAN
                    }
                })
                .collect());
        }
        let trend = self.rcompat_fit()?;
        means
            .iter()
            .copied()
            .map(|mean| {
                if mean.is_finite() && mean > 0.0 {
                    let predicted = trend
                        .predict_one(mean.ln())
                        .map_err(rcompat_local_dispersion_error)?;
                    let dispersion = predicted.exp();
                    if !dispersion.is_finite() || dispersion <= 0.0 {
                        return Err(DeseqError::InvalidDispersion {
                            reason: "local dispersion trend produced a non-finite or non-positive fitted value"
                                .to_string(),
                        });
                    }
                    Ok(dispersion)
                } else {
                    Ok(f64::NAN)
                }
            })
            .collect()
    }

    fn single_row_dispersion(&self) -> Result<f64, DeseqError> {
        let dispersion = self.log_disps[0].exp();
        if !dispersion.is_finite() || dispersion <= 0.0 {
            return Err(DeseqError::InvalidDispersion {
                reason: "local dispersion trend produced a non-finite or non-positive fitted value"
                    .to_string(),
            });
        }
        Ok(dispersion)
    }

    fn validate_fit(&self) -> Result<(), DeseqError> {
        if !self.min_disp.is_finite() || self.min_disp <= 0.0 {
            return Err(DeseqError::InvalidDispersion {
                reason: "minimum dispersion must be finite and positive".to_string(),
            });
        }
        validate_local_trend_shape(self.span, self.degree)?;
        if self.log_means.is_empty() && (!self.log_disps.is_empty() || !self.weights.is_empty()) {
            return Err(DeseqError::InvalidDispersion {
                reason: "empty local dispersion trend must not carry fit values or weights"
                    .to_string(),
            });
        }
        if !self.log_means.is_empty() {
            validate_local_trend_fit_rows(&self.log_means, &self.log_disps, &self.weights)?;
        }
        Ok(())
    }

    fn rcompat_fit(&self) -> Result<LocalFit, DeseqError> {
        let degree = self.degree.min(self.log_means.len().saturating_sub(1));
        let min_points = (degree + 1).min(self.log_means.len()).max(1);
        // R locfit predicts from an interpolated evaluation grid. Use the
        // Hermite-backed locfit approximation for ordinary DESeq2 fits once
        // there are enough points to support its local quadratic evaluations;
        // keep direct prediction for non-default options and tiny edge cases.
        let prediction_method = if self.span == LocalDispersionTrendOptions::default().span
            && self.degree == LocalDispersionTrendOptions::default().degree
            && self.log_means.len() >= 4
        {
            PredictionMethod::LocfitHermiteApprox
        } else {
            PredictionMethod::Direct
        };
        LocalFit::fit(
            &self.log_means,
            &self.log_disps,
            Some(&self.weights),
            LocalRegressionConfig {
                alpha: self.span,
                degree,
                kernel: Kernel::Tricube,
                min_points,
                allow_degree_downgrade: true,
                prediction_method,
            },
        )
        .map_err(rcompat_local_dispersion_error)
    }
}

/// Output from a local dispersion trend fit.
#[derive(Clone, Debug, PartialEq)]
pub struct LocalDispersionTrendFit {
    /// Fitted local trend.
    pub trend: LocalDispersionTrend,
    /// Fitted dispersion values for every input row; missing rows are `NaN`.
    pub disp_fit: Vec<f64>,
    /// Rows used by DESeq2's local rule, `dispGeneEst >= 10 * minDisp`.
    pub use_for_fit: Vec<bool>,
    /// Number of rows retained by `use_for_fit`.
    pub genes_used: usize,
    /// Whether the DESeq2 all-near-minimum fallback was used.
    pub used_min_disp_floor: bool,
}

/// Fit DESeq2's `fitType="local"` dispersion trend.
///
/// DESeq2 fits local regression on `log(dispersion)` versus `log(mean)` after
/// selecting rows with `dispGeneEst >= 10 * minDisp`, using base means as
/// weights. This implementation keeps that data contract and evaluates a
/// pure-Rust locfit-compatible local polynomial smoother.
pub fn fit_local_dispersion_trend(
    base_mean: &[f64],
    disp_gene_est: &[f64],
    options: LocalDispersionTrendOptions,
) -> Result<LocalDispersionTrendFit, DeseqError> {
    validate_local_trend_inputs(base_mean, disp_gene_est, options)?;
    let use_for_fit = local_trend_use_for_fit(base_mean, disp_gene_est, options.min_disp)?;
    let genes_used = use_for_fit.iter().filter(|value| **value).count();
    let used_min_disp_floor = genes_used == 0;

    let mut local_rows = Vec::with_capacity(genes_used);
    if !used_min_disp_floor {
        for ((mean, disp), used) in base_mean
            .iter()
            .copied()
            .zip(disp_gene_est.iter().copied())
            .zip(use_for_fit.iter().copied())
        {
            if used {
                local_rows.push((mean.ln(), disp.ln(), mean));
            }
        }
    }
    local_rows.sort_by(|left, right| left.0.total_cmp(&right.0));
    let mut log_means = Vec::with_capacity(local_rows.len());
    let mut log_disps = Vec::with_capacity(local_rows.len());
    let mut weights = Vec::with_capacity(local_rows.len());
    for (log_mean, log_disp, weight) in local_rows {
        log_means.push(log_mean);
        log_disps.push(log_disp);
        weights.push(weight);
    }

    let trend = LocalDispersionTrend {
        min_disp: options.min_disp,
        span: options.span,
        degree: options.degree,
        log_means,
        log_disps,
        weights,
    };

    Ok(LocalDispersionTrendFit {
        disp_fit: trend.evaluate_many_allow_missing(base_mean)?,
        trend,
        use_for_fit,
        genes_used,
        used_min_disp_floor,
    })
}

fn rcompat_local_dispersion_error(error: rcompat_locfit::LocfitError) -> DeseqError {
    DeseqError::InvalidDispersion {
        reason: format!("local dispersion locfit-compatible smoother failed: {error}"),
    }
}

fn validate_local_trend_fit_rows(
    xs: &[f64],
    ys: &[f64],
    weights: &[f64],
) -> Result<(), DeseqError> {
    if xs.len() != ys.len() {
        return Err(invalid_dimensions(
            "local dispersion trend fit rows",
            xs.len(),
            ys.len(),
        ));
    }
    if xs.len() != weights.len() {
        return Err(invalid_dimensions(
            "local dispersion trend fit weights",
            xs.len(),
            weights.len(),
        ));
    }
    validate_sorted_local_means(xs)?;
    validate_local_trend_values(ys)?;
    validate_local_trend_weights(weights)
}

fn validate_sorted_local_means(xs: &[f64]) -> Result<(), DeseqError> {
    if xs.iter().any(|value| !value.is_finite()) {
        return Err(DeseqError::InvalidDispersion {
            reason: "local dispersion trend log means must be finite".to_string(),
        });
    }
    if xs.windows(2).any(|window| window[0] > window[1]) {
        return Err(DeseqError::InvalidDispersion {
            reason: "local dispersion trend log means must be sorted".to_string(),
        });
    }
    Ok(())
}

fn validate_local_trend_values(ys: &[f64]) -> Result<(), DeseqError> {
    if ys.iter().any(|value| !value.is_finite()) {
        return Err(DeseqError::InvalidDispersion {
            reason: "local dispersion trend log dispersions must be finite".to_string(),
        });
    }
    Ok(())
}

fn validate_local_trend_weights(weights: &[f64]) -> Result<(), DeseqError> {
    if weights
        .iter()
        .any(|value| !value.is_finite() || *value <= 0.0)
    {
        return Err(DeseqError::InvalidDispersion {
            reason: "local dispersion trend weights must be finite and positive".to_string(),
        });
    }
    Ok(())
}

fn validate_local_trend_shape(span: f64, degree: usize) -> Result<(), DeseqError> {
    if !span.is_finite() || span <= 0.0 || span > 1.0 {
        return Err(DeseqError::InvalidDispersion {
            reason: "local dispersion span must be finite and in (0, 1]".to_string(),
        });
    }
    if degree > 2 {
        return Err(DeseqError::InvalidDispersion {
            reason: "local dispersion polynomial degree must be 0, 1, or 2".to_string(),
        });
    }
    Ok(())
}

fn validate_local_trend_inputs(
    base_mean: &[f64],
    disp_gene_est: &[f64],
    options: LocalDispersionTrendOptions,
) -> Result<(), DeseqError> {
    validate_trend_fit_rows(
        base_mean,
        disp_gene_est,
        options.min_disp,
        "local dispersion trend rows",
    )?;
    validate_local_trend_shape(options.span, options.degree)
}
