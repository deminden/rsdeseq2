use crate::errors::{invalid_dimensions, DeseqError};

/// Options for DESeq2-style parametric dispersion trend fitting.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ParametricDispersionTrendOptions {
    /// Minimum dispersion used to select rows for trend fitting.
    pub min_disp: f64,
    /// Initial asymptotic dispersion coefficient.
    pub initial_asympt_disp: f64,
    /// Initial extra-Poisson coefficient.
    pub initial_extra_pois: f64,
    /// Minimum residual ratio retained during DESeq2's robust outer loop.
    pub min_residual: f64,
    /// Maximum residual ratio retained during DESeq2's robust outer loop.
    pub max_residual: f64,
    /// DESeq2 outer-loop coefficient convergence threshold.
    pub coefficient_tol: f64,
    /// Inner Gamma identity-link IRLS convergence threshold.
    pub glm_tol: f64,
    /// Maximum DESeq2 robust outer-loop iterations.
    pub max_outer_iter: usize,
    /// Maximum inner Gamma identity-link IRLS iterations.
    pub max_irls_iter: usize,
}

impl Default for ParametricDispersionTrendOptions {
    fn default() -> Self {
        Self {
            min_disp: 1e-8,
            initial_asympt_disp: 0.1,
            initial_extra_pois: 1.0,
            min_residual: 1e-4,
            max_residual: 15.0,
            coefficient_tol: 1e-6,
            glm_tol: 1e-8,
            max_outer_iter: 10,
            max_irls_iter: 100,
        }
    }
}

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

/// Parametric dispersion trend `dispersion = asympt_disp + extra_pois / mean`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ParametricDispersionTrend {
    /// DESeq2 `asymptDisp` coefficient.
    pub asympt_disp: f64,
    /// DESeq2 `extraPois` coefficient.
    pub extra_pois: f64,
}

impl ParametricDispersionTrend {
    /// Evaluate the parametric dispersion trend at one positive mean.
    pub fn evaluate(&self, mean: f64) -> Result<f64, DeseqError> {
        validate_positive_mean(mean, None)?;
        if !self.asympt_disp.is_finite() || self.asympt_disp <= 0.0 {
            return Err(DeseqError::InvalidDispersion {
                reason: "asymptotic dispersion coefficient must be finite and positive".to_string(),
            });
        }
        if !self.extra_pois.is_finite() || self.extra_pois <= 0.0 {
            return Err(DeseqError::InvalidDispersion {
                reason: "extra-Poisson dispersion coefficient must be finite and positive"
                    .to_string(),
            });
        }
        Ok(self.asympt_disp + self.extra_pois / mean)
    }

    /// Evaluate the trend for all finite positive means, returning `NaN` for missing rows.
    pub fn evaluate_many_allow_missing(&self, means: &[f64]) -> Result<Vec<f64>, DeseqError> {
        means
            .iter()
            .copied()
            .map(|mean| {
                if mean.is_finite() && mean > 0.0 {
                    self.evaluate(mean)
                } else {
                    Ok(f64::NAN)
                }
            })
            .collect()
    }
}

/// Output from a parametric dispersion trend fit.
#[derive(Clone, Debug, PartialEq)]
pub struct ParametricDispersionTrendFit {
    /// Fitted parametric trend coefficients.
    pub trend: ParametricDispersionTrend,
    /// Fitted dispersion values for every input row; missing rows are `NaN`.
    pub disp_fit: Vec<f64>,
    /// Rows used by DESeq2's `dispGeneEst > 100 * minDisp` trend fit rule.
    pub use_for_fit: Vec<bool>,
    /// Number of rows retained by `use_for_fit`.
    pub genes_used: usize,
    /// Number of robust outer-loop iterations.
    pub outer_iterations: usize,
    /// Number of inner Gamma identity-link IRLS iterations from the last outer fit.
    pub irls_iterations: usize,
    /// Whether the final inner Gamma identity-link IRLS fit converged.
    pub converged: bool,
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
        validate_positive_dispersion(self.mean_disp, None)?;
        Ok(self.mean_disp)
    }

    /// Evaluate the trend for all finite positive means, returning `NaN` for missing rows.
    pub fn evaluate_many_allow_missing(&self, means: &[f64]) -> Result<Vec<f64>, DeseqError> {
        means
            .iter()
            .copied()
            .map(|mean| {
                if mean.is_finite() && mean > 0.0 {
                    self.evaluate(mean)
                } else {
                    Ok(f64::NAN)
                }
            })
            .collect()
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

#[derive(Clone, Copy, Debug, PartialEq)]
struct GammaIdentityFit {
    trend: ParametricDispersionTrend,
    iterations: usize,
    converged: bool,
}

/// Fit DESeq2's parametric dispersion trend.
///
/// This mirrors the shape of `parametricDispersionFit`: start at `(0.1, 1)`,
/// keep rows with residual ratios in `(1e-4, 15)`, fit a Gamma GLM with
/// identity link for `disp ~ 1 + 1 / mean`, and stop when the log coefficient
/// change is below `1e-6`.
pub fn fit_parametric_dispersion_trend(
    base_mean: &[f64],
    disp_gene_est: &[f64],
    options: ParametricDispersionTrendOptions,
) -> Result<ParametricDispersionTrendFit, DeseqError> {
    validate_parametric_trend_inputs(base_mean, disp_gene_est, options)?;
    let use_for_fit = parametric_trend_use_for_fit(base_mean, disp_gene_est, options.min_disp)?;
    let genes_used = use_for_fit.iter().filter(|value| **value).count();
    if genes_used == 0 {
        return Err(DeseqError::InvalidDispersion {
            reason: "all gene-wise dispersion estimates are within 2 orders of magnitude from the minimum value".to_string(),
        });
    }

    let means = base_mean
        .iter()
        .copied()
        .zip(use_for_fit.iter().copied())
        .filter_map(|(mean, used)| used.then_some(mean))
        .collect::<Vec<_>>();
    let disps = disp_gene_est
        .iter()
        .copied()
        .zip(use_for_fit.iter().copied())
        .filter_map(|(disp, used)| used.then_some(disp))
        .collect::<Vec<_>>();
    let mut trend = ParametricDispersionTrend {
        asympt_disp: options.initial_asympt_disp,
        extra_pois: options.initial_extra_pois,
    };

    for outer in 0..=options.max_outer_iter {
        let good = robust_parametric_trend_rows(&means, &disps, trend, options)?;
        let fit = gamma_identity_irls(&means, &disps, &good, trend, options)?;
        let old = trend;
        trend = fit.trend;
        if trend.asympt_disp <= 0.0 || trend.extra_pois <= 0.0 {
            return Err(DeseqError::InvalidDispersion {
                reason: "parametric dispersion fit failed".to_string(),
            });
        }
        let coefficient_change = (trend.asympt_disp / old.asympt_disp).ln().powi(2)
            + (trend.extra_pois / old.extra_pois).ln().powi(2);
        if coefficient_change < options.coefficient_tol && fit.converged {
            return Ok(ParametricDispersionTrendFit {
                trend,
                disp_fit: trend.evaluate_many_allow_missing(base_mean)?,
                use_for_fit,
                genes_used,
                outer_iterations: outer,
                irls_iterations: fit.iterations,
                converged: true,
            });
        }
    }

    Err(DeseqError::InvalidDispersion {
        reason: "dispersion fit did not converge".to_string(),
    })
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
    let use_for_fit = parametric_trend_use_for_fit(base_mean, disp_gene_est, options.min_disp)?;
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

/// DESeq2's parametric trend row-selection rule.
pub fn parametric_trend_use_for_fit(
    base_mean: &[f64],
    disp_gene_est: &[f64],
    min_disp: f64,
) -> Result<Vec<bool>, DeseqError> {
    if base_mean.len() != disp_gene_est.len() {
        return Err(invalid_dimensions(
            "parametric dispersion trend rows",
            base_mean.len(),
            disp_gene_est.len(),
        ));
    }
    if !min_disp.is_finite() || min_disp <= 0.0 {
        return Err(DeseqError::InvalidDispersion {
            reason: "min_disp must be finite and positive".to_string(),
        });
    }
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
    if base_mean.len() != disp_gene_est.len() {
        return Err(invalid_dimensions(
            "mean dispersion trend rows",
            base_mean.len(),
            disp_gene_est.len(),
        ));
    }
    if !min_disp.is_finite() || min_disp <= 0.0 {
        return Err(DeseqError::InvalidDispersion {
            reason: "min_disp must be finite and positive".to_string(),
        });
    }
    Ok(base_mean
        .iter()
        .copied()
        .zip(disp_gene_est.iter().copied())
        .map(|(mean, disp)| {
            mean.is_finite() && mean > 0.0 && disp.is_finite() && disp > 10.0 * min_disp
        })
        .collect())
}

/// Placeholder for non-parametric dispersion trend types.
pub fn fit_dispersion_trend() -> Result<(), DeseqError> {
    Err(DeseqError::UnsupportedFeature {
        feature: "non-parametric dispersion trend fitting".to_string(),
    })
}

fn robust_parametric_trend_rows(
    means: &[f64],
    disps: &[f64],
    trend: ParametricDispersionTrend,
    options: ParametricDispersionTrendOptions,
) -> Result<Vec<bool>, DeseqError> {
    let mut good = Vec::with_capacity(means.len());
    for (idx, (mean, disp)) in means.iter().copied().zip(disps.iter().copied()).enumerate() {
        validate_positive_mean(mean, Some(idx))?;
        validate_positive_dispersion(disp, Some(idx))?;
        let fitted = trend.evaluate(mean)?;
        let residual = disp / fitted;
        good.push(
            residual.is_finite()
                && residual > options.min_residual
                && residual < options.max_residual,
        );
    }
    if good.iter().filter(|value| **value).count() < 2 {
        return Err(DeseqError::InvalidDispersion {
            reason: "not enough genes remain for parametric dispersion trend fitting".to_string(),
        });
    }
    Ok(good)
}

fn gamma_identity_irls(
    means: &[f64],
    disps: &[f64],
    good: &[bool],
    start: ParametricDispersionTrend,
    options: ParametricDispersionTrendOptions,
) -> Result<GammaIdentityFit, DeseqError> {
    let mut trend = start;
    let mut old_deviance = gamma_deviance(means, disps, good, trend)?;
    for iteration in 1..=options.max_irls_iter {
        let candidate = weighted_gamma_identity_least_squares(means, disps, good, trend)?;
        let candidate = positive_step_halving(means, good, trend, candidate)?;
        let deviance = gamma_deviance(means, disps, good, candidate)?;
        let converged =
            ((deviance - old_deviance).abs() / (deviance.abs() + 0.1)) < options.glm_tol;
        trend = candidate;
        if converged {
            return Ok(GammaIdentityFit {
                trend,
                iterations: iteration,
                converged: true,
            });
        }
        old_deviance = deviance;
    }
    Ok(GammaIdentityFit {
        trend,
        iterations: options.max_irls_iter,
        converged: false,
    })
}

fn weighted_gamma_identity_least_squares(
    means: &[f64],
    disps: &[f64],
    good: &[bool],
    trend: ParametricDispersionTrend,
) -> Result<ParametricDispersionTrend, DeseqError> {
    let mut s00 = 0.0;
    let mut s01 = 0.0;
    let mut s11 = 0.0;
    let mut t0 = 0.0;
    let mut t1 = 0.0;
    for (idx, ((mean, disp), used)) in means
        .iter()
        .copied()
        .zip(disps.iter().copied())
        .zip(good.iter().copied())
        .enumerate()
    {
        if !used {
            continue;
        }
        validate_positive_mean(mean, Some(idx))?;
        validate_positive_dispersion(disp, Some(idx))?;
        let fitted = trend.evaluate(mean)?;
        let weight = fitted.powi(-2);
        let x = mean.recip();
        s00 += weight;
        s01 += weight * x;
        s11 += weight * x * x;
        t0 += weight * disp;
        t1 += weight * x * disp;
    }
    let determinant = s00 * s11 - s01 * s01;
    if !determinant.is_finite() || determinant.abs() <= f64::EPSILON {
        return Err(DeseqError::InvalidDimensions {
            context: "parametric dispersion trend weighted design".to_string(),
            expected: 2,
            actual: 0,
        });
    }
    Ok(ParametricDispersionTrend {
        asympt_disp: (t0 * s11 - t1 * s01) / determinant,
        extra_pois: (s00 * t1 - s01 * t0) / determinant,
    })
}

fn positive_step_halving(
    means: &[f64],
    good: &[bool],
    current: ParametricDispersionTrend,
    target: ParametricDispersionTrend,
) -> Result<ParametricDispersionTrend, DeseqError> {
    let mut step = 1.0;
    for _ in 0..50 {
        let candidate = ParametricDispersionTrend {
            asympt_disp: current.asympt_disp + step * (target.asympt_disp - current.asympt_disp),
            extra_pois: current.extra_pois + step * (target.extra_pois - current.extra_pois),
        };
        if candidate.asympt_disp > 0.0
            && candidate.extra_pois > 0.0
            && fitted_values_are_positive(means, good, candidate)?
        {
            return Ok(candidate);
        }
        step *= 0.5;
    }
    Err(DeseqError::InvalidDispersion {
        reason: "parametric dispersion fit failed".to_string(),
    })
}

fn fitted_values_are_positive(
    means: &[f64],
    good: &[bool],
    trend: ParametricDispersionTrend,
) -> Result<bool, DeseqError> {
    for (idx, (mean, used)) in means.iter().copied().zip(good.iter().copied()).enumerate() {
        if used {
            validate_positive_mean(mean, Some(idx))?;
            let fitted = trend.evaluate(mean)?;
            if !fitted.is_finite() || fitted <= 0.0 {
                return Ok(false);
            }
        }
    }
    Ok(true)
}

fn gamma_deviance(
    means: &[f64],
    disps: &[f64],
    good: &[bool],
    trend: ParametricDispersionTrend,
) -> Result<f64, DeseqError> {
    let mut deviance = 0.0;
    for (idx, ((mean, disp), used)) in means
        .iter()
        .copied()
        .zip(disps.iter().copied())
        .zip(good.iter().copied())
        .enumerate()
    {
        if !used {
            continue;
        }
        validate_positive_mean(mean, Some(idx))?;
        validate_positive_dispersion(disp, Some(idx))?;
        let fitted = trend.evaluate(mean)?;
        let ratio = disp / fitted;
        deviance += 2.0 * ((disp - fitted) / fitted - ratio.ln());
    }
    Ok(deviance)
}

fn validate_parametric_trend_inputs(
    base_mean: &[f64],
    disp_gene_est: &[f64],
    options: ParametricDispersionTrendOptions,
) -> Result<(), DeseqError> {
    if base_mean.len() != disp_gene_est.len() {
        return Err(invalid_dimensions(
            "parametric dispersion trend rows",
            base_mean.len(),
            disp_gene_est.len(),
        ));
    }
    if base_mean.is_empty() {
        return Err(DeseqError::InvalidDimensions {
            context: "parametric dispersion trend rows".to_string(),
            expected: 1,
            actual: 0,
        });
    }
    if !options.min_disp.is_finite() || options.min_disp <= 0.0 {
        return Err(DeseqError::InvalidDispersion {
            reason: "min_disp must be finite and positive".to_string(),
        });
    }
    if !options.initial_asympt_disp.is_finite() || options.initial_asympt_disp <= 0.0 {
        return Err(DeseqError::InvalidDispersion {
            reason: "initial asymptotic dispersion must be finite and positive".to_string(),
        });
    }
    if !options.initial_extra_pois.is_finite() || options.initial_extra_pois <= 0.0 {
        return Err(DeseqError::InvalidDispersion {
            reason: "initial extra-Poisson dispersion must be finite and positive".to_string(),
        });
    }
    if !options.min_residual.is_finite()
        || !options.max_residual.is_finite()
        || options.min_residual <= 0.0
        || options.max_residual <= options.min_residual
    {
        return Err(DeseqError::InvalidDispersion {
            reason: "parametric trend residual bounds must be finite and ordered".to_string(),
        });
    }
    if !options.coefficient_tol.is_finite() || options.coefficient_tol <= 0.0 {
        return Err(DeseqError::InvalidDispersion {
            reason: "parametric trend coefficient tolerance must be finite and positive"
                .to_string(),
        });
    }
    if !options.glm_tol.is_finite() || options.glm_tol <= 0.0 {
        return Err(DeseqError::InvalidDispersion {
            reason: "parametric trend GLM tolerance must be finite and positive".to_string(),
        });
    }
    if options.max_irls_iter == 0 {
        return Err(DeseqError::InvalidDimensions {
            context: "parametric trend max IRLS iterations".to_string(),
            expected: 1,
            actual: 0,
        });
    }
    Ok(())
}

fn validate_mean_trend_inputs(
    base_mean: &[f64],
    disp_gene_est: &[f64],
    options: MeanDispersionTrendOptions,
) -> Result<(), DeseqError> {
    if base_mean.len() != disp_gene_est.len() {
        return Err(invalid_dimensions(
            "mean dispersion trend rows",
            base_mean.len(),
            disp_gene_est.len(),
        ));
    }
    if base_mean.is_empty() {
        return Err(DeseqError::InvalidDimensions {
            context: "mean dispersion trend rows".to_string(),
            expected: 1,
            actual: 0,
        });
    }
    if !options.min_disp.is_finite() || options.min_disp <= 0.0 {
        return Err(DeseqError::InvalidDispersion {
            reason: "min_disp must be finite and positive".to_string(),
        });
    }
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
        return if sorted.len() % 2 == 0 {
            Ok((sorted[mid - 1] + sorted[mid]) / 2.0)
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
    Ok(sorted[start..end].iter().sum::<f64>() / (end - start) as f64)
}

fn validate_positive_mean(mean: f64, index: Option<usize>) -> Result<(), DeseqError> {
    if !mean.is_finite() || mean <= 0.0 {
        return Err(DeseqError::NonFiniteValue {
            context: "parametric dispersion trend mean".to_string(),
            index,
            value: mean,
        });
    }
    Ok(())
}

fn validate_positive_dispersion(dispersion: f64, index: Option<usize>) -> Result<(), DeseqError> {
    if !dispersion.is_finite() || dispersion <= 0.0 {
        return Err(DeseqError::InvalidDispersion {
            reason: format!(
                "gene-wise dispersion at index {} must be finite and positive",
                index
                    .map(|value| value.to_string())
                    .unwrap_or_else(|| "unknown".to_string())
            ),
        });
    }
    Ok(())
}
