//! Parametric trend model and Gamma identity-link fitting.

use super::numerics::{
    checked_add_assign, checked_div2, checked_product_difference, checked_product2,
    checked_product3, checked_sub2, checked_sum2, validate_positive_dispersion,
    validate_positive_mean,
};
use super::selection::{parametric_trend_use_for_fit, validate_trend_fit_rows};
use crate::errors::{DeseqError, invalid_dimensions};

const PARAMETRIC_WLS_SUM_CONTEXT: &str = "parametric dispersion weighted least-squares sums";

const PARAMETRIC_WLS_DETERMINANT_CONTEXT: &str =
    "parametric dispersion weighted least-squares determinant";

const PARAMETRIC_WLS_NUMERATOR_CONTEXT: &str =
    "parametric dispersion weighted least-squares coefficient numerator";

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
        self.validate_fit()?;
        self.evaluate_validated(mean)
    }

    /// Evaluate the trend for all finite positive means, returning `NaN` for missing rows.
    pub fn evaluate_many_allow_missing(&self, means: &[f64]) -> Result<Vec<f64>, DeseqError> {
        self.validate_fit()?;
        means
            .iter()
            .copied()
            .map(|mean| {
                if mean.is_finite() && mean > 0.0 {
                    self.evaluate_validated(mean)
                } else {
                    Ok(f64::NAN)
                }
            })
            .collect()
    }

    fn evaluate_validated(&self, mean: f64) -> Result<f64, DeseqError> {
        let fitted = self.asympt_disp + self.extra_pois / mean;
        if !fitted.is_finite() || fitted <= 0.0 {
            return Err(DeseqError::InvalidDispersion {
                reason:
                    "parametric dispersion trend produced a non-finite or non-positive fitted value"
                        .to_string(),
            });
        }
        Ok(fitted)
    }

    fn validate_fit(&self) -> Result<(), DeseqError> {
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
        Ok(())
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
        let coefficient_change = parametric_coefficient_change(old, trend)?;
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

pub(super) fn robust_parametric_trend_rows(
    means: &[f64],
    disps: &[f64],
    trend: ParametricDispersionTrend,
    options: ParametricDispersionTrendOptions,
) -> Result<Vec<bool>, DeseqError> {
    validate_parametric_work_rows(means, disps)?;
    trend.validate_fit()?;
    let mut good = Vec::with_capacity(means.len());
    for (idx, (mean, disp)) in means.iter().copied().zip(disps.iter().copied()).enumerate() {
        validate_positive_mean(mean, Some(idx))?;
        validate_positive_dispersion(disp, Some(idx))?;
        let fitted = trend.evaluate_validated(mean)?;
        let residual = checked_div2(disp, fitted, "parametric dispersion residual filter")?;
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
        let deviance_change = checked_sub2(
            deviance,
            old_deviance,
            "parametric dispersion gamma IRLS convergence",
        )?
        .abs();
        let deviance_scale = checked_sum2(
            deviance.abs(),
            0.1,
            "parametric dispersion gamma IRLS convergence",
        )?;
        let converged = checked_div2(
            deviance_change,
            deviance_scale,
            "parametric dispersion gamma IRLS convergence",
        )? < options.glm_tol;
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

pub(super) fn weighted_gamma_identity_least_squares(
    means: &[f64],
    disps: &[f64],
    good: &[bool],
    trend: ParametricDispersionTrend,
) -> Result<ParametricDispersionTrend, DeseqError> {
    validate_parametric_good_rows(means, disps, good)?;
    let mut s00 = 0.0;
    let mut s01 = 0.0;
    let mut s11 = 0.0;
    let mut t0 = 0.0;
    let mut t1 = 0.0;
    trend.validate_fit()?;
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
        let fitted = trend.evaluate_validated(mean)?;
        let inv_fitted = checked_div2(1.0, fitted, PARAMETRIC_WLS_SUM_CONTEXT)?;
        let weight = checked_product2(
            inv_fitted,
            inv_fitted,
            "parametric dispersion weighted least-squares weights",
        )?;
        let x = mean.recip();
        checked_add_assign(&mut s00, weight, PARAMETRIC_WLS_SUM_CONTEXT)?;
        checked_add_assign(
            &mut s01,
            checked_product2(weight, x, PARAMETRIC_WLS_SUM_CONTEXT)?,
            PARAMETRIC_WLS_SUM_CONTEXT,
        )?;
        checked_add_assign(
            &mut s11,
            checked_product3(weight, x, x, PARAMETRIC_WLS_SUM_CONTEXT)?,
            PARAMETRIC_WLS_SUM_CONTEXT,
        )?;
        checked_add_assign(
            &mut t0,
            checked_product2(weight, disp, PARAMETRIC_WLS_SUM_CONTEXT)?,
            PARAMETRIC_WLS_SUM_CONTEXT,
        )?;
        checked_add_assign(
            &mut t1,
            checked_product3(weight, x, disp, PARAMETRIC_WLS_SUM_CONTEXT)?,
            PARAMETRIC_WLS_SUM_CONTEXT,
        )?;
    }
    let determinant =
        checked_product_difference(s00, s11, s01, s01, PARAMETRIC_WLS_DETERMINANT_CONTEXT)?;
    if !determinant.is_finite() || determinant.abs() <= f64::EPSILON {
        return Err(DeseqError::InvalidDimensions {
            context: "parametric dispersion trend weighted design".to_string(),
            expected: 2,
            actual: 0,
        });
    }
    let asympt_numerator =
        checked_product_difference(t0, s11, t1, s01, PARAMETRIC_WLS_NUMERATOR_CONTEXT)?;
    let extra_pois_numerator =
        checked_product_difference(s00, t1, s01, t0, PARAMETRIC_WLS_NUMERATOR_CONTEXT)?;
    let target = ParametricDispersionTrend {
        asympt_disp: checked_div2(
            asympt_numerator,
            determinant,
            "parametric dispersion weighted least-squares coefficients",
        )?,
        extra_pois: checked_div2(
            extra_pois_numerator,
            determinant,
            "parametric dispersion weighted least-squares coefficients",
        )?,
    };
    if !target.asympt_disp.is_finite() || !target.extra_pois.is_finite() {
        return Err(DeseqError::InvalidDispersion {
            reason: "parametric dispersion weighted least-squares produced non-finite coefficients"
                .to_string(),
        });
    }
    Ok(target)
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
            asympt_disp: interpolate_parametric_coefficient(
                current.asympt_disp,
                target.asympt_disp,
                step,
            ),
            extra_pois: interpolate_parametric_coefficient(
                current.extra_pois,
                target.extra_pois,
                step,
            ),
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

pub(super) fn interpolate_parametric_coefficient(current: f64, target: f64, step: f64) -> f64 {
    current * (1.0 - step) + target * step
}

pub(super) fn fitted_values_are_positive(
    means: &[f64],
    good: &[bool],
    trend: ParametricDispersionTrend,
) -> Result<bool, DeseqError> {
    if means.len() != good.len() {
        return Err(invalid_dimensions(
            "parametric dispersion trend good rows",
            means.len(),
            good.len(),
        ));
    }
    trend.validate_fit()?;
    for (idx, (mean, used)) in means.iter().copied().zip(good.iter().copied()).enumerate() {
        if used {
            validate_positive_mean(mean, Some(idx))?;
            let fitted = trend.asympt_disp + trend.extra_pois / mean;
            if !fitted.is_finite() || fitted <= 0.0 {
                return Ok(false);
            }
        }
    }
    Ok(true)
}

pub(super) fn gamma_deviance(
    means: &[f64],
    disps: &[f64],
    good: &[bool],
    trend: ParametricDispersionTrend,
) -> Result<f64, DeseqError> {
    validate_parametric_good_rows(means, disps, good)?;
    trend.validate_fit()?;
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
        let fitted = trend.evaluate_validated(mean)?;
        let ratio = disp / fitted;
        let relative_delta = ratio - 1.0;
        checked_add_assign(
            &mut deviance,
            2.0 * (relative_delta - relative_delta.ln_1p()),
            "parametric dispersion gamma deviance",
        )?;
    }
    Ok(deviance)
}

fn validate_parametric_work_rows(means: &[f64], disps: &[f64]) -> Result<(), DeseqError> {
    if means.len() != disps.len() {
        return Err(invalid_dimensions(
            "parametric dispersion trend working rows",
            means.len(),
            disps.len(),
        ));
    }
    Ok(())
}

fn validate_parametric_good_rows(
    means: &[f64],
    disps: &[f64],
    good: &[bool],
) -> Result<(), DeseqError> {
    validate_parametric_work_rows(means, disps)?;
    if means.len() != good.len() {
        return Err(invalid_dimensions(
            "parametric dispersion trend good rows",
            means.len(),
            good.len(),
        ));
    }
    Ok(())
}

pub(super) fn parametric_coefficient_change(
    old: ParametricDispersionTrend,
    new: ParametricDispersionTrend,
) -> Result<f64, DeseqError> {
    old.validate_fit()?;
    new.validate_fit()?;
    let context = "parametric dispersion coefficient change";
    let asympt_change = checked_div2(new.asympt_disp, old.asympt_disp, context)?.ln();
    let extra_pois_change = checked_div2(new.extra_pois, old.extra_pois, context)?.ln();
    checked_sum2(
        checked_product2(asympt_change, asympt_change, context)?,
        checked_product2(extra_pois_change, extra_pois_change, context)?,
        context,
    )
}

fn validate_parametric_trend_inputs(
    base_mean: &[f64],
    disp_gene_est: &[f64],
    options: ParametricDispersionTrendOptions,
) -> Result<(), DeseqError> {
    validate_trend_fit_rows(
        base_mean,
        disp_gene_est,
        options.min_disp,
        "parametric dispersion trend rows",
    )?;
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
    if options.max_outer_iter == 0 {
        return Err(DeseqError::InvalidDimensions {
            context: "parametric trend max outer iterations".to_string(),
            expected: 1,
            actual: 0,
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
