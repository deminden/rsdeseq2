use crate::errors::{invalid_dimensions, DeseqError};
use crate::options::FitType;

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

/// Options for the pure-Rust local log-dispersion trend.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LocalDispersionTrendOptions {
    /// Minimum dispersion used by DESeq2's local trend floor and fit rule.
    pub min_disp: f64,
    /// Fraction of fitted rows used in the adaptive local neighborhood.
    pub span: f64,
    /// Local polynomial degree on log(mean), capped to quadratic.
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

/// Local dispersion trend fit on log mean/log dispersion scale.
///
/// DESeq2 uses the R `locfit` package for `fitType="local"`. This pure-Rust
/// representation keeps the same fitted rows and base-mean weights, then
/// evaluates a deterministic adaptive local polynomial smoother.
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
        if !self.min_disp.is_finite() || self.min_disp <= 0.0 {
            return Err(DeseqError::InvalidDispersion {
                reason: "minimum dispersion must be finite and positive".to_string(),
            });
        }
        if self.log_means.is_empty() {
            return Ok(self.min_disp);
        }
        let predicted = local_polynomial_predict(
            mean.ln(),
            &self.log_means,
            &self.log_disps,
            &self.weights,
            self.span,
            self.degree,
        )?;
        let dispersion = predicted.exp();
        if !dispersion.is_finite() || dispersion <= 0.0 {
            return Err(DeseqError::InvalidDispersion {
                reason: "local dispersion trend produced a non-positive fitted value".to_string(),
            });
        }
        Ok(dispersion)
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

/// Fitted dispersion trend selected by a DESeq2-style `fitType`.
#[derive(Clone, Debug, PartialEq)]
pub enum DispersionTrendFit {
    /// Parametric `asymptDisp + extraPois / mean` trend.
    Parametric(ParametricDispersionTrendFit),
    /// Pure-Rust local log-dispersion trend.
    Local(LocalDispersionTrendFit),
    /// Constant mean dispersion trend.
    Mean(MeanDispersionTrendFit),
}

impl DispersionTrendFit {
    /// Stable DESeq2-style fit type label for this fitted trend.
    pub fn fit_type_label(&self) -> &'static str {
        match self {
            Self::Parametric(_) => "parametric",
            Self::Local(_) => "local",
            Self::Mean(_) => "mean",
        }
    }

    /// Fitted dispersion values for every input row; missing rows are `NaN`.
    pub fn disp_fit(&self) -> &[f64] {
        match self {
            Self::Parametric(fit) => &fit.disp_fit,
            Self::Local(fit) => &fit.disp_fit,
            Self::Mean(fit) => &fit.disp_fit,
        }
    }

    /// Rows retained by the selected trend fit rule.
    pub fn use_for_fit(&self) -> &[bool] {
        match self {
            Self::Parametric(fit) => &fit.use_for_fit,
            Self::Local(fit) => &fit.use_for_fit,
            Self::Mean(fit) => &fit.use_for_fit,
        }
    }

    /// Number of rows retained by the selected trend fit rule.
    pub fn genes_used(&self) -> usize {
        match self {
            Self::Parametric(fit) => fit.genes_used,
            Self::Local(fit) => fit.genes_used,
            Self::Mean(fit) => fit.genes_used,
        }
    }
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
        let asympt_change = (trend.asympt_disp / old.asympt_disp).ln();
        let extra_pois_change = (trend.extra_pois / old.extra_pois).ln();
        let coefficient_change =
            asympt_change * asympt_change + extra_pois_change * extra_pois_change;
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

/// Fit DESeq2's `fitType="local"` dispersion trend.
///
/// DESeq2 fits local regression on `log(dispersion)` versus `log(mean)` after
/// selecting rows with `dispGeneEst >= 10 * minDisp`, using base means as
/// weights. This implementation keeps that data contract and evaluates a
/// deterministic adaptive local polynomial smoother in Rust.
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
        degree: options.degree.min(2),
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

/// DESeq2's `fitType="local"` row-selection rule for local regression.
pub fn local_trend_use_for_fit(
    base_mean: &[f64],
    disp_gene_est: &[f64],
    min_disp: f64,
) -> Result<Vec<bool>, DeseqError> {
    if base_mean.len() != disp_gene_est.len() {
        return Err(invalid_dimensions(
            "local dispersion trend rows",
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
            mean.is_finite() && mean > 0.0 && disp.is_finite() && disp >= 10.0 * min_disp
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
        let inv_fitted = fitted.recip();
        let weight = inv_fitted * inv_fitted;
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
        let relative_delta = ratio - 1.0;
        deviance += 2.0 * (relative_delta - relative_delta.ln_1p());
    }
    Ok(deviance)
}

fn local_polynomial_predict(
    x0: f64,
    xs: &[f64],
    ys: &[f64],
    weights: &[f64],
    span: f64,
    degree: usize,
) -> Result<f64, DeseqError> {
    if xs.len() != ys.len() || xs.len() != weights.len() {
        return Err(invalid_dimensions(
            "local dispersion trend fit rows",
            xs.len(),
            ys.len(),
        ));
    }
    if xs.is_empty() {
        return Err(DeseqError::InvalidDispersion {
            reason: "no local dispersion fit rows are available".to_string(),
        });
    }
    if xs.len() == 1 {
        return ys
            .first()
            .copied()
            .filter(|value| value.is_finite())
            .ok_or_else(|| DeseqError::InvalidDispersion {
                reason: "single-row local dispersion trend value must be finite".to_string(),
            });
    }

    validate_sorted_local_means(xs)?;
    let degree = degree.min(2).min(xs.len().saturating_sub(1));
    let window = adaptive_nearest_window(x0, xs, span, degree + 1);
    let bandwidth = local_window_bandwidth(x0, xs, window);
    for local_degree in (0..=degree).rev() {
        if let Some(beta0) =
            local_polynomial_intercept(x0, xs, ys, weights, window, bandwidth, local_degree)
        {
            if beta0.is_finite() {
                return Ok(beta0);
            }
        }
    }

    let mut numerator = 0.0;
    let mut denominator = 0.0;
    for idx in window.0..window.1 {
        let weight = weights[idx] * tricube_weight((xs[idx] - x0).abs(), bandwidth);
        let y = ys[idx];
        if weight.is_finite() && weight > 0.0 && y.is_finite() {
            numerator += weight * y;
            denominator += weight;
        }
    }
    if denominator > 0.0 {
        return Ok(numerator / denominator);
    }
    xs.iter()
        .copied()
        .zip(ys.iter().copied())
        .min_by(|(left_x, _), (right_x, _)| (left_x - x0).abs().total_cmp(&(right_x - x0).abs()))
        .map(|(_, y)| y)
        .filter(|value| value.is_finite())
        .ok_or_else(|| DeseqError::InvalidDispersion {
            reason: "local dispersion trend has zero usable neighborhood weight".to_string(),
        })
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

fn adaptive_nearest_window(x0: f64, xs: &[f64], span: f64, min_neighbors: usize) -> (usize, usize) {
    let neighbors =
        ((xs.len() as f64 * span).ceil() as usize).clamp(min_neighbors.max(1), xs.len());
    let insertion = xs.partition_point(|x| *x < x0);
    let mut start = insertion;
    let mut end = insertion;
    while end - start < neighbors {
        if start == 0 {
            end = (start + neighbors).min(xs.len());
            break;
        }
        if end == xs.len() {
            start = end.saturating_sub(neighbors);
            break;
        }
        let left_distance = (xs[start - 1] - x0).abs();
        let right_distance = (xs[end] - x0).abs();
        if left_distance <= right_distance {
            start -= 1;
        } else {
            end += 1;
        }
    }
    (start, end)
}

fn local_window_bandwidth(x0: f64, xs: &[f64], window: (usize, usize)) -> f64 {
    let mut bandwidth = 0.0_f64;
    for x in xs[window.0..window.1].iter().copied() {
        bandwidth = bandwidth.max((x - x0).abs());
    }
    bandwidth
}

fn tricube_weight(distance: f64, bandwidth: f64) -> f64 {
    if bandwidth <= f64::EPSILON {
        if distance <= f64::EPSILON {
            1.0
        } else {
            0.0
        }
    } else if distance <= bandwidth {
        let scaled = distance / bandwidth;
        let scaled_cubed = scaled * scaled * scaled;
        let one_minus_scaled_cubed = 1.0 - scaled_cubed;
        one_minus_scaled_cubed * one_minus_scaled_cubed * one_minus_scaled_cubed
    } else {
        0.0
    }
}

fn local_polynomial_intercept(
    x0: f64,
    xs: &[f64],
    ys: &[f64],
    weights: &[f64],
    window: (usize, usize),
    bandwidth: f64,
    degree: usize,
) -> Option<f64> {
    let size = degree + 1;
    let mut lhs = [[0.0; 3]; 3];
    let mut rhs = [0.0; 3];
    for idx in window.0..window.1 {
        let x = xs[idx];
        let y = ys[idx];
        let weight = weights[idx] * tricube_weight((x - x0).abs(), bandwidth);
        if !weight.is_finite() || weight <= 0.0 || !x.is_finite() || !y.is_finite() {
            continue;
        }
        let dx = x - x0;
        let powers = [1.0, dx, dx * dx, dx * dx * dx, dx * dx * dx * dx];
        for row in 0..size {
            rhs[row] += weight * powers[row] * y;
            for col in 0..size {
                lhs[row][col] += weight * powers[row + col];
            }
        }
    }
    solve_small_linear_system(lhs, rhs, size).map(|beta| beta[0])
}

fn solve_small_linear_system(
    mut lhs: [[f64; 3]; 3],
    mut rhs: [f64; 3],
    size: usize,
) -> Option<[f64; 3]> {
    for pivot in 0..size {
        let mut pivot_row = pivot;
        for row in (pivot + 1)..size {
            if lhs[row][pivot].abs() > lhs[pivot_row][pivot].abs() {
                pivot_row = row;
            }
        }
        if !lhs[pivot_row][pivot].is_finite() || lhs[pivot_row][pivot].abs() <= 1e-12 {
            return None;
        }
        if pivot_row != pivot {
            lhs.swap(pivot, pivot_row);
            rhs.swap(pivot, pivot_row);
        }
        let pivot_value = lhs[pivot][pivot];
        for value in lhs[pivot].iter_mut().take(size).skip(pivot) {
            *value /= pivot_value;
        }
        rhs[pivot] /= pivot_value;
        let pivot_values = lhs[pivot];
        for row in 0..size {
            if row == pivot {
                continue;
            }
            let factor = lhs[row][pivot];
            if factor == 0.0 {
                continue;
            }
            for (col, pivot_entry) in pivot_values.iter().enumerate().take(size).skip(pivot) {
                lhs[row][col] -= factor * pivot_entry;
            }
            rhs[row] -= factor * rhs[pivot];
        }
    }
    Some(rhs)
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

fn validate_local_trend_inputs(
    base_mean: &[f64],
    disp_gene_est: &[f64],
    options: LocalDispersionTrendOptions,
) -> Result<(), DeseqError> {
    if base_mean.len() != disp_gene_est.len() {
        return Err(invalid_dimensions(
            "local dispersion trend rows",
            base_mean.len(),
            disp_gene_est.len(),
        ));
    }
    if base_mean.is_empty() {
        return Err(DeseqError::InvalidDimensions {
            context: "local dispersion trend rows".to_string(),
            expected: 1,
            actual: 0,
        });
    }
    if !options.min_disp.is_finite() || options.min_disp <= 0.0 {
        return Err(DeseqError::InvalidDispersion {
            reason: "min_disp must be finite and positive".to_string(),
        });
    }
    if !options.span.is_finite() || options.span <= 0.0 || options.span > 1.0 {
        return Err(DeseqError::InvalidDispersion {
            reason: "local dispersion span must be finite and in (0, 1]".to_string(),
        });
    }
    if options.degree > 2 {
        return Err(DeseqError::InvalidDispersion {
            reason: "local dispersion polynomial degree must be 0, 1, or 2".to_string(),
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

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn gamma_deviance_keeps_near_perfect_fit_precision() {
        let trend = ParametricDispersionTrend {
            asympt_disp: 0.05,
            extra_pois: 2.0,
        };
        let means = [20.0, 40.0, 80.0, 160.0];
        let relative_delta = 1.0e-8;
        let disps = means
            .iter()
            .map(|mean| trend.evaluate(*mean).unwrap() * (1.0 + relative_delta))
            .collect::<Vec<_>>();
        let good = [true; 4];

        let deviance = gamma_deviance(&means, &disps, &good, trend).unwrap();
        let expected = means.len() as f64 * 2.0 * (relative_delta - relative_delta.ln_1p());

        assert!(deviance.is_finite());
        assert!(deviance > 0.0);
        assert_relative_eq!(deviance, expected, max_relative = 1e-12);
    }
}
