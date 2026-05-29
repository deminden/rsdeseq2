use crate::errors::{invalid_dimensions, DeseqError};
use crate::options::FitType;

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
        self.validate_fit()?;
        self.evaluate_validated(mean)
    }

    fn evaluate_validated(&self, mean: f64) -> Result<f64, DeseqError> {
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
                reason: "local dispersion trend produced a non-finite or non-positive fitted value"
                    .to_string(),
            });
        }
        Ok(dispersion)
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
    ///
    /// For `fitType="mean"`, this is the preliminary viability mask using
    /// `dispGeneEst > 100 * minDisp`; use [`Self::use_for_mean`] to inspect the
    /// trimmed-mean mask.
    pub fn use_for_fit(&self) -> &[bool] {
        match self {
            Self::Parametric(fit) => &fit.use_for_fit,
            Self::Local(fit) => &fit.use_for_fit,
            Self::Mean(fit) => &fit.use_for_fit,
        }
    }

    /// Rows used by the `fitType="mean"` trimmed-mean calculation.
    ///
    /// Returns `None` for parametric and local trends.
    pub fn use_for_mean(&self) -> Option<&[bool]> {
        match self {
            Self::Mean(fit) => Some(&fit.use_for_mean),
            Self::Parametric(_) | Self::Local(_) => None,
        }
    }

    /// Number of rows used by the `fitType="mean"` trimmed-mean calculation.
    ///
    /// Returns `None` for parametric and local trends.
    pub fn genes_used_for_mean(&self) -> Option<usize> {
        match self {
            Self::Mean(fit) => Some(fit.genes_used_for_mean),
            Self::Parametric(_) | Self::Local(_) => None,
        }
    }

    /// Whether `fitType="local"` used DESeq2's all-near-minimum dispersion floor.
    ///
    /// Returns `None` for parametric and mean trends.
    pub fn used_min_disp_floor(&self) -> Option<bool> {
        match self {
            Self::Local(fit) => Some(fit.used_min_disp_floor),
            Self::Parametric(_) | Self::Mean(_) => None,
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

fn robust_parametric_trend_rows(
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

fn weighted_gamma_identity_least_squares(
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

fn checked_add_assign(total: &mut f64, value: f64, context: &str) -> Result<(), DeseqError> {
    let next = *total + value;
    if !value.is_finite() || !next.is_finite() {
        return Err(DeseqError::InvalidDispersion {
            reason: format!("{context} must remain finite"),
        });
    }
    *total = next;
    Ok(())
}

fn checked_product2(left: f64, right: f64, context: &str) -> Result<f64, DeseqError> {
    let product = left * right;
    if !left.is_finite() || !right.is_finite() || !product.is_finite() {
        return Err(DeseqError::InvalidDispersion {
            reason: format!("{context} must remain finite"),
        });
    }
    Ok(product)
}

fn checked_product3(left: f64, middle: f64, right: f64, context: &str) -> Result<f64, DeseqError> {
    checked_product2(checked_product2(left, middle, context)?, right, context)
}

fn checked_sum2(left: f64, right: f64, context: &str) -> Result<f64, DeseqError> {
    let sum = left + right;
    if !left.is_finite() || !right.is_finite() || !sum.is_finite() {
        return Err(DeseqError::InvalidDispersion {
            reason: format!("{context} must remain finite"),
        });
    }
    Ok(sum)
}

fn checked_sub2(left: f64, right: f64, context: &str) -> Result<f64, DeseqError> {
    let difference = left - right;
    if !left.is_finite() || !right.is_finite() || !difference.is_finite() {
        return Err(DeseqError::InvalidDispersion {
            reason: format!("{context} must remain finite"),
        });
    }
    Ok(difference)
}

fn checked_div2(left: f64, right: f64, context: &str) -> Result<f64, DeseqError> {
    let quotient = left / right;
    if !left.is_finite() || !right.is_finite() || right == 0.0 || !quotient.is_finite() {
        return Err(DeseqError::InvalidDispersion {
            reason: format!("{context} must remain finite"),
        });
    }
    Ok(quotient)
}

fn checked_product_difference(
    left_a: f64,
    left_b: f64,
    right_a: f64,
    right_b: f64,
    context: &str,
) -> Result<f64, DeseqError> {
    let left = left_a * left_b;
    let right = right_a * right_b;
    let difference = left - right;
    if !left.is_finite() || !right.is_finite() || !difference.is_finite() {
        return Err(DeseqError::InvalidDispersion {
            reason: format!("{context} must remain finite"),
        });
    }
    Ok(difference)
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

fn interpolate_parametric_coefficient(current: f64, target: f64, step: f64) -> f64 {
    current * (1.0 - step) + target * step
}

fn fitted_values_are_positive(
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

fn gamma_deviance(
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

fn local_polynomial_predict(
    x0: f64,
    xs: &[f64],
    ys: &[f64],
    weights: &[f64],
    span: f64,
    degree: usize,
) -> Result<f64, DeseqError> {
    debug_assert_eq!(xs.len(), ys.len());
    debug_assert_eq!(xs.len(), weights.len());
    debug_assert!(!xs.is_empty());
    if xs.is_empty() {
        return Err(DeseqError::InvalidDispersion {
            reason: "no local dispersion fit rows are available".to_string(),
        });
    }
    if xs.len() == 1 {
        return ys
            .first()
            .copied()
            .ok_or_else(|| DeseqError::InvalidDispersion {
                reason: "single-row local dispersion trend value must be finite".to_string(),
            });
    }

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

    if let Some(beta0) = local_weighted_mean(x0, xs, ys, weights, window, bandwidth) {
        return Ok(beta0);
    }
    xs.iter()
        .copied()
        .zip(ys.iter().copied())
        .min_by(|(left_x, _), (right_x, _)| (left_x - x0).abs().total_cmp(&(right_x - x0).abs()))
        .map(|(_, y)| y)
        .ok_or_else(|| DeseqError::InvalidDispersion {
            reason: "local dispersion trend has zero usable neighborhood weight".to_string(),
        })
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
    if degree == 0 {
        return local_weighted_mean(x0, xs, ys, weights, window, bandwidth);
    }

    let size = degree + 1;
    let mut lhs = [[0.0; 3]; 3];
    let mut rhs = [0.0; 3];
    for idx in window.0..window.1 {
        let x = xs[idx];
        let y = ys[idx];
        let weight = weights[idx] * tricube_weight((x - x0).abs(), bandwidth);
        if weight <= 0.0 {
            continue;
        }
        let dx = x - x0;
        let dx2 = checked_option_product2(dx, dx)?;
        let dx3 = checked_option_product2(dx2, dx)?;
        let dx4 = checked_option_product2(dx2, dx2)?;
        let powers = [1.0, dx, dx2, dx3, dx4];
        for row in 0..size {
            rhs[row] =
                checked_option_sum2(rhs[row], checked_option_product3(weight, powers[row], y)?)?;
            for col in 0..size {
                lhs[row][col] = checked_option_sum2(
                    lhs[row][col],
                    checked_option_product2(weight, powers[row + col])?,
                )?;
            }
        }
    }
    solve_small_linear_system(lhs, rhs, size).map(|beta| beta[0])
}

fn local_weighted_mean(
    x0: f64,
    xs: &[f64],
    ys: &[f64],
    weights: &[f64],
    window: (usize, usize),
    bandwidth: f64,
) -> Option<f64> {
    let mut total_weight = 0.0;
    let mut mean = 0.0;
    for idx in window.0..window.1 {
        let weight = checked_option_product2(
            weights[idx],
            tricube_weight((xs[idx] - x0).abs(), bandwidth),
        )?;
        if weight <= 0.0 {
            continue;
        }
        let next_total = checked_option_sum2(total_weight, weight)?;
        let delta = checked_option_sub2(ys[idx], mean)?;
        let step = checked_option_product2(checked_option_div2(weight, next_total)?, delta)?;
        mean = checked_option_sum2(mean, step)?;
        total_weight = next_total;
    }
    (total_weight > 0.0 && mean.is_finite()).then_some(mean)
}

fn checked_option_sum2(left: f64, right: f64) -> Option<f64> {
    let sum = left + right;
    (left.is_finite() && right.is_finite() && sum.is_finite()).then_some(sum)
}

fn checked_option_product2(left: f64, right: f64) -> Option<f64> {
    let product = left * right;
    (left.is_finite() && right.is_finite() && product.is_finite()).then_some(product)
}

fn checked_option_sub2(left: f64, right: f64) -> Option<f64> {
    let difference = left - right;
    (left.is_finite() && right.is_finite() && difference.is_finite()).then_some(difference)
}

fn checked_option_div2(left: f64, right: f64) -> Option<f64> {
    let quotient = left / right;
    (left.is_finite() && right.is_finite() && right != 0.0 && quotient.is_finite())
        .then_some(quotient)
}

fn checked_option_product3(left: f64, middle: f64, right: f64) -> Option<f64> {
    checked_option_product2(checked_option_product2(left, middle)?, right)
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
            if !value.is_finite() {
                return None;
            }
        }
        rhs[pivot] /= pivot_value;
        if !rhs[pivot].is_finite() {
            return None;
        }
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
                let product = checked_option_product2(factor, *pivot_entry)?;
                lhs[row][col] = checked_option_sub2(lhs[row][col], product)?;
            }
            let rhs_product = checked_option_product2(factor, rhs[pivot])?;
            rhs[row] = checked_option_sub2(rhs[row], rhs_product)?;
        }
    }
    rhs.iter()
        .take(size)
        .all(|value| value.is_finite())
        .then_some(rhs)
}

fn parametric_coefficient_change(
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

fn validate_trend_fit_rows(
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

fn stable_midpoint(left: f64, right: f64) -> f64 {
    left / 2.0 + right / 2.0
}

fn validate_positive_mean(mean: f64, index: Option<usize>) -> Result<(), DeseqError> {
    if !mean.is_finite() || mean <= 0.0 {
        return Err(DeseqError::NonFiniteValue {
            context: "dispersion trend mean".to_string(),
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

    #[test]
    fn gamma_deviance_rejects_nonfinite_accumulation() {
        let trend = ParametricDispersionTrend {
            asympt_disp: f64::MIN_POSITIVE,
            extra_pois: f64::MIN_POSITIVE,
        };
        let means = [1.0, 2.0];
        let disps = [f64::MAX, f64::MAX];
        let good = [true, true];

        let err = gamma_deviance(&means, &disps, &good, trend).unwrap_err();

        assert!(err
            .to_string()
            .contains("parametric dispersion gamma deviance"));
    }

    #[test]
    fn checked_div2_rejects_zero_and_nonfinite_quotients() {
        let zero_err = checked_div2(1.0, 0.0, "test division").unwrap_err();
        let overflow_err = checked_div2(f64::MAX, f64::MIN_POSITIVE, "test division").unwrap_err();

        assert!(zero_err.to_string().contains("test division"));
        assert!(overflow_err.to_string().contains("test division"));
    }

    #[test]
    fn stable_midpoint_handles_opposite_extreme_endpoints() {
        assert_eq!(stable_midpoint(-f64::MAX, f64::MAX), 0.0);
        assert_eq!(stable_midpoint(f64::MAX, f64::MAX), f64::MAX);
    }

    #[test]
    fn parametric_work_helpers_validate_aligned_rows() {
        let trend = ParametricDispersionTrend {
            asympt_disp: 0.05,
            extra_pois: 2.0,
        };
        let means = [20.0, 40.0];
        let disps = [0.15];
        let good = [true, true];

        assert!(robust_parametric_trend_rows(
            &means,
            &disps,
            trend,
            ParametricDispersionTrendOptions::default(),
        )
        .is_err());
        assert!(gamma_deviance(&means, &disps, &good, trend).is_err());
        assert!(weighted_gamma_identity_least_squares(&means, &disps, &good, trend).is_err());
        assert!(fitted_values_are_positive(&means, &[true], trend).is_err());
    }

    #[test]
    fn weighted_gamma_identity_least_squares_rejects_nonfinite_coefficients() {
        let trend = ParametricDispersionTrend {
            asympt_disp: 0.05,
            extra_pois: 2.0,
        };
        let means = [20.0, 40.0];
        let disps = [f64::MAX, f64::MAX];
        let good = [true, true];

        assert!(weighted_gamma_identity_least_squares(&means, &disps, &good, trend).is_err());
    }

    #[test]
    fn weighted_gamma_identity_least_squares_rejects_nonfinite_weights() {
        let trend = ParametricDispersionTrend {
            asympt_disp: f64::MIN_POSITIVE,
            extra_pois: f64::MIN_POSITIVE,
        };
        let means = [1.0, 2.0];
        let disps = [0.1, 0.2];
        let good = [true, true];

        let err = weighted_gamma_identity_least_squares(&means, &disps, &good, trend).unwrap_err();

        assert!(err
            .to_string()
            .contains("parametric dispersion weighted least-squares"));
    }

    #[test]
    fn weighted_gamma_identity_least_squares_rejects_nonfinite_determinant_product() {
        let trend = ParametricDispersionTrend {
            asympt_disp: 1.0,
            extra_pois: f64::MIN_POSITIVE,
        };
        let means = [1.0e-154, 2.0e-154];
        let disps = [0.1, 0.2];
        let good = [true, true];

        let err = weighted_gamma_identity_least_squares(&means, &disps, &good, trend).unwrap_err();

        assert!(err
            .to_string()
            .contains("parametric dispersion weighted least-squares determinant"));
    }

    #[test]
    fn weighted_gamma_identity_least_squares_rejects_nonfinite_numerator_product() {
        let trend = ParametricDispersionTrend {
            asympt_disp: 1.0,
            extra_pois: f64::MIN_POSITIVE,
        };
        let means = [1.0, 0.5];
        let disps = [f64::MAX / 4.0, f64::MAX / 4.0];
        let good = [true, true];

        let err = weighted_gamma_identity_least_squares(&means, &disps, &good, trend).unwrap_err();

        assert!(err
            .to_string()
            .contains("parametric dispersion weighted least-squares coefficient numerator"));
    }

    #[test]
    fn checked_add_assign_rejects_nonfinite_sum() {
        let mut total = f64::MAX;

        let err = checked_add_assign(&mut total, f64::MAX, "test accumulation").unwrap_err();

        assert!(err.to_string().contains("test accumulation"));
        assert_eq!(total, f64::MAX);
    }

    #[test]
    fn checked_product_difference_rejects_nonfinite_product() {
        let err = checked_product_difference(f64::MAX, 2.0, 1.0, 1.0, "test product").unwrap_err();

        assert!(err.to_string().contains("test product"));
    }

    #[test]
    fn checked_product_difference_matches_finite_difference() {
        let observed = checked_product_difference(6.0, 7.0, 5.0, 3.0, "test product").unwrap();

        assert_eq!(observed, 27.0);
    }

    #[test]
    fn interpolate_parametric_coefficient_avoids_difference_overflow() {
        let interpolated = interpolate_parametric_coefficient(f64::MAX, -f64::MAX, 0.5);

        assert_eq!(interpolated, 0.0);
    }

    #[test]
    fn parametric_coefficient_change_rejects_nonfinite_ratio() {
        let old = ParametricDispersionTrend {
            asympt_disp: f64::MIN_POSITIVE,
            extra_pois: 1.0,
        };
        let new = ParametricDispersionTrend {
            asympt_disp: f64::MAX,
            extra_pois: 1.0,
        };

        let err = parametric_coefficient_change(old, new).unwrap_err();

        assert!(err
            .to_string()
            .contains("parametric dispersion coefficient change"));
    }

    #[test]
    fn parametric_coefficient_change_matches_log_squared_update() {
        let old = ParametricDispersionTrend {
            asympt_disp: 0.1,
            extra_pois: 2.0,
        };
        let new = ParametricDispersionTrend {
            asympt_disp: 0.2,
            extra_pois: 1.0,
        };

        let observed = parametric_coefficient_change(old, new).unwrap();
        let expected = 2.0 * 2.0_f64.ln().powi(2);

        assert_relative_eq!(observed, expected, epsilon = 1e-15);
    }

    #[test]
    fn small_linear_solver_rejects_nonfinite_solution() {
        let lhs = [[1.0, 0.0, 0.0], [0.0; 3], [0.0; 3]];
        let rhs = [f64::INFINITY, 0.0, 0.0];

        assert!(solve_small_linear_system(lhs, rhs, 1).is_none());
    }

    #[test]
    fn local_weighted_mean_rejects_nonfinite_online_update() {
        let xs = [0.0, 0.0];
        let ys = [f64::MAX, -f64::MAX];
        let weights = [1.0, 1.0];

        assert!(local_weighted_mean(0.0, &xs, &ys, &weights, (0, 2), 1.0).is_none());
    }

    #[test]
    fn small_linear_solver_rejects_nonfinite_pivot_row_normalization() {
        let lhs = [[2.0e-12, f64::MAX, 0.0], [0.0; 3], [0.0; 3]];
        let rhs = [1.0, 0.0, 0.0];

        assert!(solve_small_linear_system(lhs, rhs, 2).is_none());
    }

    #[test]
    fn small_linear_solver_rejects_nonfinite_elimination_update() {
        let lhs = [[1.0, f64::MAX, 0.0], [1.0, -f64::MAX, 0.0], [0.0; 3]];
        let rhs = [1.0, 2.0, 0.0];

        assert!(solve_small_linear_system(lhs, rhs, 2).is_none());
    }
}
