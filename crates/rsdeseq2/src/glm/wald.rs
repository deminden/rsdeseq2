use statrs::distribution::{Continuous, ContinuousCDF, Normal, StudentsT};
use statrs::function::erf::erfc;

use crate::errors::{invalid_dimensions, DeseqError};
use crate::glm::{NbinomGlmFit, WaldOutput};

/// Options controlling Wald p-value calculation.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct WaldTestOptions {
    /// Null distribution used to convert Wald statistics into p-values.
    pub pvalue_type: WaldPvalueType,
    /// Alternative hypothesis used by DESeq2's `results(lfcThreshold=...)` path.
    pub alternative: WaldAlternative,
    /// Non-negative log2 fold-change threshold.
    pub lfc_threshold: f64,
}

impl WaldTestOptions {
    /// DESeq2 default `useT=FALSE`: standard-Normal p-values.
    pub fn normal() -> Self {
        Self {
            pvalue_type: WaldPvalueType::Normal,
            alternative: WaldAlternative::GreaterAbs,
            lfc_threshold: 0.0,
        }
    }

    /// DESeq2 `useT=TRUE` with residual degrees of freedom.
    pub fn t_residual_degrees_of_freedom() -> Self {
        Self {
            pvalue_type: WaldPvalueType::T {
                degrees_of_freedom: WaldDegreesOfFreedom::Residual,
            },
            alternative: WaldAlternative::GreaterAbs,
            lfc_threshold: 0.0,
        }
    }

    /// DESeq2 `useT=TRUE` with one caller-supplied df recycled over genes.
    pub fn t_degrees_of_freedom(df: f64) -> Self {
        Self {
            pvalue_type: WaldPvalueType::T {
                degrees_of_freedom: WaldDegreesOfFreedom::Scalar(df),
            },
            alternative: WaldAlternative::GreaterAbs,
            lfc_threshold: 0.0,
        }
    }

    /// DESeq2 `useT=TRUE` with one caller-supplied df per input gene.
    pub fn t_per_gene_degrees_of_freedom(df: Vec<f64>) -> Self {
        Self {
            pvalue_type: WaldPvalueType::T {
                degrees_of_freedom: WaldDegreesOfFreedom::PerGene(df),
            },
            alternative: WaldAlternative::GreaterAbs,
            lfc_threshold: 0.0,
        }
    }

    /// Set a DESeq2-style log2 fold-change threshold alternative.
    pub fn with_lfc_threshold(mut self, threshold: f64, alternative: WaldAlternative) -> Self {
        self.lfc_threshold = threshold;
        self.alternative = alternative;
        self
    }
}

/// Wald p-value null distribution.
#[derive(Clone, Debug, Default, PartialEq)]
pub enum WaldPvalueType {
    /// DESeq2 default `useT=FALSE`.
    #[default]
    Normal,
    /// DESeq2 `useT=TRUE`.
    T {
        /// Source of t degrees of freedom.
        degrees_of_freedom: WaldDegreesOfFreedom,
    },
}

/// Source of degrees of freedom for t-distribution Wald p-values.
#[derive(Clone, Debug, PartialEq)]
pub enum WaldDegreesOfFreedom {
    /// Use `n_samples - n_coefficients`.
    Residual,
    /// Recycle one df value over genes.
    Scalar(f64),
    /// Use one df value per full input gene row.
    PerGene(Vec<f64>),
}

/// DESeq2 Wald alternative hypothesis for thresholded result p-values.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum WaldAlternative {
    /// `|beta| > lfcThreshold`.
    #[default]
    GreaterAbs,
    /// UP-SHOT version of `|beta| > lfcThreshold`; unsupported with t p-values, matching DESeq2.
    GreaterAbsUpshot,
    /// Older 2014 implementation of `greaterAbs`.
    GreaterAbs2014,
    /// `|beta| < lfcThreshold`.
    LessAbs,
    /// `beta > lfcThreshold`.
    Greater,
    /// `beta < -lfcThreshold`.
    Less,
}

/// Wald output for an explicit linear contrast.
#[derive(Clone, Debug, PartialEq)]
pub struct WaldContrastOutput {
    /// Log2-scale contrast estimate `c' beta`.
    pub log2_fold_change: Vec<Option<f64>>,
    /// Log2-scale contrast standard error `sqrt(c' Sigma c)`.
    pub lfc_se: Vec<Option<f64>>,
    /// Wald statistic, p-value, and optional degrees of freedom.
    pub wald: WaldOutput,
}

/// Compute DESeq2-style Wald statistics for one coefficient with explicit options.
///
/// This is a convenience alias for the currently implemented selected-
/// coefficient Wald path. Use `wald_test_contrast_with_options` for explicit
/// linear contrasts.
pub fn wald_test(
    fit: &NbinomGlmFit,
    coefficient: usize,
    options: &WaldTestOptions,
) -> Result<WaldOutput, DeseqError> {
    wald_test_coefficient_with_options(fit, coefficient, options)
}

/// Compute default DESeq2-style Wald statistics for one coefficient.
///
/// This implements the default `useT=FALSE` path:
///
/// `stat = beta / betaSE`
///
/// `pvalue = 2 * pnorm(abs(stat), lower.tail = FALSE)`
pub fn wald_test_coefficient(
    fit: &NbinomGlmFit,
    coefficient: usize,
) -> Result<WaldOutput, DeseqError> {
    wald_test_coefficient_with_options(fit, coefficient, &WaldTestOptions::normal())
}

/// Compute DESeq2-style Wald statistics for an explicit coefficient contrast.
///
/// The contrast vector must have one finite value per coefficient. The fit must
/// contain per-gene beta covariance matrices on log2 scale.
pub fn wald_test_contrast(
    fit: &NbinomGlmFit,
    contrast: &[f64],
) -> Result<WaldContrastOutput, DeseqError> {
    wald_test_contrast_with_options(fit, contrast, &WaldTestOptions::normal())
}

/// Compute DESeq2-style Wald statistics for an explicit coefficient contrast
/// with explicit p-value options.
pub fn wald_test_contrast_with_options(
    fit: &NbinomGlmFit,
    contrast: &[f64],
    options: &WaldTestOptions,
) -> Result<WaldContrastOutput, DeseqError> {
    validate_contrast_inputs(fit, contrast)?;
    validate_wald_options(options)?;
    let degrees_of_freedom = resolve_wald_degrees_of_freedom(fit, options)?;
    let covariance =
        fit.beta_covariance
            .as_ref()
            .ok_or_else(|| DeseqError::UnsupportedFeature {
                feature: "Wald contrast requires beta covariance matrices".to_string(),
            })?;

    let mut log2_fold_change = Vec::with_capacity(fit.beta.n_rows());
    let mut lfc_se = Vec::with_capacity(fit.beta.n_rows());
    let mut stat = Vec::with_capacity(fit.beta.n_rows());
    let mut pvalue = Vec::with_capacity(fit.beta.n_rows());
    let p = fit.beta.n_cols();

    for gene in 0..fit.beta.n_rows() {
        let beta = fit.beta.row(gene)?;
        let covariance = covariance.row(gene)?;
        let Some(estimate) = contrast_estimate(contrast, beta) else {
            log2_fold_change.push(None);
            lfc_se.push(None);
            stat.push(None);
            pvalue.push(None);
            continue;
        };
        let Some(variance) = contrast_variance(contrast, covariance, p) else {
            log2_fold_change.push(None);
            lfc_se.push(None);
            stat.push(None);
            pvalue.push(None);
            continue;
        };
        let Some(se) = valid_contrast_se(estimate, variance) else {
            log2_fold_change.push(None);
            lfc_se.push(None);
            stat.push(None);
            pvalue.push(None);
            continue;
        };
        match wald_stat_and_pvalue_with_options(
            estimate,
            se,
            gene,
            options,
            degrees_of_freedom.as_deref(),
        )? {
            Some((gene_stat, gene_pvalue)) => {
                log2_fold_change.push(Some(estimate));
                lfc_se.push(Some(se));
                stat.push(Some(gene_stat));
                pvalue.push(gene_pvalue);
            }
            None => {
                log2_fold_change.push(None);
                lfc_se.push(None);
                stat.push(None);
                pvalue.push(None);
            }
        }
    }

    Ok(WaldContrastOutput {
        log2_fold_change,
        lfc_se,
        wald: WaldOutput {
            stat,
            pvalue,
            degrees_of_freedom,
        },
    })
}

/// Compute DESeq2-style Wald statistics for one coefficient with explicit p-value options.
pub fn wald_test_coefficient_with_options(
    fit: &NbinomGlmFit,
    coefficient: usize,
    options: &WaldTestOptions,
) -> Result<WaldOutput, DeseqError> {
    if coefficient >= fit.beta.n_cols() {
        return Err(DeseqError::InvalidDimensions {
            context: "Wald coefficient index".to_string(),
            expected: fit.beta.n_cols().saturating_sub(1),
            actual: coefficient,
        });
    }
    if fit.beta.n_rows() != fit.beta_se.n_rows() || fit.beta.n_cols() != fit.beta_se.n_cols() {
        return Err(invalid_dimensions(
            "Wald beta/betaSE matrix values",
            fit.beta.len(),
            fit.beta_se.len(),
        ));
    }
    validate_wald_options(options)?;
    let degrees_of_freedom = resolve_wald_degrees_of_freedom(fit, options)?;

    let mut stat = Vec::with_capacity(fit.beta.n_rows());
    let mut pvalue = Vec::with_capacity(fit.beta.n_rows());
    for gene in 0..fit.beta.n_rows() {
        let beta = fit.beta.row(gene)?[coefficient];
        let se = fit.beta_se.row(gene)?[coefficient];
        match wald_stat_and_pvalue_with_options(
            beta,
            se,
            gene,
            options,
            degrees_of_freedom.as_deref(),
        )? {
            Some((gene_stat, gene_pvalue)) => {
                stat.push(Some(gene_stat));
                pvalue.push(gene_pvalue);
            }
            None => {
                stat.push(None);
                pvalue.push(None);
            }
        }
    }
    Ok(WaldOutput {
        stat,
        pvalue,
        degrees_of_freedom,
    })
}

fn validate_contrast_inputs(fit: &NbinomGlmFit, contrast: &[f64]) -> Result<(), DeseqError> {
    if contrast.len() != fit.beta.n_cols() {
        return Err(invalid_dimensions(
            "Wald contrast coefficients",
            fit.beta.n_cols(),
            contrast.len(),
        ));
    }
    if fit.beta.n_rows() != fit.beta_se.n_rows() || fit.beta.n_cols() != fit.beta_se.n_cols() {
        return Err(invalid_dimensions(
            "Wald beta/betaSE matrix values",
            fit.beta.len(),
            fit.beta_se.len(),
        ));
    }
    if let Some(covariance) = &fit.beta_covariance {
        if covariance.n_rows() != fit.beta.n_rows()
            || covariance.n_cols() != fit.beta.n_cols() * fit.beta.n_cols()
        {
            return Err(DeseqError::InvalidDimensions {
                context: "Wald beta covariance matrix".to_string(),
                expected: fit.beta.n_rows() * fit.beta.n_cols() * fit.beta.n_cols(),
                actual: covariance.len(),
            });
        }
    }
    let mut any_nonzero = false;
    for (idx, value) in contrast.iter().copied().enumerate() {
        if !value.is_finite() {
            return Err(DeseqError::NonFiniteValue {
                context: "Wald contrast".to_string(),
                index: Some(idx),
                value,
            });
        }
        any_nonzero |= value != 0.0;
    }
    if !any_nonzero {
        return Err(DeseqError::InvalidOptions {
            reason: "Wald contrast must contain at least one non-zero coefficient".to_string(),
        });
    }
    Ok(())
}

fn contrast_estimate(contrast: &[f64], beta: &[f64]) -> Option<f64> {
    let mut estimate = 0.0;
    for (contrast, beta) in contrast.iter().copied().zip(beta.iter().copied()) {
        let term = contrast * beta;
        let next = estimate + term;
        if !term.is_finite() || !next.is_finite() {
            return None;
        }
        estimate = next;
    }
    Some(estimate)
}

fn contrast_variance(contrast: &[f64], covariance: &[f64], p: usize) -> Option<f64> {
    let mut variance = 0.0;
    for row in 0..p {
        for col in 0..p {
            let term = contrast[row] * covariance[row * p + col] * contrast[col];
            let next = variance + term;
            if !term.is_finite() || !next.is_finite() {
                return None;
            }
            variance = next;
        }
    }
    Some(variance)
}

fn valid_contrast_se(estimate: f64, variance: f64) -> Option<f64> {
    if !estimate.is_finite() || !variance.is_finite() || variance <= 0.0 {
        return None;
    }
    let se = variance.sqrt();
    se.is_finite().then_some(se)
}

/// Compute one Wald statistic and two-sided standard-Normal p-value.
pub fn wald_stat_and_pvalue(beta: f64, beta_se: f64) -> Option<(f64, f64)> {
    let stat = wald_stat(beta, beta_se)?;
    Some((stat, two_sided_normal_pvalue(stat)))
}

/// Compute one Wald statistic and p-value using explicit DESeq2-style options.
pub fn wald_stat_and_pvalue_with_options(
    beta: f64,
    beta_se: f64,
    gene: usize,
    options: &WaldTestOptions,
    degrees_of_freedom: Option<&[Option<f64>]>,
) -> Result<Option<(f64, Option<f64>)>, DeseqError> {
    validate_wald_options(options)?;
    let Some(default_stat) = wald_stat(beta, beta_se) else {
        return Ok(None);
    };
    let Some(tail) = WaldTail::new(gene, options, degrees_of_freedom)? else {
        return Ok(Some((default_stat, None)));
    };
    let threshold = options.lfc_threshold;
    let abs_beta = beta.abs();
    let (stat, pvalue) = match options.alternative {
        WaldAlternative::GreaterAbs if threshold == 0.0 => {
            (default_stat, tail.two_sided(default_stat))
        }
        WaldAlternative::GreaterAbs => {
            let q1 = (-abs_beta + threshold) / beta_se;
            let q2 = (-abs_beta - threshold) / beta_se;
            (
                default_stat,
                tail.lower(q1)
                    .zip(tail.lower(q2))
                    .and_then(|(a, b)| clamp_probability(a + b)),
            )
        }
        WaldAlternative::GreaterAbsUpshot if threshold == 0.0 => {
            (default_stat, tail.two_sided(default_stat))
        }
        WaldAlternative::GreaterAbsUpshot => {
            let pvalue = match tail {
                WaldTail::Normal => greater_abs_upshot_normal_pvalue(abs_beta, beta_se, threshold),
                WaldTail::T { .. } => None,
            };
            (default_stat, pvalue)
        }
        WaldAlternative::GreaterAbs2014 => {
            let shifted = (abs_beta - threshold) / beta_se;
            let stat = beta.signum() * shifted.max(0.0);
            let pvalue = tail
                .upper(shifted)
                .and_then(|pvalue| clamp_probability(2.0 * pvalue));
            (stat, pvalue)
        }
        WaldAlternative::LessAbs => {
            let above = ((threshold - beta) / beta_se).max(0.0);
            let below = ((beta + threshold) / beta_se).max(0.0);
            let pvalue_above = tail.upper((threshold - beta) / beta_se);
            let pvalue_below = tail.upper((beta + threshold) / beta_se);
            (
                above.min(below),
                pvalue_above.zip(pvalue_below).map(|(a, b)| a.max(b)),
            )
        }
        WaldAlternative::Greater => {
            let shifted = (beta - threshold) / beta_se;
            (shifted.max(0.0), tail.upper(shifted))
        }
        WaldAlternative::Less => {
            let shifted_stat = ((beta + threshold) / beta_se).min(0.0);
            let shifted_pvalue = (-threshold - beta) / beta_se;
            (shifted_stat, tail.upper(shifted_pvalue))
        }
    };
    Ok(Some((stat, pvalue)))
}

fn wald_stat(beta: f64, beta_se: f64) -> Option<f64> {
    if !beta.is_finite() || !beta_se.is_finite() || beta_se <= 0.0 {
        return None;
    }
    let stat = beta / beta_se;
    if !stat.is_finite() {
        return None;
    }
    Some(stat)
}

/// Two-sided standard-Normal p-value.
pub fn two_sided_normal_pvalue(stat: f64) -> f64 {
    if stat.is_nan() {
        return 1.0;
    }
    erfc(stat.abs() / std::f64::consts::SQRT_2).clamp(0.0, 1.0)
}

/// Two-sided Student t p-value.
pub fn two_sided_t_pvalue(stat: f64, degrees_of_freedom: f64) -> Option<f64> {
    if !stat.is_finite() || !degrees_of_freedom.is_finite() || degrees_of_freedom <= 0.0 {
        return None;
    }
    let distribution = StudentsT::new(0.0, 1.0, degrees_of_freedom).ok()?;
    clamp_probability(2.0 * (1.0 - distribution.cdf(stat.abs())))
}

fn resolve_wald_degrees_of_freedom(
    fit: &NbinomGlmFit,
    options: &WaldTestOptions,
) -> Result<Option<Vec<Option<f64>>>, DeseqError> {
    match &options.pvalue_type {
        WaldPvalueType::Normal => Ok(None),
        WaldPvalueType::T { degrees_of_freedom } => match degrees_of_freedom {
            WaldDegreesOfFreedom::Residual => {
                let df =
                    fit.model_matrix.n_samples() as f64 - fit.model_matrix.n_coefficients() as f64;
                Ok(Some(vec![valid_t_df(df); fit.beta.n_rows()]))
            }
            WaldDegreesOfFreedom::Scalar(df) => Ok(Some(vec![valid_t_df(*df); fit.beta.n_rows()])),
            WaldDegreesOfFreedom::PerGene(df) => {
                if df.len() != fit.beta.n_rows() {
                    return Err(invalid_dimensions(
                        "Wald per-gene degrees of freedom",
                        fit.beta.n_rows(),
                        df.len(),
                    ));
                }
                Ok(Some(df.iter().copied().map(valid_t_df).collect()))
            }
        },
    }
}

fn valid_t_df(df: f64) -> Option<f64> {
    (df.is_finite() && df > 0.0).then_some(df)
}

fn validate_wald_options(options: &WaldTestOptions) -> Result<(), DeseqError> {
    if !options.lfc_threshold.is_finite() || options.lfc_threshold < 0.0 {
        return Err(DeseqError::InvalidOptions {
            reason: "Wald LFC threshold must be finite and non-negative".to_string(),
        });
    }
    if options.alternative == WaldAlternative::LessAbs && options.lfc_threshold == 0.0 {
        return Err(DeseqError::UnsupportedFeature {
            feature: "altHypothesis='lessAbs' requires a positive lfcThreshold".to_string(),
        });
    }
    if matches!(options.pvalue_type, WaldPvalueType::T { .. })
        && options.alternative == WaldAlternative::GreaterAbsUpshot
    {
        return Err(DeseqError::UnsupportedFeature {
            feature: "greaterAbsUPSHOT with useT=TRUE".to_string(),
        });
    }
    Ok(())
}

fn greater_abs_upshot_normal_pvalue(abs_beta: f64, beta_se: f64, threshold: f64) -> Option<f64> {
    if threshold == 0.0 {
        return Some(two_sided_normal_pvalue(abs_beta / beta_se));
    }
    let distribution = Normal::new(0.0, 1.0).ok()?;
    let a = (abs_beta + threshold) / beta_se;
    let b = (abs_beta - threshold) / beta_se;
    let value = 2.0 / (b - a)
        * (-a * distribution.cdf(-a) + distribution.pdf(a) + b * distribution.cdf(-b)
            - distribution.pdf(b));
    clamp_probability(value)
}

fn clamp_probability(value: f64) -> Option<f64> {
    value.is_finite().then_some(value.clamp(0.0, 1.0))
}

#[derive(Clone, Copy, Debug)]
enum WaldTail {
    Normal,
    T { degrees_of_freedom: f64 },
}

impl WaldTail {
    fn new(
        gene: usize,
        options: &WaldTestOptions,
        degrees_of_freedom: Option<&[Option<f64>]>,
    ) -> Result<Option<Self>, DeseqError> {
        match options.pvalue_type {
            WaldPvalueType::Normal => Ok(Some(Self::Normal)),
            WaldPvalueType::T { .. } => {
                let Some(degrees_of_freedom) = degrees_of_freedom
                    .and_then(|df| df.get(gene))
                    .copied()
                    .flatten()
                else {
                    return Ok(None);
                };
                Ok(Some(Self::T { degrees_of_freedom }))
            }
        }
    }

    fn lower(self, q: f64) -> Option<f64> {
        if !q.is_finite() {
            return None;
        }
        match self {
            Self::Normal => {
                let distribution = Normal::new(0.0, 1.0).ok()?;
                clamp_probability(distribution.cdf(q))
            }
            Self::T { degrees_of_freedom } => {
                let distribution = StudentsT::new(0.0, 1.0, degrees_of_freedom).ok()?;
                clamp_probability(distribution.cdf(q))
            }
        }
    }

    fn upper(self, q: f64) -> Option<f64> {
        self.lower(q).and_then(|cdf| clamp_probability(1.0 - cdf))
    }

    fn two_sided(self, stat: f64) -> Option<f64> {
        match self {
            Self::Normal => Some(two_sided_normal_pvalue(stat)),
            Self::T { degrees_of_freedom } => two_sided_t_pvalue(stat, degrees_of_freedom),
        }
    }
}
