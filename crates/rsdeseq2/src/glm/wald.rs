use crate::probability::{normal_two_sided_tail, normal_upper_tail, student_t_upper_tail};
use statrs::distribution::{Continuous, Normal};

use crate::errors::{DeseqError, invalid_dimensions};
use crate::glm::inference::GlmInferenceView;
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
/// contain per-gene beta covariance matrices on log2 scale for multi-coefficient
/// contrasts. A single active coefficient uses its supplied marginal SE.
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
    wald_test_contrast_view(&GlmInferenceView::from(fit), contrast, options)
}

pub(crate) fn wald_test_contrast_view(
    fit: &GlmInferenceView<'_>,
    contrast: &[f64],
    options: &WaldTestOptions,
) -> Result<WaldContrastOutput, DeseqError> {
    use super::contrast_projection::wald_contrast_with_context;
    let estimates = fit.estimates()?;
    let residual = fit.model_matrix.n_samples() as f64 - fit.model_matrix.n_coefficients() as f64;
    wald_contrast_with_context(&estimates, contrast, options, Some(residual))
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
    let p = fit.beta.n_cols();
    let beta = fit.beta.as_slice();
    let se = fit.beta_se.as_slice();
    let rows = (0..fit.beta.n_rows())
        .map(|gene| (beta[gene * p + coefficient], se[gene * p + coefficient]));
    let residual = fit.model_matrix.n_samples() as f64 - fit.model_matrix.n_coefficients() as f64;
    super::inference_primitives::wald_from_estimate_rows(rows, options, Some(residual))
}

pub(super) fn checked_scaled_sum(values: &[f64]) -> Option<f64> {
    let mut scale = 0.0_f64;
    for value in values.iter().copied() {
        if !value.is_finite() {
            return None;
        }
        scale = scale.max(value.abs());
    }
    if scale == 0.0 {
        return Some(0.0);
    }
    let mut normalized_sum = 0.0;
    for value in values.iter().copied() {
        let term = value / scale;
        let next = normalized_sum + term;
        if !term.is_finite() || !next.is_finite() {
            return None;
        }
        normalized_sum = next;
    }
    let sum = normalized_sum * scale;
    sum.is_finite().then_some(sum)
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
    let tail = WaldTail::new(gene, options, degrees_of_freedom)?.unwrap_or(WaldTail::Missing);
    let threshold = options.lfc_threshold;
    let abs_beta = beta.abs();
    let (stat, pvalue) = match options.alternative {
        WaldAlternative::GreaterAbs if threshold == 0.0 => {
            (default_stat, tail.two_sided(default_stat))
        }
        WaldAlternative::GreaterAbs => {
            let q1 = checked_shifted_stat(-abs_beta, threshold, beta_se);
            let q2 = checked_shifted_stat(-abs_beta, -threshold, beta_se);
            (
                default_stat,
                q1.and_then(|q| tail.lower(q))
                    .zip(q2.and_then(|q| tail.lower(q)))
                    .and_then(|(a, b)| clamp_probability(a + b)),
            )
        }
        WaldAlternative::GreaterAbsUpshot if threshold == 0.0 => {
            (default_stat, tail.two_sided(default_stat))
        }
        WaldAlternative::GreaterAbsUpshot => {
            let pvalue = match tail {
                WaldTail::Normal => greater_abs_upshot_normal_pvalue(abs_beta, beta_se, threshold),
                WaldTail::T { .. } | WaldTail::Missing => None,
            };
            (default_stat, pvalue)
        }
        WaldAlternative::GreaterAbs2014 => {
            let Some(shifted) = checked_shifted_stat(abs_beta, -threshold, beta_se) else {
                return Ok(None);
            };
            let stat = beta.signum() * shifted.max(0.0);
            let pvalue = tail
                .upper(shifted)
                .and_then(|pvalue| checked_product2(2.0, pvalue).and_then(clamp_probability));
            (stat, pvalue)
        }
        WaldAlternative::LessAbs => {
            let Some(above_shift) = checked_shifted_stat(threshold, -beta, beta_se) else {
                return Ok(None);
            };
            let Some(below_shift) = checked_shifted_stat(beta, threshold, beta_se) else {
                return Ok(None);
            };
            let above = above_shift.max(0.0);
            let below = below_shift.max(0.0);
            let pvalue_above = tail.upper(above_shift);
            let pvalue_below = tail.upper(below_shift);
            (
                above.min(below),
                pvalue_above.zip(pvalue_below).map(|(a, b)| a.max(b)),
            )
        }
        WaldAlternative::Greater => {
            let Some(shifted) = checked_shifted_stat(beta, -threshold, beta_se) else {
                return Ok(None);
            };
            (shifted.max(0.0), tail.upper(shifted))
        }
        WaldAlternative::Less => {
            let Some(shifted_stat) = checked_shifted_stat(beta, threshold, beta_se) else {
                return Ok(None);
            };
            let Some(shifted_pvalue) = checked_shifted_stat(-threshold, -beta, beta_se) else {
                return Ok(None);
            };
            (shifted_stat.min(0.0), tail.upper(shifted_pvalue))
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
    normal_two_sided_tail(stat)
}

/// Two-sided Student t p-value.
pub fn two_sided_t_pvalue(stat: f64, degrees_of_freedom: f64) -> Option<f64> {
    if !stat.is_finite() || !degrees_of_freedom.is_finite() || degrees_of_freedom <= 0.0 {
        return None;
    }
    clamp_probability(2.0 * student_t_upper_tail(stat.abs(), degrees_of_freedom)?)
}

pub(super) fn valid_t_df(df: f64) -> Option<f64> {
    (df.is_finite() && df > 0.0).then_some(df)
}

pub(crate) fn validate_wald_options(options: &WaldTestOptions) -> Result<(), DeseqError> {
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
    let a = checked_shifted_stat(abs_beta, threshold, beta_se)?;
    let b = checked_shifted_stat(abs_beta, -threshold, beta_se)?;
    let denominator = checked_sub(b, a)?;
    let scale = checked_div(2.0, denominator)?;
    let left_tail_a = normal_upper_tail(a);
    let left_tail_b = normal_upper_tail(b);
    let value = checked_product2(
        scale,
        checked_scaled_sum(&[
            checked_product2(-a, left_tail_a)?,
            distribution.pdf(a),
            checked_product2(b, left_tail_b)?,
            -distribution.pdf(b),
        ])?,
    )?;
    clamp_probability(value)
}

fn clamp_probability(value: f64) -> Option<f64> {
    value.is_finite().then_some(value.clamp(0.0, 1.0))
}

fn checked_sub(left: f64, right: f64) -> Option<f64> {
    let value = left - right;
    value.is_finite().then_some(value)
}

fn checked_add(left: f64, right: f64) -> Option<f64> {
    let value = left + right;
    (left.is_finite() && right.is_finite() && value.is_finite()).then_some(value)
}

fn checked_div(left: f64, right: f64) -> Option<f64> {
    let value = left / right;
    (left.is_finite() && right.is_finite() && right != 0.0 && value.is_finite()).then_some(value)
}

fn checked_shifted_stat(left: f64, right: f64, beta_se: f64) -> Option<f64> {
    checked_div(checked_add(left, right)?, beta_se)
}

fn checked_product2(left: f64, right: f64) -> Option<f64> {
    let value = left * right;
    (left.is_finite() && right.is_finite() && value.is_finite()).then_some(value)
}

pub(super) fn checked_product3_ordered_by_magnitude(
    left: f64,
    middle: f64,
    right: f64,
) -> Option<f64> {
    let mut factors = [left, middle, right];
    factors.sort_by(|a, b| a.abs().total_cmp(&b.abs()));
    checked_product2(checked_product2(factors[0], factors[1])?, factors[2])
}

#[derive(Clone, Copy, Debug)]
enum WaldTail {
    Missing,
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
            Self::Missing => None,
            Self::Normal => clamp_probability(normal_upper_tail(-q)),
            Self::T { degrees_of_freedom } => student_t_upper_tail(-q, degrees_of_freedom),
        }
    }

    fn upper(self, q: f64) -> Option<f64> {
        if !q.is_finite() {
            return None;
        }
        match self {
            Self::Missing => None,
            Self::Normal => clamp_probability(normal_upper_tail(q)),
            Self::T { degrees_of_freedom } => student_t_upper_tail(q, degrees_of_freedom),
        }
    }

    fn two_sided(self, stat: f64) -> Option<f64> {
        match self {
            Self::Missing => None,
            Self::Normal => Some(two_sided_normal_pvalue(stat)),
            Self::T { degrees_of_freedom } => two_sided_t_pvalue(stat, degrees_of_freedom),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn checked_shifted_stat_rejects_overflowed_threshold_sum() {
        assert_eq!(checked_shifted_stat(f64::MAX, f64::MAX, 1.0), None);
        assert_eq!(checked_shifted_stat(4.0, -1.0, 2.0), Some(1.5));
    }

    #[test]
    fn thresholded_wald_returns_missing_for_overflowed_shift() {
        let options =
            WaldTestOptions::normal().with_lfc_threshold(f64::MAX, WaldAlternative::LessAbs);
        let result = wald_stat_and_pvalue_with_options(-f64::MAX, 1.0, 0, &options, None).unwrap();

        assert!(result.is_none());
    }
}
