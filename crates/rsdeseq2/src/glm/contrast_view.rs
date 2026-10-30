//! Prepared borrowed inputs reused across independent contrast requests.
use super::contrast_projection::wald_contrast_with_context;
use super::{ContrastProjection, InferenceEstimates, WaldOutput, WaldTestOptions};
use crate::contrasts::{ContrastCounts, ContrastZeroPolicy};
use crate::errors::{DeseqError, invalid_dimensions};
use crate::matrix::MatrixView;

/// Contrast-specific zero rule, separate from coefficient weights.
#[derive(Clone, Copy, Debug)]
pub enum ContrastZeroRule<'a> {
    /// Named coefficients and inputs without count context have no override.
    None,
    /// Numeric/list contrasts use the nonzero design-column pattern.
    Numeric,
    /// Character contrasts select the two requested sample levels.
    Factor {
        /// One level per count sample.
        levels: &'a [&'a str],
        /// Numerator level or its supported alias.
        numerator: &'a str,
        /// Denominator level or its supported alias.
        denominator: &'a str,
    },
}

/// One owned result from a prepared contrast view.
#[derive(Clone, Debug, PartialEq)]
pub struct ContrastInference {
    /// Projected effect and SE columns.
    pub projection: ContrastProjection,
    /// Wald inference; absent for an LRT display-only projection.
    pub wald: Option<WaldOutput>,
}

/// Shape-checked inference and count inputs for one or many requests.
///
/// Coefficient/covariance/count matrices remain borrowed. Global-zero flags are
/// prepared once. Each evaluation owns only its gene-length results and scratch;
/// callers can consume results in order without retaining a second result batch.
#[derive(Debug)]
pub struct PreparedContrastView<'a> {
    estimates: InferenceEstimates<'a>,
    counts: Option<ContrastCounts<'a>>,
    design: Option<MatrixView<'a, f64>>,
}

impl<'a> PreparedContrastView<'a> {
    /// Validate shared shapes once, preserving primitive inputs with no design
    /// or count context. Actual requests diagnose unavailable variance/design.
    pub fn new(
        estimates: InferenceEstimates<'a>,
        counts: Option<ContrastCounts<'a>>,
        design: Option<MatrixView<'a, f64>>,
    ) -> Result<Self, DeseqError> {
        if let Some(counts) = &counts {
            if counts.n_genes() != estimates.n_genes() {
                return Err(invalid_dimensions(
                    "contrast counts genes",
                    estimates.n_genes(),
                    counts.n_genes(),
                ));
            }
            if let Some(design) = design
                && counts.n_samples() != design.n_rows()
            {
                return Err(invalid_dimensions(
                    "contrast design samples",
                    counts.n_samples(),
                    design.n_rows(),
                ));
            }
        }
        if let Some(design) = design
            && design.n_cols() != estimates.n_coefficients()
        {
            return Err(invalid_dimensions(
                "contrast design coefficients",
                estimates.n_coefficients(),
                design.n_cols(),
            ));
        }
        Ok(Self {
            estimates,
            counts,
            design,
        })
    }

    /// Evaluate one request through the shared projection, Wald and zero policy.
    /// Passing no Wald options requests only the LRT display-effect columns.
    /// No statistical decisions from earlier requests are reused.
    pub fn evaluate(
        &self,
        contrast: &[f64],
        options: Option<&WaldTestOptions>,
        zero: ContrastZeroRule<'_>,
    ) -> Result<ContrastInference, DeseqError> {
        let mut result = if let Some(options) = options {
            let residual = self.design.map(|d| d.n_rows() as f64 - d.n_cols() as f64);
            let output = wald_contrast_with_context(&self.estimates, contrast, options, residual)?;
            ContrastInference {
                projection: ContrastProjection {
                    log2_fold_change: output.log2_fold_change,
                    lfc_se: output.lfc_se,
                },
                wald: Some(output.wald),
            }
        } else {
            ContrastInference {
                projection: self.estimates.project_contrast(contrast)?,
                wald: None,
            }
        };
        let Some(counts) = &self.counts else {
            return Ok(result);
        };
        let flags = match zero {
            ContrastZeroRule::None => return Ok(result),
            ContrastZeroRule::Numeric => counts.numeric_flags(self.design, contrast)?,
            ContrastZeroRule::Factor {
                levels,
                numerator,
                denominator,
            } => counts.factor_flags(levels, numerator, denominator)?,
        };
        let policy = ContrastZeroPolicy::new(&flags, counts.all_zero())?;
        if let Some(wald) = &mut result.wald {
            policy.apply_wald(
                &mut result.projection.log2_fold_change,
                &mut wald.stat,
                &mut wald.pvalue,
            )?;
        } else {
            policy.apply_effects(result.projection.log2_fold_change.iter_mut())?;
        }
        Ok(result)
    }
}
