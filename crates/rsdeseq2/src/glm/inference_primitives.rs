//! Inference from borrowed estimates, without reconstructing a fitted model.

use super::wald::{valid_t_df, validate_wald_options, wald_stat_and_pvalue_with_options};
use super::{WaldDegreesOfFreedom, WaldOutput, WaldPvalueType, WaldTestOptions};
use crate::errors::{DeseqError, invalid_dimensions};

/// Evaluate Wald hypotheses from aligned effect estimates and standard errors.
///
/// Uses the same scalar inference as fitted coefficient and contrast queries.
/// Nonfinite estimates and nonpositive standard errors produce missing rows.
/// Missing or nonpositive t degrees of freedom produce missing p-values while
/// retaining the hypothesis statistic. Residual degrees of freedom require a
/// fitted design; primitive callers must supply scalar or per-row values.
pub fn wald_test_from_estimates(
    beta: &[f64],
    beta_se: &[f64],
    options: &WaldTestOptions,
) -> Result<WaldOutput, DeseqError> {
    if beta.len() != beta_se.len() {
        return Err(invalid_dimensions(
            "Wald estimate standard errors",
            beta.len(),
            beta_se.len(),
        ));
    }
    validate_wald_options(options)?;
    let degrees_of_freedom = match &options.pvalue_type {
        WaldPvalueType::Normal => None,
        WaldPvalueType::T { degrees_of_freedom } => Some(match degrees_of_freedom {
            WaldDegreesOfFreedom::Residual => {
                return Err(DeseqError::InvalidOptions {
                    reason: "primitive Wald inference requires explicit degrees of freedom".into(),
                });
            }
            WaldDegreesOfFreedom::Scalar(df) => vec![valid_t_df(*df); beta.len()],
            WaldDegreesOfFreedom::PerGene(df) => {
                if df.len() != beta.len() {
                    return Err(invalid_dimensions(
                        "Wald per-gene degrees of freedom",
                        beta.len(),
                        df.len(),
                    ));
                }
                df.iter().copied().map(valid_t_df).collect()
            }
        }),
    };
    let mut stat = Vec::with_capacity(beta.len());
    let mut pvalue = Vec::with_capacity(beta.len());
    for (gene, (&beta, &se)) in beta.iter().zip(beta_se).enumerate() {
        let result = wald_stat_and_pvalue_with_options(
            beta,
            se,
            gene,
            options,
            degrees_of_freedom.as_deref(),
        )?;
        stat.push(result.map(|(stat, _)| stat));
        pvalue.push(result.and_then(|(_, pvalue)| pvalue));
    }
    Ok(WaldOutput {
        stat,
        pvalue,
        degrees_of_freedom,
    })
}
