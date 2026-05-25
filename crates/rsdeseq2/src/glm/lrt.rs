use statrs::distribution::{ChiSquared, ContinuousCDF};

use crate::errors::{invalid_dimensions, DeseqError};
use crate::glm::{LrtOutput, NbinomGlmFit};

/// Compute DESeq2-style likelihood-ratio statistics from full and reduced fits.
///
/// This mirrors `2 * (fullModel$logLike - reducedModel$logLike)` with a
/// chi-square upper-tail p-value using `df = ncol(full) - ncol(reduced)`.
pub fn lrt_test(full: &NbinomGlmFit, reduced: &NbinomGlmFit) -> Result<LrtOutput, DeseqError> {
    if full.log_like.len() != reduced.log_like.len() {
        return Err(invalid_dimensions(
            "LRT log-likelihood rows",
            full.log_like.len(),
            reduced.log_like.len(),
        ));
    }
    if full.beta.n_rows() != reduced.beta.n_rows() {
        return Err(invalid_dimensions(
            "LRT beta rows",
            full.beta.n_rows(),
            reduced.beta.n_rows(),
        ));
    }
    if reduced.beta_converged.len() != reduced.beta.n_rows() {
        return Err(invalid_dimensions(
            "LRT reduced convergence flags",
            reduced.beta.n_rows(),
            reduced.beta_converged.len(),
        ));
    }
    if full.beta.n_cols() <= reduced.beta.n_cols() {
        return Err(DeseqError::InvalidDimensions {
            context: "LRT degrees of freedom".to_string(),
            expected: reduced.beta.n_cols() + 1,
            actual: full.beta.n_cols(),
        });
    }

    let degrees_of_freedom = full.beta.n_cols() - reduced.beta.n_cols();
    let distribution =
        ChiSquared::new(degrees_of_freedom as f64).map_err(|error| DeseqError::InvalidCounts {
            reason: format!("invalid LRT chi-square distribution: {error}"),
        })?;
    let mut deviance = Vec::with_capacity(full.log_like.len());
    let mut pvalue = Vec::with_capacity(full.log_like.len());
    for (full_log_like, reduced_log_like) in full
        .log_like
        .iter()
        .copied()
        .zip(reduced.log_like.iter().copied())
    {
        if let Some(statistic) = lrt_deviance_statistic(full_log_like, reduced_log_like) {
            deviance.push(Some(statistic));
            pvalue.push(lrt_pvalue(&distribution, statistic));
        } else {
            deviance.push(None);
            pvalue.push(None);
        }
    }
    Ok(LrtOutput {
        deviance,
        pvalue,
        degrees_of_freedom,
        reduced_converged: reduced.beta_converged.clone(),
    })
}

fn lrt_deviance_statistic(full_log_like: f64, reduced_log_like: f64) -> Option<f64> {
    if !full_log_like.is_finite() || !reduced_log_like.is_finite() {
        return None;
    }
    let statistic = checked_product2(2.0, checked_sub(full_log_like, reduced_log_like)?)?;
    statistic.is_finite().then_some(statistic)
}

fn lrt_pvalue(distribution: &ChiSquared, statistic: f64) -> Option<f64> {
    let pvalue = distribution.sf(statistic.max(0.0));
    pvalue.is_finite().then_some(pvalue.clamp(0.0, 1.0))
}

fn checked_sub(left: f64, right: f64) -> Option<f64> {
    let value = left - right;
    value.is_finite().then_some(value)
}

fn checked_product2(left: f64, right: f64) -> Option<f64> {
    let value = left * right;
    (left.is_finite() && right.is_finite() && value.is_finite()).then_some(value)
}
