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
        let statistic = 2.0 * (full_log_like - reduced_log_like);
        if statistic.is_finite() {
            deviance.push(Some(statistic));
            pvalue.push(Some(1.0 - distribution.cdf(statistic.max(0.0))));
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
