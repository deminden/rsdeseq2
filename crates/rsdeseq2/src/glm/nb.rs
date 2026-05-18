use statrs::function::gamma::ln_gamma;

use crate::core::CountMatrix;
use crate::errors::{invalid_dimensions, DeseqError};
use crate::matrix::RowMajorMatrix;

/// Negative-binomial log PMF using DESeq2's `mu` and dispersion parameterization.
///
/// This matches R's `dnbinom(x, mu = mu, size = 1 / dispersion, log = TRUE)`.
pub fn nbinom_log_pmf(count: u32, mu: f64, dispersion: f64) -> Result<f64, DeseqError> {
    validate_mu(mu, Some(0))?;
    validate_dispersion(dispersion, Some(0))?;
    let y = f64::from(count);
    let size = dispersion.recip();
    Ok(ln_gamma(y + size) - ln_gamma(size) - ln_gamma(y + 1.0)
        + size * (size / (size + mu)).ln()
        + y * (mu / (size + mu)).ln())
}

/// Row log likelihood for one gene.
pub fn nbinom_log_likelihood(
    counts: &[u32],
    mu: &[f64],
    dispersion: f64,
) -> Result<f64, DeseqError> {
    nbinom_log_likelihood_weighted(counts, mu, dispersion, None)
}

/// Weighted row log likelihood for one gene.
pub fn nbinom_log_likelihood_weighted(
    counts: &[u32],
    mu: &[f64],
    dispersion: f64,
    weights: Option<&[f64]>,
) -> Result<f64, DeseqError> {
    if counts.len() != mu.len() {
        return Err(invalid_dimensions(
            "NB likelihood mu",
            counts.len(),
            mu.len(),
        ));
    }
    if let Some(weights) = weights {
        if weights.len() != counts.len() {
            return Err(invalid_dimensions(
                "NB likelihood weights",
                counts.len(),
                weights.len(),
            ));
        }
    }
    validate_dispersion(dispersion, None)?;
    let mut log_like = 0.0;
    for (idx, (count, mu)) in counts.iter().copied().zip(mu.iter().copied()).enumerate() {
        validate_mu(mu, Some(idx))?;
        let term = nbinom_log_pmf_unchecked(count, mu, dispersion);
        log_like += match weights {
            Some(weights) => {
                let weight = weights[idx];
                validate_weight(weight, Some(idx))?;
                weight * term
            }
            None => term,
        };
    }
    Ok(log_like)
}

/// DESeq2-style `-2 * logLik` for one gene.
pub fn nbinom_negative_twice_log_likelihood(
    counts: &[u32],
    mu: &[f64],
    dispersion: f64,
) -> Result<f64, DeseqError> {
    Ok(-2.0 * nbinom_log_likelihood(counts, mu, dispersion)?)
}

/// Row-wise log likelihoods for a count matrix.
///
/// This mirrors DESeq2's `nbinomLogLike`: each row is summed over samples using
/// `dnbinom(..., mu=mu, size=1/disp, log=TRUE)`.
pub fn nbinom_log_likelihood_matrix(
    counts: &CountMatrix,
    mu: &RowMajorMatrix<f64>,
    dispersions: &[f64],
    weights: Option<&RowMajorMatrix<f64>>,
) -> Result<Vec<f64>, DeseqError> {
    if mu.n_rows() != counts.n_genes() || mu.n_cols() != counts.n_samples() {
        return Err(DeseqError::InvalidDimensions {
            context: "NB likelihood mu matrix".to_string(),
            expected: counts.n_genes() * counts.n_samples(),
            actual: mu.len(),
        });
    }
    if dispersions.len() != counts.n_genes() {
        return Err(invalid_dimensions(
            "NB likelihood dispersions",
            counts.n_genes(),
            dispersions.len(),
        ));
    }
    if let Some(weights) = weights {
        if weights.n_rows() != counts.n_genes() || weights.n_cols() != counts.n_samples() {
            return Err(DeseqError::InvalidDimensions {
                context: "NB likelihood weights matrix".to_string(),
                expected: counts.n_genes() * counts.n_samples(),
                actual: weights.len(),
            });
        }
    }

    let mut out = Vec::with_capacity(counts.n_genes());
    for (gene, dispersion) in dispersions.iter().copied().enumerate() {
        let count_row = counts.row(gene)?;
        let mu_row = mu.row(gene)?;
        let weight_row = weights.map(|matrix| matrix.row(gene)).transpose()?;
        out.push(nbinom_log_likelihood_weighted(
            count_row, mu_row, dispersion, weight_row,
        )?);
    }
    Ok(out)
}

fn nbinom_log_pmf_unchecked(count: u32, mu: f64, dispersion: f64) -> f64 {
    let y = f64::from(count);
    let size = dispersion.recip();
    ln_gamma(y + size) - ln_gamma(size) - ln_gamma(y + 1.0)
        + size * (size / (size + mu)).ln()
        + y * (mu / (size + mu)).ln()
}

fn validate_mu(mu: f64, index: Option<usize>) -> Result<(), DeseqError> {
    if !mu.is_finite() || mu <= 0.0 {
        return Err(DeseqError::NonFiniteValue {
            context: "negative-binomial mean".to_string(),
            index,
            value: mu,
        });
    }
    Ok(())
}

fn validate_dispersion(dispersion: f64, index: Option<usize>) -> Result<(), DeseqError> {
    if !dispersion.is_finite() || dispersion <= 0.0 {
        return Err(DeseqError::InvalidDispersion {
            reason: match index {
                Some(index) => {
                    format!("dispersion at index {index} must be finite and positive")
                }
                None => "dispersion must be finite and positive".to_string(),
            },
        });
    }
    Ok(())
}

fn validate_weight(weight: f64, index: Option<usize>) -> Result<(), DeseqError> {
    if !weight.is_finite() || weight < 0.0 {
        return Err(DeseqError::NonFiniteValue {
            context: "negative-binomial weight".to_string(),
            index,
            value: weight,
        });
    }
    Ok(())
}
