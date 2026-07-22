use statrs::function::gamma::ln_gamma;

use crate::core::CountMatrix;
use crate::errors::{DeseqError, invalid_dimensions};
use crate::matrix::RowMajorMatrix;

/// Negative-binomial log PMF using DESeq2's `mu` and dispersion parameterization.
///
/// This matches R's `dnbinom(x, mu = mu, size = 1 / dispersion, log = TRUE)`.
pub fn nbinom_log_pmf(count: u32, mu: f64, dispersion: f64) -> Result<f64, DeseqError> {
    validate_mu(mu, Some(0))?;
    validate_dispersion(dispersion, Some(0))?;
    let term = nbinom_log_pmf_unchecked(count, mu, dispersion);
    validate_log_pmf_term(term, Some(0))?;
    Ok(term)
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
    if let Some(weights) = weights
        && weights.len() != counts.len()
    {
        return Err(invalid_dimensions(
            "NB likelihood weights",
            counts.len(),
            weights.len(),
        ));
    }
    validate_dispersion(dispersion, None)?;
    let mut log_like = 0.0;
    for (idx, (count, mu)) in counts.iter().copied().zip(mu.iter().copied()).enumerate() {
        validate_mu(mu, Some(idx))?;
        let term = nbinom_log_pmf_unchecked(count, mu, dispersion);
        validate_log_pmf_term(term, Some(idx))?;
        let weighted_term = match weights {
            Some(weights) => {
                let weight = weights[idx];
                validate_weight(weight, Some(idx))?;
                checked_weighted_log_likelihood_term(weight, term).ok_or_else(|| {
                    DeseqError::InvalidDispersion {
                        reason: format!(
                            "negative-binomial weighted log-likelihood term at sample {idx} must be finite"
                        ),
                    }
                })?
            }
            None => term,
        };
        validate_weighted_log_likelihood_term(weighted_term, Some(idx))?;
        checked_add_log_likelihood(&mut log_like, weighted_term)?;
    }
    Ok(log_like)
}

/// DESeq2-style `-2 * logLik` for one gene.
pub fn nbinom_negative_twice_log_likelihood(
    counts: &[u32],
    mu: &[f64],
    dispersion: f64,
) -> Result<f64, DeseqError> {
    negative_twice_log_likelihood_from_log_like(nbinom_log_likelihood(counts, mu, dispersion)?)
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
    if let Some(weights) = weights
        && (weights.n_rows() != counts.n_genes() || weights.n_cols() != counts.n_samples())
    {
        return Err(DeseqError::InvalidDimensions {
            context: "NB likelihood weights matrix".to_string(),
            expected: counts.n_genes() * counts.n_samples(),
            actual: weights.len(),
        });
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
    let mu_dispersion = mu * dispersion;
    ln_gamma(y + size)
        - ln_gamma(size)
        - ln_gamma(y + 1.0)
        - size * log1p_mu_dispersion(mu, dispersion, mu_dispersion)
        + count_log_term(y, mu_dispersion)
}

fn count_log_term(y: f64, mu_dispersion: f64) -> f64 {
    if y == 0.0 || (mu_dispersion.is_infinite() && mu_dispersion.is_sign_positive()) {
        0.0
    } else {
        -y * mu_dispersion.recip().ln_1p()
    }
}

fn log1p_mu_dispersion(mu: f64, dispersion: f64, mu_dispersion: f64) -> f64 {
    if mu_dispersion.is_finite() {
        mu_dispersion.ln_1p()
    } else {
        mu.ln() + dispersion.ln()
    }
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

fn validate_log_pmf_term(term: f64, index: Option<usize>) -> Result<(), DeseqError> {
    if !term.is_finite() {
        return Err(DeseqError::InvalidDispersion {
            reason: match index {
                Some(index) => {
                    format!("negative-binomial log PMF at sample {index} must be finite")
                }
                None => "negative-binomial log PMF must be finite".to_string(),
            },
        });
    }
    Ok(())
}

fn validate_weighted_log_likelihood_term(
    term: f64,
    index: Option<usize>,
) -> Result<(), DeseqError> {
    if !term.is_finite() {
        return Err(DeseqError::InvalidDispersion {
            reason: match index {
                Some(index) => {
                    format!(
                        "negative-binomial weighted log-likelihood term at sample {index} must be finite"
                    )
                }
                None => "negative-binomial weighted log-likelihood term must be finite".to_string(),
            },
        });
    }
    Ok(())
}

fn checked_weighted_log_likelihood_term(weight: f64, term: f64) -> Option<f64> {
    let product = weight * term;
    (weight.is_finite() && term.is_finite() && product.is_finite()).then_some(product)
}

fn checked_add_log_likelihood(total: &mut f64, term: f64) -> Result<(), DeseqError> {
    let next = *total + term;
    if !next.is_finite() {
        return Err(DeseqError::InvalidDispersion {
            reason: "negative-binomial row log-likelihood sum must remain finite".to_string(),
        });
    }
    *total = next;
    Ok(())
}

fn negative_twice_log_likelihood_from_log_like(log_like: f64) -> Result<f64, DeseqError> {
    let deviance = -2.0 * log_like;
    if log_like.is_finite() && deviance.is_finite() {
        Ok(deviance)
    } else {
        Err(DeseqError::InvalidDispersion {
            reason: "negative-binomial deviance must remain finite".to_string(),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::negative_twice_log_likelihood_from_log_like;

    #[test]
    fn negative_twice_log_likelihood_rejects_nonfinite_scaling() {
        assert_eq!(
            negative_twice_log_likelihood_from_log_like(-2.0).unwrap(),
            4.0
        );
        assert!(negative_twice_log_likelihood_from_log_like(f64::NAN).is_err());
        assert!(negative_twice_log_likelihood_from_log_like(f64::MAX).is_err());
    }
}
