//! DESeq2-parameterized negative-binomial likelihoods.

use crate::errors::{DeseqError, NegativeBinomialStage, invalid_dimensions};
use crate::matrix::{CountMatrix, RowMajorMatrix, try_alloc_matrix};
use crate::numeric::checked;
use crate::numeric::math::ProfileMath;
use crate::numeric::special::ln_gamma;

/// Negative-binomial log PMF using DESeq2's `mu` and dispersion parameterization.
///
/// This matches R's `dnbinom(x, mu = mu, size = 1 / dispersion, log = TRUE)`.
pub fn nbinom_log_pmf(count: u32, mu: f64, dispersion: f64) -> Result<f64, DeseqError> {
    validate_mu(mu, None, 0)?;
    validate_dispersion(dispersion, None, Some(0))?;
    let term = nbinom_log_pmf_unchecked(count, mu, dispersion);
    validate_log_pmf_term(term, None, 0)?;
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
    nbinom_log_likelihood_with_cached_gamma(counts, mu, dispersion, weights, None, None)
}

/// Invariant NB terms for repeated fitted means at fixed counts/dispersion.
/// All likelihood validation and the sample reduction retain the public path.
pub(crate) struct FixedDispersionLikelihood<'a> {
    counts: &'a [u32],
    dispersion: f64,
    gene: usize,
    size: f64,
    gamma_terms: Vec<f64>,
}

impl<'a> FixedDispersionLikelihood<'a> {
    pub(crate) fn new(counts: &'a [u32], dispersion: f64, gene: usize) -> Self {
        let size = dispersion.recip();
        let gamma_size = ln_gamma(size);
        let gamma_one = ln_gamma(1.0);
        let gamma_terms = counts
            .iter()
            .copied()
            .map(|count| {
                let y = f64::from(count);
                // Preserve both subtractions before the mean-dependent terms.
                let gamma_y_size = if count == 0 {
                    gamma_size
                } else {
                    ln_gamma(y + size)
                };
                let gamma_y_one = if count == 0 {
                    gamma_one
                } else {
                    ln_gamma(y + 1.0)
                };
                gamma_y_size - gamma_size - gamma_y_one
            })
            .collect();
        Self {
            counts,
            dispersion,
            gene,
            size,
            gamma_terms,
        }
    }

    pub(crate) fn log_likelihood(
        &self,
        mu: &[f64],
        weights: Option<&[f64]>,
    ) -> Result<f64, DeseqError> {
        nbinom_log_likelihood_with_cached_gamma(
            self.counts,
            mu,
            self.dispersion,
            weights,
            Some(self.gene),
            Some((&self.gamma_terms, self.size)),
        )
    }
}

fn nbinom_log_likelihood_with_cached_gamma(
    counts: &[u32],
    mu: &[f64],
    dispersion: f64,
    weights: Option<&[f64]>,
    gene: Option<usize>,
    cached_gamma: Option<(&[f64], f64)>,
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
    validate_dispersion(dispersion, gene, None)?;
    let size = dispersion.recip();
    let gamma_size = cached_gamma.is_none().then(|| ln_gamma(size));
    let gamma_one = cached_gamma.is_none().then(|| ln_gamma(1.0));
    let mut log_like = 0.0;
    for (idx, (count, mu)) in counts.iter().copied().zip(mu.iter().copied()).enumerate() {
        validate_mu(mu, gene, idx)?;
        let term = match cached_gamma {
            Some((gamma_terms, size)) => {
                let mu_dispersion = mu * dispersion;
                gamma_terms[idx] - size * log1p_mu_dispersion(mu, dispersion, mu_dispersion)
                    + count_log_term(f64::from(count), mu_dispersion)
            }
            None => {
                let y = f64::from(count);
                let mu_dispersion = mu * dispersion;
                let gamma_size = gamma_size.unwrap_or_else(|| ln_gamma(size));
                let gamma_y_size = if count == 0 {
                    gamma_size
                } else {
                    ln_gamma(y + size)
                };
                let gamma_y_one = if count == 0 {
                    gamma_one.unwrap_or_else(|| ln_gamma(1.0))
                } else {
                    ln_gamma(y + 1.0)
                };
                gamma_y_size
                    - gamma_size
                    - gamma_y_one
                    - size * log1p_mu_dispersion(mu, dispersion, mu_dispersion)
                    + count_log_term(y, mu_dispersion)
            }
        };
        validate_log_pmf_term(term, gene, idx)?;
        let weighted_term = match weights {
            Some(weights) => {
                let weight = weights[idx];
                validate_weight(weight, gene, idx)?;
                checked::try_mul(weight, term).map_err(|failure| {
                    invalid_value(
                        NegativeBinomialStage::WeightedTerm,
                        gene,
                        Some(idx),
                        failure.value,
                    )
                })?
            }
            None => term,
        };
        log_like = checked::try_add(log_like, weighted_term).map_err(|failure| {
            invalid_value(
                NegativeBinomialStage::RowSum,
                gene,
                Some(idx),
                failure.value,
            )
        })?;
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
        out.push(nbinom_log_likelihood_with_cached_gamma(
            count_row,
            mu_row,
            dispersion,
            weight_row,
            Some(gene),
            None,
        )?);
    }
    Ok(out)
}

/// Log likelihoods for selected source rows, indexed in compact fit order.
pub(crate) fn nbinom_log_likelihood_selected_rows(
    counts: &CountMatrix,
    mu: &RowMajorMatrix<f64>,
    dispersions: &[f64],
    weights: Option<&RowMajorMatrix<f64>>,
    genes: &[usize],
) -> Result<Vec<f64>, DeseqError> {
    let mut result = try_alloc_matrix(genes.len(), 1, "selected NB log likelihoods")?;
    for (gene, &source_gene) in genes.iter().enumerate() {
        let count_row = counts.row(source_gene)?;
        let weight_row = weights.map(|matrix| matrix.row(source_gene)).transpose()?;
        result.push(nbinom_log_likelihood_with_cached_gamma(
            count_row,
            mu.row(gene)?,
            dispersions[source_gene],
            weight_row,
            Some(gene),
            None,
        )?);
    }
    Ok(result)
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
        -y * mu_dispersion.recip().ln_1p_p()
    }
}

fn log1p_mu_dispersion(mu: f64, dispersion: f64, mu_dispersion: f64) -> f64 {
    if mu_dispersion.is_finite() {
        mu_dispersion.ln_1p_p()
    } else {
        mu.ln_p() + dispersion.ln_p()
    }
}

fn invalid_value(
    stage: NegativeBinomialStage,
    gene: Option<usize>,
    sample: Option<usize>,
    value: f64,
) -> DeseqError {
    DeseqError::InvalidNegativeBinomialValue {
        stage,
        gene,
        sample,
        value,
    }
}

fn validate_mu(mu: f64, gene: Option<usize>, sample: usize) -> Result<(), DeseqError> {
    if !mu.is_finite() || mu <= 0.0 {
        return Err(invalid_value(
            NegativeBinomialStage::Mean,
            gene,
            Some(sample),
            mu,
        ));
    }
    Ok(())
}

fn validate_dispersion(
    dispersion: f64,
    gene: Option<usize>,
    sample: Option<usize>,
) -> Result<(), DeseqError> {
    if !dispersion.is_finite() || dispersion <= 0.0 {
        return Err(invalid_value(
            NegativeBinomialStage::Dispersion,
            gene,
            sample,
            dispersion,
        ));
    }
    Ok(())
}

fn validate_weight(weight: f64, gene: Option<usize>, sample: usize) -> Result<(), DeseqError> {
    if !weight.is_finite() || weight < 0.0 {
        return Err(invalid_value(
            NegativeBinomialStage::Weight,
            gene,
            Some(sample),
            weight,
        ));
    }
    Ok(())
}

fn validate_log_pmf_term(term: f64, gene: Option<usize>, sample: usize) -> Result<(), DeseqError> {
    if !term.is_finite() {
        return Err(invalid_value(
            NegativeBinomialStage::LogPmf,
            gene,
            Some(sample),
            term,
        ));
    }
    Ok(())
}

fn negative_twice_log_likelihood_from_log_like(log_like: f64) -> Result<f64, DeseqError> {
    checked::try_mul(-2.0, log_like).map_err(|failure| {
        invalid_value(NegativeBinomialStage::Deviance, None, None, failure.value)
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn negative_twice_log_likelihood_rejects_nonfinite_scaling() {
        assert_eq!(
            negative_twice_log_likelihood_from_log_like(-2.0).unwrap(),
            4.0
        );
        std::assert_matches!(
            negative_twice_log_likelihood_from_log_like(f64::NAN).unwrap_err(),
            DeseqError::InvalidNegativeBinomialValue {
                stage: NegativeBinomialStage::Deviance, gene: None, sample: None, value,
            } if value.is_nan()
        );
        std::assert_matches!(
            negative_twice_log_likelihood_from_log_like(f64::MAX).unwrap_err(),
            DeseqError::InvalidNegativeBinomialValue {
                stage: NegativeBinomialStage::Deviance,
                gene: None,
                sample: None,
                value: f64::NEG_INFINITY,
            }
        );
    }
}
