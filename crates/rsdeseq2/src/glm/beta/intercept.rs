use super::*;
use crate::matrix::try_alloc_matrix;
use crate::numeric::math::ProfileMath;
// Intercept-only fits, normalization and working weights.

/// Fit DESeq2's intercept-only fixed-dispersion shortcut with size factors.
///
/// This mirrors the `justIntercept` branch in DESeq2 `fitNbinomGLMs` for the
/// common unweighted size-factor case.
#[deprecated(
    since = "0.6.0",
    note = "use fit_intercept_only_glm(GlmInput::new(counts, Offsets::PerSample(size_factors), dispersions))"
)]
pub fn fit_intercept_only_fixed_dispersion(
    counts: &CountMatrix,
    size_factors: &[f64],
    dispersions: &[f64],
) -> Result<NbinomGlmFit, DeseqError> {
    fit_intercept_only_glm(GlmInput::new(
        counts,
        Offsets::PerSample(size_factors),
        dispersions,
    ))
}

/// Fit DESeq2's intercept-only fixed-dispersion shortcut with optional weights.
#[deprecated(
    since = "0.6.0",
    note = "use fit_intercept_only_glm with Offsets::PerSample and GlmInput::weights"
)]
pub fn fit_intercept_only_fixed_dispersion_with_weights(
    counts: &CountMatrix,
    size_factors: &[f64],
    dispersions: &[f64],
    weights: Option<&RowMajorMatrix<f64>>,
) -> Result<NbinomGlmFit, DeseqError> {
    fit_intercept_only_glm(GlmInput {
        counts,
        offsets: Offsets::PerSample(size_factors),
        dispersions,
        weights,
    })
}

/// Fit the intercept-only fixed-dispersion GLM in closed form (DESeq2's
/// `fitNbinomGLMs` intercept shortcut) with size-factor or
/// normalization-factor offsets and optional observation weights.
pub fn fit_intercept_only_glm(input: GlmInput<'_>) -> Result<NbinomGlmFit, DeseqError> {
    validate_intercept_fit_inputs(input)?;
    let GlmInput {
        counts,
        offsets,
        dispersions,
        weights,
    } = input;

    let mut normalized_row = Vec::with_capacity(counts.n_samples());
    let mut weighted_values = Vec::new();
    let mut beta_values = Vec::with_capacity(counts.n_genes());
    let mut mu_values =
        try_alloc_matrix(counts.n_genes(), counts.n_samples(), "intercept-only fit")?;
    let mut hat_values =
        try_alloc_matrix(counts.n_genes(), counts.n_samples(), "intercept-only fit")?;
    let mut beta_se_values = Vec::with_capacity(counts.n_genes());
    let mut beta_covariance_values = Vec::with_capacity(counts.n_genes());

    for (gene, dispersion) in dispersions.iter().copied().enumerate() {
        let factor_row = offsets.row(gene)?;
        normalized_row.clear();
        normalized_row.extend(
            counts
                .row_values(gene)
                .iter()
                .zip(factor_row)
                .map(|(&count, &factor)| f64::from(count) / factor),
        );
        let weight_row = weights.map(|matrix| matrix.row(gene)).transpose()?;
        let mean_normalized = match weight_row {
            Some(weights) => weighted_mean(&normalized_row, weights, gene, &mut weighted_values)?,
            None => {
                if counts.is_all_zero_gene(gene)? {
                    return Err(DeseqError::InvalidCounts {
                        reason: format!(
                            "gene {gene} is all zero; DESeq2 GLM fitting excludes allZero rows"
                        ),
                    });
                }
                checked::mean(&normalized_row).ok_or_else(|| DeseqError::InvalidCounts {
                    reason: format!("gene {gene} has non-finite normalized intercept mean"),
                })?
            }
        };
        if !mean_normalized.is_finite() || mean_normalized <= 0.0 {
            return Err(DeseqError::InvalidCounts {
                reason: format!("gene {gene} has non-positive normalized intercept mean"),
            });
        }

        let beta_log2 = mean_normalized.log2_p();
        beta_values.push(beta_log2);

        for factor in factor_row {
            let mu = *factor * 2.0_f64.powf_p(beta_log2);
            if !mu.is_finite() || mu <= 0.0 {
                return Err(DeseqError::InvalidCounts {
                    reason: format!("gene {gene} has non-finite fitted intercept mean"),
                });
            }
            mu_values.push(mu);
        }

        let mu_start = gene * counts.n_samples();
        let mu_row = &mu_values[mu_start..mu_start + counts.n_samples()];
        let working_weights = intercept_working_weights(mu_row, dispersion, weight_row)?;
        let xtwx = checked::sum(working_weights.iter().copied()).ok_or_else(|| {
            DeseqError::InvalidDispersion {
                reason: format!("gene {gene} has non-finite intercept information"),
            }
        })?;
        if !xtwx.is_finite() || xtwx <= 0.0 {
            return Err(DeseqError::InvalidDispersion {
                reason: format!("gene {gene} has non-positive intercept information"),
            });
        }
        let sigma = xtwx.recip();
        if !sigma.is_finite() {
            return Err(DeseqError::InvalidDispersion {
                reason: format!("gene {gene} has non-finite intercept covariance"),
            });
        }
        beta_se_values.push(checked_intercept_beta_se(sigma, gene)?);
        beta_covariance_values.push(checked_intercept_beta_covariance(sigma, gene)?);
        for (sample, value) in working_weights.into_iter().enumerate() {
            let Some(hat) = checked::mul(value, sigma) else {
                return Err(DeseqError::InvalidDispersion {
                    reason: format!(
                        "gene {gene} sample {sample} has non-finite intercept hat diagonal"
                    ),
                });
            };
            hat_values.push(hat);
        }
    }

    let beta = RowMajorMatrix::from_row_major(counts.n_genes(), 1, beta_values)?;
    let beta_se = RowMajorMatrix::from_row_major(counts.n_genes(), 1, beta_se_values)?;
    let beta_covariance =
        RowMajorMatrix::from_row_major(counts.n_genes(), 1, beta_covariance_values)?;
    let mu = RowMajorMatrix::from_row_major(counts.n_genes(), counts.n_samples(), mu_values)?;
    let hat_diagonal =
        RowMajorMatrix::from_row_major(counts.n_genes(), counts.n_samples(), hat_values)?;
    let log_like = nbinom_log_likelihood_matrix(counts, &mu, dispersions, weights)?;
    let model_matrix = DesignMatrix::from_row_major(
        counts.n_samples(),
        1,
        vec![1.0; counts.n_samples()],
        Some(vec!["Intercept".to_string()]),
    )?;

    Ok(NbinomGlmFit {
        log_like,
        beta_converged: vec![true; counts.n_genes()],
        beta,
        beta_se,
        beta_optim_start: RowMajorMatrix::from_elem(counts.n_genes(), 1, f64::NAN)?,
        beta_covariance: Some(beta_covariance),
        mu,
        beta_iter: vec![1; counts.n_genes()],
        beta_optim_iter: vec![f64::NAN; counts.n_genes()],
        beta_optim_start_objective: vec![f64::NAN; counts.n_genes()],
        beta_optim_objective: vec![f64::NAN; counts.n_genes()],
        beta_optim_gradient_norm: vec![f64::NAN; counts.n_genes()],
        model_matrix,
        n_terms: 1,
        hat_diagonal,
    })
}

/// Fit DESeq2's intercept-only fixed-dispersion shortcut with normalization factors.
///
/// DESeq2's low-level code reconstructs `mu` as
/// `normalizationFactors * 2^betaMatrix`; this function exposes that matrix
/// form directly for future support of gene/sample normalization factors.
#[deprecated(
    since = "0.6.0",
    note = "use fit_intercept_only_glm with Offsets::PerCell and GlmInput::weights"
)]
pub fn fit_intercept_only_fixed_dispersion_with_normalization_factors(
    counts: &CountMatrix,
    normalization_factors: &RowMajorMatrix<f64>,
    dispersions: &[f64],
    weights: Option<&RowMajorMatrix<f64>>,
) -> Result<NbinomGlmFit, DeseqError> {
    fit_intercept_only_glm(GlmInput {
        counts,
        offsets: Offsets::PerCell(normalization_factors),
        dispersions,
        weights,
    })
}

pub(super) fn weighted_mean(
    values: &[f64],
    weights: &[f64],
    gene: usize,
    weighted_values: &mut Vec<f64>,
) -> Result<f64, DeseqError> {
    weighted_values.clear();
    for (sample, (value, weight)) in values
        .iter()
        .copied()
        .zip(weights.iter().copied())
        .enumerate()
    {
        validate_nonnegative_finite(weight, "weight", sample)?;
        let Some(value_term) = checked::mul(weight, value) else {
            return Err(DeseqError::InvalidCounts {
                reason: format!("gene {gene} has non-finite weighted normalized mean"),
            });
        };
        weighted_values.push(value_term);
    }
    let Some(numerator_mean) = checked::mean(weighted_values) else {
        return Err(DeseqError::InvalidCounts {
            reason: format!("gene {gene} has non-finite weighted normalized mean"),
        });
    };
    let Some(denominator) = checked::scaled_sum(weights) else {
        return Err(DeseqError::InvalidCounts {
            reason: format!("gene {gene} has non-finite total weight"),
        });
    };
    if !denominator.is_finite() {
        return Err(DeseqError::InvalidCounts {
            reason: format!("gene {gene} has non-finite total weight"),
        });
    }
    if denominator <= 0.0 {
        return Err(DeseqError::InvalidCounts {
            reason: format!("gene {gene} has zero total weight"),
        });
    }
    let mean_scale = checked::div(weighted_values.len() as f64, denominator).ok_or_else(|| {
        DeseqError::InvalidCounts {
            reason: format!("gene {gene} has non-finite weighted normalized mean"),
        }
    })?;
    checked::mul(numerator_mean, mean_scale).ok_or_else(|| DeseqError::InvalidCounts {
        reason: format!("gene {gene} has non-finite weighted normalized mean"),
    })
}

pub(super) fn intercept_working_weights(
    mu: &[f64],
    dispersion: f64,
    weights: Option<&[f64]>,
) -> Result<Vec<f64>, DeseqError> {
    validate_positive_finite(dispersion, "dispersion", 0)?;
    let mut out = Vec::with_capacity(mu.len());
    for (sample, value) in mu.iter().copied().enumerate() {
        validate_positive_finite(value, "mu", sample)?;
        let working_weight = checked_mean_dispersion_working_weight(
            value,
            dispersion,
            sample,
            "intercept working weight",
        )?;
        out.push(match weights {
            Some(weights) => {
                let weight = weights[sample];
                validate_nonnegative_finite(weight, "weight", sample)?;
                checked::try_mul(weight, working_weight)
                    .nonfinite_at("intercept working weight", Some(sample))?
            }
            None => working_weight,
        });
    }
    Ok(out)
}
