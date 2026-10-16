// Intercept-only fits, normalization and working weights.

/// Fit DESeq2's intercept-only fixed-dispersion shortcut with size factors.
///
/// This mirrors the `justIntercept` branch in DESeq2 `fitNbinomGLMs` for the
/// common unweighted size-factor case.
pub fn fit_intercept_only_fixed_dispersion(
    counts: &CountMatrix,
    size_factors: &[f64],
    dispersions: &[f64],
) -> Result<NbinomGlmFit, DeseqError> {
    fit_intercept_only_fixed_dispersion_with_weights(counts, size_factors, dispersions, None)
}

/// Fit DESeq2's intercept-only fixed-dispersion shortcut with optional weights.
pub fn fit_intercept_only_fixed_dispersion_with_weights(
    counts: &CountMatrix,
    size_factors: &[f64],
    dispersions: &[f64],
    weights: Option<&RowMajorMatrix<f64>>,
) -> Result<NbinomGlmFit, DeseqError> {
    validate_fit_inputs(counts, size_factors, dispersions, weights)?;
    let normalization_factors = normalization_factors_from_size_factors(counts, size_factors)?;
    fit_intercept_only_fixed_dispersion_with_normalization_factors(
        counts,
        &normalization_factors,
        dispersions,
        weights,
    )
}

/// Fit DESeq2's intercept-only fixed-dispersion shortcut with normalization factors.
///
/// DESeq2's low-level code reconstructs `mu` as
/// `normalizationFactors * 2^betaMatrix`; this function exposes that matrix
/// form directly for future support of gene/sample normalization factors.
pub fn fit_intercept_only_fixed_dispersion_with_normalization_factors(
    counts: &CountMatrix,
    normalization_factors: &RowMajorMatrix<f64>,
    dispersions: &[f64],
    weights: Option<&RowMajorMatrix<f64>>,
) -> Result<NbinomGlmFit, DeseqError> {
    validate_nf_fit_inputs(counts, normalization_factors, dispersions, weights)?;

    let normalized = normalized_counts_with_factors(counts, normalization_factors)?;
    let mut beta_values = Vec::with_capacity(counts.n_genes());
    let mut mu_values = Vec::with_capacity(counts.n_genes() * counts.n_samples());
    let mut hat_values = Vec::with_capacity(counts.n_genes() * counts.n_samples());
    let mut beta_se_values = Vec::with_capacity(counts.n_genes());
    let mut beta_covariance_values = Vec::with_capacity(counts.n_genes());

    for (gene, dispersion) in dispersions.iter().copied().enumerate() {
        let normalized_row = normalized.row(gene)?;
        let weight_row = weights.map(|matrix| matrix.row(gene)).transpose()?;
        let mean_normalized = match weight_row {
            Some(weights) => weighted_mean(normalized_row, weights, gene)?,
            None => {
                if counts.is_all_zero_gene(gene)? {
                    return Err(DeseqError::InvalidCounts {
                        reason: format!(
                            "gene {gene} is all zero; DESeq2 GLM fitting excludes allZero rows"
                        ),
                    });
                }
                checked_mean(normalized_row).ok_or_else(|| DeseqError::InvalidCounts {
                    reason: format!("gene {gene} has non-finite normalized intercept mean"),
                })?
            }
        };
        if !mean_normalized.is_finite() || mean_normalized <= 0.0 {
            return Err(DeseqError::InvalidCounts {
                reason: format!("gene {gene} has non-positive normalized intercept mean"),
            });
        }

        let beta_log2 = mean_normalized.log2();
        beta_values.push(beta_log2);

        for factor in normalization_factors.row(gene)? {
            let mu = *factor * 2.0_f64.powf(beta_log2);
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
        let xtwx = checked_sum(working_weights.iter().copied()).ok_or_else(|| {
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
            let Some(hat) = checked_product2(value, sigma) else {
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

fn normalization_factors_from_size_factors(
    counts: &CountMatrix,
    size_factors: &[f64],
) -> Result<RowMajorMatrix<f64>, DeseqError> {
    if size_factors.len() != counts.n_samples() {
        return Err(invalid_dimensions(
            "size factors",
            counts.n_samples(),
            size_factors.len(),
        ));
    }
    let mut values = Vec::with_capacity(counts.n_genes() * counts.n_samples());
    for _ in 0..counts.n_genes() {
        for (idx, factor) in size_factors.iter().copied().enumerate() {
            validate_positive_finite(factor, "size factor", idx)?;
            values.push(factor);
        }
    }
    RowMajorMatrix::from_row_major(counts.n_genes(), counts.n_samples(), values)
}

fn normalized_counts_with_factors(
    counts: &CountMatrix,
    normalization_factors: &RowMajorMatrix<f64>,
) -> Result<RowMajorMatrix<f64>, DeseqError> {
    let mut values = Vec::with_capacity(counts.n_genes() * counts.n_samples());
    for gene in 0..counts.n_genes() {
        let count_row = counts.row(gene)?;
        let factor_row = normalization_factors.row(gene)?;
        for (sample, (count, factor)) in count_row
            .iter()
            .copied()
            .zip(factor_row.iter().copied())
            .enumerate()
        {
            validate_positive_finite(factor, "normalization factor", sample)?;
            values.push(f64::from(count) / factor);
        }
    }
    RowMajorMatrix::from_row_major(counts.n_genes(), counts.n_samples(), values)
}

fn weighted_mean(values: &[f64], weights: &[f64], gene: usize) -> Result<f64, DeseqError> {
    let mut weighted_values = Vec::with_capacity(values.len());
    for (sample, (value, weight)) in values
        .iter()
        .copied()
        .zip(weights.iter().copied())
        .enumerate()
    {
        validate_nonnegative_finite(weight, "weight", sample)?;
        let Some(value_term) = checked_product2(weight, value) else {
            return Err(DeseqError::InvalidCounts {
                reason: format!("gene {gene} has non-finite weighted normalized mean"),
            });
        };
        weighted_values.push(value_term);
    }
    let Some(numerator_mean) = checked_scaled_mean(&weighted_values) else {
        return Err(DeseqError::InvalidCounts {
            reason: format!("gene {gene} has non-finite weighted normalized mean"),
        });
    };
    let Some(denominator) = checked_scaled_sum(weights.iter().copied()) else {
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
    let mean_scale = checked_div2(weighted_values.len() as f64, denominator).ok_or_else(|| {
        DeseqError::InvalidCounts {
            reason: format!("gene {gene} has non-finite weighted normalized mean"),
        }
    })?;
    checked_product2(numerator_mean, mean_scale).ok_or_else(|| DeseqError::InvalidCounts {
        reason: format!("gene {gene} has non-finite weighted normalized mean"),
    })
}

fn intercept_working_weights(
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
                checked_product2(weight, working_weight).ok_or_else(|| {
                    DeseqError::NonFiniteValue {
                        context: "intercept working weight".to_string(),
                        index: Some(sample),
                        value: f64::NAN,
                    }
                })?
            }
            None => working_weight,
        });
    }
    Ok(out)
}
