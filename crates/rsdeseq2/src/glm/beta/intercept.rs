use super::*;
use crate::errors::ErrorKind;
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
    fit_intercept_only_glm_rows(input, None)
}

pub(crate) fn fit_intercept_only_glm_rows(
    input: GlmInput<'_>,
    genes: Option<&[usize]>,
) -> Result<NbinomGlmFit, DeseqError> {
    validate_intercept_fit_inputs(input, genes)?;
    let GlmInput {
        counts,
        offsets,
        dispersions,
        weights,
    } = input;

    let n_genes = genes.map_or(counts.n_genes(), <[usize]>::len);
    let n_samples = counts.n_samples();
    let mut normalized_row = try_alloc_matrix(n_samples, 1, "intercept normalized row")?;
    let mut working_weights = try_alloc_matrix(n_samples, 1, "intercept working weights")?;
    let mut weighted_values = Vec::new();
    let mut beta_values = try_alloc_matrix(n_genes, 1, "intercept coefficients")?;
    let mut mu_values = try_alloc_matrix(n_genes, counts.n_samples(), "intercept-only fit")?;
    let mut hat_values = try_alloc_matrix(n_genes, counts.n_samples(), "intercept-only fit")?;
    let mut beta_se_values = try_alloc_matrix(n_genes, 1, "intercept standard errors")?;
    let mut beta_covariance_values = try_alloc_matrix(n_genes, 1, "intercept covariance")?;

    for gene in 0..n_genes {
        let source_gene = genes.map_or(gene, |rows| rows[gene]);
        let dispersion = dispersions[source_gene];
        let factor_row = offsets.row(source_gene)?;
        normalized_row.clear();
        normalized_row.extend(
            counts
                .row_values(source_gene)
                .iter()
                .zip(factor_row)
                .map(|(&count, &factor)| f64::from(count) / factor),
        );
        let weight_row = weights.map(|matrix| matrix.row(source_gene)).transpose()?;
        let mean_normalized = match weight_row {
            Some(weights) => weighted_mean(&normalized_row, weights, gene, &mut weighted_values)?,
            None => {
                if counts.is_all_zero_gene(source_gene)? {
                    return Err(DeseqError::AllZeroGene { gene, source_gene });
                }
                checked::try_mean(&normalized_row).map_err(|failure| {
                    intercept_value_error(
                        ErrorKind::Counts,
                        gene,
                        failure.index,
                        "normalized intercept mean",
                        "non-finite",
                        failure.value,
                    )
                })?
            }
        };
        if !mean_normalized.is_finite() || mean_normalized <= 0.0 {
            return Err(intercept_value_error(
                ErrorKind::Counts,
                gene,
                None,
                "normalized intercept mean",
                "non-positive",
                mean_normalized,
            ));
        }

        let beta_log2 = mean_normalized.log2_p();
        beta_values.push(beta_log2);

        for (sample, factor) in factor_row.iter().enumerate() {
            let mu = *factor * 2.0_f64.powf_p(beta_log2);
            if !mu.is_finite() || mu <= 0.0 {
                return Err(intercept_value_error(
                    ErrorKind::Counts,
                    gene,
                    Some(sample),
                    "fitted intercept mean",
                    "non-finite",
                    mu,
                ));
            }
            mu_values.push(mu);
        }

        let mu_start = gene * counts.n_samples();
        let mu_row = &mu_values[mu_start..mu_start + counts.n_samples()];
        fill_intercept_working_weights(mu_row, dispersion, weight_row, &mut working_weights)?;
        let xtwx = checked::try_sum(working_weights.iter().copied()).map_err(|failure| {
            intercept_value_error(
                ErrorKind::Dispersion,
                gene,
                failure.index,
                "intercept information",
                "non-finite",
                failure.value,
            )
        })?;
        if !xtwx.is_finite() || xtwx <= 0.0 {
            return Err(intercept_value_error(
                ErrorKind::Dispersion,
                gene,
                None,
                "intercept information",
                "non-positive",
                xtwx,
            ));
        }
        let sigma = xtwx.recip();
        if !sigma.is_finite() {
            return Err(intercept_value_error(
                ErrorKind::Dispersion,
                gene,
                None,
                "intercept covariance",
                "non-finite",
                sigma,
            ));
        }
        beta_se_values.push(checked_intercept_beta_se(sigma, gene)?);
        beta_covariance_values.push(checked_intercept_beta_covariance(sigma, gene)?);
        for (sample, value) in working_weights.iter().copied().enumerate() {
            let hat =
                checked::try_mul(value, sigma).map_err(|failure| DeseqError::InvalidGeneValue {
                    domain: ErrorKind::Dispersion,
                    operation: "intercept hat diagonal",
                    gene,
                    sample: Some(sample),
                    value: failure.value,
                    reason: format!(
                        "gene {gene} sample {sample} has non-finite intercept hat diagonal"
                    ),
                })?;
            hat_values.push(hat);
        }
    }

    let beta = RowMajorMatrix::from_row_major(n_genes, 1, beta_values)?;
    let beta_se = RowMajorMatrix::from_row_major(n_genes, 1, beta_se_values)?;
    let beta_covariance = RowMajorMatrix::from_row_major(n_genes, 1, beta_covariance_values)?;
    let mu = RowMajorMatrix::from_row_major(n_genes, counts.n_samples(), mu_values)?;
    let hat_diagonal = RowMajorMatrix::from_row_major(n_genes, counts.n_samples(), hat_values)?;
    let log_like = match genes {
        Some(genes) => {
            nbinom_log_likelihood_selected_rows(counts, &mu, dispersions, weights, genes)?
        }
        None => nbinom_log_likelihood_matrix(counts, &mu, dispersions, weights)?,
    };
    let model_matrix = DesignMatrix::from_row_major(
        counts.n_samples(),
        1,
        vec![1.0; counts.n_samples()],
        Some(vec!["Intercept".to_string()]),
    )?;

    Ok(NbinomGlmFit {
        log_like,
        beta_converged: vec![true; n_genes],
        beta,
        beta_se,
        beta_optim_start: RowMajorMatrix::from_elem(n_genes, 1, f64::NAN)?,
        beta_covariance: Some(beta_covariance),
        mu,
        beta_iter: vec![1; n_genes],
        beta_optim_iter: vec![f64::NAN; n_genes],
        beta_optim_start_objective: vec![f64::NAN; n_genes],
        beta_optim_objective: vec![f64::NAN; n_genes],
        beta_optim_gradient_norm: vec![f64::NAN; n_genes],
        model_matrix,
        n_terms: 1,
        hat_diagonal,
    })
}

pub(super) fn intercept_value_error(
    domain: ErrorKind,
    gene: usize,
    sample: Option<usize>,
    operation: &'static str,
    issue: &str,
    value: f64,
) -> DeseqError {
    DeseqError::InvalidGeneValue {
        domain,
        operation,
        gene,
        sample,
        value,
        reason: format!("gene {gene} has {issue} {operation}"),
    }
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
    let invalid_mean = |sample, value| {
        intercept_value_error(
            ErrorKind::Counts,
            gene,
            sample,
            "weighted normalized mean",
            "non-finite",
            value,
        )
    };
    for (sample, (value, weight)) in values
        .iter()
        .copied()
        .zip(weights.iter().copied())
        .enumerate()
    {
        validate_nonnegative_finite(weight, "weight", sample)?;
        let value_term = checked::try_mul(weight, value)
            .map_err(|failure| invalid_mean(Some(sample), failure.value))?;
        weighted_values.push(value_term);
    }
    let numerator_mean = checked::try_mean(weighted_values)
        .map_err(|failure| invalid_mean(failure.index, failure.value))?;
    let denominator = checked::try_scaled_sum(weights).map_err(|failure| {
        intercept_value_error(
            ErrorKind::Counts,
            gene,
            failure.index,
            "total weight",
            "non-finite",
            failure.value,
        )
    })?;
    if denominator <= 0.0 {
        return Err(intercept_value_error(
            ErrorKind::Counts,
            gene,
            None,
            "total weight",
            "zero",
            denominator,
        ));
    }
    let mean_scale = checked::try_div(weighted_values.len() as f64, denominator)
        .map_err(|failure| invalid_mean(None, failure.value))?;
    checked::try_mul(numerator_mean, mean_scale)
        .map_err(|failure| invalid_mean(None, failure.value))
}

#[cfg(test)]
pub(super) fn intercept_working_weights(
    mu: &[f64],
    dispersion: f64,
    weights: Option<&[f64]>,
) -> Result<Vec<f64>, DeseqError> {
    let mut out = Vec::new();
    fill_intercept_working_weights(mu, dispersion, weights, &mut out)?;
    Ok(out)
}

fn fill_intercept_working_weights(
    mu: &[f64],
    dispersion: f64,
    weights: Option<&[f64]>,
    out: &mut Vec<f64>,
) -> Result<(), DeseqError> {
    validate_positive_finite(dispersion, "dispersion", 0)?;
    out.clear();
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
    Ok(())
}
