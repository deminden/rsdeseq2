use crate::core::CountMatrix;
use crate::design::DesignMatrix;
use crate::errors::{invalid_dimensions, DeseqError};
use crate::glm::irls::{
    fit_fixed_dispersion_irls, fit_fixed_dispersion_irls_with_normalization_factors_and_weights,
    fit_irls, IrlsOptions,
};
use crate::glm::nb::nbinom_log_likelihood_matrix;
use crate::glm::NbinomGlmFit;
use crate::matrix::RowMajorMatrix;
use statrs::distribution::{ContinuousCDF, Normal};

/// Method used to estimate DESeq2's beta prior variance.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum BetaPriorVarianceMethod {
    /// Match the upper absolute-beta quantile without mean/dispersion weights.
    Quantile,
    /// Match the upper absolute-beta quantile using DESeq2's
    /// `1 / (1 / baseMean + dispFit)` weights.
    Weighted,
}

/// Options for DESeq2-style beta prior variance estimation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BetaPriorVarianceOptions {
    /// Estimation method. DESeq2 defaults to `weighted`.
    pub method: BetaPriorVarianceMethod,
    /// Upper tail probability matched against a zero-centered Normal.
    pub upper_quantile: f64,
    /// Wide prior variance used for intercepts and columns with no finite betas.
    pub wide_prior_variance: f64,
    /// Absolute beta cutoff used to discard near-infinite MLEs.
    pub finite_beta_cutoff: f64,
}

impl Default for BetaPriorVarianceOptions {
    fn default() -> Self {
        Self {
            method: BetaPriorVarianceMethod::Weighted,
            upper_quantile: 0.05,
            wide_prior_variance: 1e6,
            finite_beta_cutoff: 10.0,
        }
    }
}

/// Two-stage fixed-dispersion GLM output for a DESeq2-style beta prior refit.
#[derive(Clone, Debug, PartialEq)]
pub struct BetaPriorGlmFit {
    /// First pass MLE GLM fit used to estimate beta prior variances.
    pub mle_fit: NbinomGlmFit,
    /// Refit GLM using beta-prior variance as a per-coefficient ridge.
    pub prior_fit: NbinomGlmFit,
    /// Log2-scale beta prior variances, one per design coefficient.
    pub beta_prior_variance: Vec<f64>,
}

/// Options for the two-stage beta-prior fixed-dispersion refit helper.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct BetaPriorRefitOptions {
    /// Options used for the first MLE GLM fit and the final prior refit.
    ///
    /// The prior refit replaces any ridge settings with values derived from
    /// the estimated beta prior variances.
    pub fit_options: IrlsOptions,
    /// Options used when estimating beta prior variances from the MLE fit.
    pub variance_options: BetaPriorVarianceOptions,
}

/// Size-factor offsets and optional observation weights for beta-prior refits.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BetaPriorSizeFactorWeightInput<'a> {
    /// Per-sample size factors.
    pub size_factors: &'a [f64],
    /// Optional normalized observation weights.
    pub weights: Option<&'a RowMajorMatrix<f64>>,
}

/// Normalization-factor offsets and optional observation weights for beta-prior refits.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BetaPriorNormalizationFactorWeightInput<'a> {
    /// Gene x sample normalization-factor matrix.
    pub normalization_factors: &'a RowMajorMatrix<f64>,
    /// Optional normalized observation weights.
    pub weights: Option<&'a RowMajorMatrix<f64>>,
}

/// Estimate fixed-dispersion beta coefficients with DESeq2-style GLM dispatch.
///
/// This is a public beta-estimation convenience entry point for callers that
/// already have size factors and per-gene dispersions. Intercept-only designs
/// use the closed-form DESeq2 shortcut through `fit_irls`; other designs use
/// the general fixed-dispersion IRLS implementation.
pub fn estimate_beta(
    counts: &CountMatrix,
    design: &DesignMatrix,
    size_factors: &[f64],
    dispersions: &[f64],
    options: IrlsOptions,
) -> Result<NbinomGlmFit, DeseqError> {
    fit_irls(counts, design, size_factors, dispersions, options)
}

/// Convert log2-scale beta prior variances to natural-log IRLS ridge values.
///
/// DESeq2 computes `lambda = 1 / betaPriorVar` on the log2 beta scale and then
/// divides by `log(2)^2` before fitting on the natural-log scale. This helper
/// exposes that conversion for primitive Rust GLM refits.
pub fn beta_prior_variance_to_ridge_lambda(
    beta_prior_variance: &[f64],
) -> Result<Vec<f64>, DeseqError> {
    beta_prior_variance
        .iter()
        .copied()
        .enumerate()
        .map(|(idx, variance)| {
            validate_positive_finite(variance, "beta prior variance", idx)?;
            let inv_ln2 = std::f64::consts::LOG2_E;
            Ok(variance.recip() * inv_ln2 * inv_ln2)
        })
        .collect()
}

/// Refit a fixed-dispersion GLM with supplied DESeq2-style beta prior variance.
pub fn fit_glms_with_beta_prior_variance(
    counts: &CountMatrix,
    design: &DesignMatrix,
    size_factors: &[f64],
    dispersions: &[f64],
    beta_prior_variance: &[f64],
    options: IrlsOptions,
) -> Result<NbinomGlmFit, DeseqError> {
    let options = options_with_beta_prior_variance(design, beta_prior_variance, options)?;
    fit_fixed_dispersion_irls(counts, design, size_factors, dispersions, options)
}

/// Refit a fixed-dispersion GLM with supplied beta prior variance, size factors, and weights.
pub fn fit_glms_with_beta_prior_variance_and_weights(
    counts: &CountMatrix,
    design: &DesignMatrix,
    size_factors: &[f64],
    dispersions: &[f64],
    weights: Option<&RowMajorMatrix<f64>>,
    beta_prior_variance: &[f64],
    options: IrlsOptions,
) -> Result<NbinomGlmFit, DeseqError> {
    let normalization_factors = normalization_factors_from_size_factors(counts, size_factors)?;
    fit_glms_with_beta_prior_variance_and_normalization_factors_and_weights(
        counts,
        design,
        &normalization_factors,
        dispersions,
        weights,
        beta_prior_variance,
        options,
    )
}

/// Refit a fixed-dispersion GLM with supplied beta prior variance and offsets.
pub fn fit_glms_with_beta_prior_variance_and_normalization_factors(
    counts: &CountMatrix,
    design: &DesignMatrix,
    normalization_factors: &RowMajorMatrix<f64>,
    dispersions: &[f64],
    beta_prior_variance: &[f64],
    options: IrlsOptions,
) -> Result<NbinomGlmFit, DeseqError> {
    fit_glms_with_beta_prior_variance_and_normalization_factors_and_weights(
        counts,
        design,
        normalization_factors,
        dispersions,
        None,
        beta_prior_variance,
        options,
    )
}

/// Refit a fixed-dispersion GLM with supplied beta prior variance, offsets, and weights.
pub fn fit_glms_with_beta_prior_variance_and_normalization_factors_and_weights(
    counts: &CountMatrix,
    design: &DesignMatrix,
    normalization_factors: &RowMajorMatrix<f64>,
    dispersions: &[f64],
    weights: Option<&RowMajorMatrix<f64>>,
    beta_prior_variance: &[f64],
    options: IrlsOptions,
) -> Result<NbinomGlmFit, DeseqError> {
    let options = options_with_beta_prior_variance(design, beta_prior_variance, options)?;
    fit_fixed_dispersion_irls_with_normalization_factors_and_weights(
        counts,
        design,
        normalization_factors,
        dispersions,
        weights,
        options,
    )
}

/// Run an MLE fixed-dispersion fit, estimate beta prior variance, then refit.
pub fn fit_glms_with_estimated_beta_prior_variance(
    counts: &CountMatrix,
    design: &DesignMatrix,
    size_factors: &[f64],
    dispersions: &[f64],
    base_mean: &[f64],
    disp_fit: &[f64],
    options: BetaPriorRefitOptions,
) -> Result<BetaPriorGlmFit, DeseqError> {
    let fit_options = options.fit_options;
    let mle_fit = fit_irls(
        counts,
        design,
        size_factors,
        dispersions,
        fit_options.clone(),
    )?;
    let beta_prior_variance = estimate_beta_prior_variance(
        &mle_fit.beta,
        base_mean,
        disp_fit,
        design.coefficient_names(),
        options.variance_options,
    )?;
    let prior_fit = fit_glms_with_beta_prior_variance(
        counts,
        design,
        size_factors,
        dispersions,
        &beta_prior_variance,
        fit_options,
    )?;

    Ok(BetaPriorGlmFit {
        mle_fit,
        prior_fit,
        beta_prior_variance,
    })
}

/// Run an MLE fixed-dispersion fit with size factors and weights, estimate beta prior variance, then refit.
pub fn fit_glms_with_estimated_beta_prior_variance_and_weights(
    counts: &CountMatrix,
    design: &DesignMatrix,
    input: BetaPriorSizeFactorWeightInput<'_>,
    dispersions: &[f64],
    base_mean: &[f64],
    disp_fit: &[f64],
    options: BetaPriorRefitOptions,
) -> Result<BetaPriorGlmFit, DeseqError> {
    let normalization_factors =
        normalization_factors_from_size_factors(counts, input.size_factors)?;
    fit_glms_with_estimated_beta_prior_variance_and_normalization_factors_and_weights(
        counts,
        design,
        BetaPriorNormalizationFactorWeightInput {
            normalization_factors: &normalization_factors,
            weights: input.weights,
        },
        dispersions,
        base_mean,
        disp_fit,
        options,
    )
}

/// Run an MLE fixed-dispersion fit with offsets, estimate beta prior variance, then refit.
pub fn fit_glms_with_estimated_beta_prior_variance_and_normalization_factors(
    counts: &CountMatrix,
    design: &DesignMatrix,
    normalization_factors: &RowMajorMatrix<f64>,
    dispersions: &[f64],
    base_mean: &[f64],
    disp_fit: &[f64],
    options: BetaPriorRefitOptions,
) -> Result<BetaPriorGlmFit, DeseqError> {
    fit_glms_with_estimated_beta_prior_variance_and_normalization_factors_and_weights(
        counts,
        design,
        BetaPriorNormalizationFactorWeightInput {
            normalization_factors,
            weights: None,
        },
        dispersions,
        base_mean,
        disp_fit,
        options,
    )
}

/// Run an MLE fixed-dispersion fit with offsets and weights, estimate beta prior variance, then refit.
pub fn fit_glms_with_estimated_beta_prior_variance_and_normalization_factors_and_weights(
    counts: &CountMatrix,
    design: &DesignMatrix,
    input: BetaPriorNormalizationFactorWeightInput<'_>,
    dispersions: &[f64],
    base_mean: &[f64],
    disp_fit: &[f64],
    options: BetaPriorRefitOptions,
) -> Result<BetaPriorGlmFit, DeseqError> {
    let fit_options = options.fit_options;
    let mle_fit = fit_fixed_dispersion_irls_with_normalization_factors_and_weights(
        counts,
        design,
        input.normalization_factors,
        dispersions,
        input.weights,
        fit_options.clone(),
    )?;
    let beta_prior_variance = estimate_beta_prior_variance(
        &mle_fit.beta,
        base_mean,
        disp_fit,
        design.coefficient_names(),
        options.variance_options,
    )?;
    let prior_fit = fit_glms_with_beta_prior_variance_and_normalization_factors_and_weights(
        counts,
        design,
        input.normalization_factors,
        dispersions,
        input.weights,
        &beta_prior_variance,
        fit_options,
    )?;

    Ok(BetaPriorGlmFit {
        mle_fit,
        prior_fit,
        beta_prior_variance,
    })
}

/// Estimate DESeq2-style beta prior variances from unshrunken MLE betas.
///
/// This mirrors the computational core of DESeq2 `estimateBetaPriorVar` for
/// already-built primitive matrices: each coefficient uses MLE betas whose
/// absolute value is below the finite-beta cutoff, then matches an upper
/// absolute-beta quantile to a zero-centered Normal. The weighted method uses
/// DESeq2's `1 / (1 / baseMean + dispFit)` row weights. Intercept columns are
/// assigned the configured wide prior variance.
pub fn estimate_beta_prior_variance(
    beta_matrix: &RowMajorMatrix<f64>,
    base_mean: &[f64],
    disp_fit: &[f64],
    coefficient_names: Option<&[String]>,
    options: BetaPriorVarianceOptions,
) -> Result<Vec<f64>, DeseqError> {
    validate_beta_prior_inputs(beta_matrix, base_mean, disp_fit, coefficient_names, options)?;
    let weights = match options.method {
        BetaPriorVarianceMethod::Weighted => Some(beta_prior_weights(
            base_mean,
            disp_fit,
            beta_matrix.n_rows(),
        )?),
        BetaPriorVarianceMethod::Quantile => None,
    };

    let mut prior_variance = Vec::with_capacity(beta_matrix.n_cols());
    for coefficient in 0..beta_matrix.n_cols() {
        let value = if beta_matrix.n_rows() == 1 {
            let beta = beta_matrix.row(0)?[coefficient];
            if beta.is_finite() {
                beta * beta
            } else {
                options.wide_prior_variance
            }
        } else {
            let mut betas = Vec::new();
            let mut selected_weights = Vec::new();
            for row in 0..beta_matrix.n_rows() {
                let beta = beta_matrix.row(row)?[coefficient];
                if beta.is_finite() && beta.abs() < options.finite_beta_cutoff {
                    betas.push(beta);
                    if let Some(weights) = &weights {
                        selected_weights.push(weights[row]);
                    }
                }
            }
            if betas.is_empty() {
                options.wide_prior_variance
            } else {
                match options.method {
                    BetaPriorVarianceMethod::Quantile => {
                        match_upper_quantile_for_variance(&betas, options.upper_quantile)?
                    }
                    BetaPriorVarianceMethod::Weighted => {
                        match_weighted_upper_quantile_for_variance(
                            &betas,
                            &selected_weights,
                            options.upper_quantile,
                        )?
                    }
                }
            }
        };
        prior_variance.push(value);
    }

    if let Some(names) = coefficient_names {
        for (idx, name) in names.iter().enumerate() {
            if name == "Intercept" || name == "(Intercept)" {
                prior_variance[idx] = options.wide_prior_variance;
            }
        }
    }
    Ok(prior_variance)
}

/// Match an upper absolute-beta quantile to a zero-centered Normal variance.
pub fn match_upper_quantile_for_variance(
    betas: &[f64],
    upper_quantile: f64,
) -> Result<f64, DeseqError> {
    validate_upper_quantile(upper_quantile)?;
    let abs_betas = finite_abs_values(betas)?;
    let quantile = quantile_type7(abs_betas, 1.0 - upper_quantile)?;
    let normal = Normal::new(0.0, 1.0).map_err(|error| DeseqError::InvalidOptions {
        reason: format!("normal quantile construction failed: {error}"),
    })?;
    let normal_quantile = normal.inverse_cdf(1.0 - upper_quantile / 2.0);
    let scale = quantile / normal_quantile;
    Ok(scale * scale)
}

/// Weighted version of [`match_upper_quantile_for_variance`].
pub fn match_weighted_upper_quantile_for_variance(
    betas: &[f64],
    weights: &[f64],
    upper_quantile: f64,
) -> Result<f64, DeseqError> {
    validate_upper_quantile(upper_quantile)?;
    if betas.len() != weights.len() {
        return Err(invalid_dimensions(
            "beta prior variance weights",
            betas.len(),
            weights.len(),
        ));
    }
    let weighted_quantile = weighted_abs_quantile(betas, weights, 1.0 - upper_quantile)?;
    let normal = Normal::new(0.0, 1.0).map_err(|error| DeseqError::InvalidOptions {
        reason: format!("normal quantile construction failed: {error}"),
    })?;
    let normal_quantile = normal.inverse_cdf(1.0 - upper_quantile / 2.0);
    let scale = weighted_quantile / normal_quantile;
    Ok(scale * scale)
}

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
                normalized_row.iter().sum::<f64>() / counts.n_samples() as f64
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
            mu_values.push(*factor * 2.0_f64.powf(beta_log2));
        }

        let mu_start = gene * counts.n_samples();
        let mu_row = &mu_values[mu_start..mu_start + counts.n_samples()];
        let working_weights = intercept_working_weights(mu_row, dispersion, weight_row)?;
        let xtwx = working_weights.iter().sum::<f64>();
        if !xtwx.is_finite() || xtwx <= 0.0 {
            return Err(DeseqError::InvalidDispersion {
                reason: format!("gene {gene} has non-positive intercept information"),
            });
        }
        let sigma = xtwx.recip();
        beta_se_values.push(std::f64::consts::LOG2_E * sigma.sqrt());
        let log2_e = std::f64::consts::LOG2_E;
        beta_covariance_values.push(log2_e * log2_e * sigma);
        hat_values.extend(working_weights.into_iter().map(|value| value * sigma));
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
        beta_covariance: Some(beta_covariance),
        mu,
        beta_iter: vec![1; counts.n_genes()],
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
    let mut numerator = 0.0;
    let mut denominator = 0.0;
    for (sample, (value, weight)) in values
        .iter()
        .copied()
        .zip(weights.iter().copied())
        .enumerate()
    {
        validate_nonnegative_finite(weight, "weight", sample)?;
        numerator += weight * value;
        denominator += weight;
    }
    if denominator <= 0.0 {
        return Err(DeseqError::InvalidCounts {
            reason: format!("gene {gene} has zero total weight"),
        });
    }
    Ok(numerator / denominator)
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
        let working_weight = mean_dispersion_working_weight(value, dispersion);
        out.push(match weights {
            Some(weights) => {
                let weight = weights[sample];
                validate_nonnegative_finite(weight, "weight", sample)?;
                weight * working_weight
            }
            None => working_weight,
        });
    }
    Ok(out)
}

fn options_with_beta_prior_variance(
    design: &DesignMatrix,
    beta_prior_variance: &[f64],
    options: IrlsOptions,
) -> Result<IrlsOptions, DeseqError> {
    if beta_prior_variance.len() != design.n_coefficients() {
        return Err(invalid_dimensions(
            "beta prior variance coefficients",
            design.n_coefficients(),
            beta_prior_variance.len(),
        ));
    }
    let ridge_lambda = beta_prior_variance_to_ridge_lambda(beta_prior_variance)?;
    Ok(options.ridge_lambda_by_coefficient(ridge_lambda))
}

fn validate_beta_prior_inputs(
    beta_matrix: &RowMajorMatrix<f64>,
    base_mean: &[f64],
    disp_fit: &[f64],
    coefficient_names: Option<&[String]>,
    options: BetaPriorVarianceOptions,
) -> Result<(), DeseqError> {
    if base_mean.len() != beta_matrix.n_rows() {
        return Err(invalid_dimensions(
            "beta prior variance baseMean rows",
            beta_matrix.n_rows(),
            base_mean.len(),
        ));
    }
    if disp_fit.len() != beta_matrix.n_rows() {
        return Err(invalid_dimensions(
            "beta prior variance dispFit rows",
            beta_matrix.n_rows(),
            disp_fit.len(),
        ));
    }
    if let Some(names) = coefficient_names {
        if names.len() != beta_matrix.n_cols() {
            return Err(invalid_dimensions(
                "beta prior variance coefficient names",
                beta_matrix.n_cols(),
                names.len(),
            ));
        }
    }
    validate_upper_quantile(options.upper_quantile)?;
    if !options.wide_prior_variance.is_finite() || options.wide_prior_variance <= 0.0 {
        return Err(DeseqError::InvalidOptions {
            reason: "wide beta prior variance must be finite and positive".to_string(),
        });
    }
    if !options.finite_beta_cutoff.is_finite() || options.finite_beta_cutoff <= 0.0 {
        return Err(DeseqError::InvalidOptions {
            reason: "finite beta cutoff must be finite and positive".to_string(),
        });
    }
    Ok(())
}

fn beta_prior_weights(
    base_mean: &[f64],
    disp_fit: &[f64],
    n_rows: usize,
) -> Result<Vec<f64>, DeseqError> {
    let mut weights = Vec::with_capacity(n_rows);
    for row in 0..n_rows {
        validate_positive_finite(base_mean[row], "beta prior baseMean", row)?;
        validate_positive_finite(disp_fit[row], "beta prior dispFit", row)?;
        weights.push(mean_dispersion_working_weight(
            base_mean[row],
            disp_fit[row],
        ));
    }
    Ok(weights)
}

fn mean_dispersion_working_weight(mean: f64, dispersion: f64) -> f64 {
    mean / (1.0 + mean * dispersion)
}

fn validate_upper_quantile(upper_quantile: f64) -> Result<(), DeseqError> {
    if !upper_quantile.is_finite() || upper_quantile <= 0.0 || upper_quantile >= 1.0 {
        return Err(DeseqError::InvalidOptions {
            reason: "upper quantile must be finite and between 0 and 1".to_string(),
        });
    }
    Ok(())
}

fn finite_abs_values(values: &[f64]) -> Result<Vec<f64>, DeseqError> {
    let out = values
        .iter()
        .copied()
        .filter(|value| value.is_finite())
        .map(f64::abs)
        .collect::<Vec<_>>();
    if out.is_empty() {
        return Err(DeseqError::InvalidOptions {
            reason: "beta prior variance quantile needs at least one finite beta".to_string(),
        });
    }
    Ok(out)
}

fn quantile_type7(mut values: Vec<f64>, probability: f64) -> Result<f64, DeseqError> {
    if !probability.is_finite() || !(0.0..=1.0).contains(&probability) {
        return Err(DeseqError::InvalidOptions {
            reason: "quantile probability must be finite and between 0 and 1".to_string(),
        });
    }
    values.sort_by(|a, b| a.total_cmp(b));
    if values.len() == 1 {
        return Ok(values[0]);
    }
    let h = (values.len() as f64 - 1.0) * probability + 1.0;
    let lower = h.floor() as usize;
    let fraction = h - lower as f64;
    if lower == 0 {
        Ok(values[0])
    } else if lower >= values.len() {
        Ok(values[values.len() - 1])
    } else {
        Ok(values[lower - 1] + fraction * (values[lower] - values[lower - 1]))
    }
}

fn weighted_abs_quantile(
    betas: &[f64],
    weights: &[f64],
    probability: f64,
) -> Result<f64, DeseqError> {
    if !probability.is_finite() || !(0.0..=1.0).contains(&probability) {
        return Err(DeseqError::InvalidOptions {
            reason: "weighted quantile probability must be finite and between 0 and 1".to_string(),
        });
    }
    let mut pairs = Vec::with_capacity(betas.len());
    for (idx, (beta, weight)) in betas
        .iter()
        .copied()
        .zip(weights.iter().copied())
        .enumerate()
    {
        if beta.is_finite() && weight.is_finite() && weight > 0.0 {
            pairs.push((beta.abs(), weight));
        } else if !weight.is_finite() || weight < 0.0 {
            return Err(DeseqError::NonFiniteValue {
                context: "beta prior variance weight".to_string(),
                index: Some(idx),
                value: weight,
            });
        }
    }
    if pairs.is_empty() {
        return Err(DeseqError::InvalidOptions {
            reason: "weighted beta prior variance quantile needs positive finite weights"
                .to_string(),
        });
    }
    pairs.sort_by(|a, b| a.0.total_cmp(&b.0));

    let weight_sum = pairs.iter().map(|(_, weight)| *weight).sum::<f64>();
    if !weight_sum.is_finite() || weight_sum <= 0.0 {
        return Err(DeseqError::InvalidOptions {
            reason: "weighted beta prior variance quantile needs positive finite total weight"
                .to_string(),
        });
    }

    let norm_scale = pairs.len() as f64 / weight_sum;
    let mut unique = Vec::<(f64, f64)>::with_capacity(pairs.len());
    for (value, weight) in pairs {
        let normalized_weight = weight * norm_scale;
        if let Some((last_value, last_weight)) = unique.last_mut() {
            if *last_value == value {
                *last_weight += normalized_weight;
                continue;
            }
        }
        unique.push((value, normalized_weight));
    }

    if unique.len() == 1 {
        return Ok(unique[0].0);
    }

    let n = unique.iter().map(|(_, weight)| *weight).sum::<f64>();
    let order = 1.0 + (n - 1.0) * probability;
    let low = order.floor().max(1.0);
    let high = (low + 1.0).min(n);
    let fraction = order.fract();

    let mut cumulative_weights = Vec::with_capacity(unique.len());
    let mut cumulative = 0.0;
    for (_, weight) in &unique {
        cumulative += *weight;
        cumulative_weights.push(cumulative);
    }

    let low_quantile = weighted_order_statistic(&unique, &cumulative_weights, low);
    let high_quantile = weighted_order_statistic(&unique, &cumulative_weights, high);
    Ok((1.0 - fraction) * low_quantile + fraction * high_quantile)
}

fn weighted_order_statistic(
    values_and_weights: &[(f64, f64)],
    cumulative_weights: &[f64],
    position: f64,
) -> f64 {
    for (idx, cumulative) in cumulative_weights.iter().copied().enumerate() {
        if position <= cumulative {
            return values_and_weights[idx].0;
        }
    }
    values_and_weights[values_and_weights.len() - 1].0
}

fn validate_fit_inputs(
    counts: &CountMatrix,
    size_factors: &[f64],
    dispersions: &[f64],
    weights: Option<&RowMajorMatrix<f64>>,
) -> Result<(), DeseqError> {
    if size_factors.len() != counts.n_samples() {
        return Err(invalid_dimensions(
            "size factors",
            counts.n_samples(),
            size_factors.len(),
        ));
    }
    for (idx, factor) in size_factors.iter().copied().enumerate() {
        validate_positive_finite(factor, "size factor", idx)?;
    }
    let normalization_factors = normalization_factors_from_size_factors(counts, size_factors)?;
    validate_nf_fit_inputs(counts, &normalization_factors, dispersions, weights)
}

fn validate_nf_fit_inputs(
    counts: &CountMatrix,
    normalization_factors: &RowMajorMatrix<f64>,
    dispersions: &[f64],
    weights: Option<&RowMajorMatrix<f64>>,
) -> Result<(), DeseqError> {
    if normalization_factors.n_rows() != counts.n_genes()
        || normalization_factors.n_cols() != counts.n_samples()
    {
        return Err(DeseqError::InvalidDimensions {
            context: "normalization factors".to_string(),
            expected: counts.n_genes() * counts.n_samples(),
            actual: normalization_factors.len(),
        });
    }
    if dispersions.len() != counts.n_genes() {
        return Err(invalid_dimensions(
            "dispersions",
            counts.n_genes(),
            dispersions.len(),
        ));
    }
    for (idx, dispersion) in dispersions.iter().copied().enumerate() {
        validate_positive_finite(dispersion, "dispersion", idx)?;
    }
    if let Some(weights) = weights {
        if weights.n_rows() != counts.n_genes() || weights.n_cols() != counts.n_samples() {
            return Err(DeseqError::InvalidDimensions {
                context: "weights".to_string(),
                expected: counts.n_genes() * counts.n_samples(),
                actual: weights.len(),
            });
        }
        weights.validate_finite("weights")?;
    }
    normalization_factors.validate_finite("normalization factors")?;
    Ok(())
}

fn validate_positive_finite(value: f64, context: &str, index: usize) -> Result<(), DeseqError> {
    if !value.is_finite() || value <= 0.0 {
        return Err(DeseqError::NonFiniteValue {
            context: context.to_string(),
            index: Some(index),
            value,
        });
    }
    Ok(())
}

fn validate_nonnegative_finite(value: f64, context: &str, index: usize) -> Result<(), DeseqError> {
    if !value.is_finite() || value < 0.0 {
        return Err(DeseqError::NonFiniteValue {
            context: context.to_string(),
            index: Some(index),
            value,
        });
    }
    Ok(())
}
