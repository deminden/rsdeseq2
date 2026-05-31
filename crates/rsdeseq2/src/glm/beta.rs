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

/// Expanded-model fixed-dispersion GLM output for a DESeq2-style beta-prior refit.
#[derive(Clone, Debug, PartialEq)]
pub struct ExpandedModelBetaPriorGlmFit {
    /// First pass MLE GLM fit on the expanded design.
    pub expanded_mle_fit: NbinomGlmFit,
    /// Refit GLM on the expanded design using beta-prior variance as ridge.
    pub expanded_prior_fit: NbinomGlmFit,
    /// Prior fit collapsed onto the caller-supplied standard design surface.
    pub prior_fit: NbinomGlmFit,
    /// Log2-scale beta prior variances, one per expanded-design coefficient.
    pub beta_prior_variance: Vec<f64>,
}

/// Expanded and reported design surfaces for a beta-prior expanded-model refit.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ExpandedModelBetaPriorDesignInput<'a> {
    /// Design used for expanded-model fitting.
    pub expanded_design: &'a DesignMatrix,
    /// Reported standard design used after grouped coefficient collapse.
    pub standard_design: &'a DesignMatrix,
    /// Expanded coefficient columns averaged into each standard coefficient.
    pub coefficient_groups: &'a [Vec<usize>],
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

/// Average expanded-model coefficients into standard model-term coefficients.
///
/// DESeq2's beta-prior workflow can fit an expanded model matrix and then
/// collapse groups of expanded columns back to the reported standard
/// coefficients. This primitive performs that numeric collapse for an existing
/// log2-scale beta matrix: each output column is the per-gene average of the
/// input columns listed in the corresponding group.
pub fn average_expanded_model_coefficients(
    expanded_beta: &RowMajorMatrix<f64>,
    coefficient_groups: &[Vec<usize>],
) -> Result<RowMajorMatrix<f64>, DeseqError> {
    validate_expanded_coefficient_groups(
        expanded_beta.n_cols(),
        coefficient_groups,
        "expanded beta coefficient groups",
    )?;

    let mut values = Vec::with_capacity(expanded_beta.n_rows() * coefficient_groups.len());
    for gene in 0..expanded_beta.n_rows() {
        let row = expanded_beta.row(gene)?;
        for group in coefficient_groups {
            let grouped = group.iter().map(|&column| row[column]).collect::<Vec<_>>();
            values.push(checked_scaled_mean(&grouped).ok_or_else(|| {
                DeseqError::NonFiniteValue {
                    context: "expanded beta coefficient average".to_string(),
                    index: Some(gene),
                    value: f64::NAN,
                }
            })?);
        }
    }

    RowMajorMatrix::from_row_major(expanded_beta.n_rows(), coefficient_groups.len(), values)
}

/// Propagate expanded-model covariance through grouped coefficient averaging.
///
/// The input covariance matrix stores each gene's expanded-model coefficient
/// covariance as one row-major `p x p` row. The output stores each gene's
/// collapsed `q x q` covariance matrix, where `q = coefficient_groups.len()`
/// and every output cell is `A Sigma A'` for the averaging matrix implied by
/// `coefficient_groups`.
pub fn average_expanded_model_covariances(
    expanded_covariance: &RowMajorMatrix<f64>,
    n_expanded_coefficients: usize,
    coefficient_groups: &[Vec<usize>],
) -> Result<RowMajorMatrix<f64>, DeseqError> {
    validate_expanded_covariance_inputs(
        expanded_covariance,
        n_expanded_coefficients,
        coefficient_groups,
    )?;

    let n_groups = coefficient_groups.len();
    let mut values = Vec::with_capacity(expanded_covariance.n_rows() * n_groups * n_groups);
    for gene in 0..expanded_covariance.n_rows() {
        let covariance = expanded_covariance.row(gene)?;
        for left_group in coefficient_groups {
            for right_group in coefficient_groups {
                values.push(average_covariance_block(
                    covariance,
                    n_expanded_coefficients,
                    left_group,
                    right_group,
                    gene,
                )?);
            }
        }
    }

    RowMajorMatrix::from_row_major(expanded_covariance.n_rows(), n_groups * n_groups, values)
}

/// Collapse an expanded-model GLM fit to a standard coefficient surface.
///
/// This keeps gene/sample diagnostics and fitted means from the expanded fit,
/// replaces the reported model matrix with `standard_design`, averages beta
/// coefficients by group, propagates covariance through the same averaging
/// matrix, and recomputes standard errors from the collapsed covariance
/// diagonal. It is a primitive building block for DESeq2's beta-prior expanded
/// model workflow; callers remain responsible for constructing the expanded
/// and standard design matrices from formula metadata.
pub fn collapse_expanded_model_fit(
    expanded_fit: &NbinomGlmFit,
    standard_design: &DesignMatrix,
    coefficient_groups: &[Vec<usize>],
) -> Result<NbinomGlmFit, DeseqError> {
    validate_expanded_fit_collapse_inputs(expanded_fit, standard_design, coefficient_groups)?;
    let covariance =
        expanded_fit
            .beta_covariance
            .as_ref()
            .ok_or_else(|| DeseqError::UnsupportedFeature {
                feature: "expanded model fit collapse requires beta covariance matrices"
                    .to_string(),
            })?;
    let beta = average_expanded_model_coefficients(&expanded_fit.beta, coefficient_groups)?;
    let beta_covariance = average_expanded_model_covariances(
        covariance,
        expanded_fit.beta.n_cols(),
        coefficient_groups,
    )?;
    let beta_se = beta_se_from_covariance(&beta_covariance, beta.n_cols())?;

    Ok(NbinomGlmFit {
        log_like: expanded_fit.log_like.clone(),
        beta_converged: expanded_fit.beta_converged.clone(),
        beta,
        beta_se,
        beta_optim_start: expanded_fit.beta_optim_start.clone(),
        beta_covariance: Some(beta_covariance),
        mu: expanded_fit.mu.clone(),
        beta_iter: expanded_fit.beta_iter.clone(),
        beta_optim_iter: expanded_fit.beta_optim_iter.clone(),
        beta_optim_start_objective: expanded_fit.beta_optim_start_objective.clone(),
        beta_optim_objective: expanded_fit.beta_optim_objective.clone(),
        beta_optim_gradient_norm: expanded_fit.beta_optim_gradient_norm.clone(),
        model_matrix: standard_design.clone(),
        n_terms: standard_design.n_coefficients(),
        hat_diagonal: expanded_fit.hat_diagonal.clone(),
    })
}

/// Build a numeric contrast vector from two expanded-model coefficient groups.
///
/// The returned vector has `+1 / n` weights across numerator columns and
/// `-1 / m` weights across denominator columns, matching the averaging shape
/// needed to compare collapsed expanded-model terms.
pub fn expanded_model_group_contrast(
    n_coefficients: usize,
    numerator_columns: &[usize],
    denominator_columns: &[usize],
) -> Result<Vec<f64>, DeseqError> {
    validate_expanded_contrast_group(
        n_coefficients,
        numerator_columns,
        "expanded contrast numerator columns",
    )?;
    validate_expanded_contrast_group(
        n_coefficients,
        denominator_columns,
        "expanded contrast denominator columns",
    )?;
    for &column in numerator_columns {
        if denominator_columns.contains(&column) {
            return Err(DeseqError::InvalidOptions {
                reason: "expanded contrast numerator and denominator columns must be disjoint"
                    .to_string(),
            });
        }
    }

    let numerator_weight = checked_div2(1.0, numerator_columns.len() as f64).ok_or_else(|| {
        DeseqError::InvalidOptions {
            reason: "expanded contrast numerator weight is non-finite".to_string(),
        }
    })?;
    let denominator_weight =
        checked_div2(-1.0, denominator_columns.len() as f64).ok_or_else(|| {
            DeseqError::InvalidOptions {
                reason: "expanded contrast denominator weight is non-finite".to_string(),
            }
        })?;

    let mut contrast = vec![0.0; n_coefficients];
    for &column in numerator_columns {
        contrast[column] = numerator_weight;
    }
    for &column in denominator_columns {
        contrast[column] = denominator_weight;
    }
    Ok(contrast)
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

/// Run an expanded-design MLE fit, beta-prior refit, and standard-design collapse.
///
/// This helper covers the numeric core of DESeq2's expanded-model beta-prior
/// path for callers that already constructed both model matrices and the
/// mapping from expanded columns to reported standard coefficients.
pub fn fit_expanded_glms_with_estimated_beta_prior_variance(
    counts: &CountMatrix,
    design: ExpandedModelBetaPriorDesignInput<'_>,
    size_factors: &[f64],
    dispersions: &[f64],
    base_mean: &[f64],
    disp_fit: &[f64],
    options: BetaPriorRefitOptions,
) -> Result<ExpandedModelBetaPriorGlmFit, DeseqError> {
    let normalization_factors = normalization_factors_from_size_factors(counts, size_factors)?;
    fit_expanded_glms_with_estimated_beta_prior_variance_and_normalization_factors_and_weights(
        counts,
        design,
        BetaPriorNormalizationFactorWeightInput {
            normalization_factors: &normalization_factors,
            weights: None,
        },
        dispersions,
        base_mean,
        disp_fit,
        options,
    )
}

/// Run an expanded-design MLE fit, beta-prior refit, weights, and standard-design collapse.
pub fn fit_expanded_glms_with_estimated_beta_prior_variance_and_weights(
    counts: &CountMatrix,
    design: ExpandedModelBetaPriorDesignInput<'_>,
    input: BetaPriorSizeFactorWeightInput<'_>,
    dispersions: &[f64],
    base_mean: &[f64],
    disp_fit: &[f64],
    options: BetaPriorRefitOptions,
) -> Result<ExpandedModelBetaPriorGlmFit, DeseqError> {
    let normalization_factors =
        normalization_factors_from_size_factors(counts, input.size_factors)?;
    fit_expanded_glms_with_estimated_beta_prior_variance_and_normalization_factors_and_weights(
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

/// Run an expanded-design MLE fit, beta-prior refit with offsets, and standard-design collapse.
pub fn fit_expanded_glms_with_estimated_beta_prior_variance_and_normalization_factors(
    counts: &CountMatrix,
    design: ExpandedModelBetaPriorDesignInput<'_>,
    normalization_factors: &RowMajorMatrix<f64>,
    dispersions: &[f64],
    base_mean: &[f64],
    disp_fit: &[f64],
    options: BetaPriorRefitOptions,
) -> Result<ExpandedModelBetaPriorGlmFit, DeseqError> {
    fit_expanded_glms_with_estimated_beta_prior_variance_and_normalization_factors_and_weights(
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

/// Run an expanded-design MLE fit, beta-prior refit with offsets/weights, and collapse.
pub fn fit_expanded_glms_with_estimated_beta_prior_variance_and_normalization_factors_and_weights(
    counts: &CountMatrix,
    design: ExpandedModelBetaPriorDesignInput<'_>,
    input: BetaPriorNormalizationFactorWeightInput<'_>,
    dispersions: &[f64],
    base_mean: &[f64],
    disp_fit: &[f64],
    options: BetaPriorRefitOptions,
) -> Result<ExpandedModelBetaPriorGlmFit, DeseqError> {
    validate_expanded_coefficient_groups(
        design.expanded_design.n_coefficients(),
        design.coefficient_groups,
        "expanded beta-prior coefficient groups",
    )?;
    if design.coefficient_groups.len() != design.standard_design.n_coefficients() {
        return Err(invalid_dimensions(
            "expanded beta-prior collapsed coefficient groups",
            design.standard_design.n_coefficients(),
            design.coefficient_groups.len(),
        ));
    }

    let fit_options = options.fit_options;
    let expanded_mle_fit = fit_fixed_dispersion_irls_with_normalization_factors_and_weights(
        counts,
        design.expanded_design,
        input.normalization_factors,
        dispersions,
        input.weights,
        fit_options.clone(),
    )?;
    let beta_prior_variance = estimate_beta_prior_variance(
        &expanded_mle_fit.beta,
        base_mean,
        disp_fit,
        design.expanded_design.coefficient_names(),
        options.variance_options,
    )?;
    let expanded_prior_fit =
        fit_glms_with_beta_prior_variance_and_normalization_factors_and_weights(
            counts,
            design.expanded_design,
            input.normalization_factors,
            dispersions,
            input.weights,
            &beta_prior_variance,
            fit_options,
        )?;
    let prior_fit = collapse_expanded_model_fit(
        &expanded_prior_fit,
        design.standard_design,
        design.coefficient_groups,
    )?;

    Ok(ExpandedModelBetaPriorGlmFit {
        expanded_mle_fit,
        expanded_prior_fit,
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
    checked_square(scale).ok_or_else(|| DeseqError::InvalidOptions {
        reason: "beta prior variance quantile produced non-finite variance".to_string(),
    })
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
    checked_square(scale).ok_or_else(|| DeseqError::InvalidOptions {
        reason: "weighted beta prior variance quantile produced non-finite variance".to_string(),
    })
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

fn validate_expanded_coefficient_groups(
    n_coefficients: usize,
    coefficient_groups: &[Vec<usize>],
    context: &str,
) -> Result<(), DeseqError> {
    if coefficient_groups.is_empty() {
        return Err(DeseqError::InvalidDimensions {
            context: context.to_string(),
            expected: 1,
            actual: 0,
        });
    }
    for (group_idx, group) in coefficient_groups.iter().enumerate() {
        validate_expanded_contrast_group(
            n_coefficients,
            group,
            &format!("{context} group {group_idx}"),
        )?;
    }
    Ok(())
}

fn validate_expanded_covariance_inputs(
    expanded_covariance: &RowMajorMatrix<f64>,
    n_expanded_coefficients: usize,
    coefficient_groups: &[Vec<usize>],
) -> Result<(), DeseqError> {
    if n_expanded_coefficients == 0 {
        return Err(DeseqError::InvalidDimensions {
            context: "expanded covariance coefficient count".to_string(),
            expected: 1,
            actual: 0,
        });
    }
    let expected = n_expanded_coefficients
        .checked_mul(n_expanded_coefficients)
        .ok_or_else(|| DeseqError::InvalidDimensions {
            context: "expanded covariance columns".to_string(),
            expected: usize::MAX,
            actual: expanded_covariance.n_cols(),
        })?;
    if expanded_covariance.n_cols() != expected {
        return Err(invalid_dimensions(
            "expanded covariance columns",
            expected,
            expanded_covariance.n_cols(),
        ));
    }
    validate_expanded_coefficient_groups(
        n_expanded_coefficients,
        coefficient_groups,
        "expanded covariance coefficient groups",
    )
}

fn validate_expanded_fit_collapse_inputs(
    expanded_fit: &NbinomGlmFit,
    standard_design: &DesignMatrix,
    coefficient_groups: &[Vec<usize>],
) -> Result<(), DeseqError> {
    validate_expanded_coefficient_groups(
        expanded_fit.beta.n_cols(),
        coefficient_groups,
        "expanded fit coefficient groups",
    )?;
    if standard_design.n_coefficients() != coefficient_groups.len() {
        return Err(invalid_dimensions(
            "standard design coefficients",
            coefficient_groups.len(),
            standard_design.n_coefficients(),
        ));
    }
    if standard_design.n_samples() != expanded_fit.mu.n_cols() {
        return Err(invalid_dimensions(
            "standard design samples",
            expanded_fit.mu.n_cols(),
            standard_design.n_samples(),
        ));
    }
    validate_expanded_fit_row_vector(
        "expanded fit log-likelihood rows",
        expanded_fit.beta.n_rows(),
        expanded_fit.log_like.len(),
    )?;
    validate_expanded_fit_row_vector(
        "expanded fit convergence rows",
        expanded_fit.beta.n_rows(),
        expanded_fit.beta_converged.len(),
    )?;
    validate_expanded_fit_row_vector(
        "expanded fit beta-iteration rows",
        expanded_fit.beta.n_rows(),
        expanded_fit.beta_iter.len(),
    )?;
    if expanded_fit.beta_se.n_rows() != expanded_fit.beta.n_rows()
        || expanded_fit.beta_se.n_cols() != expanded_fit.beta.n_cols()
    {
        return Err(DeseqError::InvalidDimensions {
            context: "expanded fit beta SE matrix".to_string(),
            expected: expanded_fit.beta.len(),
            actual: expanded_fit.beta_se.len(),
        });
    }
    if expanded_fit.mu.n_rows() != expanded_fit.beta.n_rows() {
        return Err(invalid_dimensions(
            "expanded fit fitted mean rows",
            expanded_fit.beta.n_rows(),
            expanded_fit.mu.n_rows(),
        ));
    }
    if expanded_fit.hat_diagonal.n_rows() != expanded_fit.beta.n_rows()
        || expanded_fit.hat_diagonal.n_cols() != expanded_fit.mu.n_cols()
    {
        return Err(DeseqError::InvalidDimensions {
            context: "expanded fit hat diagonal matrix".to_string(),
            expected: expanded_fit.beta.n_rows() * expanded_fit.mu.n_cols(),
            actual: expanded_fit.hat_diagonal.len(),
        });
    }
    Ok(())
}

fn validate_expanded_fit_row_vector(
    context: &str,
    expected: usize,
    actual: usize,
) -> Result<(), DeseqError> {
    if actual != expected {
        return Err(invalid_dimensions(context, expected, actual));
    }
    Ok(())
}

fn validate_expanded_contrast_group(
    n_coefficients: usize,
    columns: &[usize],
    context: &str,
) -> Result<(), DeseqError> {
    if columns.is_empty() {
        return Err(DeseqError::InvalidDimensions {
            context: context.to_string(),
            expected: 1,
            actual: 0,
        });
    }
    let mut seen = vec![false; n_coefficients];
    for &column in columns {
        if column >= n_coefficients {
            return Err(DeseqError::InvalidDimensions {
                context: context.to_string(),
                expected: n_coefficients.saturating_sub(1),
                actual: column,
            });
        }
        if seen[column] {
            return Err(DeseqError::InvalidOptions {
                reason: format!("{context} contains duplicate column {column}"),
            });
        }
        seen[column] = true;
    }
    Ok(())
}

fn average_covariance_block(
    covariance: &[f64],
    n_coefficients: usize,
    left_group: &[usize],
    right_group: &[usize],
    gene: usize,
) -> Result<f64, DeseqError> {
    let terms = left_group
        .iter()
        .flat_map(|&left| {
            right_group
                .iter()
                .map(move |&right| covariance[left * n_coefficients + right])
        })
        .collect::<Vec<_>>();
    let sum = checked_scaled_sum(terms).ok_or_else(|| DeseqError::NonFiniteValue {
        context: "expanded beta covariance average".to_string(),
        index: Some(gene),
        value: f64::NAN,
    })?;
    let denominator = left_group
        .len()
        .checked_mul(right_group.len())
        .ok_or_else(|| DeseqError::InvalidDimensions {
            context: "expanded covariance average denominator".to_string(),
            expected: usize::MAX,
            actual: 0,
        })?;
    checked_div2(sum, denominator as f64).ok_or_else(|| DeseqError::NonFiniteValue {
        context: "expanded beta covariance average".to_string(),
        index: Some(gene),
        value: f64::NAN,
    })
}

fn beta_se_from_covariance(
    covariance: &RowMajorMatrix<f64>,
    n_coefficients: usize,
) -> Result<RowMajorMatrix<f64>, DeseqError> {
    let expected_cols = n_coefficients.checked_mul(n_coefficients).ok_or_else(|| {
        DeseqError::InvalidDimensions {
            context: "collapsed covariance columns".to_string(),
            expected: usize::MAX,
            actual: covariance.n_cols(),
        }
    })?;
    if covariance.n_cols() != expected_cols {
        return Err(invalid_dimensions(
            "collapsed covariance columns",
            expected_cols,
            covariance.n_cols(),
        ));
    }
    let mut values = Vec::with_capacity(covariance.n_rows() * n_coefficients);
    for gene in 0..covariance.n_rows() {
        let row = covariance.row(gene)?;
        for coefficient in 0..n_coefficients {
            let variance = row[coefficient * n_coefficients + coefficient];
            if !variance.is_finite() {
                return Err(DeseqError::NonFiniteValue {
                    context: "collapsed beta standard error".to_string(),
                    index: Some(gene),
                    value: variance,
                });
            }
            values.push(variance.max(0.0).sqrt());
        }
    }
    RowMajorMatrix::from_row_major(covariance.n_rows(), n_coefficients, values)
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
        weights.push(checked_mean_dispersion_working_weight(
            base_mean[row],
            disp_fit[row],
            row,
            "beta prior row weight",
        )?);
    }
    Ok(weights)
}

fn checked_mean_dispersion_working_weight(
    mean: f64,
    dispersion: f64,
    index: usize,
    context: &str,
) -> Result<f64, DeseqError> {
    let Some(inv_mean) = checked_div2(1.0, mean) else {
        return Err(DeseqError::NonFiniteValue {
            context: context.to_string(),
            index: Some(index),
            value: f64::NAN,
        });
    };
    let Some(denominator) = checked_sum2(inv_mean, dispersion) else {
        return Err(DeseqError::NonFiniteValue {
            context: context.to_string(),
            index: Some(index),
            value: f64::NAN,
        });
    };
    checked_div2(1.0, denominator).ok_or_else(|| DeseqError::NonFiniteValue {
        context: context.to_string(),
        index: Some(index),
        value: f64::NAN,
    })
}

fn checked_mean(values: &[f64]) -> Option<f64> {
    checked_scaled_mean(values)
}

fn checked_scaled_sum(values: impl IntoIterator<Item = f64>) -> Option<f64> {
    let values = values.into_iter().collect::<Vec<_>>();
    let mut scale = 0.0_f64;
    for value in values.iter().copied() {
        if !value.is_finite() {
            return None;
        }
        scale = scale.max(value.abs());
    }
    if scale == 0.0 {
        return Some(0.0);
    }
    let normalized_sum = checked_sum(
        values
            .iter()
            .copied()
            .map(|value| checked_div2(value, scale))
            .collect::<Option<Vec<_>>>()?,
    )?;
    checked_product2(normalized_sum, scale)
}

fn checked_scaled_mean(values: &[f64]) -> Option<f64> {
    let mut scale = 0.0_f64;
    for value in values.iter().copied() {
        if !value.is_finite() {
            return None;
        }
        scale = scale.max(value.abs());
    }
    if scale == 0.0 {
        return Some(0.0);
    }
    let normalized_sum = checked_sum(
        values
            .iter()
            .copied()
            .map(|value| checked_div2(value, scale))
            .collect::<Option<Vec<_>>>()?,
    )?;
    checked_product2(checked_div2(normalized_sum, values.len() as f64)?, scale)
}

fn checked_sum(values: impl IntoIterator<Item = f64>) -> Option<f64> {
    values.into_iter().try_fold(0.0, checked_sum2)
}

fn checked_sum2(left: f64, right: f64) -> Option<f64> {
    let sum = left + right;
    (left.is_finite() && right.is_finite() && sum.is_finite()).then_some(sum)
}

fn checked_square(value: f64) -> Option<f64> {
    let square = value * value;
    (value.is_finite() && square.is_finite()).then_some(square)
}

fn checked_product2(left: f64, right: f64) -> Option<f64> {
    let product = left * right;
    (left.is_finite() && right.is_finite() && product.is_finite()).then_some(product)
}

fn checked_div2(left: f64, right: f64) -> Option<f64> {
    let quotient = left / right;
    (left.is_finite() && right.is_finite() && right != 0.0 && quotient.is_finite())
        .then_some(quotient)
}

fn checked_intercept_beta_se(sigma: f64, gene: usize) -> Result<f64, DeseqError> {
    let se = sigma.sqrt();
    let Some(value) = checked_product2(std::f64::consts::LOG2_E, se) else {
        return Err(DeseqError::InvalidDispersion {
            reason: format!("gene {gene} has non-finite intercept beta standard error"),
        });
    };
    Ok(value)
}

fn checked_intercept_beta_covariance(sigma: f64, gene: usize) -> Result<f64, DeseqError> {
    let log2_e = std::f64::consts::LOG2_E;
    let Some(log2_e_squared) = checked_product2(log2_e, log2_e) else {
        return Err(DeseqError::InvalidDispersion {
            reason: "non-finite log2 covariance scaling factor".to_string(),
        });
    };
    let Some(value) = checked_product2(log2_e_squared, sigma) else {
        return Err(DeseqError::InvalidDispersion {
            reason: format!("gene {gene} has non-finite intercept beta covariance"),
        });
    };
    Ok(value)
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
        let quantile = stable_interpolate(values[lower - 1], values[lower], fraction);
        if quantile.is_finite() {
            Ok(quantile)
        } else {
            Err(DeseqError::InvalidOptions {
                reason: "beta prior variance quantile produced non-finite quantile".to_string(),
            })
        }
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

    let Some(weight_sum) = checked_sum(pairs.iter().map(|(_, weight)| *weight)) else {
        return Err(DeseqError::InvalidOptions {
            reason: "weighted beta prior variance quantile needs finite total weight".to_string(),
        });
    };
    if weight_sum <= 0.0 {
        return Err(DeseqError::InvalidOptions {
            reason: "weighted beta prior variance quantile needs positive finite total weight"
                .to_string(),
        });
    }

    let Some(norm_scale) = checked_div2(pairs.len() as f64, weight_sum) else {
        return Err(DeseqError::InvalidOptions {
            reason: "weighted beta prior variance quantile produced non-finite normalization scale"
                .to_string(),
        });
    };
    let mut unique = Vec::<(f64, f64)>::with_capacity(pairs.len());
    for (value, weight) in pairs {
        let Some(normalized_weight) = checked_product2(weight, norm_scale) else {
            return Err(DeseqError::InvalidOptions {
                reason:
                    "weighted beta prior variance quantile produced non-finite normalized weight"
                        .to_string(),
            });
        };
        if normalized_weight <= 0.0 {
            return Err(DeseqError::InvalidOptions {
                reason: "weighted beta prior variance quantile needs positive normalized weight"
                    .to_string(),
            });
        }
        if let Some((last_value, last_weight)) = unique.last_mut() {
            if *last_value == value {
                let Some(next_weight) = checked_sum2(*last_weight, normalized_weight) else {
                    return Err(DeseqError::InvalidOptions {
                        reason: "weighted beta prior variance quantile produced non-finite merged weight"
                            .to_string(),
                    });
                };
                *last_weight = next_weight;
                continue;
            }
        }
        unique.push((value, normalized_weight));
    }

    if unique.len() == 1 {
        return Ok(unique[0].0);
    }

    let Some(n) = checked_sum(unique.iter().map(|(_, weight)| *weight)) else {
        return Err(DeseqError::InvalidOptions {
            reason:
                "weighted beta prior variance quantile produced non-finite total normalized weight"
                    .to_string(),
        });
    };
    let order = 1.0 + (n - 1.0) * probability;
    if !order.is_finite() {
        return Err(DeseqError::InvalidOptions {
            reason: "weighted beta prior variance quantile produced non-finite order".to_string(),
        });
    }
    let low = order.floor().max(1.0);
    let high = (low + 1.0).min(n);
    let fraction = order.fract();

    let mut cumulative_weights = Vec::with_capacity(unique.len());
    let mut cumulative = 0.0;
    for (_, weight) in &unique {
        let Some(next_cumulative) = checked_sum2(cumulative, *weight) else {
            return Err(DeseqError::InvalidOptions {
                reason:
                    "weighted beta prior variance quantile produced non-finite cumulative weight"
                        .to_string(),
            });
        };
        cumulative = next_cumulative;
        cumulative_weights.push(cumulative);
    }

    let low_quantile = weighted_order_statistic(&unique, &cumulative_weights, low);
    let high_quantile = weighted_order_statistic(&unique, &cumulative_weights, high);
    let quantile = stable_interpolate(low_quantile, high_quantile, fraction);
    if !quantile.is_finite() {
        return Err(DeseqError::InvalidOptions {
            reason: "weighted beta prior variance quantile produced non-finite quantile"
                .to_string(),
        });
    }
    Ok(quantile)
}

fn stable_interpolate(left: f64, right: f64, fraction: f64) -> f64 {
    if left == right {
        return left;
    }
    let Some(delta) = checked_sum2(right, -left) else {
        return f64::NAN;
    };
    let Some(offset) = checked_product2(fraction, delta) else {
        return f64::NAN;
    };
    checked_sum2(left, offset).unwrap_or(f64::NAN)
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

#[cfg(test)]
mod tests {
    use super::{
        checked_div2, checked_mean_dispersion_working_weight, checked_scaled_sum,
        intercept_working_weights, stable_interpolate, weighted_mean,
    };
    use crate::errors::DeseqError;

    #[test]
    fn stable_interpolate_preserves_equal_extreme_endpoints() {
        let interpolated = stable_interpolate(f64::MAX, f64::MAX, 0.5);

        assert_eq!(interpolated, f64::MAX);
    }

    #[test]
    fn stable_interpolate_rejects_overflowed_delta() {
        let interpolated = stable_interpolate(-f64::MAX, f64::MAX, 0.5);

        assert!(interpolated.is_nan());
    }

    #[test]
    fn quantile_type7_rejects_overflowed_interpolation() {
        let err = super::quantile_type7(vec![-f64::MAX, f64::MAX], 0.5).unwrap_err();

        assert!(matches!(
            err,
            crate::errors::DeseqError::InvalidOptions { reason }
                if reason == "beta prior variance quantile produced non-finite quantile"
        ));
    }

    #[test]
    fn checked_scaled_sum_rejects_overflowed_rescale() {
        assert_eq!(checked_scaled_sum([f64::MAX, f64::MAX]), None);
    }

    #[test]
    fn checked_div2_rejects_zero_and_nonfinite_inputs() {
        assert_eq!(checked_div2(1.0, 0.0), None);
        assert_eq!(checked_div2(f64::NAN, 1.0), None);
        assert_eq!(checked_div2(4.0, 2.0), Some(2.0));
    }

    #[test]
    fn mean_dispersion_weight_rejects_nonfinite_arithmetic() {
        let err = checked_mean_dispersion_working_weight(
            f64::MIN_POSITIVE,
            f64::MAX,
            2,
            "test mean-dispersion weight",
        )
        .unwrap_err();

        assert!(matches!(
            err,
            DeseqError::NonFiniteValue { context, index, .. }
                if context == "test mean-dispersion weight" && index == Some(2)
        ));
    }

    #[test]
    fn intercept_working_weights_reject_nonfinite_weight_scaling() {
        let err =
            intercept_working_weights(&[f64::MAX], f64::MIN_POSITIVE, Some(&[10.0])).unwrap_err();

        assert!(matches!(
            err,
            DeseqError::NonFiniteValue { context, index, .. }
                if context == "intercept working weight" && index == Some(0)
        ));
    }

    #[test]
    fn weighted_mean_rejects_nonfinite_products() {
        let err = weighted_mean(&[f64::MAX], &[2.0], 0).unwrap_err();

        assert!(err
            .to_string()
            .contains("gene 0 has non-finite weighted normalized mean"));
    }
}
