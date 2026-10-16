// Fixed and estimated beta-prior fitting workflows.

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
