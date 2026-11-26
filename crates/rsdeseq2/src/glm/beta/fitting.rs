use super::*;
use crate::glm::irls::fit_selected_irls;
// Fixed and estimated beta-prior fitting workflows.

/// Counts, normalization offsets, dispersions and optional observation
/// weights for the beta-prior GLM fits.
///
/// Build with [`GlmInput::new`] and [`GlmInput::weights`].
#[derive(Clone, Copy, Debug)]
#[non_exhaustive]
pub struct GlmInput<'a> {
    /// Raw counts.
    pub counts: &'a CountMatrix,
    /// Size factors or normalization factors.
    pub offsets: Offsets<'a>,
    /// One dispersion per gene.
    pub dispersions: &'a [f64],
    /// Optional normalized observation weights.
    pub weights: Option<&'a RowMajorMatrix<f64>>,
}

impl<'a> GlmInput<'a> {
    /// Unweighted input.
    pub fn new(counts: &'a CountMatrix, offsets: Offsets<'a>, dispersions: &'a [f64]) -> Self {
        Self {
            counts,
            offsets,
            dispersions,
            weights: None,
        }
    }

    /// Add observation weights.
    pub fn weights(mut self, weights: &'a RowMajorMatrix<f64>) -> Self {
        self.weights = Some(weights);
        self
    }
}

/// Refit a fixed-dispersion GLM with a supplied DESeq2-style beta prior
/// variance (one value per coefficient).
pub fn fit_glms_with_beta_prior(
    input: GlmInput<'_>,
    design: &DesignMatrix,
    beta_prior_variance: &[f64],
    options: IrlsOptions,
) -> Result<NbinomGlmFit, DeseqError> {
    fit_glms_with_beta_prior_rows(input, design, beta_prior_variance, options, None)
}

pub(crate) fn fit_glms_with_beta_prior_rows(
    input: GlmInput<'_>,
    design: &DesignMatrix,
    beta_prior_variance: &[f64],
    options: IrlsOptions,
    genes: Option<&[usize]>,
) -> Result<NbinomGlmFit, DeseqError> {
    // Preserve the validation order of the weighted size-factor fit.
    if let (Offsets::PerSample(size_factors), Some(_)) = (input.offsets, input.weights) {
        validate_size_factor_offsets(input.counts, size_factors)?;
    }
    let options = options_with_beta_prior_variance(design, beta_prior_variance, options)?;
    match genes {
        Some(genes) => fit_selected_irls(input, design, options, genes),
        None => fit_fixed_dispersion_glm(input, design, options),
    }
}

/// Run an MLE fixed-dispersion fit, estimate the beta prior variance
/// (`estimateBetaPriorVar`), then refit with that prior.
pub fn fit_glms_with_estimated_beta_prior(
    input: GlmInput<'_>,
    design: &DesignMatrix,
    base_mean: &[f64],
    disp_fit: &[f64],
    options: BetaPriorRefitOptions,
) -> Result<BetaPriorGlmFit, DeseqError> {
    let fit_options = options.fit_options;
    let mle_fit = match (input.offsets, input.weights) {
        // Unweighted size factors keep DESeq2's intercept-only shortcut.
        (Offsets::PerSample(size_factors), None) => fit_irls(
            input.counts,
            design,
            size_factors,
            input.dispersions,
            fit_options.clone(),
        )?,
        _ => fit_fixed_dispersion_glm(input, design, fit_options.clone())?,
    };
    let beta_prior_variance = estimate_beta_prior_variance(
        &mle_fit.beta,
        base_mean,
        disp_fit,
        design.coefficient_names(),
        options.variance_options,
    )?;
    let prior_fit = fit_glms_with_beta_prior(input, design, &beta_prior_variance, fit_options)?;
    Ok(BetaPriorGlmFit {
        mle_fit,
        prior_fit,
        beta_prior_variance,
    })
}

/// Run an expanded-design MLE fit, estimate the beta prior variance, refit,
/// and collapse to the standard design (DESeq2's expanded-model beta prior).
pub fn fit_expanded_glms_with_estimated_beta_prior(
    input: GlmInput<'_>,
    design: ExpandedModelBetaPriorDesignInput<'_>,
    base_mean: &[f64],
    disp_fit: &[f64],
    options: BetaPriorRefitOptions,
) -> Result<ExpandedModelBetaPriorGlmFit, DeseqError> {
    if let Offsets::PerSample(size_factors) = input.offsets {
        validate_size_factor_offsets(input.counts, size_factors)?;
    }
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
    let expanded_mle_fit =
        fit_fixed_dispersion_glm(input, design.expanded_design, fit_options.clone())?;
    let beta_prior_variance = estimate_beta_prior_variance(
        &expanded_mle_fit.beta,
        base_mean,
        disp_fit,
        design.expanded_design.coefficient_names(),
        options.variance_options,
    )?;
    let expanded_prior_fit = fit_glms_with_beta_prior(
        input,
        design.expanded_design,
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

pub(crate) fn fit_prior_for_output(
    input: GlmInput<'_>,
    design: &DesignMatrix,
    beta_prior_variance: &[f64],
    options: IrlsOptions,
    genes: Option<&[usize]>,
    offset_scales: Option<&[f64]>,
    requirements: crate::glm::irls::OutputRequirements,
) -> Result<crate::glm::irls::IrlsOutput, DeseqError> {
    use crate::glm::irls::{
        fit_irls_with_offset_scales, fit_nf_irls_for_output, fit_selected_irls_output,
        fit_sf_irls_for_output,
    };
    let options = options_with_beta_prior_variance(design, beta_prior_variance, options)?;
    match (genes, offset_scales) {
        (_, Some(scales)) => {
            fit_irls_with_offset_scales(input, design, options, genes, scales, requirements)
        }
        (Some(genes), None) => {
            fit_selected_irls_output(input, design, options, genes, requirements)
        }
        (None, None) => match input.offsets {
            Offsets::PerSample(factors) => fit_sf_irls_for_output(
                input.counts,
                design,
                factors,
                input.dispersions,
                input.weights,
                options,
                requirements,
            ),
            Offsets::PerCell(factors) => fit_nf_irls_for_output(
                input.counts,
                design,
                factors,
                input.dispersions,
                input.weights,
                options,
                requirements,
            ),
        },
    }
}
