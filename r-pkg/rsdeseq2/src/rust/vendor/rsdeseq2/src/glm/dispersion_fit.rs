use crate::core::CountMatrix;
use crate::design::DesignMatrix;
use crate::errors::DeseqError;
use crate::glm::NbinomGlmFit;
use crate::glm::irls::{IrlsOptions, fit_fixed_dispersion_irls};

/// Fit a supplied-dispersion negative-binomial GLM.
///
/// This is a DESeq2-shaped convenience wrapper around the implemented
/// fixed-dispersion IRLS path. Use the more specific `fit_fixed_dispersion_*`
/// helpers when normalization factors or observation weights are needed.
pub fn fit_with_dispersion(
    counts: &CountMatrix,
    design: &DesignMatrix,
    size_factors: &[f64],
    dispersions: &[f64],
    options: IrlsOptions,
) -> Result<NbinomGlmFit, DeseqError> {
    fit_fixed_dispersion_irls(counts, design, size_factors, dispersions, options)
}
