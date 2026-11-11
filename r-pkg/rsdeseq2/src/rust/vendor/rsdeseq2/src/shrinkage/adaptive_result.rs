//! Shared result population and threshold policy for the ashr-style interface.
use super::{
    AdaptiveOutputSelection, LocationCdfOptions, LocationShrinkage, LocationShrinkageOptions,
    location_batch::shrink_location_batch_with_outputs, location_cdf_batch, threshold_svalues,
};
use crate::DeseqError;
/// Original result population, retained separately from forwarded fitting data.
#[derive(Clone, Copy, Debug)]
pub struct AdaptiveShrinkageQuery<'a> {
    pub estimates: &'a [f64],
    pub standard_errors: &'a [f64],
    pub lfc_threshold: f64,
    pub svalues: bool,
}
/// Numerical fields ready for either frontend's result-table presentation.
#[derive(Clone, Debug)]
pub struct AdaptiveShrinkageResult {
    pub fit: LocationShrinkage,
    pub posterior_mean: Vec<f64>,
    pub posterior_sd: Vec<f64>,
    pub svalues: Option<Vec<f64>>,
}
fn invalid(reason: &str) -> DeseqError {
    DeseqError::InvalidOptions {
        reason: format!("adaptive shrinkage result: {reason}"),
    }
}
/// Fit on the supplied data and present the original result population. A
/// threshold query always uses the original estimates/SEs and normal alpha=0
/// likelihood, retaining the established wrapper policy under fit overrides.
pub fn adaptive_shrinkage(
    estimates: &[f64],
    standard_errors: &[f64],
    options: &LocationShrinkageOptions,
    query: AdaptiveShrinkageQuery<'_>,
) -> Result<AdaptiveShrinkageResult, DeseqError> {
    adaptive_shrinkage_impl(
        estimates,
        standard_errors,
        options,
        query,
        AdaptiveOutputSelection::level(2.0)?,
    )
}
pub(super) fn adaptive_shrinkage_impl(
    estimates: &[f64],
    standard_errors: &[f64],
    options: &LocationShrinkageOptions,
    query: AdaptiveShrinkageQuery<'_>,
    output: AdaptiveOutputSelection,
) -> Result<AdaptiveShrinkageResult, DeseqError> {
    if !query.lfc_threshold.is_finite() || query.lfc_threshold < 0.0 {
        return Err(invalid("LFC threshold must be finite and nonnegative"));
    }
    let fit = shrink_location_batch_with_outputs(estimates, standard_errors, options, output)?;
    let n = query.estimates.len();
    let fitted = fit.posterior_mean.len();
    if n == 0 || fitted == 0 || n < fitted || !n.is_multiple_of(fitted) {
        return Err(invalid(
            "posterior length cannot be assigned to the original result population",
        ));
    }
    let align = |v: &[f64]| (0..n).map(|i| v[i % v.len()]).collect::<Vec<_>>();
    let posterior_mean = align(&fit.posterior_mean);
    let posterior_sd = align(&fit.posterior_sd);
    let svalues = if query.lfc_threshold > 0.0 {
        let cdf = location_cdf_batch(
            &fit.fitted_grid,
            &fit.fitted_weights,
            query.estimates,
            query.standard_errors,
            &LocationCdfOptions {
                thresholds: vec![query.lfc_threshold, -query.lfc_threshold],
                ..Default::default()
            },
        )?;
        let positive = cdf
            .as_chunks::<2>()
            .0
            .iter()
            .map(|r| r[0])
            .collect::<Vec<_>>();
        let negative = cdf
            .as_chunks::<2>()
            .0
            .iter()
            .map(|r| r[1])
            .collect::<Vec<_>>();
        Some(threshold_svalues(&posterior_mean, &positive, &negative)?)
    } else if query.svalues {
        Some(align(&fit.svalues))
    } else {
        None
    };
    Ok(AdaptiveShrinkageResult {
        fit,
        posterior_mean,
        posterior_sd,
        svalues,
    })
}
