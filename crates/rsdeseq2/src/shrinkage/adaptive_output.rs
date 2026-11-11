//! Selection of fields consumed by the DESeq-style adaptive shrinkage wrapper.
use super::{
    AdaptiveShrinkageQuery, LocationShrinkage, LocationShrinkageOptions,
    adaptive_result::adaptive_shrinkage_impl, location_batch::shrink_location_batch_with_outputs,
};
use crate::DeseqError;

/// Original ashr output selection, projected onto fields used by lfcShrink.
/// Additional original summary selections preserve their evaluation diagnostics.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct AdaptiveOutputSelection {
    pub posterior_mean: bool,
    pub posterior_sd: bool,
    pub svalues: bool,
    pub fitted_prior: bool,
    /// Q-value selection also controls its missing-value recycling warning.
    pub qvalues: bool,
    pub negative_probability: bool,
    pub positive_probability: bool,
    pub local_false_sign_rate: bool,
    pub local_false_discovery_rate: bool,
    pub log_likelihood: bool,
    pub log_likelihood_ratio: bool,
    pub flash_data: bool,
}
impl AdaptiveOutputSelection {
    /// Numeric levels use inequalities, including fractional and infinite levels.
    /// Levels four and five have the original special field selections.
    pub fn level(level: f64) -> Result<Self, DeseqError> {
        if level.is_nan() {
            return Err(invalid("outputlevel must not be missing"));
        }
        Ok(Self {
            posterior_mean: level > 0.0 && level != 5.0,
            posterior_sd: level > 0.0 && level != 5.0,
            svalues: level > 1.0 && level != 4.0 && level != 5.0,
            fitted_prior: level != 5.0,
            qvalues: level > 1.0 && level != 4.0 && level != 5.0,
            negative_probability: level > 1.0 && level != 4.0 && level != 5.0,
            positive_probability: level > 1.0 && level != 4.0 && level != 5.0,
            local_false_sign_rate: level > 1.0 && level != 4.0 && level != 5.0,
            local_false_discovery_rate: level > 1.0 && level != 4.0 && level != 5.0,
            log_likelihood: level > 0.0 && level != 4.0 && level != 5.0,
            log_likelihood_ratio: level > 0.0 && level != 4.0 && level != 5.0,
            flash_data: level == 4.0 || level == 5.0,
        })
    }
    /// Named output subsets use exact names; duplicates/unknown names are inert.
    pub fn names<'a>(names: impl IntoIterator<Item = &'a str>) -> Self {
        let mut selected = Self::default();
        for name in names {
            match name {
                "PosteriorMean" => selected.posterior_mean = true,
                "PosteriorSD" => selected.posterior_sd = true,
                "svalue" => selected.svalues = true,
                "fitted_g" => selected.fitted_prior = true,
                "qvalue" => selected.qvalues = true,
                "NegativeProb" => selected.negative_probability = true,
                "PositiveProb" => selected.positive_probability = true,
                "lfsr" => selected.local_false_sign_rate = true,
                "lfdr" => selected.local_false_discovery_rate = true,
                "loglik" => selected.log_likelihood = true,
                "logLR" => selected.log_likelihood_ratio = true,
                "flash_data" => selected.flash_data = true,
                _ => {}
            }
        }
        selected
    }
}
/// Selected numerical results. Absent fields remain absent in both frontends.
#[derive(Clone, Debug)]
pub struct SelectedAdaptiveShrinkageResult {
    pub fit: LocationShrinkage,
    pub posterior_mean: Option<Vec<f64>>,
    pub posterior_sd: Option<Vec<f64>>,
    pub svalues: Option<Vec<f64>>,
    pub fitted_prior: bool,
}
fn invalid(reason: &str) -> DeseqError {
    DeseqError::InvalidOptions {
        reason: format!("adaptive shrinkage output: {reason}"),
    }
}
/// Fit and select the wrapper's numerical fields in the shared Rust engine.
/// Threshold queries require a fitted prior and posterior means even when the
/// ordinary s-value output is deselected. Without any selected summary there is
/// no assignment/recycling into the query population.
pub fn adaptive_shrinkage_with_output(
    estimates: &[f64],
    standard_errors: &[f64],
    options: &LocationShrinkageOptions,
    query: AdaptiveShrinkageQuery<'_>,
    output: AdaptiveOutputSelection,
) -> Result<SelectedAdaptiveShrinkageResult, DeseqError> {
    if !query.lfc_threshold.is_finite() || query.lfc_threshold < 0.0 {
        return Err(invalid("LFC threshold must be finite and nonnegative"));
    }
    if query.lfc_threshold > 0.0 && (!output.posterior_mean || !output.fitted_prior) {
        return Err(invalid(
            "threshold queries require PosteriorMean and fitted_g outputs",
        ));
    }
    let use_svalues = output.svalues && query.svalues;
    if !output.posterior_mean && !output.posterior_sd && !use_svalues {
        return Ok(SelectedAdaptiveShrinkageResult {
            fit: shrink_location_batch_with_outputs(estimates, standard_errors, options, output)?,
            posterior_mean: None,
            posterior_sd: None,
            svalues: None,
            fitted_prior: output.fitted_prior,
        });
    }
    let result = adaptive_shrinkage_impl(
        estimates,
        standard_errors,
        options,
        AdaptiveShrinkageQuery {
            svalues: use_svalues,
            ..query
        },
        output,
    )?;
    Ok(SelectedAdaptiveShrinkageResult {
        fit: result.fit,
        posterior_mean: output.posterior_mean.then_some(result.posterior_mean),
        posterior_sd: output.posterior_sd.then_some(result.posterior_sd),
        svalues: result.svalues,
        fitted_prior: output.fitted_prior,
    })
}
