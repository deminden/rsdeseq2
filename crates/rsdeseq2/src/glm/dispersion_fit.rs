use crate::errors::DeseqError;

/// Placeholder for future fixed-dispersion GLM fitting.
pub fn fit_with_dispersion() -> Result<(), DeseqError> {
    Err(DeseqError::UnsupportedFeature {
        feature: "fixed-dispersion GLM fitting".to_string(),
    })
}
