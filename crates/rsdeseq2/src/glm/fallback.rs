use crate::errors::DeseqError;

/// Placeholder for future fallback fitting paths.
pub fn fallback_fit() -> Result<(), DeseqError> {
    Err(DeseqError::UnsupportedFeature {
        feature: "fallback GLM fitting".to_string(),
    })
}
