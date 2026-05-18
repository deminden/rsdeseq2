use crate::errors::DeseqError;

/// Placeholder for future variance-stabilizing transformation.
pub fn vst() -> Result<(), DeseqError> {
    Err(DeseqError::UnsupportedFeature {
        feature: "variance-stabilizing transformation".to_string(),
    })
}
