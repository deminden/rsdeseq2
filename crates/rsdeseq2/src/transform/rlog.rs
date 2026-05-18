use crate::errors::DeseqError;

/// Placeholder for future regularized-log transformation.
pub fn rlog() -> Result<(), DeseqError> {
    Err(DeseqError::UnsupportedFeature {
        feature: "regularized-log transformation".to_string(),
    })
}
