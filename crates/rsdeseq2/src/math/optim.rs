use crate::errors::DeseqError;

/// Placeholder for future optimization routines.
pub fn optimization_routines() -> Result<(), DeseqError> {
    Err(DeseqError::UnsupportedFeature {
        feature: "optimization routines".to_string(),
    })
}
