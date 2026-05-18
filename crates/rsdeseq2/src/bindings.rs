use crate::errors::DeseqError;

/// Placeholder for future language bindings.
pub fn language_bindings() -> Result<(), DeseqError> {
    Err(DeseqError::UnsupportedFeature {
        feature: "language bindings".to_string(),
    })
}
