use crate::errors::DeseqError;

/// Placeholder for future negative-binomial distribution helpers.
pub fn negative_binomial_helpers() -> Result<(), DeseqError> {
    Err(DeseqError::UnsupportedFeature {
        feature: "negative-binomial distribution helpers".to_string(),
    })
}
