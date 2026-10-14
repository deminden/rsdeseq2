use thiserror::Error;

/// Error type used by the `rsdeseq2` Rust core.
#[derive(Debug, Error)]
pub enum DeseqError {
    /// Matrix or vector dimensions do not match the requested shape.
    #[error("invalid dimensions for {context}: expected {expected}, got {actual}")]
    InvalidDimensions {
        context: String,
        expected: usize,
        actual: usize,
    },

    /// Count data failed validation.
    #[error("invalid counts: {reason}")]
    InvalidCounts { reason: String },

    /// Size factors failed validation or could not be estimated.
    #[error("invalid size factors: {reason}")]
    InvalidSizeFactors { reason: String },

    /// Dispersion values failed validation.
    #[error("invalid dispersion: {reason}")]
    InvalidDispersion { reason: String },

    /// Analysis options failed validation.
    #[error("invalid options: {reason}")]
    InvalidOptions { reason: String },

    /// No rows can contribute to median-ratio size-factor estimation.
    #[error("no usable genes for size-factor estimation")]
    NoUsableGenesForSizeFactors,

    /// A value was expected to be finite.
    #[error("non-finite value in {context} at index {index:?}: {value}")]
    NonFiniteValue {
        context: String,
        index: Option<usize>,
        value: f64,
    },

    /// The requested behavior is not implemented yet.
    #[error("unsupported feature: {feature}")]
    UnsupportedFeature { feature: String },

    /// I/O error.
    #[error(transparent)]
    Io(#[from] std::io::Error),

    /// CSV/TSV parser error.
    #[error(transparent)]
    Csv(#[from] csv::Error),

    /// Integer parsing error in an input table.
    #[error("failed to parse integer value '{value}' in {context}")]
    ParseInt { context: String, value: String },

    /// Floating-point parsing error in an input table.
    #[error("failed to parse floating-point value '{value}' in {context}")]
    ParseFloat { context: String, value: String },
}

pub(crate) fn invalid_dimensions(
    context: impl Into<String>,
    expected: usize,
    actual: usize,
) -> DeseqError {
    DeseqError::InvalidDimensions {
        context: context.into(),
        expected,
        actual,
    }
}
