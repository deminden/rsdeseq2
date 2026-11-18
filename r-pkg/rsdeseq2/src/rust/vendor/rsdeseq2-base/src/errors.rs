//! The error type shared by every rsdeseq2 stage.

use thiserror::Error;

/// Error type used by the `rsdeseq2` Rust core.
///
/// The enum is `#[non_exhaustive]`: match on the variants you handle and keep a
/// wildcard arm. [`DeseqError::kind`] gives a stable coarse classification.
#[derive(Debug, Error)]
#[non_exhaustive]
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

    /// A value was expected to be finite, or finite and within a range
    /// (for example positive). The message names which case occurred.
    #[error(
        "{} value in {context} at index {index:?}: {value}",
        if value.is_finite() { "out-of-range" } else { "non-finite" }
    )]
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

/// Coarse classification of a [`DeseqError`], stable across message changes.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
#[non_exhaustive]
pub enum ErrorKind {
    /// Matrix or vector shapes disagree.
    Shape,
    /// Count data are invalid.
    Counts,
    /// Size or normalization factors are invalid or cannot be estimated.
    Normalization,
    /// Dispersion inputs or estimates are invalid.
    Dispersion,
    /// Options, designs, contrasts or other analysis inputs are invalid.
    Options,
    /// A computation produced or received a non-finite value.
    Numeric,
    /// The request is outside the supported contract.
    Unsupported,
    /// Reading or writing files failed.
    Io,
    /// An input table could not be parsed.
    Parse,
}

impl DeseqError {
    /// Coarse classification of this error.
    pub fn kind(&self) -> ErrorKind {
        match self {
            Self::InvalidDimensions { .. } => ErrorKind::Shape,
            Self::InvalidCounts { .. } => ErrorKind::Counts,
            Self::InvalidSizeFactors { .. } | Self::NoUsableGenesForSizeFactors => {
                ErrorKind::Normalization
            }
            Self::InvalidDispersion { .. } => ErrorKind::Dispersion,
            Self::InvalidOptions { .. } => ErrorKind::Options,
            Self::NonFiniteValue { .. } => ErrorKind::Numeric,
            Self::UnsupportedFeature { .. } => ErrorKind::Unsupported,
            Self::Io(_) => ErrorKind::Io,
            Self::Csv(_) | Self::ParseInt { .. } | Self::ParseFloat { .. } => ErrorKind::Parse,
        }
    }
}

/// `InvalidOptions` with the reason as given.
#[doc(hidden)]
pub fn invalid_options(reason: impl Into<String>) -> DeseqError {
    DeseqError::InvalidOptions {
        reason: reason.into(),
    }
}

/// `InvalidOptions` whose reason is prefixed with the component that raised
/// it: `"{scope}: {reason}"`.
#[doc(hidden)]
pub fn scoped_invalid_options(scope: &str, reason: impl std::fmt::Display) -> DeseqError {
    DeseqError::InvalidOptions {
        reason: format!("{scope}: {reason}"),
    }
}

/// An internal consistency check failed: a state that earlier validation
/// rules out was reached. Reported as `InvalidOptions` so callers never panic.
#[doc(hidden)]
pub fn invariant(what: &str) -> DeseqError {
    DeseqError::InvalidOptions {
        reason: format!("internal invariant violated: {what}"),
    }
}

/// `InvalidDispersion` with the reason as given.
#[doc(hidden)]
pub fn invalid_dispersion(reason: impl Into<String>) -> DeseqError {
    DeseqError::InvalidDispersion {
        reason: reason.into(),
    }
}

#[doc(hidden)]
pub fn invalid_dimensions(
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn kinds_classify_every_variant() {
        assert_eq!(invalid_dimensions("x", 1, 2).kind(), ErrorKind::Shape);
        assert_eq!(invalid_options("x").kind(), ErrorKind::Options);
        assert_eq!(invalid_dispersion("x").kind(), ErrorKind::Dispersion);
        assert_eq!(invariant("x").kind(), ErrorKind::Options);
        assert_eq!(
            DeseqError::NoUsableGenesForSizeFactors.kind(),
            ErrorKind::Normalization
        );
        let io = DeseqError::from(std::io::Error::other("x"));
        assert_eq!(io.kind(), ErrorKind::Io);
    }

    #[test]
    fn finite_out_of_range_values_are_not_called_non_finite() {
        let negative = DeseqError::NonFiniteValue {
            context: "dispersion".into(),
            index: Some(1),
            value: -0.1,
        };
        assert_eq!(
            negative.to_string(),
            "out-of-range value in dispersion at index Some(1): -0.1"
        );
        let nan = DeseqError::NonFiniteValue {
            context: "dispersion".into(),
            index: Some(1),
            value: f64::NAN,
        };
        assert_eq!(
            nan.to_string(),
            "non-finite value in dispersion at index Some(1): NaN"
        );
    }
}
