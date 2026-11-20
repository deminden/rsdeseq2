//! The error type of the rcompat crates: input validation, numeric domain
//! and shape errors, with stable messages and a coarse [`ErrorKind`].

/// Error type of the rcompat crates (re-exported by `rsdeseq2` as
/// `Error`).
///
/// The enum is `#[non_exhaustive]`: match on the variants you handle and keep a
/// wildcard arm. [`Error::kind`] gives a stable coarse classification.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum Error {
    /// Matrix or vector dimensions do not match the requested shape.
    #[error("invalid dimensions for {context}: expected {expected}, got {actual}")]
    InvalidDimensions {
        /// What was being checked.
        context: String,
        /// The required length or count.
        expected: usize,
        /// The length or count supplied.
        actual: usize,
    },

    /// Count data failed validation.
    #[error("invalid counts: {reason}")]
    InvalidCounts {
        /// What is wrong with the counts.
        reason: String,
    },

    /// A fixed-dispersion GLM received an all-zero count row.
    #[error("invalid counts: gene {gene} is all zero; DESeq2 GLM fitting excludes allZero rows")]
    AllZeroGene {
        /// Row in the fit, after any row selection.
        gene: usize,
        /// Row in the count matrix supplied to the fit.
        source_gene: usize,
    },

    /// A gene-level calculation produced a value outside its domain.
    #[error("{}: {reason}", domain_prefix(*domain))]
    InvalidGeneValue {
        /// The existing classification of the failing stage.
        domain: ErrorKind,
        /// The calculation that failed.
        operation: &'static str,
        /// Row in the input to the stage.
        gene: usize,
        /// Sample involved in the calculation, when known.
        sample: Option<usize>,
        /// The value the stage rejected.
        value: f64,
        /// The diagnostic presented to the caller.
        reason: String,
    },

    /// Size factors failed validation or could not be estimated.
    #[error("invalid size factors: {reason}")]
    InvalidSizeFactors {
        /// What is wrong with the size factors.
        reason: String,
    },

    /// Dispersion values failed validation.
    #[error("invalid dispersion: {reason}")]
    InvalidDispersion {
        /// What is wrong with the dispersions.
        reason: String,
    },

    /// Analysis options failed validation.
    #[error("invalid options: {reason}")]
    InvalidOptions {
        /// What is wrong with the options.
        reason: String,
    },

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
        /// Where the value occurred.
        context: String,
        /// Its index (gene, sample or element), when known.
        index: Option<usize>,
        /// The offending value.
        value: f64,
    },

    /// The requested behavior is not implemented yet.
    #[error("unsupported feature: {feature}")]
    UnsupportedFeature {
        /// The unsupported request.
        feature: String,
    },

    /// I/O error.
    #[error(transparent)]
    Io(#[from] std::io::Error),

    /// CSV/TSV parser error (feature `csv`).
    #[cfg(feature = "csv")]
    #[error(transparent)]
    Csv(#[from] csv::Error),

    /// Integer parsing error in an input table.
    #[error("failed to parse integer value '{value}' in {context}")]
    ParseInt {
        /// Where the value occurred.
        context: String,
        /// The text that failed to parse.
        value: String,
    },

    /// Floating-point parsing error in an input table.
    #[error("failed to parse floating-point value '{value}' in {context}")]
    ParseFloat {
        /// Where the value occurred.
        context: String,
        /// The text that failed to parse.
        value: String,
    },

    /// A named parameter is outside its domain, for example a tolerance that
    /// must be finite and positive.
    #[error("{}: {parameter} must be {requirement}", domain_prefix(*domain))]
    InvalidParameter {
        /// Which input family the parameter belongs to: `Options`,
        /// `Dispersion`, `Normalization` or `Counts`. It is also
        /// [`Error::kind`].
        domain: ErrorKind,
        /// The parameter as the message names it.
        parameter: String,
        /// The condition it failed.
        requirement: Requirement,
    },

    /// An operation was asked for before something it needs: an earlier
    /// stage's output, an input or a companion option.
    #[error("{}: {operation} requires {prerequisite}", domain_prefix(*domain))]
    MissingPrerequisite {
        /// Which input family the operation belongs to: `Options`,
        /// `Dispersion`, `Normalization` or `Counts`. It is also
        /// [`Error::kind`].
        domain: ErrorKind,
        /// What was asked for, for example "Wald fitting".
        operation: String,
        /// What it needs, for example "MAP dispersions".
        prerequisite: String,
    },

    /// A design formula could not be parsed or evaluated.
    #[error("invalid options: {reason}")]
    InvalidFormula {
        /// What is wrong with the formula.
        reason: String,
    },

    /// A requested coefficient is not among the fit's coefficient names.
    #[error("invalid options: {reason}")]
    CoefficientNotFound {
        /// The full message.
        reason: String,
        /// The names that were looked for, in search order.
        candidates: Vec<String>,
    },

    /// A coefficient request admits more than one match.
    #[error("invalid options: {reason}")]
    AmbiguousCoefficient {
        /// The full message.
        reason: String,
        /// The coefficient names or levels requested, in search order.
        candidates: Vec<String>,
        /// The matching coefficients in model-matrix column order.
        matches: Vec<String>,
    },

    /// The model matrix has dependent or zero columns.
    #[error(
        "invalid options: {context} model matrix is not full rank: rank {rank}, columns {columns}{}",
        zero_columns_suffix(zero_columns)
    )]
    DesignRankDeficient {
        /// The model or fitting stage that requires full rank.
        context: String,
        /// The computed rank.
        rank: usize,
        /// The number of columns in the design.
        columns: usize,
        /// Names or zero-based positions of columns found to be zero.
        zero_columns: Vec<String>,
    },

    /// The reduced LRT design contains directions outside the full design.
    #[error("invalid options: LRT reduced design must be nested in the full design column space")]
    LrtDesignNotNested {
        /// Rank of the full design.
        full_rank: usize,
        /// Rank of the full and reduced designs combined.
        combined_rank: usize,
    },

    /// An internal consistency check failed: a state that earlier validation
    /// rules out was reached. This is a bug; please report it.
    #[error("invalid options: internal invariant violated: {what}")]
    Internal {
        /// The violated condition.
        what: String,
    },
}

fn zero_columns_suffix(columns: &[String]) -> String {
    if columns.is_empty() {
        String::new()
    } else {
        format!("; zero columns: {}", columns.join(", "))
    }
}

/// The condition an [`Error::InvalidParameter`] failed. `Display` gives
/// the phrase used in the message.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
#[non_exhaustive]
pub enum Requirement {
    /// "finite"
    Finite,
    /// "positive"
    Positive,
    /// "non-negative"
    NonNegative,
    /// "finite and positive"
    FiniteAndPositive,
    /// "finite and non-negative"
    FiniteAndNonNegative,
    /// "finite and within [0, 1]"
    FiniteInUnitInterval,
    /// "non-empty"
    NonEmpty,
}

impl Requirement {
    /// The phrase used in messages.
    pub fn phrase(self) -> &'static str {
        match self {
            Self::Finite => "finite",
            Self::Positive => "positive",
            Self::NonNegative => "non-negative",
            Self::FiniteAndPositive => "finite and positive",
            Self::FiniteAndNonNegative => "finite and non-negative",
            Self::FiniteInUnitInterval => "finite and within [0, 1]",
            Self::NonEmpty => "non-empty",
        }
    }
}

impl std::fmt::Display for Requirement {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.phrase())
    }
}

/// The message prefix of the free-text variant for `domain`.
fn domain_prefix(domain: ErrorKind) -> &'static str {
    match domain {
        ErrorKind::Dispersion => "invalid dispersion",
        ErrorKind::Normalization => "invalid size factors",
        ErrorKind::Counts => "invalid counts",
        _ => "invalid options",
    }
}

/// Coarse classification of a [`Error`], stable across message changes.
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
    /// A design formula could not be parsed or evaluated.
    Formula,
    /// An internal consistency check failed (a bug).
    Internal,
}

impl Error {
    /// Coarse classification of this error.
    pub fn kind(&self) -> ErrorKind {
        match self {
            Self::InvalidDimensions { .. } => ErrorKind::Shape,
            Self::InvalidCounts { .. } | Self::AllZeroGene { .. } => ErrorKind::Counts,
            Self::InvalidSizeFactors { .. } | Self::NoUsableGenesForSizeFactors => {
                ErrorKind::Normalization
            }
            Self::InvalidDispersion { .. } => ErrorKind::Dispersion,
            Self::InvalidOptions { .. } => ErrorKind::Options,
            Self::NonFiniteValue { .. } => ErrorKind::Numeric,
            Self::UnsupportedFeature { .. } => ErrorKind::Unsupported,
            Self::Io(_) => ErrorKind::Io,
            #[cfg(feature = "csv")]
            Self::Csv(_) => ErrorKind::Parse,
            Self::ParseInt { .. } | Self::ParseFloat { .. } => ErrorKind::Parse,
            Self::InvalidParameter { domain, .. }
            | Self::MissingPrerequisite { domain, .. }
            | Self::InvalidGeneValue { domain, .. } => *domain,
            Self::InvalidFormula { .. } => ErrorKind::Formula,
            Self::CoefficientNotFound { .. }
            | Self::AmbiguousCoefficient { .. }
            | Self::DesignRankDeficient { .. }
            | Self::LrtDesignNotNested { .. } => ErrorKind::Options,
            Self::Internal { .. } => ErrorKind::Internal,
        }
    }
}

/// `InvalidOptions` with the reason as given.
#[doc(hidden)]
pub fn invalid_options(reason: impl Into<String>) -> Error {
    Error::InvalidOptions {
        reason: reason.into(),
    }
}

/// `InvalidOptions` whose reason is prefixed with the component that raised
/// it: `"{scope}: {reason}"`.
#[doc(hidden)]
pub fn scoped_invalid_options(scope: &str, reason: impl std::fmt::Display) -> Error {
    Error::InvalidOptions {
        reason: format!("{scope}: {reason}"),
    }
}

/// An internal consistency check failed: a state that earlier validation
/// rules out was reached. Reported as `InvalidOptions` so callers never panic.
#[doc(hidden)]
pub fn invariant(what: &str) -> Error {
    Error::Internal { what: what.into() }
}

/// `InvalidParameter` in `domain`.
#[doc(hidden)]
pub fn invalid_parameter(
    domain: ErrorKind,
    parameter: impl Into<String>,
    requirement: Requirement,
) -> Error {
    Error::InvalidParameter {
        domain,
        parameter: parameter.into(),
        requirement,
    }
}

/// `MissingPrerequisite` in `domain`.
#[doc(hidden)]
pub fn missing_prerequisite(
    domain: ErrorKind,
    operation: impl Into<String>,
    prerequisite: impl Into<String>,
) -> Error {
    Error::MissingPrerequisite {
        domain,
        operation: operation.into(),
        prerequisite: prerequisite.into(),
    }
}

/// `InvalidFormula` with the reason as given.
#[doc(hidden)]
pub fn invalid_formula(reason: impl Into<String>) -> Error {
    Error::InvalidFormula {
        reason: reason.into(),
    }
}

/// `InvalidDispersion` with the reason as given.
#[doc(hidden)]
pub fn invalid_dispersion(reason: impl Into<String>) -> Error {
    Error::InvalidDispersion {
        reason: reason.into(),
    }
}

#[doc(hidden)]
pub fn invalid_dimensions(context: impl Into<String>, expected: usize, actual: usize) -> Error {
    Error::InvalidDimensions {
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
        assert_eq!(invariant("x").kind(), ErrorKind::Internal);
        assert_eq!(invalid_formula("x").kind(), ErrorKind::Formula);
        let missing = missing_prerequisite(ErrorKind::Options, "VST", "a fitted dispersion trend");
        assert_eq!(missing.kind(), ErrorKind::Options);
        assert_eq!(
            missing.to_string(),
            "invalid options: VST requires a fitted dispersion trend"
        );
        assert_eq!(
            invalid_parameter(ErrorKind::Dispersion, "x", Requirement::Finite).kind(),
            ErrorKind::Dispersion
        );
        assert_eq!(
            Error::NoUsableGenesForSizeFactors.kind(),
            ErrorKind::Normalization
        );
        let io = Error::from(std::io::Error::other("x"));
        assert_eq!(io.kind(), ErrorKind::Io);
    }

    #[test]
    fn structured_variants_keep_the_free_text_messages() {
        let free = |domain, reason: String| match domain {
            ErrorKind::Dispersion => Error::InvalidDispersion { reason },
            ErrorKind::Normalization => Error::InvalidSizeFactors { reason },
            ErrorKind::Counts => Error::InvalidCounts { reason },
            _ => Error::InvalidOptions { reason },
        };
        for domain in [
            ErrorKind::Options,
            ErrorKind::Dispersion,
            ErrorKind::Normalization,
            ErrorKind::Counts,
        ] {
            for requirement in [
                Requirement::Finite,
                Requirement::Positive,
                Requirement::NonNegative,
                Requirement::FiniteAndPositive,
                Requirement::FiniteAndNonNegative,
                Requirement::FiniteInUnitInterval,
                Requirement::NonEmpty,
            ] {
                let structured = invalid_parameter(domain, "gene 3 tolerance", requirement);
                let reason = format!("gene 3 tolerance must be {}", requirement.phrase());
                assert_eq!(structured.to_string(), free(domain, reason).to_string());
            }
        }
        assert_eq!(
            invalid_formula("formula is empty").to_string(),
            invalid_options("formula is empty").to_string()
        );
        assert_eq!(
            invariant("x").to_string(),
            invalid_options("internal invariant violated: x").to_string()
        );
    }

    #[test]
    fn finite_out_of_range_values_are_not_called_non_finite() {
        let negative = Error::NonFiniteValue {
            context: "dispersion".into(),
            index: Some(1),
            value: -0.1,
        };
        assert_eq!(
            negative.to_string(),
            "out-of-range value in dispersion at index Some(1): -0.1"
        );
        let nan = Error::NonFiniteValue {
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
