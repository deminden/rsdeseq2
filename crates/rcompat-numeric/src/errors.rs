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

    /// A calculation produced a value outside its domain.
    #[error("{message}")]
    InvalidCalculation {
        /// Classification of the failing calculation.
        domain: ErrorKind,
        /// The calculation that failed.
        operation: &'static str,
        /// The input position associated with the failure, when known.
        location: Option<CalculationLocation>,
        /// The value rejected by the calculation.
        value: f64,
        /// The original diagnostic, including its domain prefix.
        message: String,
    },

    /// No gene passed the dispersion trend's eligibility rule.
    #[error(
        "invalid dispersion: all gene-wise dispersion estimates are within 2 orders of magnitude from the minimum value"
    )]
    NoUsableDispersionTrendRows {
        /// The configured minimum dispersion.
        min_disp: f64,
        /// Number of input rows examined.
        genes_checked: usize,
    },

    /// A negative-binomial likelihood input or intermediate failed validation.
    #[error("{}", stage.error_message(*sample, *value))]
    InvalidNegativeBinomialValue {
        /// Input or intermediate that failed.
        stage: NegativeBinomialStage,
        /// Zero-based count-matrix row, when the caller supplied a matrix.
        gene: Option<usize>,
        /// Zero-based sample involved, when known.
        sample: Option<usize>,
        /// The rejected input or arithmetic result.
        value: f64,
    },

    /// A column sum overflowed while checking the weighted Cox-Reid design.
    #[error(
        "out-of-range value in observation weight Cox-Reid column sum at index Some({column}): {term}"
    )]
    CoxReidColumnSumOverflow {
        /// Zero-based row of the observation weights.
        gene: usize,
        /// Zero-based sample that caused the partial sum to overflow.
        sample: usize,
        /// Zero-based column of the design matrix.
        column: usize,
        /// The finite absolute design entry added at this sample.
        term: f64,
        /// The overflowing partial sum.
        value: f64,
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
        numeric_value_label(*value)
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

    /// A name lookup or declaration has more than one matching entry.
    #[error("invalid options: {reason}")]
    AmbiguousName {
        /// Classification of the lookup or validation stage.
        domain: ErrorKind,
        /// The collection and role of the name.
        scope: NameScope,
        /// The requested or conflicting declared name.
        name: String,
        /// Names of matching entries in input order.
        matches: Vec<String>,
        /// The original diagnostic text.
        reason: String,
    },

    /// A requested name is absent from the matching collection.
    #[error("invalid options: {reason}")]
    NameNotFound {
        /// Classification of the lookup stage.
        domain: ErrorKind,
        /// The collection and role of the name.
        scope: NameScope,
        /// The requested name.
        name: String,
        /// The original diagnostic text.
        reason: String,
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

/// The position associated with a calculation failure.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[non_exhaustive]
pub enum CalculationLocation {
    /// Zero-based row in an input of genes.
    Gene(usize),
    /// Zero-based column in an input of coefficients.
    Coefficient(usize),
}

/// Input and arithmetic stages of a negative-binomial likelihood.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
#[non_exhaustive]
pub enum NegativeBinomialStage {
    /// Input mean.
    Mean,
    /// Input dispersion.
    Dispersion,
    /// Input observation weight.
    Weight,
    /// Per-sample log probability mass.
    LogPmf,
    /// Log probability mass scaled by the observation weight.
    WeightedTerm,
    /// Partial sum of the per-sample terms.
    RowSum,
    /// Negative twice the row log likelihood.
    Deviance,
}

/// Names searched or checked during formula design and contrast resolution.
#[derive(Clone, Debug, Eq, PartialEq)]
#[non_exhaustive]
pub enum NameScope {
    /// Variables in the combined factor and numeric model frame.
    FormulaVariable,
    /// Factor columns.
    Factor,
    /// Numeric covariate columns.
    NumericCovariate,
    /// Levels of a single factor.
    FactorLevel {
        /// The factor name when known.
        factor: Option<String>,
        /// How the requested level is used.
        role: FactorLevelRole,
    },
}

/// The use of a level in a design or contrast request.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[non_exhaustive]
pub enum FactorLevelRole {
    /// The contrast numerator.
    Numerator,
    /// The contrast denominator.
    Denominator,
    /// The treatment reference.
    Reference,
    /// An explicitly declared factor level.
    Declared,
    /// An observed sample's level.
    Sample {
        /// The zero-based sample index.
        index: usize,
    },
}

impl std::fmt::Display for FactorLevelRole {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Self::Numerator => "numerator",
            Self::Denominator => "denominator",
            Self::Reference => "reference",
            Self::Declared => "declared",
            Self::Sample { .. } => "sample",
        })
    }
}

impl NegativeBinomialStage {
    fn kind(self) -> ErrorKind {
        match self {
            Self::Mean | Self::Weight => ErrorKind::Numeric,
            _ => ErrorKind::Dispersion,
        }
    }

    fn error_message(self, sample: Option<usize>, value: f64) -> String {
        match self {
            Self::Mean | Self::Weight => {
                let quantity = if self == Self::Mean { "mean" } else { "weight" };
                format!(
                    "{} value in negative-binomial {quantity} at index {sample:?}: {value}",
                    numeric_value_label(value)
                )
            }
            Self::Dispersion => match sample {
                Some(sample) => format!(
                    "invalid dispersion: dispersion at index {sample} must be finite and positive"
                ),
                None => "invalid dispersion: dispersion must be finite and positive".into(),
            },
            Self::LogPmf | Self::WeightedTerm => {
                let quantity = if self == Self::LogPmf {
                    "log PMF"
                } else {
                    "weighted log-likelihood term"
                };
                match sample {
                    Some(sample) => format!(
                        "invalid dispersion: negative-binomial {quantity} at sample {sample} must be finite"
                    ),
                    None => {
                        format!("invalid dispersion: negative-binomial {quantity} must be finite")
                    }
                }
            }
            Self::RowSum => {
                "invalid dispersion: negative-binomial row log-likelihood sum must remain finite"
                    .into()
            }
            Self::Deviance => {
                "invalid dispersion: negative-binomial deviance must remain finite".into()
            }
        }
    }
}

fn numeric_value_label(value: f64) -> &'static str {
    if value.is_finite() {
        "out-of-range"
    } else {
        "non-finite"
    }
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
            Self::InvalidDispersion { .. } | Self::NoUsableDispersionTrendRows { .. } => {
                ErrorKind::Dispersion
            }
            Self::InvalidOptions { .. } => ErrorKind::Options,
            Self::NonFiniteValue { .. } | Self::CoxReidColumnSumOverflow { .. } => {
                ErrorKind::Numeric
            }
            Self::InvalidNegativeBinomialValue { stage, .. } => stage.kind(),
            Self::UnsupportedFeature { .. } => ErrorKind::Unsupported,
            Self::Io(_) => ErrorKind::Io,
            #[cfg(feature = "csv")]
            Self::Csv(_) => ErrorKind::Parse,
            Self::ParseInt { .. } | Self::ParseFloat { .. } => ErrorKind::Parse,
            Self::InvalidParameter { domain, .. }
            | Self::MissingPrerequisite { domain, .. }
            | Self::InvalidGeneValue { domain, .. }
            | Self::InvalidCalculation { domain, .. }
            | Self::AmbiguousName { domain, .. }
            | Self::NameNotFound { domain, .. } => *domain,
            Self::InvalidFormula { .. } => ErrorKind::Formula,
            Self::CoefficientNotFound { .. }
            | Self::AmbiguousCoefficient { .. }
            | Self::DesignRankDeficient { .. }
            | Self::LrtDesignNotNested { .. } => ErrorKind::Options,
            Self::Internal { .. } => ErrorKind::Internal,
        }
    }
}

#[doc(hidden)]
pub fn ambiguous_name(
    domain: ErrorKind,
    scope: NameScope,
    name: &str,
    matches: impl IntoIterator<Item = impl Into<String>>,
    reason: String,
) -> Error {
    Error::AmbiguousName {
        domain,
        scope,
        name: name.into(),
        matches: matches.into_iter().map(Into::into).collect(),
        reason,
    }
}

#[doc(hidden)]
pub fn invalid_calculation(
    domain: ErrorKind,
    operation: &'static str,
    location: Option<CalculationLocation>,
    value: f64,
    reason: impl std::fmt::Display,
) -> Error {
    Error::InvalidCalculation {
        domain,
        operation,
        location,
        value,
        message: format!("{}: {reason}", domain_prefix(domain)),
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
