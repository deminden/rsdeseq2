//! Public contrast requests and resolved result metadata.

use super::aliases::{factor_level_result_name, list_contrast_comparison};

/// Primitive contrast specification for already-built design matrices.
///
/// This intentionally covers only the parts of DESeq2 contrast handling that
/// can be resolved from coefficient names and optional reference metadata.
/// Formula-aware callers can use [`factor_level_contrast_from_model_frame`](crate::contrasts::factor_level_contrast_from_model_frame) to
/// supply the matching factor levels for character-contrast all-zero handling.
#[derive(Clone, Debug, PartialEq)]
pub enum ContrastSpec {
    /// Explicit numeric contrast with one value per coefficient.
    Numeric(Vec<f64>),
    /// One named coefficient, equivalent to `results(name=...)` in DESeq2.
    CoefficientName(String),
    /// Positive and negative coefficient-name lists with DESeq2-style list values.
    List {
        /// Coefficients receiving `positive_weight`.
        positive: Vec<String>,
        /// Coefficients receiving `negative_weight`.
        negative: Vec<String>,
        /// Weight for coefficients in `positive`.
        positive_weight: f64,
        /// Weight for coefficients in `negative`.
        negative_weight: f64,
    },
    /// Factor-level contrast resolved from DESeq2-shaped coefficient names.
    ///
    /// This covers common coefficient shapes such as `condition_B_vs_A` and
    /// expanded/no-intercept shapes such as `conditionB`. Supplying a reference
    /// level allows non-reference comparisons such as `C` vs `B` to resolve to
    /// `condition_C_vs_A - condition_B_vs_A`.
    FactorLevel {
        /// Factor or variable name.
        factor: String,
        /// Numerator level.
        numerator: String,
        /// Denominator level.
        denominator: String,
        /// Optional reference/base level.
        reference: Option<String>,
    },
}

/// Factor-level contrast request with sample labels for character-style
/// `contrastAllZero` handling.
///
/// This request keeps the primitive values together after a caller has built
/// the model matrix and extracted sample levels. Formula-model-frame callers
/// can construct it with [`factor_level_contrast_from_model_frame`](crate::contrasts::factor_level_contrast_from_model_frame).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FactorLevelContrast<'a> {
    /// Factor or variable name.
    pub factor: &'a str,
    /// Numerator level.
    pub numerator: &'a str,
    /// Denominator level.
    pub denominator: &'a str,
    /// Optional reference/base level.
    pub reference: Option<&'a str>,
    /// One factor level per sample/column in the count matrix.
    pub sample_levels: &'a [String],
}

/// DESeq2 `results(contrast=...)` request shape.
///
/// This keeps the user-facing contrast form separate from the already-resolved
/// primitive [`ContrastSpec`]. In DESeq2, character triplet contrasts use
/// factor-level `contrastAllZero` handling, while list and numeric contrasts
/// use numeric-model-matrix handling. Preserving the source shape lets pipeline
/// callers choose the matching all-zero rule after resolving coefficient names.
#[derive(Clone, Debug, PartialEq)]
pub enum ResultsContrast {
    /// `contrast = c(factor, numerator, denominator)`.
    Character {
        /// Factor or variable name.
        factor: String,
        /// Numerator level.
        numerator: String,
        /// Denominator level.
        denominator: String,
        /// Optional reference/base level, supplied by formula-aware callers.
        reference: Option<String>,
    },
    /// `contrast = list(numeratorNames, denominatorNames)`.
    List {
        /// Coefficients receiving `list_values[0]`.
        positive: Vec<String>,
        /// Coefficients receiving `list_values[1]`.
        negative: Vec<String>,
        /// DESeq2 `listValues`, conventionally `c(1, -1)`.
        list_values: [f64; 2],
    },
    /// Numeric contrast with one value per `resultsNames(object)`.
    Numeric(Vec<f64>),
}

/// Which DESeq2 all-zero contrast rule applies after resolving a contrast.
#[derive(Clone, Debug, PartialEq)]
pub enum ResultsContrastAllZero {
    /// Use `contrastAllZeroCharacter`: inspect samples in the two factor levels.
    Character {
        /// Factor or variable name supplied by the character contrast.
        factor: String,
        /// Numerator level.
        numerator: String,
        /// Denominator level.
        denominator: String,
    },
    /// Use `contrastAllZeroNumeric`: inspect samples selected by model matrix
    /// columns with non-zero numeric contrast coefficients.
    Numeric,
}

/// Resolved `results(contrast=...)` information.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedResultsContrast {
    /// Primitive numeric contrast vector.
    pub numeric: Vec<f64>,
    /// Metadata-compatible result name.
    pub result_name: String,
    /// Metadata-compatible comparison label.
    pub comparison: String,
    /// Matching DESeq2 all-zero handling rule.
    pub all_zero: ResultsContrastAllZero,
}

impl<'a> FactorLevelContrast<'a> {
    /// Build a factor-level contrast without an explicit reference level.
    pub fn new(
        factor: &'a str,
        numerator: &'a str,
        denominator: &'a str,
        sample_levels: &'a [String],
    ) -> Self {
        Self {
            factor,
            numerator,
            denominator,
            reference: None,
            sample_levels,
        }
    }

    /// Build a factor-level contrast with an explicit reference/base level.
    pub fn with_reference(
        factor: &'a str,
        numerator: &'a str,
        denominator: &'a str,
        reference: &'a str,
        sample_levels: &'a [String],
    ) -> Self {
        Self {
            factor,
            numerator,
            denominator,
            reference: Some(reference),
            sample_levels,
        }
    }
}

impl ResultsContrast {
    /// Build `contrast = c(factor, numerator, denominator)`.
    pub fn character(
        factor: impl Into<String>,
        numerator: impl Into<String>,
        denominator: impl Into<String>,
    ) -> Self {
        Self::Character {
            factor: factor.into(),
            numerator: numerator.into(),
            denominator: denominator.into(),
            reference: None,
        }
    }

    /// Build a character triplet contrast with an explicit reference/base level.
    pub fn character_with_reference(
        factor: impl Into<String>,
        numerator: impl Into<String>,
        denominator: impl Into<String>,
        reference: impl Into<String>,
    ) -> Self {
        Self::Character {
            factor: factor.into(),
            numerator: numerator.into(),
            denominator: denominator.into(),
            reference: Some(reference.into()),
        }
    }

    /// Build `contrast = list(positive, negative)` with `listValues = c(1, -1)`.
    pub fn list(positive: Vec<String>, negative: Vec<String>) -> Self {
        Self::list_with_values(positive, negative, 1.0, -1.0)
    }

    /// Build a list contrast with explicit DESeq2 `listValues`.
    pub fn list_with_values(
        positive: Vec<String>,
        negative: Vec<String>,
        positive_weight: f64,
        negative_weight: f64,
    ) -> Self {
        Self::List {
            positive,
            negative,
            list_values: [positive_weight, negative_weight],
        }
    }

    /// Build a numeric contrast vector.
    pub fn numeric(values: Vec<f64>) -> Self {
        Self::Numeric(values)
    }

    /// Convert to the lower-level primitive contrast specification.
    pub fn as_contrast_spec(&self) -> ContrastSpec {
        match self {
            Self::Character {
                factor,
                numerator,
                denominator,
                reference,
            } => ContrastSpec::FactorLevel {
                factor: factor.clone(),
                numerator: numerator.clone(),
                denominator: denominator.clone(),
                reference: reference.clone(),
            },
            Self::List {
                positive,
                negative,
                list_values,
            } => ContrastSpec::List {
                positive: positive.clone(),
                negative: negative.clone(),
                positive_weight: list_values[0],
                negative_weight: list_values[1],
            },
            Self::Numeric(values) => ContrastSpec::Numeric(values.clone()),
        }
    }
}

impl ContrastSpec {
    /// Build a coefficient-name contrast.
    pub fn coefficient_name(name: impl Into<String>) -> Self {
        Self::CoefficientName(name.into())
    }

    /// Build a DESeq2-style list contrast with default `listValues=c(1, -1)`.
    pub fn list(positive: Vec<String>, negative: Vec<String>) -> Self {
        Self::List {
            positive,
            negative,
            positive_weight: 1.0,
            negative_weight: -1.0,
        }
    }

    /// Build a list contrast with explicit list values.
    pub fn list_with_values(
        positive: Vec<String>,
        negative: Vec<String>,
        positive_weight: f64,
        negative_weight: f64,
    ) -> Self {
        Self::List {
            positive,
            negative,
            positive_weight,
            negative_weight,
        }
    }

    /// Build a factor-level contrast without an explicit reference level.
    pub fn factor_level(
        factor: impl Into<String>,
        numerator: impl Into<String>,
        denominator: impl Into<String>,
    ) -> Self {
        Self::FactorLevel {
            factor: factor.into(),
            numerator: numerator.into(),
            denominator: denominator.into(),
            reference: None,
        }
    }

    /// Build a factor-level contrast with an explicit reference/base level.
    pub fn factor_level_with_reference(
        factor: impl Into<String>,
        numerator: impl Into<String>,
        denominator: impl Into<String>,
        reference: impl Into<String>,
    ) -> Self {
        Self::FactorLevel {
            factor: factor.into(),
            numerator: numerator.into(),
            denominator: denominator.into(),
            reference: Some(reference.into()),
        }
    }

    /// Stable result name for metadata after resolving this contrast.
    pub fn result_name(&self) -> String {
        match self {
            Self::Numeric(_) => "contrast".to_string(),
            Self::CoefficientName(name) => name.clone(),
            Self::List { .. } => "contrast".to_string(),
            Self::FactorLevel {
                factor,
                numerator,
                denominator,
                ..
            } => factor_level_result_name(factor, numerator, denominator),
        }
    }

    /// Stable comparison label for result-table metadata.
    pub fn comparison(&self) -> String {
        match self {
            Self::Numeric(_) => "primitive numeric contrast".to_string(),
            Self::CoefficientName(name) => format!("coefficient {name}"),
            Self::List {
                positive,
                negative,
                positive_weight,
                negative_weight,
            } => list_contrast_comparison(positive, negative, *positive_weight, *negative_weight),
            Self::FactorLevel {
                factor,
                numerator,
                denominator,
                ..
            } => format!("factor-level contrast: {factor} {numerator} vs {denominator}"),
        }
    }
}
