use std::collections::HashSet;

use crate::core::CountMatrix;
use crate::design::DesignMatrix;
use crate::errors::{invalid_dimensions, DeseqError};

/// Primitive contrast specification for already-built design matrices.
///
/// This intentionally covers only the parts of DESeq2 contrast handling that
/// can be resolved from coefficient names alone. Formula parsing and factor
/// level lookup remain caller or future-wrapper work.
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
/// The Rust core does not parse R formulas or own colData. This request keeps
/// the necessary primitive values together after a caller or future wrapper has
/// already built the model matrix and extracted sample levels.
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
            } => format!("{factor}_{numerator}_vs_{denominator}"),
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

/// Resolve a primitive contrast specification into a numeric contrast vector.
pub fn resolve_contrast(
    design: &DesignMatrix,
    contrast: &ContrastSpec,
) -> Result<Vec<f64>, DeseqError> {
    match contrast {
        ContrastSpec::Numeric(values) => validate_numeric_contrast(values, design.n_coefficients()),
        ContrastSpec::CoefficientName(name) => {
            let names = coefficient_names(design)?;
            let index = resolve_coefficient_name(names, name)?;
            let mut values = vec![0.0; design.n_coefficients()];
            values[index] = 1.0;
            Ok(values)
        }
        ContrastSpec::List {
            positive,
            negative,
            positive_weight,
            negative_weight,
        } => resolve_list_contrast(
            design,
            positive,
            negative,
            *positive_weight,
            *negative_weight,
        ),
        ContrastSpec::FactorLevel {
            factor,
            numerator,
            denominator,
            reference,
        } => resolve_factor_level_contrast(design, factor, numerator, denominator, reference),
    }
}

/// Resolve a design coefficient name with the same DESeq2-style aliases used by
/// named primitive contrasts.
///
/// Exact coefficient names win first. If no exact name exists, the resolver
/// accepts R-cleaned aliases, including the `Intercept`/`(Intercept)` spelling
/// pair, and reports ambiguous aliases instead of choosing arbitrarily.
pub fn resolve_coefficient_index(
    design: &DesignMatrix,
    coefficient_name: &str,
) -> Result<usize, DeseqError> {
    resolve_coefficient_name(coefficient_names(design)?, coefficient_name)
}

/// Identify rows where every sample involved in a numeric contrast has zero counts.
///
/// This mirrors DESeq2's `contrastAllZeroNumeric` helper for primitive
/// matrices: if a numeric contrast contains both positive and negative
/// coefficients, non-zero contrast coefficients are converted to one, samples
/// with non-zero `modelMatrix %*% contrastBinary` are selected, and rows whose
/// selected raw counts sum to zero are flagged. One-sided numeric contrasts
/// return all `false`, matching DESeq2.
pub fn contrast_all_zero_numeric(
    counts: &CountMatrix,
    design: &DesignMatrix,
    contrast: &[f64],
) -> Result<Vec<bool>, DeseqError> {
    if design.n_samples() != counts.n_samples() {
        return Err(invalid_dimensions(
            "contrastAllZero design samples",
            counts.n_samples(),
            design.n_samples(),
        ));
    }
    let contrast = validate_numeric_contrast(contrast, design.n_coefficients())?;
    if contrast.iter().all(|value| *value >= 0.0) || contrast.iter().all(|value| *value <= 0.0) {
        return Ok(vec![false; counts.n_genes()]);
    }

    let contrast_binary = contrast
        .iter()
        .map(|value| if *value == 0.0 { 0.0 } else { 1.0 })
        .collect::<Vec<_>>();
    let mut selected_samples = Vec::with_capacity(design.n_samples());
    for sample in 0..design.n_samples() {
        let row = design.matrix().row(sample)?;
        let mut terms = Vec::with_capacity(row.len());
        for (design_value, contrast_value) in row.iter().zip(contrast_binary.iter()) {
            let term = design_value * contrast_value;
            if !term.is_finite() {
                return Err(DeseqError::NonFiniteValue {
                    context: "contrastAllZero design score".to_string(),
                    index: Some(sample),
                    value: f64::NAN,
                });
            }
            terms.push(term);
        }
        let Some(score) = checked_scaled_sum(&terms) else {
            return Err(DeseqError::NonFiniteValue {
                context: "contrastAllZero design score".to_string(),
                index: Some(sample),
                value: f64::NAN,
            });
        };
        selected_samples.push(score != 0.0);
    }

    let mut flags = Vec::with_capacity(counts.n_genes());
    for gene in 0..counts.n_genes() {
        let selected_sum = counts
            .row(gene)?
            .iter()
            .zip(selected_samples.iter())
            .filter_map(|(count, selected)| selected.then_some(*count as u64))
            .sum::<u64>();
        flags.push(selected_sum == 0);
    }
    Ok(flags)
}

fn checked_scaled_sum(values: &[f64]) -> Option<f64> {
    let mut scale = 0.0_f64;
    for value in values.iter().copied() {
        if !value.is_finite() {
            return None;
        }
        scale = scale.max(value.abs());
    }
    if scale == 0.0 {
        return Some(0.0);
    }
    let mut normalized_sum = 0.0;
    for value in values.iter().copied() {
        let term = value / scale;
        let next = normalized_sum + term;
        if !term.is_finite() || !next.is_finite() {
            return None;
        }
        normalized_sum = next;
    }
    let sum = normalized_sum * scale;
    sum.is_finite().then_some(sum)
}

/// Identify rows where both requested factor levels have zero counts.
///
/// This is the primitive-matrix analogue of DESeq2's
/// `contrastAllZeroCharacter`: select samples whose supplied level is either
/// `numerator` or `denominator`, then flag genes for which all selected raw
/// counts are zero. Full factor validation remains caller or future-wrapper
/// responsibility.
pub fn contrast_all_zero_factor_levels<S: AsRef<str>>(
    counts: &CountMatrix,
    sample_levels: &[S],
    numerator: &str,
    denominator: &str,
) -> Result<Vec<bool>, DeseqError> {
    validate_factor_level_contrast("factor", numerator, denominator)?;
    if sample_levels.len() != counts.n_samples() {
        return Err(invalid_dimensions(
            "contrastAllZero sample levels",
            counts.n_samples(),
            sample_levels.len(),
        ));
    }
    let selected_samples = sample_levels
        .iter()
        .map(|level| {
            let level = level.as_ref();
            level == numerator || level == denominator
        })
        .collect::<Vec<_>>();
    let selected_count = selected_samples
        .iter()
        .filter(|selected| **selected)
        .count();

    let mut flags = Vec::with_capacity(counts.n_genes());
    for gene in 0..counts.n_genes() {
        let zero_count = counts
            .row(gene)?
            .iter()
            .zip(selected_samples.iter())
            .filter(|(_, selected)| **selected)
            .filter(|(count, _)| **count == 0)
            .count();
        flags.push(zero_count == selected_count);
    }
    Ok(flags)
}

fn resolve_factor_level_contrast(
    design: &DesignMatrix,
    factor: &str,
    numerator: &str,
    denominator: &str,
    reference: &Option<String>,
) -> Result<Vec<f64>, DeseqError> {
    validate_factor_level_contrast(factor, numerator, denominator)?;
    let names = coefficient_names(design)?;

    if let Some(reference) = reference {
        if numerator == reference {
            let denominator_index = find_first_coefficient(
                names,
                &standard_coefficient_names(factor, denominator, reference),
            )
            .or_else(|_| {
                find_first_coefficient(names, &expanded_coefficient_names(factor, denominator))
            })?;
            let mut values = vec![0.0; design.n_coefficients()];
            values[denominator_index] = -1.0;
            return Ok(values);
        }
        if denominator == reference {
            let numerator_index = find_first_coefficient(
                names,
                &standard_coefficient_names(factor, numerator, reference),
            )
            .or_else(|_| {
                find_first_coefficient(names, &expanded_coefficient_names(factor, numerator))
            })?;
            let mut values = vec![0.0; design.n_coefficients()];
            values[numerator_index] = 1.0;
            return Ok(values);
        }

        let numerator_index = find_first_coefficient(
            names,
            &standard_coefficient_names(factor, numerator, reference),
        )?;
        let denominator_index = find_first_coefficient(
            names,
            &standard_coefficient_names(factor, denominator, reference),
        )?;
        let mut values = vec![0.0; design.n_coefficients()];
        values[numerator_index] = 1.0;
        values[denominator_index] = -1.0;
        return Ok(values);
    }

    if let Ok(index) = find_first_coefficient(
        names,
        &standard_coefficient_names(factor, numerator, denominator),
    ) {
        let mut values = vec![0.0; design.n_coefficients()];
        values[index] = 1.0;
        return Ok(values);
    }
    if let Ok(index) = find_first_coefficient(
        names,
        &standard_coefficient_names(factor, denominator, numerator),
    ) {
        let mut values = vec![0.0; design.n_coefficients()];
        values[index] = -1.0;
        return Ok(values);
    }

    if let Some((numerator_index, denominator_index)) =
        find_shared_reference_standard_coefficients(names, factor, numerator, denominator)
    {
        let mut values = vec![0.0; design.n_coefficients()];
        values[numerator_index] = 1.0;
        values[denominator_index] = -1.0;
        return Ok(values);
    }

    let numerator_index =
        find_first_coefficient(names, &expanded_coefficient_names(factor, numerator));
    let denominator_index =
        find_first_coefficient(names, &expanded_coefficient_names(factor, denominator));
    match (numerator_index, denominator_index) {
        (Ok(numerator_index), Ok(denominator_index)) => {
            let mut values = vec![0.0; design.n_coefficients()];
            values[numerator_index] = 1.0;
            values[denominator_index] = -1.0;
            Ok(values)
        }
        _ => Err(DeseqError::InvalidOptions {
            reason: format!(
                "factor-level contrast {factor}: {numerator} vs {denominator} could not be resolved from coefficient names; provide a reference level or a numeric contrast"
            ),
        }),
    }
}

fn resolve_list_contrast(
    design: &DesignMatrix,
    positive: &[String],
    negative: &[String],
    positive_weight: f64,
    negative_weight: f64,
) -> Result<Vec<f64>, DeseqError> {
    validate_list_values(positive, negative, positive_weight, negative_weight)?;
    let names = coefficient_names(design)?;
    let positive_indices = resolve_coefficient_name_list(names, positive)?;
    let negative_indices = resolve_coefficient_name_list(names, negative)?;
    if positive_indices
        .iter()
        .any(|index| negative_indices.contains(index))
    {
        return Err(DeseqError::InvalidOptions {
            reason: "contrast list entries must not appear in both numerator and denominator"
                .to_string(),
        });
    }

    let mut values = vec![0.0; design.n_coefficients()];
    for index in positive_indices {
        values[index] = positive_weight;
    }
    for index in negative_indices {
        values[index] = negative_weight;
    }
    validate_numeric_contrast(&values, design.n_coefficients())
}

fn validate_numeric_contrast(
    values: &[f64],
    n_coefficients: usize,
) -> Result<Vec<f64>, DeseqError> {
    if values.len() != n_coefficients {
        return Err(invalid_dimensions(
            "contrast coefficients",
            n_coefficients,
            values.len(),
        ));
    }
    let mut any_nonzero = false;
    for (idx, value) in values.iter().copied().enumerate() {
        if !value.is_finite() {
            return Err(DeseqError::NonFiniteValue {
                context: "contrast coefficient".to_string(),
                index: Some(idx),
                value,
            });
        }
        any_nonzero |= value != 0.0;
    }
    if !any_nonzero {
        return Err(DeseqError::InvalidOptions {
            reason: "contrast vector cannot be all zero".to_string(),
        });
    }
    Ok(values.to_vec())
}

fn validate_list_values(
    positive: &[String],
    negative: &[String],
    positive_weight: f64,
    negative_weight: f64,
) -> Result<(), DeseqError> {
    if positive.is_empty() && negative.is_empty() {
        return Err(DeseqError::InvalidOptions {
            reason: "contrast list must contain at least one coefficient name".to_string(),
        });
    }
    for (label, value) in [
        ("positive contrast weight", positive_weight),
        ("negative contrast weight", negative_weight),
    ] {
        if !value.is_finite() {
            return Err(DeseqError::NonFiniteValue {
                context: label.to_string(),
                index: None,
                value,
            });
        }
    }
    if positive_weight <= 0.0 {
        return Err(DeseqError::InvalidOptions {
            reason: "positive contrast weight must be greater than zero".to_string(),
        });
    }
    if negative_weight >= 0.0 {
        return Err(DeseqError::InvalidOptions {
            reason: "negative contrast weight must be less than zero".to_string(),
        });
    }
    Ok(())
}

fn validate_factor_level_contrast(
    factor: &str,
    numerator: &str,
    denominator: &str,
) -> Result<(), DeseqError> {
    if factor.is_empty() || numerator.is_empty() || denominator.is_empty() {
        return Err(DeseqError::InvalidOptions {
            reason: "factor-level contrast factor and levels must be non-empty".to_string(),
        });
    }
    if numerator == denominator {
        return Err(DeseqError::InvalidOptions {
            reason: "factor-level contrast numerator and denominator must differ".to_string(),
        });
    }
    Ok(())
}

fn coefficient_names(design: &DesignMatrix) -> Result<&[String], DeseqError> {
    design
        .coefficient_names()
        .ok_or_else(|| DeseqError::InvalidOptions {
            reason: "coefficient names are required to resolve named contrasts".to_string(),
        })
}

fn find_coefficient(names: &[String], wanted: &str) -> Result<usize, DeseqError> {
    names
        .iter()
        .position(|name| name == wanted)
        .ok_or_else(|| DeseqError::InvalidOptions {
            reason: format!("coefficient '{wanted}' is not present in coefficient names"),
        })
}

fn resolve_coefficient_name(names: &[String], wanted: &str) -> Result<usize, DeseqError> {
    if let Ok(index) = find_coefficient(names, wanted) {
        return Ok(index);
    }
    let candidates = coefficient_name_candidates(wanted)
        .into_iter()
        .filter(|candidate| candidate != wanted)
        .collect::<Vec<_>>();
    let matches = candidates
        .iter()
        .filter_map(|candidate| names.iter().position(|name| name == candidate))
        .collect::<HashSet<_>>();
    match matches.len() {
        0 => Err(DeseqError::InvalidOptions {
            reason: format!(
                "coefficient '{wanted}' is not present in coefficient names or R-cleaned aliases"
            ),
        }),
        1 => Ok(*matches.iter().next().unwrap()),
        _ => Err(DeseqError::InvalidOptions {
            reason: format!("coefficient '{wanted}' resolves ambiguously after R-style cleanup"),
        }),
    }
}

fn resolve_coefficient_name_list(
    names: &[String],
    wanted: &[String],
) -> Result<Vec<usize>, DeseqError> {
    let mut indices = Vec::with_capacity(wanted.len());
    let mut seen = HashSet::with_capacity(wanted.len());
    for name in wanted {
        let index = resolve_coefficient_name(names, name)?;
        if !seen.insert(index) {
            return Err(DeseqError::InvalidOptions {
                reason: format!(
                    "contrast list contains duplicate coefficient '{name}' after R-style cleanup"
                ),
            });
        }
        indices.push(index);
    }
    Ok(indices)
}

fn coefficient_name_candidates(name: &str) -> Vec<String> {
    let mut candidates = candidate_names(name);
    if name == "(Intercept)" {
        push_unique_candidate(&mut candidates, "Intercept".to_string());
    } else if name == "Intercept" {
        push_unique_candidate(&mut candidates, "(Intercept)".to_string());
    }
    candidates
}

fn list_contrast_comparison(
    positive: &[String],
    negative: &[String],
    positive_weight: f64,
    negative_weight: f64,
) -> String {
    let positive_label = weighted_list_label(positive, positive_weight.abs());
    let negative_label = weighted_list_label(negative, negative_weight.abs());
    if !positive.is_empty() && !negative.is_empty() {
        format!("coefficient list contrast: {positive_label} vs {negative_label}")
    } else if !positive.is_empty() {
        format!("coefficient list contrast: {positive_label} effect")
    } else {
        format!("coefficient list contrast: -{negative_label} effect")
    }
}

fn weighted_list_label(names: &[String], weight: f64) -> String {
    let names = names.join("+");
    if (weight - 1.0).abs() <= f64::EPSILON {
        names
    } else {
        format!("{} {names}", format_rounded_weight(weight))
    }
}

fn format_rounded_weight(weight: f64) -> String {
    let rounded = (weight * 1000.0).round() / 1000.0;
    let formatted = format!("{rounded:.3}");
    formatted
        .trim_end_matches('0')
        .trim_end_matches('.')
        .to_string()
}

fn find_first_coefficient(names: &[String], candidates: &[String]) -> Result<usize, DeseqError> {
    candidates
        .iter()
        .find_map(|candidate| names.iter().position(|name| name == candidate))
        .ok_or_else(|| DeseqError::InvalidOptions {
            reason: format!(
                "none of the candidate coefficients are present: {}",
                candidates.join(", ")
            ),
        })
}

fn find_shared_reference_standard_coefficients(
    names: &[String],
    factor: &str,
    numerator: &str,
    denominator: &str,
) -> Option<(usize, usize)> {
    let numerator_pairs = standard_coefficients_for_level(names, factor, numerator);
    let denominator_pairs = standard_coefficients_for_level(names, factor, denominator);
    numerator_pairs
        .iter()
        .find_map(|(numerator_index, numerator_reference)| {
            denominator_pairs
                .iter()
                .find(|(_, denominator_reference)| denominator_reference == numerator_reference)
                .map(|(denominator_index, _)| (*numerator_index, *denominator_index))
        })
}

fn standard_coefficients_for_level(
    names: &[String],
    factor: &str,
    level: &str,
) -> Vec<(usize, String)> {
    standard_coefficient_prefixes(factor, level)
        .into_iter()
        .flat_map(|prefix| {
            names.iter().enumerate().filter_map(move |(index, name)| {
                name.strip_prefix(&prefix)
                    .map(|reference| (index, reference.to_string()))
            })
        })
        .collect()
}

fn standard_coefficient_prefixes(factor: &str, level: &str) -> Vec<String> {
    let raw = format!("{factor}_{level}_vs_");
    let mut candidates = candidate_names(&raw);
    push_unique_candidate(
        &mut candidates,
        format!(
            "{}_{}_vs_",
            r_like_make_name(factor),
            r_like_make_name(level)
        ),
    );
    candidates
}

fn standard_coefficient_names(factor: &str, level: &str, reference: &str) -> Vec<String> {
    standard_coefficient_candidates(factor, level, reference)
}

fn expanded_coefficient_names(factor: &str, level: &str) -> Vec<String> {
    let mut candidates = candidate_names(&format!("{factor}{level}"));
    push_unique_candidate(
        &mut candidates,
        format!("{}{}", r_like_make_name(factor), r_like_make_name(level)),
    );
    candidates
}

fn standard_coefficient_candidates(factor: &str, level: &str, reference: &str) -> Vec<String> {
    let raw = format!("{factor}_{level}_vs_{reference}");
    let mut candidates = candidate_names(&raw);
    push_unique_candidate(
        &mut candidates,
        format!(
            "{}_{}_vs_{}",
            r_like_make_name(factor),
            r_like_make_name(level),
            r_like_make_name(reference)
        ),
    );
    candidates
}

fn candidate_names(raw: &str) -> Vec<String> {
    let made = r_like_make_name(raw);
    if made == raw {
        vec![raw.to_string()]
    } else {
        vec![raw.to_string(), made]
    }
}

fn push_unique_candidate(candidates: &mut Vec<String>, candidate: String) {
    if !candidates.iter().any(|existing| existing == &candidate) {
        candidates.push(candidate);
    }
}

fn r_like_make_name(raw: &str) -> String {
    let mut out = String::with_capacity(raw.len().max(1));
    for ch in raw.chars() {
        if ch.is_ascii_alphanumeric() || ch == '.' || ch == '_' {
            out.push(ch);
        } else {
            out.push('.');
        }
    }
    if out.is_empty() {
        return "X".to_string();
    }
    let mut chars = out.chars();
    let first = chars.next().unwrap_or('X');
    let second = chars.next();
    let invalid_start = !first.is_ascii_alphabetic() && first != '.'
        || (first == '.' && second.is_some_and(|ch| ch.is_ascii_digit()));
    if invalid_start {
        out.insert(0, 'X');
    }
    if is_r_reserved_word(&out) {
        out.push('.');
    }
    out
}

fn is_r_reserved_word(name: &str) -> bool {
    matches!(
        name,
        "if" | "else"
            | "repeat"
            | "while"
            | "function"
            | "for"
            | "in"
            | "next"
            | "break"
            | "TRUE"
            | "FALSE"
            | "NULL"
            | "Inf"
            | "NaN"
            | "NA"
            | "NA_integer_"
            | "NA_real_"
            | "NA_complex_"
            | "NA_character_"
    )
}
