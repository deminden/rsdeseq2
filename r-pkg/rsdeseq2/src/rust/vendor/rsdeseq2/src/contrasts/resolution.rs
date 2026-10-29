//! Validation and resolution of numeric, named, list, and factor contrasts.

use super::aliases::{
    coefficient_names, expanded_coefficient_names, find_first_coefficient,
    find_optional_coefficient, find_shared_reference_standard_coefficients,
    find_treatment_or_expanded_coefficient, resolve_coefficient_name,
    resolve_coefficient_name_list, standard_coefficient_names,
};
use super::types::{
    ContrastSpec, ResolvedResultsContrast, ResultsContrast, ResultsContrastAllZero,
};
use crate::design::DesignMatrix;
use crate::errors::{DeseqError, invalid_dimensions};

/// Resolve a DESeq2 `results(contrast=...)` request into the primitive numeric
/// vector plus the all-zero rule implied by the original contrast form.
pub fn resolve_results_contrast(
    design: &DesignMatrix,
    contrast: &ResultsContrast,
) -> Result<ResolvedResultsContrast, DeseqError> {
    let spec = contrast.as_contrast_spec();
    let numeric = resolve_contrast(design, &spec)?;
    let all_zero = match contrast {
        ResultsContrast::Character {
            factor,
            numerator,
            denominator,
            ..
        } => ResultsContrastAllZero::Character {
            factor: factor.clone(),
            numerator: numerator.clone(),
            denominator: denominator.clone(),
        },
        ResultsContrast::List { .. } | ResultsContrast::Numeric(_) => {
            ResultsContrastAllZero::Numeric
        }
    };
    Ok(ResolvedResultsContrast {
        numeric,
        result_name: spec.result_name(),
        comparison: spec.comparison(),
        all_zero,
    })
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
            let denominator_index = find_treatment_or_expanded_coefficient(
                names,
                &standard_coefficient_names(factor, denominator, reference),
                &expanded_coefficient_names(factor, denominator),
            )?;
            let mut values = vec![0.0; design.n_coefficients()];
            values[denominator_index] = -1.0;
            return Ok(values);
        }
        if denominator == reference {
            let numerator_index = find_treatment_or_expanded_coefficient(
                names,
                &standard_coefficient_names(factor, numerator, reference),
                &expanded_coefficient_names(factor, numerator),
            )?;
            let mut values = vec![0.0; design.n_coefficients()];
            values[numerator_index] = 1.0;
            return Ok(values);
        }

        let numerator_index = find_treatment_or_expanded_coefficient(
            names,
            &standard_coefficient_names(factor, numerator, reference),
            &expanded_coefficient_names(factor, numerator),
        )?;
        let denominator_index = find_treatment_or_expanded_coefficient(
            names,
            &standard_coefficient_names(factor, denominator, reference),
            &expanded_coefficient_names(factor, denominator),
        )?;
        let mut values = vec![0.0; design.n_coefficients()];
        values[numerator_index] = 1.0;
        values[denominator_index] = -1.0;
        return Ok(values);
    }

    if let Some(index) = find_optional_coefficient(
        names,
        &standard_coefficient_names(factor, numerator, denominator),
    )? {
        let mut values = vec![0.0; design.n_coefficients()];
        values[index] = 1.0;
        return Ok(values);
    }
    if let Some(index) = find_optional_coefficient(
        names,
        &standard_coefficient_names(factor, denominator, numerator),
    )? {
        let mut values = vec![0.0; design.n_coefficients()];
        values[index] = -1.0;
        return Ok(values);
    }

    if let Some((numerator_index, denominator_index)) =
        find_shared_reference_standard_coefficients(names, factor, numerator, denominator)?
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

pub(super) fn validate_numeric_contrast(
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

pub(super) fn validate_factor_level_contrast(
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
