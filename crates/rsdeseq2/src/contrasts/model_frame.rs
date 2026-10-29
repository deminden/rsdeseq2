//! Factor-level contrast lookup in formula model frames.

use super::resolution::validate_factor_level_contrast;
use super::types::{FactorLevelContrast, ResultsContrast};
use crate::design::{
    FormulaFactorColumn, FormulaModelFrame, r_like_name_candidates as candidate_names,
};
use crate::errors::DeseqError;
use std::collections::HashSet;

/// Resolve a character-style `results(contrast=...)` request against owned
/// formula model-frame metadata.
///
/// List and numeric contrasts return `Ok(None)` because they do not need
/// factor sample levels or reference metadata. Character contrasts return a
/// [`FactorLevelContrast`] borrowing the matching model-frame factor column.
/// When the contrast omits an explicit reference, the factor column reference
/// is used; if the column has no explicit reference, the first declared factor
/// level is used when present, otherwise the first observed sample level is
/// used, matching [`crate::design::expanded_formula_design_from_model_frame`].
pub fn factor_level_contrast_from_model_frame<'a>(
    contrast: &'a ResultsContrast,
    model_frame: &'a FormulaModelFrame,
) -> Result<Option<FactorLevelContrast<'a>>, DeseqError> {
    let ResultsContrast::Character {
        factor,
        numerator,
        denominator,
        reference,
    } = contrast
    else {
        return Ok(None);
    };
    validate_factor_level_contrast(factor, numerator, denominator)?;
    model_frame.validate()?;
    let column = resolve_model_frame_factor_column(model_frame, factor)?;
    if column.sample_levels.is_empty() {
        return Err(DeseqError::InvalidOptions {
            reason: format!("factor '{factor}' has no sample levels"),
        });
    }
    let resolved_numerator =
        resolve_observed_model_frame_factor_level(column, factor, numerator, "numerator")?;
    let resolved_denominator =
        resolve_observed_model_frame_factor_level(column, factor, denominator, "denominator")?;
    if resolved_numerator == resolved_denominator {
        return Err(DeseqError::InvalidOptions {
            reason: format!(
                "factor '{factor}' numerator and denominator resolve to the same level '{resolved_numerator}'"
            ),
        });
    }
    let inferred_reference = if let Some(reference) = reference.as_deref() {
        Some(resolve_any_model_frame_factor_level(
            column,
            factor,
            reference,
            "reference",
        )?)
    } else if let Some(reference) = column.reference.as_deref() {
        Some(resolve_any_model_frame_factor_level(
            column,
            factor,
            reference,
            "reference",
        )?)
    } else {
        column
            .levels
            .as_ref()
            .and_then(|levels| levels.first().map(String::as_str))
            .or_else(|| column.sample_levels.first().map(String::as_str))
    };
    Ok(Some(FactorLevelContrast {
        factor: &column.name,
        numerator: resolved_numerator,
        denominator: resolved_denominator,
        reference: inferred_reference,
        sample_levels: &column.sample_levels,
    }))
}

fn resolve_observed_model_frame_factor_level<'a>(
    column: &'a FormulaFactorColumn,
    factor: &str,
    requested: &str,
    role: &str,
) -> Result<&'a str, DeseqError> {
    let levels = unique_level_candidates(column.sample_levels.iter().map(String::as_str));
    resolve_model_frame_factor_level_from_candidates(&levels, factor, requested, role)
}

fn resolve_any_model_frame_factor_level<'a>(
    column: &'a FormulaFactorColumn,
    factor: &str,
    requested: &str,
    role: &str,
) -> Result<&'a str, DeseqError> {
    let levels = unique_level_candidates(
        column
            .levels
            .iter()
            .flatten()
            .map(String::as_str)
            .chain(column.sample_levels.iter().map(String::as_str)),
    );
    resolve_model_frame_factor_level_from_candidates(&levels, factor, requested, role)
}

fn unique_level_candidates<'a>(levels: impl Iterator<Item = &'a str>) -> Vec<&'a str> {
    let mut seen = HashSet::new();
    levels.filter(|level| seen.insert(*level)).collect()
}

fn resolve_model_frame_factor_level_from_candidates<'a>(
    levels: &[&'a str],
    factor: &str,
    requested: &str,
    role: &str,
) -> Result<&'a str, DeseqError> {
    let exact = levels
        .iter()
        .copied()
        .filter(|level| *level == requested)
        .collect::<Vec<_>>();
    match exact.as_slice() {
        [level] => return Ok(*level),
        [] => {}
        _ => {
            return Err(DeseqError::InvalidOptions {
                reason: format!(
                    "factor '{factor}' {role} level '{requested}' appears more than once"
                ),
            });
        }
    }

    let matches = levels
        .iter()
        .copied()
        .filter(|level| {
            candidate_names(level)
                .into_iter()
                .any(|candidate| candidate == requested)
        })
        .collect::<Vec<_>>();
    match matches.as_slice() {
        [level] => Ok(*level),
        [] => Err(DeseqError::InvalidOptions {
            reason: format!("factor '{factor}' does not contain {role} level '{requested}'"),
        }),
        _ => Err(DeseqError::InvalidOptions {
            reason: format!(
                "factor '{factor}' {role} level '{requested}' resolves ambiguously after R-style cleanup"
            ),
        }),
    }
}

fn resolve_model_frame_factor_column<'a>(
    model_frame: &'a FormulaModelFrame,
    factor: &str,
) -> Result<&'a FormulaFactorColumn, DeseqError> {
    let exact = model_frame
        .factors
        .iter()
        .filter(|column| column.name == factor)
        .collect::<Vec<_>>();
    match exact.as_slice() {
        [column] => return Ok(*column),
        [] => {}
        _ => {
            return Err(DeseqError::InvalidOptions {
                reason: format!("factor '{factor}' appears more than once in formula model frame"),
            });
        }
    }

    let matches = model_frame
        .factors
        .iter()
        .filter(|column| {
            candidate_names(&column.name)
                .into_iter()
                .any(|candidate| candidate == factor)
        })
        .collect::<Vec<_>>();
    match matches.as_slice() {
        [column] => Ok(*column),
        [] => Err(DeseqError::InvalidOptions {
            reason: format!("factor '{factor}' is not present in formula model frame"),
        }),
        _ => Err(DeseqError::InvalidOptions {
            reason: format!("factor '{factor}' resolves ambiguously after R-style cleanup"),
        }),
    }
}
