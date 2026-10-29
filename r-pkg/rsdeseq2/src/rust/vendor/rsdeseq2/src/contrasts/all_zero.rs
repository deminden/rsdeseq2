//! Contrast-specific all-zero rules and sample-level aliases.

use super::resolution::{validate_factor_level_contrast, validate_numeric_contrast};
use crate::core::CountMatrix;
use crate::design::{DesignMatrix, r_like_name_candidates as candidate_names};
use crate::errors::{DeseqError, invalid_dimensions};
use std::collections::HashSet;

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
/// counts are zero. Use [`factor_level_contrast_from_model_frame`](crate::contrasts::factor_level_contrast_from_model_frame) when the
/// levels originate from owned formula metadata.
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
    let observed_levels = unique_sample_levels(sample_levels);
    let numerator = resolve_sample_level_alias(&observed_levels, numerator, "numerator")?;
    let denominator = resolve_sample_level_alias(&observed_levels, denominator, "denominator")?;
    if numerator == denominator {
        return Err(DeseqError::InvalidOptions {
            reason: format!(
                "contrastAllZero numerator and denominator resolve to the same level '{numerator}'"
            ),
        });
    }
    let mut numerator_count = 0usize;
    let mut denominator_count = 0usize;
    let selected_samples = sample_levels
        .iter()
        .map(|level| {
            let level = level.as_ref();
            if level == numerator {
                numerator_count += 1;
            }
            if level == denominator {
                denominator_count += 1;
            }
            level == numerator || level == denominator
        })
        .collect::<Vec<_>>();
    if numerator_count == 0 || denominator_count == 0 {
        return Err(DeseqError::InvalidOptions {
            reason: format!(
                "contrastAllZero sample levels must contain numerator level '{numerator}' and denominator level '{denominator}'"
            ),
        });
    }
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

fn unique_sample_levels<S: AsRef<str>>(sample_levels: &[S]) -> Vec<&str> {
    let mut seen = HashSet::new();
    sample_levels
        .iter()
        .filter_map(|level| {
            let level = level.as_ref();
            if seen.insert(level) {
                Some(level)
            } else {
                None
            }
        })
        .collect()
}

fn resolve_sample_level_alias<'a>(
    observed_levels: &[&'a str],
    requested: &str,
    role: &str,
) -> Result<&'a str, DeseqError> {
    let exact = observed_levels
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
                    "contrastAllZero {role} level '{requested}' appears more than once"
                ),
            });
        }
    }

    let matches = observed_levels
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
            reason: format!(
                "contrastAllZero sample levels do not contain {role} level '{requested}'"
            ),
        }),
        _ => Err(DeseqError::InvalidOptions {
            reason: format!(
                "contrastAllZero {role} level '{requested}' resolves ambiguously after R-style cleanup"
            ),
        }),
    }
}
