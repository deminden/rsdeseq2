//! Contrast-specific all-zero rules and sample-level aliases.

use super::counts_view::ContrastCounts;
use super::resolution::{validate_factor_level_contrast, validate_numeric_contrast};
use crate::core::CountMatrix;
use crate::design::{DesignMatrix, r_like_name_candidates as candidate_names};
use crate::errors::{DeseqError, invalid_dimensions};
use crate::math::sum::checked_scaled_sum;
use crate::matrix::MatrixView;
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
    match numeric_samples(counts.n_samples(), Some(design.matrix().into()), contrast)? {
        Some(selected) => ContrastCounts::from(counts).selected_flags(&selected),
        None => Ok(vec![false; counts.n_genes()]),
    }
}

pub(super) fn numeric_samples(
    samples: usize,
    design: Option<MatrixView<'_, f64>>,
    contrast: &[f64],
) -> Result<Option<Vec<bool>>, DeseqError> {
    if let Some(design) = design {
        if design.n_rows() != samples {
            return Err(invalid_dimensions(
                "contrastAllZero design samples",
                samples,
                design.n_rows(),
            ));
        }
        if design.n_cols() != contrast.len() {
            return Err(invalid_dimensions(
                "contrastAllZero design coefficients",
                design.n_cols(),
                contrast.len(),
            ));
        }
    }
    let contrast = validate_numeric_contrast(contrast, contrast.len())?;
    if contrast.iter().all(|value| *value >= 0.0) || contrast.iter().all(|value| *value <= 0.0) {
        return Ok(None);
    }
    let design = design.ok_or_else(|| DeseqError::InvalidOptions {
        reason: "numeric contrast all-zero handling requires modelMatrix".into(),
    })?;
    let mut terms = Vec::with_capacity(contrast.len());
    let mut selected = Vec::with_capacity(samples);
    for sample in 0..samples {
        terms.clear();
        for (column, &weight) in contrast.iter().enumerate() {
            terms.push(design.at(sample, column) * if weight == 0.0 { 0.0 } else { 1.0 });
        }
        let score = checked_scaled_sum(&terms).ok_or_else(|| DeseqError::NonFiniteValue {
            context: "contrastAllZero design score".into(),
            index: Some(sample),
            value: f64::NAN,
        })?;
        selected.push(score != 0.0);
    }
    Ok(Some(selected))
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
    let selected = factor_level_samples(counts.n_samples(), sample_levels, numerator, denominator)?;
    ContrastCounts::from(counts).selected_flags(&selected)
}

pub(super) fn factor_level_samples<S: AsRef<str>>(
    samples: usize,
    sample_levels: &[S],
    numerator: &str,
    denominator: &str,
) -> Result<Vec<bool>, DeseqError> {
    validate_factor_level_contrast("factor", numerator, denominator)?;
    if sample_levels.len() != samples {
        return Err(invalid_dimensions(
            "contrastAllZero sample levels",
            samples,
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
    Ok(selected_samples)
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
