//! Coefficient aliases and result labels.

use crate::design::{DesignMatrix, r_like_make_name, r_like_name_candidates as candidate_names};
use crate::errors::DeseqError;
use std::collections::HashSet;

pub(super) fn coefficient_names(design: &DesignMatrix) -> Result<&[String], DeseqError> {
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

pub(super) fn resolve_coefficient_name(
    names: &[String],
    wanted: &str,
) -> Result<usize, DeseqError> {
    if let Ok(index) = find_coefficient(names, wanted) {
        return Ok(index);
    }
    let candidates = coefficient_name_candidates(wanted);
    let matches = candidates
        .iter()
        .flat_map(|candidate| {
            names.iter().enumerate().filter_map(|(idx, name)| {
                coefficient_name_candidates(name)
                    .into_iter()
                    .any(|name_candidate| name_candidate == *candidate)
                    .then_some(idx)
            })
        })
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

pub(super) fn resolve_coefficient_name_list(
    names: &[String],
    wanted: &[String],
) -> Result<Vec<usize>, DeseqError> {
    let mut indices = Vec::with_capacity(wanted.len());
    let mut seen = HashSet::with_capacity(wanted.len());
    for name in wanted {
        let index = resolve_coefficient_name(names, name)?;
        if seen.insert(index) {
            indices.push(index);
        }
    }
    Ok(indices)
}

fn coefficient_name_candidates(name: &str) -> Vec<String> {
    let mut candidates = candidate_names(name);
    if name.contains(':') {
        for candidate in interaction_coefficient_name_candidates(name) {
            push_unique_candidate(&mut candidates, candidate);
        }
    }
    if name == "(Intercept)" {
        push_unique_candidate(&mut candidates, "Intercept".to_string());
    } else if name == "Intercept" {
        push_unique_candidate(&mut candidates, "(Intercept)".to_string());
    }
    if let Some((factor_level, reference)) = name.split_once("_vs_")
        && let Some((factor, level)) = factor_level.rsplit_once('_')
    {
        push_unique_candidate(
            &mut candidates,
            format!(
                "{}_{}_vs_{}",
                r_like_make_name(factor),
                r_like_make_name(level),
                r_like_make_name(reference)
            ),
        );
        push_unique_candidate(
            &mut candidates,
            format!(
                "{}{}_vs_{}",
                r_like_make_name(factor),
                r_like_make_name(level),
                r_like_make_name(reference)
            ),
        );
    }
    candidates
}

fn interaction_coefficient_name_candidates(name: &str) -> Vec<String> {
    let parts = name.split(':').collect::<Vec<_>>();
    if parts.len() < 2 {
        return Vec::new();
    }
    let mut products = vec![String::new()];
    for part in parts {
        let part_candidates = candidate_names(part);
        let mut next = Vec::new();
        for prefix in &products {
            for candidate in &part_candidates {
                if prefix.is_empty() {
                    next.push(candidate.clone());
                } else {
                    next.push(format!("{prefix}:{candidate}"));
                }
            }
        }
        products = next;
    }
    products
}

pub(super) fn list_contrast_comparison(
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

pub(super) fn find_first_coefficient(
    names: &[String],
    candidates: &[String],
) -> Result<usize, DeseqError> {
    let matches = candidates
        .iter()
        .flat_map(|candidate| {
            names.iter().enumerate().filter_map(move |(index, name)| {
                if name == candidate {
                    return Some(index);
                }
                coefficient_name_candidates(name)
                    .into_iter()
                    .any(|name_candidate| name_candidate == *candidate)
                    .then_some(index)
            })
        })
        .collect::<HashSet<_>>();
    match matches.len() {
        0 => Err(DeseqError::InvalidOptions {
            reason: format!(
                "none of the candidate coefficients are present: {}",
                candidates.join(", ")
            ),
        }),
        1 => Ok(*matches.iter().next().unwrap()),
        _ => Err(DeseqError::InvalidOptions {
            reason: format!(
                "candidate coefficients resolve ambiguously after R-style cleanup: {}",
                candidates.join(", ")
            ),
        }),
    }
}

pub(super) fn find_optional_coefficient(
    names: &[String],
    candidates: &[String],
) -> Result<Option<usize>, DeseqError> {
    match find_first_coefficient(names, candidates) {
        Ok(index) => Ok(Some(index)),
        Err(DeseqError::InvalidOptions { reason })
            if reason.starts_with("none of the candidate coefficients are present") =>
        {
            Ok(None)
        }
        Err(err) => Err(err),
    }
}

pub(super) fn find_treatment_or_expanded_coefficient(
    names: &[String],
    treatment_candidates: &[String],
    expanded_candidates: &[String],
) -> Result<usize, DeseqError> {
    if let Some(index) = find_optional_coefficient(names, treatment_candidates)? {
        return Ok(index);
    }
    find_first_coefficient(names, expanded_candidates)
}

pub(super) fn find_shared_reference_standard_coefficients(
    names: &[String],
    factor: &str,
    numerator: &str,
    denominator: &str,
) -> Result<Option<(usize, usize)>, DeseqError> {
    let numerator_pairs = standard_coefficients_for_level(names, factor, numerator);
    let denominator_pairs = standard_coefficients_for_level(names, factor, denominator);
    let mut shared = Vec::new();
    for (numerator_index, numerator_reference) in numerator_pairs.iter() {
        for (denominator_index, denominator_reference) in denominator_pairs.iter() {
            if numerator_reference == denominator_reference {
                shared.push((*numerator_index, *denominator_index));
            }
        }
    }
    shared.sort_unstable();
    shared.dedup();
    match shared.as_slice() {
        [] => Ok(None),
        [(numerator_index, denominator_index)] => Ok(Some((*numerator_index, *denominator_index))),
        _ => Err(DeseqError::InvalidOptions {
            reason: format!(
                "factor-level contrast {factor}: {numerator} vs {denominator} resolves ambiguously through shared-reference coefficient names"
            ),
        }),
    }
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
                if let Some(reference) = name.strip_prefix(&prefix) {
                    return Some((index, reference.to_string()));
                }
                coefficient_name_candidates(name)
                    .into_iter()
                    .find_map(|candidate| {
                        candidate
                            .strip_prefix(&prefix)
                            .map(|reference| (index, reference.to_string()))
                    })
            })
        })
        .collect()
}

fn standard_coefficient_prefixes(factor: &str, level: &str) -> Vec<String> {
    let raw = format!("{factor}_{level}_vs_");
    let mut candidates = candidate_names(&raw);
    candidates.extend(candidate_names(&format!("{factor}{level}_vs_")));
    push_unique_candidate(
        &mut candidates,
        format!(
            "{}_{}_vs_",
            r_like_make_name(factor),
            r_like_make_name(level)
        ),
    );
    push_unique_candidate(
        &mut candidates,
        format!(
            "{}{}_vs_",
            r_like_make_name(factor),
            r_like_make_name(level)
        ),
    );
    candidates
}

pub(super) fn standard_coefficient_names(
    factor: &str,
    level: &str,
    reference: &str,
) -> Vec<String> {
    standard_coefficient_candidates(factor, level, reference)
}

pub(super) fn expanded_coefficient_names(factor: &str, level: &str) -> Vec<String> {
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
    candidates.extend(candidate_names(&format!("{factor}{level}_vs_{reference}")));
    push_unique_candidate(
        &mut candidates,
        format!(
            "{}_{}_vs_{}",
            r_like_make_name(factor),
            r_like_make_name(level),
            r_like_make_name(reference)
        ),
    );
    push_unique_candidate(
        &mut candidates,
        format!(
            "{}{}_vs_{}",
            r_like_make_name(factor),
            r_like_make_name(level),
            r_like_make_name(reference)
        ),
    );
    candidates
}

fn push_unique_candidate(candidates: &mut Vec<String>, candidate: String) {
    if !candidates.iter().any(|existing| existing == &candidate) {
        candidates.push(candidate);
    }
}

pub(super) fn factor_level_result_name(factor: &str, numerator: &str, denominator: &str) -> String {
    format!(
        "{}_{}_vs_{}",
        r_like_make_name(factor),
        r_like_make_name(numerator),
        r_like_make_name(denominator)
    )
}
