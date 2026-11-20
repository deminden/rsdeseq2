use super::*;
use crate::errors::{ErrorKind, FactorLevelRole, NameScope, ambiguous_name, invalid_formula};
// Model-frame construction, validation and alias resolution.

pub(super) fn validate_formula_factor_name_collisions(
    factors: &[ExpandedFactorSpec<'_>],
) -> Result<(), Error> {
    for (idx, factor) in factors.iter().enumerate() {
        for previous in &factors[..idx] {
            validate_formula_variable_name_pair(previous.factor, factor.factor)?;
        }
    }
    Ok(())
}

pub(super) fn validate_formula_numeric_covariate_name_collisions(
    numeric_covariates: &[ExpandedNumericSpec<'_>],
) -> Result<(), Error> {
    for (idx, covariate) in numeric_covariates.iter().enumerate() {
        for previous in &numeric_covariates[..idx] {
            validate_formula_variable_name_pair(previous.name, covariate.name)?;
        }
    }
    Ok(())
}

pub(super) fn validate_formula_factor_numeric_name_collisions(
    factors: &[ExpandedFactorSpec<'_>],
    numeric_covariates: &[ExpandedNumericSpec<'_>],
) -> Result<(), Error> {
    for factor in factors {
        for covariate in numeric_covariates {
            validate_formula_variable_name_pair(factor.factor, covariate.name)?;
        }
    }
    Ok(())
}

fn validate_formula_variable_name_pair(previous: &str, name: &str) -> Result<(), Error> {
    let reason = if previous == name {
        format!("formula variable '{name}' appears more than once")
    } else if name_aliases_overlap(previous, name) {
        format!("formula variable '{name}' resolves ambiguously after R-style cleanup")
    } else {
        return Ok(());
    };
    Err(ambiguous_name(
        ErrorKind::Formula,
        NameScope::FormulaVariable,
        name,
        [previous, name],
        reason,
    ))
}

fn name_aliases_overlap(left: &str, right: &str) -> bool {
    let right_aliases = r_like_name_candidates(right);
    r_like_name_candidates(left)
        .iter()
        .any(|alias| right_aliases.contains(alias))
}

pub(super) fn formula_model_frame_from_specs(
    factors: &[ExpandedFactorSpec<'_>],
    numeric_covariates: &[ExpandedNumericSpec<'_>],
) -> FormulaModelFrame {
    FormulaModelFrame {
        factors: factors
            .iter()
            .map(|factor| FormulaFactorColumn {
                name: factor.factor.to_string(),
                sample_levels: factor.sample_levels.to_vec(),
                levels: factor.levels.map(<[String]>::to_vec),
                reference: Some(factor.reference.to_string()),
            })
            .collect(),
        numeric_covariates: numeric_covariates
            .iter()
            .map(|covariate| FormulaNumericColumn {
                name: covariate.name.to_string(),
                values: covariate.values.to_vec(),
            })
            .collect(),
    }
}

/// Return whether a primitive formula contains at least one `offset(...)` term.
///
/// This syntax-level helper validates the formula shape enough to distinguish
/// real offset terms from ordinary variable names before any model-frame lookup
/// is required.
pub fn formula_has_offset_terms(formula: &str) -> Result<bool, Error> {
    let rhs = formula_rhs(formula)?;
    for (_, term) in split_formula_signed_terms(rhs)? {
        if formula_term_has_offset(&term)? {
            return Ok(true);
        }
    }
    Ok(false)
}

/// Build an expanded formula design from owned model-frame metadata.
///
/// This is the wrapper/object-facing companion to [`expanded_formula_design`].
/// It derives borrowed factor/numeric specs from an owned
/// [`FormulaModelFrame`], inferring each factor reference from the first
/// declared factor level when available, otherwise the first observed sample
/// level, when the caller did not supply one.
pub fn expanded_formula_design_from_model_frame(
    formula: &str,
    model_frame: &FormulaModelFrame,
) -> Result<ExpandedAdditiveFactorDesign, Error> {
    Ok(expanded_formula_design_with_offsets_from_model_frame(formula, model_frame)?.design)
}

/// Build an expanded formula design plus offsets from owned model-frame metadata.
pub fn expanded_formula_design_with_offsets_from_model_frame(
    formula: &str,
    model_frame: &FormulaModelFrame,
) -> Result<ExpandedFormulaDesignWithOffsets, Error> {
    let resolved = ResolvedFormulaModelFrame::new(model_frame)?;
    expanded_formula_design_with_offsets(formula, &resolved.factor_specs, &resolved.numeric_specs)
}

pub(super) struct ResolvedFormulaModelFrame<'a> {
    pub(super) factor_specs: Vec<ExpandedFactorSpec<'a>>,
    pub(super) numeric_specs: Vec<ExpandedNumericSpec<'a>>,
}

impl<'a> ResolvedFormulaModelFrame<'a> {
    // Invariant: factor references are validated by FormulaModelFrame::validate before specs are borrowed.
    #[allow(clippy::unwrap_used, clippy::expect_used)]
    pub(super) fn new(model_frame: &'a FormulaModelFrame) -> Result<Self, Error> {
        validate_formula_model_frame(model_frame)?;
        let factor_specs = model_frame
            .factors
            .iter()
            .map(|factor| {
                let reference = formula_model_frame_factor_reference(factor)
                    .expect("model-frame factor reference validated before borrowing specs");
                ExpandedFactorSpec {
                    factor: factor.name.as_str(),
                    sample_levels: factor.sample_levels.as_slice(),
                    reference,
                    levels: factor.levels.as_deref(),
                }
            })
            .collect::<Vec<_>>();
        let numeric_specs = model_frame
            .numeric_covariates
            .iter()
            .map(|covariate| ExpandedNumericSpec {
                name: covariate.name.as_str(),
                values: covariate.values.as_slice(),
            })
            .collect::<Vec<_>>();
        Ok(Self {
            factor_specs,
            numeric_specs,
        })
    }
}

pub(super) fn validate_formula_model_frame(model_frame: &FormulaModelFrame) -> Result<(), Error> {
    let n_samples = formula_model_frame_sample_count(model_frame)?;
    for (idx, factor) in model_frame.factors.iter().enumerate() {
        validate_formula_model_frame_column_name(&factor.name)?;
        if model_frame.factors[..idx]
            .iter()
            .any(|previous| previous.name == factor.name)
        {
            return Err(ambiguous_name(
                ErrorKind::Formula,
                NameScope::Factor,
                &factor.name,
                [factor.name.as_str(); 2],
                format!("formula factor '{}' appears more than once", factor.name),
            ));
        }
        if factor.sample_levels.len() != n_samples {
            return Err(invalid_dimensions(
                "formula factor sample levels",
                n_samples,
                factor.sample_levels.len(),
            ));
        }
        let reference = formula_model_frame_factor_reference(factor)?;
        validate_factor_design_inputs_with_levels(
            &factor.name,
            &factor.sample_levels,
            reference,
            factor.levels.as_deref(),
        )?;
    }
    for (idx, covariate) in model_frame.numeric_covariates.iter().enumerate() {
        validate_formula_model_frame_column_name(&covariate.name)?;
        if model_frame.numeric_covariates[..idx]
            .iter()
            .any(|previous| previous.name == covariate.name)
        {
            return Err(ambiguous_name(
                ErrorKind::Formula,
                NameScope::NumericCovariate,
                &covariate.name,
                [covariate.name.as_str(); 2],
                format!(
                    "formula numeric covariate '{}' appears more than once",
                    covariate.name
                ),
            ));
        }
        if model_frame
            .factors
            .iter()
            .any(|factor| factor.name == covariate.name)
        {
            return Err(ambiguous_name(
                ErrorKind::Formula,
                NameScope::FormulaVariable,
                &covariate.name,
                [covariate.name.as_str(); 2],
                format!(
                    "formula column '{}' cannot be both a factor and numeric covariate",
                    covariate.name
                ),
            ));
        }
        if let Some(factor) = model_frame
            .factors
            .iter()
            .find(|factor| name_aliases_overlap(&factor.name, &covariate.name))
        {
            return Err(ambiguous_name(
                ErrorKind::Formula,
                NameScope::FormulaVariable,
                &covariate.name,
                [factor.name.as_str(), covariate.name.as_str()],
                format!(
                    "formula column '{}' resolves ambiguously after R-style cleanup",
                    covariate.name
                ),
            ));
        }
        if covariate.values.len() != n_samples {
            return Err(invalid_dimensions(
                "formula numeric covariate values",
                n_samples,
                covariate.values.len(),
            ));
        }
        if !covariate.values.iter().all(|value| value.is_finite()) {
            return Err(invalid_formula(format!(
                "formula numeric covariate '{}' contains non-finite values",
                covariate.name
            )));
        }
    }
    Ok(())
}

pub(super) fn formula_model_frame_factor_reference(
    factor: &FormulaFactorColumn,
) -> Result<&str, Error> {
    if let Some(reference) = factor.reference.as_deref() {
        return resolve_formula_model_frame_reference_alias(factor, reference);
    }
    factor
        .levels
        .as_deref()
        .and_then(|levels| levels.first().map(String::as_str))
        .or_else(|| factor.sample_levels.first().map(String::as_str))
        .ok_or_else(|| {
            invalid_formula(format!(
                "formula factor '{}' requires at least one sample level",
                factor.name
            ))
        })
}

pub(super) fn resolve_formula_model_frame_reference_alias<'a>(
    factor: &'a FormulaFactorColumn,
    requested: &str,
) -> Result<&'a str, Error> {
    let levels = formula_model_frame_factor_level_candidates(factor);
    resolve_formula_reference_alias_from_levels(&factor.name, &levels, requested)
}

pub(super) fn formula_model_frame_factor_level_candidates(
    factor: &FormulaFactorColumn,
) -> Vec<&str> {
    let mut levels = Vec::new();
    if let Some(declared) = &factor.levels {
        for level in declared {
            if !levels.iter().any(|candidate| candidate == &level.as_str()) {
                levels.push(level.as_str());
            }
        }
    }
    for level in &factor.sample_levels {
        if !levels.iter().any(|candidate| candidate == &level.as_str()) {
            levels.push(level.as_str());
        }
    }
    levels
}

pub(super) fn resolve_formula_reference_alias_from_levels<'a>(
    factor_name: &str,
    levels: &[&'a str],
    requested: &str,
) -> Result<&'a str, Error> {
    let scope = || NameScope::FactorLevel {
        factor: Some(factor_name.into()),
        role: FactorLevelRole::Reference,
    };
    let exact = levels
        .iter()
        .copied()
        .filter(|level| *level == requested)
        .collect::<Vec<_>>();
    match exact.as_slice() {
        [level] => return Ok(*level),
        [] => {}
        _ => {
            return Err(ambiguous_name(
                ErrorKind::Formula,
                scope(),
                requested,
                exact,
                format!(
                    "formula factor '{factor_name}' reference level '{requested}' appears more than once"
                ),
            ));
        }
    }

    let aliases = levels
        .iter()
        .copied()
        .filter(|level| {
            r_like_name_candidates(level)
                .into_iter()
                .any(|candidate| candidate == requested)
        })
        .collect::<Vec<_>>();
    match aliases.as_slice() {
        [level] => Ok(*level),
        [] => Err(Error::NameNotFound {
            domain: ErrorKind::Formula,
            scope: scope(),
            name: requested.into(),
            reason: format!(
                "formula factor '{factor_name}' does not contain reference level '{requested}'"
            ),
        }),
        _ => Err(ambiguous_name(
            ErrorKind::Formula,
            scope(),
            requested,
            aliases,
            format!(
                "formula factor '{factor_name}' reference level '{requested}' resolves ambiguously after R-style cleanup"
            ),
        )),
    }
}

pub(super) fn formula_model_frame_sample_count(
    model_frame: &FormulaModelFrame,
) -> Result<usize, Error> {
    if let Some(factor) = model_frame.factors.first() {
        if factor.sample_levels.is_empty() {
            return Err(invalid_formula(format!(
                "formula factor '{}' requires at least one sample level",
                factor.name
            )));
        }
        return Ok(factor.sample_levels.len());
    }
    if let Some(covariate) = model_frame.numeric_covariates.first() {
        if covariate.values.is_empty() {
            return Err(invalid_formula(format!(
                "formula numeric covariate '{}' requires at least one sample value",
                covariate.name
            )));
        }
        return Ok(covariate.values.len());
    }
    Err(invalid_formula(
        "formula model frame requires at least one column",
    ))
}
