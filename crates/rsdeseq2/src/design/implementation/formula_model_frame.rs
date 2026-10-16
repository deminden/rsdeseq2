// Model-frame construction, validation and alias resolution.

fn validate_formula_factor_name_collisions(
    factors: &[ExpandedFactorSpec<'_>],
) -> Result<(), DeseqError> {
    for (idx, factor) in factors.iter().enumerate() {
        for previous in &factors[..idx] {
            if previous.factor == factor.factor {
                return Err(DeseqError::InvalidOptions {
                    reason: format!("formula variable '{}' appears more than once", factor.factor),
                });
            }
            let previous_aliases = r_like_name_candidates(previous.factor);
            let current_aliases = r_like_name_candidates(factor.factor);
            if previous_aliases
                .iter()
                .any(|previous_alias| current_aliases.contains(previous_alias))
            {
                return Err(DeseqError::InvalidOptions {
                    reason: format!(
                        "formula variable '{}' resolves ambiguously after R-style cleanup",
                        factor.factor
                    ),
                });
            }
        }
    }
    Ok(())
}

fn validate_formula_numeric_covariate_name_collisions(
    numeric_covariates: &[ExpandedNumericSpec<'_>],
) -> Result<(), DeseqError> {
    for (idx, covariate) in numeric_covariates.iter().enumerate() {
        for previous in &numeric_covariates[..idx] {
            if previous.name == covariate.name {
                return Err(DeseqError::InvalidOptions {
                    reason: format!("formula variable '{}' appears more than once", covariate.name),
                });
            }
            let previous_aliases = r_like_name_candidates(previous.name);
            let current_aliases = r_like_name_candidates(covariate.name);
            if previous_aliases
                .iter()
                .any(|previous_alias| current_aliases.contains(previous_alias))
            {
                return Err(DeseqError::InvalidOptions {
                    reason: format!(
                        "formula variable '{}' resolves ambiguously after R-style cleanup",
                        covariate.name
                    ),
                });
            }
        }
    }
    Ok(())
}

fn validate_formula_factor_numeric_name_collisions(
    factors: &[ExpandedFactorSpec<'_>],
    numeric_covariates: &[ExpandedNumericSpec<'_>],
) -> Result<(), DeseqError> {
    for factor in factors {
        for covariate in numeric_covariates {
            if factor.factor == covariate.name {
                return Err(DeseqError::InvalidOptions {
                    reason: format!("formula variable '{}' appears more than once", factor.factor),
                });
            }
            let factor_aliases = r_like_name_candidates(factor.factor);
            let covariate_aliases = r_like_name_candidates(covariate.name);
            if factor_aliases
                .iter()
                .any(|factor_alias| covariate_aliases.contains(factor_alias))
            {
                return Err(DeseqError::InvalidOptions {
                    reason: format!(
                        "formula variable '{}' resolves ambiguously after R-style cleanup",
                        covariate.name
                    ),
                });
            }
        }
    }
    Ok(())
}

fn formula_model_frame_from_specs(
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
pub fn formula_has_offset_terms(formula: &str) -> Result<bool, DeseqError> {
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
) -> Result<ExpandedAdditiveFactorDesign, DeseqError> {
    Ok(expanded_formula_design_with_offsets_from_model_frame(formula, model_frame)?.design)
}

/// Build an expanded formula design plus offsets from owned model-frame metadata.
pub fn expanded_formula_design_with_offsets_from_model_frame(
    formula: &str,
    model_frame: &FormulaModelFrame,
) -> Result<ExpandedFormulaDesignWithOffsets, DeseqError> {
    let resolved = ResolvedFormulaModelFrame::new(model_frame)?;
    expanded_formula_design_with_offsets(
        formula,
        &resolved.factor_specs,
        &resolved.numeric_specs,
    )
}

struct ResolvedFormulaModelFrame<'a> {
    factor_specs: Vec<ExpandedFactorSpec<'a>>,
    numeric_specs: Vec<ExpandedNumericSpec<'a>>,
}

impl<'a> ResolvedFormulaModelFrame<'a> {
    fn new(model_frame: &'a FormulaModelFrame) -> Result<Self, DeseqError> {
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

fn validate_formula_model_frame(model_frame: &FormulaModelFrame) -> Result<(), DeseqError> {
    let n_samples = formula_model_frame_sample_count(model_frame)?;
    for (idx, factor) in model_frame.factors.iter().enumerate() {
        validate_formula_model_frame_column_name(&factor.name)?;
        if model_frame.factors[..idx]
            .iter()
            .any(|previous| previous.name == factor.name)
        {
            return Err(DeseqError::InvalidOptions {
                reason: format!("formula factor '{}' appears more than once", factor.name),
            });
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
            return Err(DeseqError::InvalidOptions {
                reason: format!(
                    "formula numeric covariate '{}' appears more than once",
                    covariate.name
                ),
            });
        }
        if model_frame
            .factors
            .iter()
            .any(|factor| factor.name == covariate.name)
        {
            return Err(DeseqError::InvalidOptions {
                reason: format!(
                    "formula column '{}' cannot be both a factor and numeric covariate",
                    covariate.name
                ),
            });
        }
        if model_frame.factors.iter().any(|factor| {
            let factor_aliases = r_like_name_candidates(&factor.name);
            let covariate_aliases = r_like_name_candidates(&covariate.name);
            factor_aliases
                .iter()
                .any(|factor_alias| covariate_aliases.contains(factor_alias))
        }) {
            return Err(DeseqError::InvalidOptions {
                reason: format!(
                    "formula column '{}' resolves ambiguously after R-style cleanup",
                    covariate.name
                ),
            });
        }
        if covariate.values.len() != n_samples {
            return Err(invalid_dimensions(
                "formula numeric covariate values",
                n_samples,
                covariate.values.len(),
            ));
        }
        if !covariate.values.iter().all(|value| value.is_finite()) {
            return Err(DeseqError::InvalidOptions {
                reason: format!(
                    "formula numeric covariate '{}' contains non-finite values",
                    covariate.name
                ),
            });
        }
    }
    Ok(())
}

fn formula_model_frame_factor_reference(
    factor: &FormulaFactorColumn,
) -> Result<&str, DeseqError> {
    if let Some(reference) = factor.reference.as_deref() {
        return resolve_formula_model_frame_reference_alias(factor, reference);
    }
    factor
        .levels
        .as_deref()
        .and_then(|levels| levels.first().map(String::as_str))
        .or_else(|| factor.sample_levels.first().map(String::as_str))
        .ok_or_else(|| DeseqError::InvalidOptions {
            reason: format!(
                "formula factor '{}' requires at least one sample level",
                factor.name
            ),
        })
}

fn resolve_formula_model_frame_reference_alias<'a>(
    factor: &'a FormulaFactorColumn,
    requested: &str,
) -> Result<&'a str, DeseqError> {
    let levels = formula_model_frame_factor_level_candidates(factor);
    resolve_formula_reference_alias_from_levels(&factor.name, &levels, requested)
}

fn formula_model_frame_factor_level_candidates(factor: &FormulaFactorColumn) -> Vec<&str> {
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

fn resolve_formula_reference_alias_from_levels<'a>(
    factor_name: &str,
    levels: &[&'a str],
    requested: &str,
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
                    "formula factor '{factor_name}' reference level '{requested}' appears more than once"
                ),
            });
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
        [] => Err(DeseqError::InvalidOptions {
            reason: format!(
                "formula factor '{factor_name}' does not contain reference level '{requested}'"
            ),
        }),
        _ => Err(DeseqError::InvalidOptions {
            reason: format!(
                "formula factor '{factor_name}' reference level '{requested}' resolves ambiguously after R-style cleanup"
            ),
        }),
    }
}

fn formula_model_frame_sample_count(model_frame: &FormulaModelFrame) -> Result<usize, DeseqError> {
    if let Some(factor) = model_frame.factors.first() {
        if factor.sample_levels.is_empty() {
            return Err(DeseqError::InvalidOptions {
                reason: format!(
                    "formula factor '{}' requires at least one sample level",
                    factor.name
                ),
            });
        }
        return Ok(factor.sample_levels.len());
    }
    if let Some(covariate) = model_frame.numeric_covariates.first() {
        if covariate.values.is_empty() {
            return Err(DeseqError::InvalidOptions {
                reason: format!(
                    "formula numeric covariate '{}' requires at least one sample value",
                    covariate.name
                ),
            });
        }
        return Ok(covariate.values.len());
    }
    Err(DeseqError::InvalidOptions {
        reason: "formula model frame requires at least one column".to_string(),
    })
}
