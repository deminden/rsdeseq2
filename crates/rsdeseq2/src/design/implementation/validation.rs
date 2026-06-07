fn validate_additive_factor_specs(factors: &[ExpandedFactorSpec<'_>]) -> Result<(), DeseqError> {
    if factors.is_empty() {
        return Err(DeseqError::InvalidOptions {
            reason: "at least one factor is required".to_string(),
        });
    }
    let n_samples = factors[0].sample_levels.len();
    for (idx, factor) in factors.iter().enumerate() {
        validate_factor_design_inputs_with_levels(
            factor.factor,
            factor.sample_levels,
            factor.reference,
            factor.levels,
        )?;
        if factor.sample_levels.len() != n_samples {
            return Err(invalid_dimensions(
                "additive factor sample levels",
                n_samples,
                factor.sample_levels.len(),
            ));
        }
        if factors[..idx]
            .iter()
            .any(|previous| previous.factor == factor.factor)
        {
            return Err(DeseqError::InvalidOptions {
                reason: format!("factor '{}' appears more than once", factor.factor),
            });
        }
    }
    Ok(())
}

fn additive_design_sample_count(
    factors: &[ExpandedFactorSpec<'_>],
    numeric_covariates: &[ExpandedNumericSpec<'_>],
) -> Result<usize, DeseqError> {
    if let Some(factor) = factors.first() {
        return Ok(factor.sample_levels.len());
    }
    let Some(covariate) = numeric_covariates.first() else {
        return Err(DeseqError::InvalidOptions {
            reason: "at least one factor or numeric covariate is required".to_string(),
        });
    };
    if covariate.values.is_empty() {
        return Err(DeseqError::InvalidOptions {
            reason: "numeric covariate values must be non-empty".to_string(),
        });
    }
    Ok(covariate.values.len())
}

fn validate_numeric_covariate_specs(
    numeric_covariates: &[ExpandedNumericSpec<'_>],
    n_samples: usize,
) -> Result<(), DeseqError> {
    for (idx, covariate) in numeric_covariates.iter().enumerate() {
        if covariate.name.is_empty() {
            return Err(DeseqError::InvalidOptions {
                reason: format!("numeric covariate {idx} name must be non-empty"),
            });
        }
        if covariate.values.len() != n_samples {
            return Err(invalid_dimensions(
                "additive numeric covariate values",
                n_samples,
                covariate.values.len(),
            ));
        }
        if numeric_covariates[..idx]
            .iter()
            .any(|previous| previous.name == covariate.name)
        {
            return Err(DeseqError::InvalidOptions {
                reason: format!(
                    "numeric covariate '{}' appears more than once",
                    covariate.name
                ),
            });
        }
        if !covariate.values.iter().all(|value| value.is_finite()) {
            return Err(DeseqError::InvalidOptions {
                reason: format!(
                    "numeric covariate '{}' contains non-finite values",
                    covariate.name
                ),
            });
        }
    }
    Ok(())
}

fn resolve_interaction_indices(
    factors: &[ExpandedFactorSpec<'_>],
    interactions: &[ExpandedFactorInteractionSpec<'_>],
) -> Result<Vec<(usize, usize)>, DeseqError> {
    let mut indices = Vec::with_capacity(interactions.len());
    for (idx, interaction) in interactions.iter().enumerate() {
        if interaction.left_factor.is_empty() || interaction.right_factor.is_empty() {
            return Err(DeseqError::InvalidOptions {
                reason: format!("interaction {idx} factor names must be non-empty"),
            });
        }
        if interaction.left_factor == interaction.right_factor {
            return Err(DeseqError::InvalidOptions {
                reason: format!(
                    "interaction {idx} cannot use factor '{}' twice",
                    interaction.left_factor
                ),
            });
        }
        let left = factor_index(factors, interaction.left_factor)?;
        let right = factor_index(factors, interaction.right_factor)?;
        let ordered = if left < right {
            (left, right)
        } else {
            (right, left)
        };
        if indices.contains(&ordered) {
            return Err(DeseqError::InvalidOptions {
                reason: format!(
                    "interaction '{}:{}' appears more than once",
                    interaction.left_factor, interaction.right_factor
                ),
            });
        }
        indices.push(ordered);
    }
    Ok(indices)
}

fn resolve_factor_numeric_interaction_indices(
    factors: &[ExpandedFactorSpec<'_>],
    numeric_covariates: &[ExpandedNumericSpec<'_>],
    interactions: &[ExpandedFactorNumericInteractionSpec<'_>],
) -> Result<Vec<(usize, usize)>, DeseqError> {
    let mut indices = Vec::with_capacity(interactions.len());
    for (idx, interaction) in interactions.iter().enumerate() {
        if interaction.factor.is_empty() || interaction.numeric.is_empty() {
            return Err(DeseqError::InvalidOptions {
                reason: format!("factor-numeric interaction {idx} names must be non-empty"),
            });
        }
        let factor_idx = factor_index(factors, interaction.factor)?;
        let numeric_idx = numeric_index(numeric_covariates, interaction.numeric)?;
        let ordered = (factor_idx, numeric_idx);
        if indices.contains(&ordered) {
            return Err(DeseqError::InvalidOptions {
                reason: format!(
                    "factor-numeric interaction '{}:{}' appears more than once",
                    interaction.factor, interaction.numeric
                ),
            });
        }
        indices.push(ordered);
    }
    Ok(indices)
}

fn resolve_numeric_interaction_indices(
    numeric_covariates: &[ExpandedNumericSpec<'_>],
    interactions: &[ExpandedNumericInteractionSpec<'_>],
) -> Result<Vec<(usize, usize)>, DeseqError> {
    let mut indices = Vec::with_capacity(interactions.len());
    for (idx, interaction) in interactions.iter().enumerate() {
        if interaction.left_numeric.is_empty() || interaction.right_numeric.is_empty() {
            return Err(DeseqError::InvalidOptions {
                reason: format!("numeric interaction {idx} covariate names must be non-empty"),
            });
        }
        if interaction.left_numeric == interaction.right_numeric {
            return Err(DeseqError::InvalidOptions {
                reason: format!(
                    "numeric interaction {idx} cannot use covariate '{}' twice",
                    interaction.left_numeric
                ),
            });
        }
        let left = numeric_index(numeric_covariates, interaction.left_numeric)?;
        let right = numeric_index(numeric_covariates, interaction.right_numeric)?;
        let ordered = if left < right {
            (left, right)
        } else {
            (right, left)
        };
        if indices.contains(&ordered) {
            return Err(DeseqError::InvalidOptions {
                reason: format!(
                    "numeric interaction '{}:{}' appears more than once",
                    interaction.left_numeric, interaction.right_numeric
                ),
            });
        }
        indices.push(ordered);
    }
    Ok(indices)
}

fn factor_index(factors: &[ExpandedFactorSpec<'_>], factor: &str) -> Result<usize, DeseqError> {
    factors
        .iter()
        .position(|candidate| candidate.factor == factor)
        .ok_or_else(|| DeseqError::InvalidOptions {
            reason: format!("interaction factor '{factor}' is not present in additive factors"),
        })
}

fn numeric_index(
    numeric_covariates: &[ExpandedNumericSpec<'_>],
    numeric: &str,
) -> Result<usize, DeseqError> {
    numeric_covariates
        .iter()
        .position(|candidate| candidate.name == numeric)
        .ok_or_else(|| DeseqError::InvalidOptions {
            reason: format!(
                "interaction numeric covariate '{numeric}' is not present in numeric covariates"
            ),
        })
}

fn validate_unique_coefficient_names(names: &[String], context: &str) -> Result<(), DeseqError> {
    for (idx, name) in names.iter().enumerate() {
        if names[..idx].iter().any(|previous| previous == name) {
            return Err(DeseqError::InvalidOptions {
                reason: format!("{context} has duplicate coefficient name '{name}'"),
            });
        }
    }
    Ok(())
}

fn validate_factor_design_inputs<S: AsRef<str>>(
    factor: &str,
    sample_levels: &[S],
    reference: &str,
) -> Result<(), DeseqError> {
    if factor.is_empty() || reference.is_empty() {
        return Err(DeseqError::InvalidOptions {
            reason: "factor and reference level must be non-empty".to_string(),
        });
    }
    if sample_levels.is_empty() {
        return Err(DeseqError::InvalidOptions {
            reason: "sample levels must be non-empty".to_string(),
        });
    }
    let mut has_reference = false;
    for (idx, level) in sample_levels.iter().enumerate() {
        let level = level.as_ref();
        if level.is_empty() {
            return Err(DeseqError::InvalidOptions {
                reason: format!("sample level {idx} must be non-empty"),
            });
        }
        has_reference |= level == reference;
    }
    if !has_reference {
        return Err(DeseqError::InvalidOptions {
            reason: format!("reference level '{reference}' is not present in sample levels"),
        });
    }
    Ok(())
}

fn validate_factor_design_inputs_with_levels<S: AsRef<str>>(
    factor: &str,
    sample_levels: &[S],
    reference: &str,
    levels: Option<&[String]>,
) -> Result<(), DeseqError> {
    let Some(levels) = levels else {
        validate_factor_design_inputs(factor, sample_levels, reference)?;
        return Ok(());
    };
    if factor.is_empty() || reference.is_empty() {
        return Err(DeseqError::InvalidOptions {
            reason: "factor and reference level must be non-empty".to_string(),
        });
    }
    if sample_levels.is_empty() {
        return Err(DeseqError::InvalidOptions {
            reason: "sample levels must be non-empty".to_string(),
        });
    }
    if levels.is_empty() {
        return Err(DeseqError::InvalidOptions {
            reason: format!("factor '{factor}' declared levels must be non-empty"),
        });
    }
    let mut has_reference = false;
    for (idx, level) in levels.iter().enumerate() {
        if level.is_empty() {
            return Err(DeseqError::InvalidOptions {
                reason: format!("factor '{factor}' declared level {idx} must be non-empty"),
            });
        }
        if levels[..idx].iter().any(|previous| previous == level) {
            return Err(DeseqError::InvalidOptions {
                reason: format!("factor '{factor}' has duplicate declared level '{level}'"),
            });
        }
        has_reference |= level == reference;
    }
    if !has_reference {
        return Err(DeseqError::InvalidOptions {
            reason: format!(
                "reference level '{reference}' is not present in factor '{factor}' declared levels"
            ),
        });
    }
    for (idx, level) in sample_levels.iter().enumerate() {
        let level = level.as_ref();
        if !levels.iter().any(|candidate| candidate == level) {
            return Err(DeseqError::InvalidOptions {
                reason: format!(
                    "sample level {idx} value '{level}' is not present in factor '{factor}' declared levels"
                ),
            });
        }
    }
    Ok(())
}

fn ordered_levels<S: AsRef<str>>(sample_levels: &[S], reference: &str) -> Vec<String> {
    let mut levels = vec![reference.to_string()];
    for level in sample_levels {
        let level = level.as_ref();
        if !levels.iter().any(|candidate| candidate == level) {
            levels.push(level.to_string());
        }
    }
    levels
}

fn ordered_levels_with_declared<S: AsRef<str>>(
    sample_levels: &[S],
    reference: &str,
    levels: Option<&[String]>,
) -> Vec<String> {
    let Some(levels) = levels else {
        return ordered_levels(sample_levels, reference);
    };
    let mut ordered = vec![reference.to_string()];
    for level in levels {
        if level != reference {
            ordered.push(level.clone());
        }
    }
    ordered
}

fn expanded_factor_coefficient_name(factor: &str, level: &str) -> String {
    format!("{factor}{level}")
}

fn standard_factor_coefficient_name(factor: &str, level: &str, reference: &str) -> String {
    format!("{factor}_{level}_vs_{reference}")
}

fn interaction_coefficient_name(
    left_factor: &str,
    left_level: &str,
    right_factor: &str,
    right_level: &str,
) -> String {
    format!("{left_factor}{left_level}:{right_factor}{right_level}")
}

fn standard_interaction_coefficient_name(
    left_factor: &str,
    left_level: &str,
    left_reference: &str,
    right_factor: &str,
    right_level: &str,
    right_reference: &str,
) -> String {
    format!(
        "{left_factor}_{left_level}_vs_{left_reference}:{right_factor}_{right_level}_vs_{right_reference}"
    )
}

fn factor_numeric_interaction_coefficient_name(factor: &str, level: &str, numeric: &str) -> String {
    format!("{factor}{level}:{numeric}")
}

fn standard_factor_numeric_interaction_coefficient_name(
    factor: &str,
    level: &str,
    reference: &str,
    numeric: &str,
) -> String {
    format!("{factor}_{level}_vs_{reference}:{numeric}")
}

fn numeric_interaction_coefficient_name(left_numeric: &str, right_numeric: &str) -> String {
    format!("{left_numeric}:{right_numeric}")
}

pub(crate) fn r_like_name_candidates(raw: &str) -> Vec<String> {
    let made = r_like_make_name(raw);
    if made == raw {
        vec![raw.to_string()]
    } else {
        vec![raw.to_string(), made]
    }
}

pub(crate) fn r_like_make_name(raw: &str) -> String {
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
