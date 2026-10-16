// Formula right-hand sides and offset extraction.

fn formula_rhs(formula: &str) -> Result<&str, DeseqError> {
    let trimmed = formula.trim();
    let Some(rhs) = trimmed.strip_prefix('~') else {
        return Err(DeseqError::InvalidOptions {
            reason: format!("formula '{formula}' must start with '~'"),
        });
    };
    let rhs = rhs.trim();
    if rhs.is_empty() {
        return Err(DeseqError::InvalidOptions {
            reason: "formula right-hand side must be non-empty".to_string(),
        });
    }
    Ok(rhs)
}

fn extract_formula_offsets(
    rhs: &str,
    factors: &[ExpandedFactorSpec<'_>],
    numeric_covariates: &[ExpandedNumericSpec<'_>],
) -> Result<(String, Vec<f64>), DeseqError> {
    let signed_terms = split_formula_signed_terms(rhs)?;
    let n_samples = numeric_covariates
        .first()
        .map(|covariate| covariate.values.len())
        .unwrap_or(0);
    let mut offsets = vec![0.0; n_samples];
    let mut saw_offset = false;
    let mut remaining_terms = Vec::new();
    for (sign, term) in signed_terms {
        extract_formula_offset_term(
            sign,
            &term,
            factors,
            numeric_covariates,
            &mut offsets,
            &mut saw_offset,
            &mut remaining_terms,
        )?;
    }
    if remaining_terms.is_empty() {
        remaining_terms.push("1".to_string());
    }
    if !saw_offset {
        offsets.clear();
    }
    Ok((join_formula_terms(&remaining_terms), offsets))
}

fn formula_term_has_offset(term: &str) -> Result<bool, DeseqError> {
    let term = term.trim();
    if formula_offset_expression(term)?.is_some() {
        return Ok(true);
    }
    let Some(stripped) = formula_outer_group_for_offset(term)? else {
        return Ok(false);
    };
    for (_, expanded_term) in split_formula_signed_terms(stripped)? {
        if formula_term_has_offset(&expanded_term)? {
            return Ok(true);
        }
    }
    Ok(false)
}

fn extract_formula_offset_term(
    sign: i8,
    term: &str,
    factors: &[ExpandedFactorSpec<'_>],
    numeric_covariates: &[ExpandedNumericSpec<'_>],
    offsets: &mut Vec<f64>,
    saw_offset: &mut bool,
    remaining_terms: &mut Vec<String>,
) -> Result<(), DeseqError> {
    let term = term.trim();
    if let Some(stripped) = formula_outer_group_for_offset(term)? {
        for (inner_sign, expanded_term) in split_formula_signed_terms(stripped)? {
            extract_formula_offset_term(
                sign * inner_sign,
                &expanded_term,
                factors,
                numeric_covariates,
                offsets,
                saw_offset,
                remaining_terms,
            )?;
        }
        return Ok(());
    }

    let Some(offset_values) = formula_offset_values(term, factors, numeric_covariates)? else {
        if sign < 0 {
            remaining_terms.push(format!("- {term}"));
        } else {
            remaining_terms.push(term.to_string());
        }
        return Ok(());
    };
    if sign < 0 {
        return Err(DeseqError::InvalidOptions {
            reason: "formula offset terms cannot be subtracted".to_string(),
        });
    }
    *saw_offset = true;
    if offsets.is_empty() {
        offsets.resize(offset_values.len(), 0.0);
    }
    if offset_values.len() != offsets.len() {
        return Err(invalid_dimensions(
            "formula offset values",
            offsets.len(),
            offset_values.len(),
        ));
    }
    for (idx, value) in offset_values.iter().copied().enumerate() {
        if !value.is_finite() {
            return Err(DeseqError::InvalidOptions {
                reason: format!("formula offset '{term}' is non-finite at sample {idx}"),
            });
        }
        offsets[idx] += value;
        if !offsets[idx].is_finite() {
            return Err(DeseqError::InvalidOptions {
                reason: format!("formula offsets sum to a non-finite value at sample {idx}"),
            });
        }
    }
    Ok(())
}

fn formula_outer_group_for_offset(term: &str) -> Result<Option<&str>, DeseqError> {
    if formula_offset_expression(term)?.is_some() {
        return Ok(None);
    }
    if !term.contains("offset") {
        return Ok(None);
    }
    let stripped = strip_formula_outer_parentheses(term)?;
    Ok((stripped != term).then_some(stripped))
}

fn formula_offset_values(
    term: &str,
    factors: &[ExpandedFactorSpec<'_>],
    numeric_covariates: &[ExpandedNumericSpec<'_>],
) -> Result<Option<Vec<f64>>, DeseqError> {
    let Some(expression) = formula_offset_expression(term)? else {
        return Ok(None);
    };
    if let Ok(offset_name) = formula_variable_name(expression) {
        let covariate = numeric_covariates
            .iter()
            .find(|candidate| candidate.name == offset_name)
            .ok_or_else(|| DeseqError::InvalidOptions {
                reason: format!(
                    "formula offset numeric covariate '{offset_name}' is not present in supplied design metadata"
                ),
            })?;
        return Ok(Some(covariate.values.to_vec()));
    }
    for transform in FORMULA_NUMERIC_TRANSFORMS {
        let Some(inner) = formula_numeric_transform_inner_expression(expression, transform) else {
            continue;
        };
        let (_replacement, transformed_covariates) =
            formula_numeric_transform_term(transform, inner.trim(), factors, numeric_covariates)?;
        let [(_name, values)] = transformed_covariates.as_slice() else {
            return Err(DeseqError::InvalidOptions {
                reason: format!("formula offset term '{term}' must produce one numeric vector"),
            });
        };
        return Ok(Some(values.clone()));
    }
    Err(DeseqError::InvalidOptions {
        reason: format!(
            "formula offset expression '{expression}' is not a supported numeric covariate or transform"
        ),
    })
}

fn formula_offset_expression(term: &str) -> Result<Option<&str>, DeseqError> {
    let term = term.trim();
    let Some(after_open) = formula_call_after_open(term, 0, "offset") else {
        if term == "offset" {
            return Err(DeseqError::InvalidOptions {
                reason: "formula offset term 'offset' must be offset(numeric)".to_string(),
            });
        }
        return Ok(None);
    };
    let Some(inner) = term[after_open..].strip_suffix(')') else {
        return Err(DeseqError::InvalidOptions {
            reason: format!("formula offset term '{term}' must be offset(numeric)"),
        });
    };
    let expression = inner.trim();
    if expression.is_empty() {
        return Err(DeseqError::InvalidOptions {
            reason: format!("formula offset term '{term}' must be offset(numeric)"),
        });
    }
    Ok(Some(expression))
}
