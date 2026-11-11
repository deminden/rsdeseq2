// Numeric transform discovery, parsing and covariate resolution.

fn expand_formula_numeric_transform_terms(
    rhs: &str,
    factors: &[ExpandedFactorSpec<'_>],
    numeric_covariates: &[ExpandedNumericSpec<'_>],
    fitting_order: bool,
) -> Result<(String, Vec<FormulaDerivedNumericCovariate>), DeseqError> {
    let mut expanded = String::with_capacity(rhs.len());
    let mut derived = Vec::new();
    let mut groups = Vec::new();
    let mut remainder = rhs;
    while let Some((start, after_open_idx, transform)) = next_formula_numeric_transform(remainder) {
        expanded.push_str(&remainder[..start]);
        let after_open = &remainder[after_open_idx..];
        let close = formula_transform_close_index(after_open)?;
        let expression = after_open[..close].trim();
        let (replacement, transformed_covariates) =
            formula_numeric_transform_term(transform, expression, factors, numeric_covariates)?;
        if fitting_order {
            let key = formula_offset_key(&format!("{}({expression})", transform.name))?;
            let group = groups.iter().position(|g| *g == key).unwrap_or_else(|| {
                groups.push(key);
                groups.len() - 1
            });
            let mut names = derived
                .iter()
                .filter(|d: &&FormulaDerivedNumericCovariate| d.group == Some(group))
                .map(|d| d.name.clone())
                .collect::<Vec<_>>();
            if names.is_empty() {
                // Different language expressions may have the same compatibility
                // alias. Keep both variables; only the completed design can
                // decide rank (subtraction and interactions may make it valid).
                for (alias, values) in transformed_covariates {
                    let mut name = alias.clone();
                    let mut suffix = 0;
                    while derived.iter().any(|d| d.name == name)
                        || numeric_covariates.iter().any(|v| v.name == name)
                        || factors.iter().any(|f| f.factor == name)
                    {
                        suffix += 1;
                        name = format!("{alias}.{suffix}");
                    }
                    names.push(name.clone());
                    derived.push(FormulaDerivedNumericCovariate {
                        name,
                        values,
                        group: Some(group),
                    });
                }
            }
            let terms = names
                .iter()
                .map(|name| formula_variable_term(name))
                .collect::<Vec<_>>();
            if terms.len() == 1 {
                expanded.push_str(&terms[0]);
            } else {
                expanded.push_str(&format!("({})", terms.join(" + ")));
            }
        } else {
            for (name, values) in transformed_covariates {
                if !derived.iter().any(|d| d.name == name) {
                    derived.push(FormulaDerivedNumericCovariate {
                        name,
                        values,
                        group: None,
                    });
                }
            }
            expanded.push_str(&replacement);
        }
        remainder = &after_open[close + 1..];
    }
    expanded.push_str(remainder);
    Ok((expanded, derived))
}

fn formula_transform_close_index(after_open: &str) -> Result<usize, DeseqError> {
    let mut depth = 0_i32;
    let mut in_backticks = false;
    let mut quote = None;
    let mut escaped_quote = false;
    for (idx, character) in after_open.char_indices() {
        if let Some(active_quote) = quote {
            if escaped_quote {
                escaped_quote = false;
            } else if character == '\\' {
                escaped_quote = true;
            } else if character == active_quote {
                quote = None;
            }
            continue;
        }
        if matches!(character, '"' | '\'') && !in_backticks {
            quote = Some(character);
            continue;
        }
        if character == '`' {
            in_backticks = !in_backticks;
            continue;
        }
        if in_backticks {
            continue;
        }
        match character {
            '(' => depth += 1,
            ')' if depth == 0 => return Ok(idx),
            ')' => depth -= 1,
            _ => {}
        }
    }
    if in_backticks {
        return Err(DeseqError::InvalidOptions {
            reason: "formula transform has unbalanced backtick-quoted variable name".to_string(),
        });
    }
    if quote.is_some() {
        return Err(DeseqError::InvalidOptions {
            reason: "formula transform has unbalanced quoted string".to_string(),
        });
    }
    Err(DeseqError::InvalidOptions {
        reason: "formula transform has unbalanced parentheses".to_string(),
    })
}

#[derive(Clone, Copy, Debug)]
struct FormulaNumericTransform {
    name: &'static str,
    label: &'static str,
}

const FORMULA_NUMERIC_TRANSFORMS: [FormulaNumericTransform; 11] = [
    FormulaNumericTransform {
        name: "poly",
        label: "poly",
    },
    FormulaNumericTransform {
        name: "scale",
        label: "scale",
    },
    FormulaNumericTransform {
        name: "log10",
        label: "log10",
    },
    FormulaNumericTransform {
        name: "log1p",
        label: "log1p",
    },
    FormulaNumericTransform {
        name: "log2",
        label: "log2",
    },
    FormulaNumericTransform {
        name: "sqrt",
        label: "sqrt",
    },
    FormulaNumericTransform {
        name: "as.numeric",
        label: "as.numeric",
    },
    FormulaNumericTransform {
        name: "as.double",
        label: "as.double",
    },
    FormulaNumericTransform {
        name: "as.integer",
        label: "as.integer",
    },
    FormulaNumericTransform {
        name: "log",
        label: "log",
    },
    FormulaNumericTransform {
        name: "I",
        label: "I",
    },
];

fn next_formula_numeric_transform(rhs: &str) -> Option<(usize, usize, FormulaNumericTransform)> {
    let mut in_backticks = false;
    let mut quote = None;
    let mut escaped_quote = false;
    let mut idx = 0_usize;
    while idx < rhs.len() {
        let character = rhs[idx..]
            .chars()
            .next()
            .expect("idx is inside rhs char boundary");
        if let Some(active_quote) = quote {
            if escaped_quote {
                escaped_quote = false;
            } else if character == '\\' {
                escaped_quote = true;
            } else if character == active_quote {
                quote = None;
            }
            idx += character.len_utf8();
            continue;
        }
        if matches!(character, '"' | '\'') && !in_backticks {
            quote = Some(character);
            idx += character.len_utf8();
            continue;
        }
        if character == '`' {
            in_backticks = !in_backticks;
            idx += character.len_utf8();
            continue;
        }
        if !in_backticks {
            for transform in FORMULA_NUMERIC_TRANSFORMS {
                if let Some(after_open_idx) = formula_call_after_open(rhs, idx, transform.name) {
                    return Some((idx, after_open_idx, transform));
                }
            }
        }
        idx += character.len_utf8();
    }
    None
}

fn formula_numeric_transform_inner_expression(
    expression: &str,
    transform: FormulaNumericTransform,
) -> Option<&str> {
    let trimmed = expression.trim();
    let after_open = formula_call_after_open(trimmed, 0, transform.name)?;
    trimmed[after_open..].strip_suffix(')')
}

fn formula_call_after_open(text: &str, idx: usize, name: &str) -> Option<usize> {
    if !formula_call_has_left_boundary(text, idx) || !text[idx..].starts_with(name) {
        return None;
    }
    let mut after_name = idx + name.len();
    while after_name < text.len() {
        let character = text[after_name..]
            .chars()
            .next()
            .expect("after_name is inside text char boundary");
        if !character.is_whitespace() {
            break;
        }
        after_name += character.len_utf8();
    }
    if text[after_name..].starts_with('(') {
        Some(after_name + 1)
    } else {
        None
    }
}

fn formula_call_has_left_boundary(text: &str, idx: usize) -> bool {
    if idx == 0 {
        return true;
    }
    let Some(previous) = text[..idx].chars().next_back() else {
        return true;
    };
    !(previous.is_alphanumeric() || previous == '_' || previous == '.')
}

fn formula_numeric_transform_term(
    transform: FormulaNumericTransform,
    expression: &str,
    factors: &[ExpandedFactorSpec<'_>],
    numeric_covariates: &[ExpandedNumericSpec<'_>],
) -> Result<FormulaNumericTransformExpansion, DeseqError> {
    match transform.label {
        "I" => {
            let (name, values) =
                formula_numeric_identity_or_power_term(expression, numeric_covariates)?;
            Ok((formula_variable_term(&name), vec![(name, values)]))
        }
        "poly" => formula_numeric_poly_term(expression, numeric_covariates),
        "scale" => {
            let (name, values) = formula_numeric_scale_term(expression, numeric_covariates)?;
            Ok((formula_variable_term(&name), vec![(name, values)]))
        }
        _ => {
            let (name, values) =
                formula_numeric_function_term(transform, expression, factors, numeric_covariates)?;
            Ok((formula_variable_term(&name), vec![(name, values)]))
        }
    }
}

fn formula_numeric_function_term(
    transform: FormulaNumericTransform,
    expression: &str,
    factors: &[ExpandedFactorSpec<'_>],
    numeric_covariates: &[ExpandedNumericSpec<'_>],
) -> Result<(String, Vec<f64>), DeseqError> {
    let (numeric_text, base) = formula_numeric_function_arguments(transform, expression)?;
    let numeric_name = formula_variable_name(numeric_text.trim())?;
    if matches!(transform.label, "as.numeric" | "as.double" | "as.integer") {
        if let Some(covariate) =
            maybe_formula_numeric_covariate(Some(numeric_name), numeric_covariates)?
        {
            return formula_numeric_covariate_function_values(
                transform, expression, covariate, base,
            );
        }
        return formula_factor_numeric_coercion_term(transform, expression, numeric_name, factors);
    }
    let covariate = find_formula_numeric_covariate(numeric_name, numeric_covariates)?;
    formula_numeric_covariate_function_values(transform, expression, covariate, base)
}

fn formula_numeric_covariate_function_values(
    transform: FormulaNumericTransform,
    expression: &str,
    covariate: &ExpandedNumericSpec<'_>,
    base: Option<f64>,
) -> Result<(String, Vec<f64>), DeseqError> {
    let numeric_name = covariate.name;
    let mut arguments = vec![(Some("x"), formula_numeric_value(covariate.values))];
    if let Some(base) = base {
        arguments.push((Some("base"), super::expression::Value::scalar(base)));
    }
    let values = formula_call_values(
        transform.name,
        arguments,
        &format!("{}({expression})", transform.label),
    )?;
    let name = if let Some(base) = base {
        format!("{numeric_name}_log_base_{}", formula_numeric_label(base))
    } else {
        format!("{}_{}", numeric_name, transform.label)
    };
    Ok((name, values))
}

fn formula_factor_numeric_coercion_term(
    transform: FormulaNumericTransform,
    expression: &str,
    factor_name: &str,
    factors: &[ExpandedFactorSpec<'_>],
) -> Result<(String, Vec<f64>), DeseqError> {
    let factor = find_formula_factor(factor_name, factors)?;
    let levels = formula_factor_code_levels(factor);
    use super::expression::{Argument, Data, Expression, Value};
    let text = |v: &[String]| {
        Expression::Literal(Value::new(Data::Text(
            v.iter().cloned().map(Some).collect(),
        )))
    };
    let input = Expression::Call {
        function: "factor".into(),
        arguments: vec![
            Argument {
                name: Some("x".into()),
                expression: text(factor.sample_levels),
            },
            Argument {
                name: Some("levels".into()),
                expression: text(&levels),
            },
        ],
    };
    let expression_ast = Expression::Call {
        function: transform.name.into(),
        arguments: vec![Argument {
            name: None,
            expression: input,
        }],
    };
    let values = formula_expression_values(
        &expression_ast,
        &format!("{}({expression})", transform.label),
    )?;
    Ok((format!("{}_{}", factor.factor, transform.label), values))
}

fn formula_factor_code_levels(factor: &ExpandedFactorSpec<'_>) -> Vec<String> {
    if let Some(levels) = factor.levels {
        return levels.to_vec();
    }
    let mut levels = Vec::new();
    for level in factor.sample_levels {
        if !levels.iter().any(|candidate| candidate == level) {
            levels.push(level.clone());
        }
    }
    levels
}

fn formula_numeric_function_arguments(
    transform: FormulaNumericTransform,
    expression: &str,
) -> Result<(String, Option<f64>), DeseqError> {
    let arguments = split_formula_transform_arguments(expression)?;
    let mut numeric_text = None;
    let mut base = None;
    let mut positional_idx = 0_usize;
    for argument in &arguments {
        let trimmed = argument.trim();
        if let Some((key, value)) = split_formula_named_argument(trimmed) {
            match key.as_str() {
                "x" => {
                    reject_duplicate_formula_transform_argument(
                        numeric_text.is_some(),
                        transform.label,
                        "x",
                        expression,
                    )?;
                    numeric_text = Some(value.to_string());
                    continue;
                }
                "base" if transform.label == "log" => {
                    reject_duplicate_formula_transform_argument(
                        base.is_some(),
                        transform.label,
                        "base",
                        expression,
                    )?;
                    base = Some(parse_formula_log_base(value, expression)?);
                    continue;
                }
                _ => {
                    return Err(DeseqError::InvalidOptions {
                        reason: format!(
                            "formula transform '{}({expression})' has unsupported argument '{argument}'",
                            transform.label
                        ),
                    });
                }
            }
        }
        match positional_idx {
            0 if numeric_text.is_none() => {
                numeric_text = Some(trimmed.to_string());
                positional_idx += 1;
                continue;
            }
            1 if transform.label == "log" && base.is_none() => {
                base = Some(parse_formula_log_base(trimmed, expression)?);
                positional_idx += 1;
                continue;
            }
            _ => {}
        }
        return Err(DeseqError::InvalidOptions {
            reason: format!(
                "formula transform '{}({expression})' has unsupported argument '{argument}'",
                transform.label
            ),
        });
    }
    let Some(numeric_text) = numeric_text else {
        return Err(DeseqError::InvalidOptions {
            reason: format!(
                "formula transform '{}({expression})' must provide a numeric covariate",
                transform.label
            ),
        });
    };
    Ok((numeric_text, base))
}

fn parse_formula_log_base(value: &str, expression: &str) -> Result<f64, DeseqError> {
    let value = strip_formula_outer_parentheses(value.trim())?;
    let parsed = value
        .parse::<f64>()
        .map_err(|_| DeseqError::InvalidOptions {
            reason: format!(
                "formula transform 'log({expression})' base must be finite and positive"
            ),
        })?;
    if !parsed.is_finite() || parsed <= 0.0 || parsed == 1.0 {
        return Err(DeseqError::InvalidOptions {
            reason: format!(
                "formula transform 'log({expression})' base must be finite and positive"
            ),
        });
    }
    Ok(parsed)
}

fn formula_numeric_label(value: f64) -> String {
    value.to_string().replace('.', "_")
}

fn find_formula_numeric_covariate<'a>(
    numeric_name: &str,
    numeric_covariates: &'a [ExpandedNumericSpec<'a>],
) -> Result<&'a ExpandedNumericSpec<'a>, DeseqError> {
    let exact = numeric_covariates
        .iter()
        .filter(|candidate| candidate.name == numeric_name)
        .collect::<Vec<_>>();
    match exact.as_slice() {
        [covariate] => return Ok(*covariate),
        [] => {}
        _ => {
            return Err(DeseqError::InvalidOptions {
                reason: format!(
                    "formula numeric covariate '{numeric_name}' appears more than once"
                ),
            });
        }
    }

    let aliases = numeric_covariates
        .iter()
        .filter(|candidate| {
            r_like_name_candidates(candidate.name)
                .into_iter()
                .any(|candidate| candidate == numeric_name)
        })
        .collect::<Vec<_>>();
    match aliases.as_slice() {
        [covariate] => Ok(*covariate),
        [] => Err(DeseqError::InvalidOptions {
            reason: format!(
                "formula numeric covariate '{numeric_name}' is not present in supplied design metadata"
            ),
        }),
        _ => Err(DeseqError::InvalidOptions {
            reason: format!(
                "formula numeric covariate '{numeric_name}' resolves ambiguously after R-style cleanup"
            ),
        }),
    }
}

fn maybe_formula_numeric_covariate<'a>(
    numeric_name: Option<&str>,
    numeric_covariates: &'a [ExpandedNumericSpec<'a>],
) -> Result<Option<&'a ExpandedNumericSpec<'a>>, DeseqError> {
    let Some(numeric_name) = numeric_name else {
        return Ok(None);
    };
    let matches_any_covariate = numeric_covariates
        .iter()
        .any(|candidate| candidate.name == numeric_name)
        || numeric_covariates.iter().any(|candidate| {
            r_like_name_candidates(candidate.name)
                .into_iter()
                .any(|candidate| candidate == numeric_name)
        });
    if !matches_any_covariate {
        return Ok(None);
    }
    find_formula_numeric_covariate(numeric_name, numeric_covariates).map(Some)
}
