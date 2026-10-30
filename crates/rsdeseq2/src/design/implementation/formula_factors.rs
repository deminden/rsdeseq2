fn expand_formula_factor_transform_terms(
    rhs: &str,
    factors: &[ExpandedFactorSpec<'_>],
) -> Result<(String, Vec<FormulaDerivedFactorCovariate>), DeseqError> {
    let mut expanded = String::with_capacity(rhs.len());
    let mut derived = Vec::new();
    let mut remainder = rhs;
    while let Some((start, after_open_idx, transform)) = next_formula_factor_transform(remainder) {
        expanded.push_str(&remainder[..start]);
        let after_open = &remainder[after_open_idx..];
        let close = formula_transform_close_index(after_open)?;
        let expression = after_open[..close].trim();
        let factor = formula_factor_transform_term(transform, expression, factors)?;
        let replacement = formula_variable_term(&factor.name);
        if !derived
            .iter()
            .any(|candidate: &FormulaDerivedFactorCovariate| candidate.name == factor.name)
        {
            derived.push(factor);
        }
        expanded.push_str(&replacement);
        remainder = &after_open[close + 1..];
    }
    expanded.push_str(remainder);
    Ok((expanded, derived))
}

#[derive(Clone, Copy, Debug)]
struct FormulaFactorTransform {
    name: &'static str,
    label: &'static str,
}

const FORMULA_FACTOR_TRANSFORMS: [FormulaFactorTransform; 6] = [
    FormulaFactorTransform {
        name: "as.factor",
        label: "as.factor",
    },
    FormulaFactorTransform {
        name: "as.ordered",
        label: "as.ordered",
    },
    FormulaFactorTransform {
        name: "droplevels",
        label: "droplevels",
    },
    FormulaFactorTransform {
        name: "relevel",
        label: "relevel",
    },
    FormulaFactorTransform {
        name: "factor",
        label: "factor",
    },
    FormulaFactorTransform {
        name: "ordered",
        label: "ordered",
    },
];

fn next_formula_factor_transform(rhs: &str) -> Option<(usize, usize, FormulaFactorTransform)> {
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
            for transform in FORMULA_FACTOR_TRANSFORMS {
                if let Some(after_open_idx) = formula_call_after_open(rhs, idx, transform.name) {
                    return Some((idx, after_open_idx, transform));
                }
            }
        }
        idx += character.len_utf8();
    }
    None
}

fn formula_factor_transform_term(
    transform: FormulaFactorTransform,
    expression: &str,
    factors: &[ExpandedFactorSpec<'_>],
) -> Result<FormulaDerivedFactorCovariate, DeseqError> {
    match transform.label {
        "relevel" => formula_relevel_transform_term(expression, factors),
        "droplevels" => formula_factor_identity_transform_term(transform, expression, factors),
        "factor" | "as.factor" | "ordered" | "as.ordered" => {
            formula_factor_identity_transform_term(transform, expression, factors)
        }
        _ => unreachable!("validated formula factor transform"),
    }
}

fn formula_factor_transform_accepts_levels(transform: FormulaFactorTransform) -> bool {
    matches!(transform.label, "factor" | "ordered")
}

fn formula_factor_identity_transform_term(
    transform: FormulaFactorTransform,
    expression: &str,
    factors: &[ExpandedFactorSpec<'_>],
) -> Result<FormulaDerivedFactorCovariate, DeseqError> {
    let arguments = split_formula_transform_arguments(expression)?;
    if arguments.is_empty() {
        return Err(DeseqError::InvalidOptions {
            reason: format!(
                "formula transform '{}({expression})' must provide a factor",
                transform.label
            ),
        });
    }
    let mut factor_text = None;
    let mut levels = None;
    let mut labels = None;
    let mut positional_idx = 0_usize;
    for argument in &arguments {
        let trimmed = argument.trim();
        if let Some((key, value)) = split_formula_named_argument(trimmed) {
            match key.as_str() {
                "x" => {
                    reject_duplicate_formula_transform_argument(
                        factor_text.is_some(),
                        transform.label,
                        "x",
                        expression,
                    )?;
                    factor_text = Some(value.to_string());
                    continue;
                }
                "levels" if formula_factor_transform_accepts_levels(transform) => {
                    reject_duplicate_formula_transform_argument(
                        levels.is_some(),
                        transform.label,
                        "levels",
                        expression,
                    )?;
                    levels = Some(parse_formula_string_vector_argument(
                        value,
                        transform.label,
                        "levels",
                        expression,
                    )?);
                    continue;
                }
                "labels" if formula_factor_transform_accepts_levels(transform) => {
                    reject_duplicate_formula_transform_argument(
                        labels.is_some(),
                        transform.label,
                        "labels",
                        expression,
                    )?;
                    labels = Some(parse_formula_string_vector_argument(
                        value,
                        transform.label,
                        "labels",
                        expression,
                    )?);
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
        match positional_idx {
            0 if factor_text.is_none() => {
                factor_text = Some(trimmed.to_string());
                positional_idx += 1;
                continue;
            }
            1 if formula_factor_transform_accepts_levels(transform) && levels.is_none() => {
                levels = Some(parse_formula_string_vector_argument(
                    trimmed,
                    transform.label,
                    "levels",
                    expression,
                )?);
                positional_idx += 1;
                continue;
            }
            2 if formula_factor_transform_accepts_levels(transform) && labels.is_none() => {
                labels = Some(parse_formula_string_vector_argument(
                    trimmed,
                    transform.label,
                    "labels",
                    expression,
                )?);
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
    let Some(factor_text) = factor_text else {
        return Err(DeseqError::InvalidOptions {
            reason: format!(
                "formula transform '{}({expression})' must provide a factor",
                transform.label
            ),
        });
    };
    let factor = formula_factor_argument_covariate(factor_text.trim(), factors)?;
    let mut arguments = vec![(Some("x"), factor.primitive.clone())];
    let levels_label = levels
        .as_ref()
        .map(|v| format!(", levels = {}", formula_string_vector_label(v)))
        .unwrap_or_default();
    let labels_label = labels
        .as_ref()
        .map(|v| format!(", labels = {}", formula_string_vector_label(v)))
        .unwrap_or_default();
    if let Some(levels) = levels {
        arguments.push((Some("levels"), formula_text_value(&levels)));
    }
    if let Some(labels) = labels {
        arguments.push((Some("labels"), formula_text_value(&labels)));
    }
    let primitive = formula_call_value(
        transform.name,
        arguments,
        &format!("{}({expression})", transform.label),
    )?;
    let reference = (transform.name == "as.factor").then_some(factor.reference.as_str());
    formula_factor_result(
        format!(
            "{}({}{levels_label}{labels_label})",
            transform.label, factor.name
        ),
        primitive,
        reference,
    )
}

fn formula_relevel_transform_term(
    expression: &str,
    factors: &[ExpandedFactorSpec<'_>],
) -> Result<FormulaDerivedFactorCovariate, DeseqError> {
    let arguments = split_formula_transform_arguments(expression)?;
    if arguments.len() < 2 {
        return Err(DeseqError::InvalidOptions {
            reason: format!(
                "formula transform 'relevel({expression})' must provide factor and ref"
            ),
        });
    }
    let mut factor_text = None;
    let mut reference = None;
    let mut positional_idx = 0_usize;
    for argument in &arguments {
        let trimmed = argument.trim();
        if let Some((key, value)) = split_formula_named_argument(trimmed) {
            match key.as_str() {
                "x" => {
                    reject_duplicate_formula_transform_argument(
                        factor_text.is_some(),
                        "relevel",
                        "x",
                        expression,
                    )?;
                    factor_text = Some(value.to_string());
                    continue;
                }
                "ref" => {
                    reject_duplicate_formula_transform_argument(
                        reference.is_some(),
                        "relevel",
                        "ref",
                        expression,
                    )?;
                    reference = Some(parse_formula_string_argument(
                        value, "relevel", "ref", expression,
                    )?);
                    continue;
                }
                _ => {
                    return Err(DeseqError::InvalidOptions {
                        reason: format!(
                            "formula transform 'relevel({expression})' has unsupported argument '{argument}'"
                        ),
                    });
                }
            }
        }
        match positional_idx {
            0 if factor_text.is_none() => {
                factor_text = Some(trimmed.to_string());
                positional_idx += 1;
                continue;
            }
            1 if reference.is_none() => {
                reference = Some(parse_formula_string_argument(
                    trimmed, "relevel", "ref", expression,
                )?);
                positional_idx += 1;
                continue;
            }
            _ => {}
        }
        return Err(DeseqError::InvalidOptions {
            reason: format!(
                "formula transform 'relevel({expression})' has unsupported argument '{argument}'"
            ),
        });
    }
    let Some(factor_text) = factor_text else {
        return Err(DeseqError::InvalidOptions {
            reason: format!("formula transform 'relevel({expression})' must provide a factor"),
        });
    };
    let Some(reference) = reference else {
        return Err(DeseqError::InvalidOptions {
            reason: format!("formula transform 'relevel({expression})' must provide ref"),
        });
    };
    let factor = formula_factor_argument_covariate(factor_text.trim(), factors)?;
    let reference = resolve_formula_factor_reference_alias(&factor, &reference)?;
    let primitive = formula_call_value(
        "relevel",
        vec![
            (Some("x"), factor.primitive.clone()),
            (
                Some("ref"),
                formula_text_value(std::slice::from_ref(&reference)),
            ),
        ],
        &format!("relevel({expression})"),
    )?;
    formula_factor_result(
        format!("relevel({}, ref = \"{reference}\")", factor.name),
        primitive,
        Some(&reference),
    )
}

fn resolve_formula_factor_reference_alias(
    factor: &FormulaDerivedFactorCovariate,
    requested: &str,
) -> Result<String, DeseqError> {
    let levels = formula_factor_level_candidates(factor);
    let borrowed_levels = levels.iter().map(String::as_str).collect::<Vec<_>>();
    resolve_formula_reference_alias_from_levels(&factor.name, &borrowed_levels, requested)
        .map(str::to_string)
}

fn formula_factor_level_candidates(factor: &FormulaDerivedFactorCovariate) -> Vec<String> {
    let mut levels = Vec::new();
    if let Some(declared) = &factor.levels {
        for level in declared {
            if !levels.iter().any(|candidate| candidate == level) {
                levels.push(level.clone());
            }
        }
    }
    for level in &factor.sample_levels {
        if !levels.iter().any(|candidate| candidate == level) {
            levels.push(level.clone());
        }
    }
    levels
}

fn formula_factor_argument_covariate(
    expression: &str,
    factors: &[ExpandedFactorSpec<'_>],
) -> Result<FormulaDerivedFactorCovariate, DeseqError> {
    if let Ok(factor_name) = formula_variable_name(expression) {
        let factor = find_formula_factor(factor_name, factors)?;
        let levels = formula_factor_code_levels(factor);
        let primitive = formula_call_value(
            "factor",
            vec![
                (Some("x"), formula_text_value(factor.sample_levels)),
                (Some("levels"), formula_text_value(&levels)),
            ],
            factor.factor,
        )?;
        return formula_factor_result(factor.factor.to_string(), primitive, Some(factor.reference));
    }
    for transform in FORMULA_FACTOR_TRANSFORMS {
        let Some(inner) = formula_transform_inner_expression(expression, transform) else {
            continue;
        };
        return formula_factor_transform_term(transform, inner.trim(), factors);
    }
    Err(DeseqError::InvalidOptions {
        reason: format!(
            "formula factor expression '{expression}' is not a supported factor column or transform"
        ),
    })
}

fn formula_transform_inner_expression(
    expression: &str,
    transform: FormulaFactorTransform,
) -> Option<&str> {
    let trimmed = expression.trim();
    let after_open = formula_call_after_open(trimmed, 0, transform.name)?;
    trimmed[after_open..].strip_suffix(')')
}

fn find_formula_factor<'a>(
    factor_name: &str,
    factors: &'a [ExpandedFactorSpec<'a>],
) -> Result<&'a ExpandedFactorSpec<'a>, DeseqError> {
    let exact = factors
        .iter()
        .filter(|candidate| candidate.factor == factor_name)
        .collect::<Vec<_>>();
    match exact.as_slice() {
        [factor] => return Ok(*factor),
        [] => {}
        _ => {
            return Err(DeseqError::InvalidOptions {
                reason: format!("formula factor '{factor_name}' appears more than once"),
            });
        }
    }

    let aliases = factors
        .iter()
        .filter(|candidate| {
            r_like_name_candidates(candidate.factor)
                .into_iter()
                .any(|candidate| candidate == factor_name)
        })
        .collect::<Vec<_>>();
    match aliases.as_slice() {
        [factor] => Ok(*factor),
        [] => Err(DeseqError::InvalidOptions {
            reason: format!(
                "formula factor '{factor_name}' is not present in supplied design metadata"
            ),
        }),
        _ => Err(DeseqError::InvalidOptions {
            reason: format!(
                "formula factor '{factor_name}' resolves ambiguously after R-style cleanup"
            ),
        }),
    }
}
