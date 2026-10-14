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
        "droplevels" => formula_droplevels_transform_term(transform, expression, factors),
        "factor" | "as.factor" | "ordered" | "as.ordered" => {
            formula_factor_identity_transform_term(transform, expression, factors)
        }
        _ => unreachable!("validated formula factor transform"),
    }
}

fn formula_factor_transform_accepts_levels(transform: FormulaFactorTransform) -> bool {
    matches!(transform.label, "factor" | "ordered")
}

fn formula_droplevels_transform_term(
    transform: FormulaFactorTransform,
    expression: &str,
    factors: &[ExpandedFactorSpec<'_>],
) -> Result<FormulaDerivedFactorCovariate, DeseqError> {
    let mut factor = formula_factor_identity_transform_term(transform, expression, factors)?;
    let dropped_levels = formula_observed_levels_for_factor(&factor);
    let reference = if dropped_levels.iter().any(|level| level == &factor.reference) {
        factor.reference
    } else {
        dropped_levels
            .first()
            .cloned()
            .ok_or_else(|| DeseqError::InvalidOptions {
                reason: format!(
                    "formula transform 'droplevels({expression})' requires at least one observed factor level"
                ),
            })?
    };
    validate_factor_design_inputs_with_levels(
        &factor.name,
        &factor.sample_levels,
        &reference,
        Some(&dropped_levels),
    )?;
    factor.reference = reference;
    factor.levels = Some(dropped_levels);
    Ok(factor)
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
    let mut factor = formula_factor_argument_covariate(factor_text.trim(), factors)?;
    let (levels_label, labels_label) = if let Some(levels) = levels {
        if let Some(labels) = labels {
            if labels.len() != levels.len() {
                return Err(DeseqError::InvalidOptions {
                    reason: format!(
                        "formula transform '{}({expression})' arguments 'levels' and 'labels' must have the same length",
                        transform.label
                    ),
                });
            }
            validate_factor_design_inputs_with_levels(
                &factor.name,
                &factor.sample_levels,
                levels.first().expect("parsed levels are non-empty"),
                Some(&levels),
            )?;
            let label = formula_string_vector_label(&labels);
            factor.sample_levels = relabel_formula_factor_sample_levels(
                &factor.name,
                &factor.sample_levels,
                &levels,
                &labels,
            )?;
            factor.reference = labels
                .first()
                .expect("parsed labels are non-empty")
                .to_string();
            factor.levels = Some(labels);
            (
                format!(", levels = {}", formula_string_vector_label(&levels)),
                format!(", labels = {label}"),
            )
        } else {
            let reference = levels
                .first()
                .ok_or_else(|| DeseqError::InvalidOptions {
                    reason: format!(
                        "formula transform '{}({expression})' argument 'levels' must contain at least one level",
                        transform.label
                    ),
                })?
                .clone();
            validate_factor_design_inputs_with_levels(
                &factor.name,
                &factor.sample_levels,
                &reference,
                Some(&levels),
            )?;
            let label = formula_string_vector_label(&levels);
            factor.reference = reference;
            factor.levels = Some(levels);
            (format!(", levels = {label}"), String::new())
        }
    } else {
        if labels.is_some() {
            return Err(DeseqError::InvalidOptions {
                reason: format!(
                    "formula transform '{}({expression})' argument 'labels' requires explicit levels",
                    transform.label
                ),
            });
        }
        (String::new(), String::new())
    };
    Ok(FormulaDerivedFactorCovariate {
        name: format!(
            "{}({}{levels_label}{labels_label})",
            transform.label, factor.name
        ),
        ..factor
    })
}

fn relabel_formula_factor_sample_levels(
    factor_name: &str,
    sample_levels: &[String],
    levels: &[String],
    labels: &[String],
) -> Result<Vec<String>, DeseqError> {
    sample_levels
        .iter()
        .map(|sample| {
            levels
                .iter()
                .position(|level| level == sample)
                .map(|idx| labels[idx].clone())
                .ok_or_else(|| DeseqError::InvalidOptions {
                    reason: format!(
                        "formula factor '{factor_name}' has sample level '{sample}' outside declared levels"
                    ),
                })
        })
        .collect()
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
            reason: format!(
                "formula transform 'relevel({expression})' must provide a factor"
            ),
        });
    };
    let Some(reference) = reference else {
        return Err(DeseqError::InvalidOptions {
            reason: format!("formula transform 'relevel({expression})' must provide ref"),
        });
    };
    let factor = formula_factor_argument_covariate(factor_text.trim(), factors)?;
    let reference = resolve_formula_factor_reference_alias(&factor, &reference)?;
    validate_factor_design_inputs_with_levels(
        &factor.name,
        &factor.sample_levels,
        &reference,
        factor.levels.as_deref(),
    )?;
    Ok(FormulaDerivedFactorCovariate {
        name: format!("relevel({}, ref = \"{reference}\")", factor.name),
        sample_levels: factor.sample_levels,
        reference,
        levels: factor.levels,
    })
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
        return Ok(FormulaDerivedFactorCovariate {
            name: factor.factor.to_string(),
            sample_levels: factor.sample_levels.to_vec(),
            reference: factor.reference.to_string(),
            levels: factor.levels.map(|levels| levels.to_vec()),
        });
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

fn formula_observed_levels_for_factor(factor: &FormulaDerivedFactorCovariate) -> Vec<String> {
    let mut levels = Vec::new();
    if let Some(declared) = factor.levels.as_ref() {
        for level in declared {
            if factor.sample_levels.iter().any(|sample| sample == level) {
                levels.push(level.clone());
            }
        }
        return levels;
    }
    for level in &factor.sample_levels {
        if !levels.iter().any(|candidate| candidate == level) {
            levels.push(level.clone());
        }
    }
    levels
}
