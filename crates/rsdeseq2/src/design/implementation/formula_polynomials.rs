// Polynomial transform argument handling.

fn formula_numeric_poly_term(
    expression: &str,
    numeric_covariates: &[ExpandedNumericSpec<'_>],
) -> Result<FormulaNumericTransformExpansion, DeseqError> {
    let arguments = split_formula_transform_arguments(expression)?;
    if arguments.len() < 2 {
        return Err(DeseqError::InvalidOptions {
            reason: format!(
                "formula transform 'poly({expression})' must provide numeric and degree"
            ),
        });
    }
    let mut numeric_text = None;
    let mut degree_text = None;
    let mut raw = false;
    let mut raw_seen = false;
    let mut simple_seen = false;
    let mut positional_idx = 0_usize;
    for argument in &arguments {
        let trimmed = argument.trim();
        let normalized = argument
            .split_whitespace()
            .collect::<String>()
            .to_ascii_lowercase();
        if normalized == "raw=true" || normalized == "raw=t" {
            reject_duplicate_formula_poly_argument(raw_seen, "raw", expression)?;
            raw_seen = true;
            raw = true;
            continue;
        }
        if normalized == "raw=false" || normalized == "raw=f" {
            reject_duplicate_formula_poly_argument(raw_seen, "raw", expression)?;
            raw_seen = true;
            raw = false;
            continue;
        }
        if normalized == "simple=true"
            || normalized == "simple=t"
            || normalized == "simple=false"
            || normalized == "simple=f"
        {
            reject_duplicate_formula_poly_argument(simple_seen, "simple", expression)?;
            simple_seen = true;
            continue;
        }
        if let Some(named_degree) = normalized.strip_prefix("degree=") {
            reject_duplicate_formula_poly_argument(degree_text.is_some(), "degree", expression)?;
            degree_text = Some(named_degree.to_string());
            continue;
        }
        if let Some((key, value)) = split_formula_named_argument(trimmed) {
            match key.as_str() {
                "x" => {
                    reject_duplicate_formula_poly_argument(numeric_text.is_some(), "x", expression)?;
                    numeric_text = Some(value.to_string());
                    continue;
                }
                "degree" => {
                    reject_duplicate_formula_poly_argument(
                        degree_text.is_some(),
                        "degree",
                        expression,
                    )?;
                    degree_text = Some(value.to_string());
                    continue;
                }
                "raw" => {
                    reject_duplicate_formula_poly_argument(raw_seen, "raw", expression)?;
                    raw_seen = true;
                    raw = parse_formula_bool_option(value, "raw", expression)?;
                    continue;
                }
                "simple" => {
                    reject_duplicate_formula_poly_argument(simple_seen, "simple", expression)?;
                    simple_seen = true;
                    parse_formula_bool_option(value, "simple", expression)?;
                    continue;
                }
                _ => {
                    return Err(DeseqError::InvalidOptions {
                        reason: format!(
                            "formula transform 'poly({expression})' has unsupported argument '{argument}'"
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
            1 if degree_text.is_none() => {
                degree_text = Some(trimmed.to_string());
                positional_idx += 1;
                continue;
            }
            _ => {}
        }
        if positional_idx < 2 {
            positional_idx += 1;
            continue;
        }
        return Err(DeseqError::InvalidOptions {
            reason: format!(
                "formula transform 'poly({expression})' has unsupported argument '{argument}'"
            ),
        });
    }
    let Some(degree_text) = degree_text else {
        return Err(DeseqError::InvalidOptions {
            reason: format!("formula transform 'poly({expression})' must provide a degree"),
        });
    };
    let Some(numeric_text) = numeric_text else {
        return Err(DeseqError::InvalidOptions {
            reason: format!(
                "formula transform 'poly({expression})' must provide a numeric covariate"
            ),
        });
    };
    let numeric_name = formula_variable_name(numeric_text.trim())?;
    let degree = degree_text
        .parse::<i32>()
        .map_err(|_| DeseqError::InvalidOptions {
            reason: format!("formula polynomial degree '{degree_text}' must be an integer"),
        })?;
    if !(1..=16).contains(&degree) {
        return Err(DeseqError::InvalidOptions {
            reason: "formula polynomial degrees must be integers from 1 through 16".to_string(),
        });
    }
    let covariate = find_formula_numeric_covariate(numeric_name, numeric_covariates)?;
    let numeric_name = covariate.name;
    if raw {
        return formula_numeric_raw_poly_columns(expression, numeric_name, covariate.values, degree);
    }
    formula_numeric_orthogonal_poly_columns(expression, numeric_name, covariate.values, degree)
}
