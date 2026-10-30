// Numeric scaling transforms and their options.

fn formula_numeric_scale_term(
    expression: &str,
    numeric_covariates: &[ExpandedNumericSpec<'_>],
) -> Result<(String, Vec<f64>), DeseqError> {
    let arguments = split_formula_transform_arguments(expression)?;
    let mut numeric_text = None;
    let mut center = FormulaScaleOption::Auto;
    let mut scale = FormulaScaleOption::Auto;
    let mut center_seen = false;
    let mut scale_seen = false;
    let mut positional_idx = 0_usize;
    for argument in &arguments {
        let trimmed = argument.trim();
        let normalized = argument
            .split_whitespace()
            .collect::<String>()
            .to_ascii_lowercase();
        if let Some(value) = normalized.strip_prefix("center=") {
            reject_duplicate_formula_transform_argument(
                center_seen,
                "scale",
                "center",
                expression,
            )?;
            center_seen = true;
            center = parse_formula_scale_option(value, "center", expression)?;
            continue;
        }
        if let Some(value) = normalized.strip_prefix("scale=") {
            reject_duplicate_formula_transform_argument(
                scale_seen,
                "scale",
                "scale",
                expression,
            )?;
            scale_seen = true;
            scale = parse_formula_scale_option(value, "scale", expression)?;
            continue;
        }
        if let Some((key, value)) = split_formula_named_argument(trimmed) {
            match key.as_str() {
                "x" => {
                    reject_duplicate_formula_transform_argument(
                        numeric_text.is_some(),
                        "scale",
                        "x",
                        expression,
                    )?;
                    numeric_text = Some(value.to_string());
                    continue;
                }
                "center" => {
                    reject_duplicate_formula_transform_argument(
                        center_seen,
                        "scale",
                        "center",
                        expression,
                    )?;
                    center_seen = true;
                    center = parse_formula_scale_option(value, "center", expression)?;
                    continue;
                }
                "scale" => {
                    reject_duplicate_formula_transform_argument(
                        scale_seen,
                        "scale",
                        "scale",
                        expression,
                    )?;
                    scale_seen = true;
                    scale = parse_formula_scale_option(value, "scale", expression)?;
                    continue;
                }
                _ => {
                    return Err(DeseqError::InvalidOptions {
                        reason: format!(
                            "formula transform 'scale({expression})' has unsupported argument '{argument}'"
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
            1 if !center_seen => {
                center_seen = true;
                center = parse_formula_scale_option(&normalized, "center", expression)?;
                positional_idx += 1;
                continue;
            }
            2 if !scale_seen => {
                scale_seen = true;
                scale = parse_formula_scale_option(&normalized, "scale", expression)?;
                positional_idx += 1;
                continue;
            }
            _ => {}
        }
        return Err(DeseqError::InvalidOptions {
            reason: format!(
                "formula transform 'scale({expression})' has unsupported argument '{argument}'"
            ),
        });
    }
    let Some(numeric_text) = numeric_text else {
        return Err(DeseqError::InvalidOptions {
            reason: format!(
                "formula transform 'scale({expression})' must provide a numeric covariate"
            ),
        });
    };
    let numeric_name = formula_variable_name(numeric_text.trim())?;
    let covariate = find_formula_numeric_covariate(numeric_name, numeric_covariates)?;
    let numeric_name = covariate.name;
    let centers;
    let divisors;
    let center_policy = match center {
        FormulaScaleOption::Auto => super::numeric_transform::ScalePolicy::Automatic,
        FormulaScaleOption::Disabled => super::numeric_transform::ScalePolicy::Disabled,
        FormulaScaleOption::Explicit(value) => { centers = [value]; super::numeric_transform::ScalePolicy::Supplied(&centers) },
    };
    let scale_policy = match scale {
        FormulaScaleOption::Auto => super::numeric_transform::ScalePolicy::Automatic,
        FormulaScaleOption::Disabled => super::numeric_transform::ScalePolicy::Disabled,
        FormulaScaleOption::Explicit(value) => { divisors = [value]; super::numeric_transform::ScalePolicy::Supplied(&divisors) },
    };
    let matrix = crate::matrix::MatrixView::new(covariate.values.len(), 1, covariate.values, crate::matrix::MatrixOrder::ColumnMajor)?;
    let values = super::numeric_transform::scale_covariates(matrix, center_policy, scale_policy)
        .map_err(|error| DeseqError::InvalidOptions { reason: format!("formula transform 'scale({expression})': {error}") })?
        .values.into_values();
    let label = match (center, scale) {
        (FormulaScaleOption::Auto, FormulaScaleOption::Auto) => "scale",
        (FormulaScaleOption::Disabled, FormulaScaleOption::Auto) => "scale_uncentered",
        (FormulaScaleOption::Auto, FormulaScaleOption::Disabled) => "center",
        (FormulaScaleOption::Disabled, FormulaScaleOption::Disabled) => "identity",
        (FormulaScaleOption::Explicit(_), FormulaScaleOption::Disabled) => "centered",
        (FormulaScaleOption::Disabled, FormulaScaleOption::Explicit(_)) => "scaled",
        _ => "centered_scaled",
    };
    Ok((format!("{numeric_name}_{label}"), values))
}

fn parse_formula_scale_option(
    value: &str,
    argument: &str,
    expression: &str,
) -> Result<FormulaScaleOption, DeseqError> {
    let value = strip_formula_outer_parentheses(value.trim())?;
    match value.to_ascii_lowercase().as_str() {
        "true" | "t" => Ok(FormulaScaleOption::Auto),
        "false" | "f" => Ok(FormulaScaleOption::Disabled),
        _ => {
            let parsed = value
                .parse::<f64>()
                .map_err(|_| DeseqError::InvalidOptions {
                    reason: format!(
                        "formula transform 'scale({expression})' argument '{argument}' must be TRUE, FALSE, or a finite number"
                    ),
                })?;
            if !parsed.is_finite() || (argument == "scale" && parsed == 0.0) {
                return Err(DeseqError::InvalidOptions {
                    reason: format!(
                        "formula transform 'scale({expression})' argument '{argument}' must be TRUE, FALSE, or a finite non-zero scale"
                    ),
                });
            }
            Ok(FormulaScaleOption::Explicit(parsed))
        }
    }
}
