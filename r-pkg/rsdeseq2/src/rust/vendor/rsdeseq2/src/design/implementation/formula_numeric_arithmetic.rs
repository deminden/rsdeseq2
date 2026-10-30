// Identity, scalar, binary and power arithmetic transforms.

fn split_formula_transform_arguments(expression: &str) -> Result<Vec<String>, DeseqError> {
    split_formula_top_level(expression, ',')
}

fn formula_numeric_identity_or_power_term(
    expression: &str,
    numeric_covariates: &[ExpandedNumericSpec<'_>],
) -> Result<(String, Vec<f64>), DeseqError> {
    if let Some((key, value)) = split_formula_named_argument(expression.trim()) {
        if key == "x" {
            return formula_numeric_identity_or_power_term(value, numeric_covariates);
        }
        return Err(DeseqError::InvalidOptions {
            reason: format!("formula transform 'I({expression})' has unsupported argument '{key}'"),
        });
    }
    let Some((numeric_name, exponent_text)) = expression.split_once('^') else {
        if let Some((left, op, right)) = split_formula_numeric_scalar_expression(expression) {
            return formula_numeric_scalar_expression_term(
                expression,
                left,
                op,
                right,
                numeric_covariates,
            );
        }
        return formula_numeric_identity_term(expression, numeric_covariates);
    };
    formula_numeric_power_term(expression, numeric_name, exponent_text, numeric_covariates)
}

fn split_formula_numeric_scalar_expression(expression: &str) -> Option<(&str, char, &str)> {
    let mut in_backticks = false;
    for (idx, op) in expression.char_indices() {
        if op == '`' {
            in_backticks = !in_backticks;
            continue;
        }
        if in_backticks {
            continue;
        }
        if idx == 0 || !matches!(op, '+' | '-' | '*' | '/') {
            continue;
        }
        let left = expression[..idx].trim();
        let right = expression[idx + op.len_utf8()..].trim();
        if left.is_empty() || right.is_empty() {
            continue;
        }
        return Some((left, op, right));
    }
    None
}

fn formula_numeric_scalar_expression_term(
    expression: &str,
    left: &str,
    op: char,
    right: &str,
    numeric_covariates: &[ExpandedNumericSpec<'_>],
) -> Result<(String, Vec<f64>), DeseqError> {
    let left_name = formula_variable_name(left).ok();
    let right_name = formula_variable_name(right).ok();
    let left_covariate = maybe_formula_numeric_covariate(left_name, numeric_covariates)?;
    let right_covariate = maybe_formula_numeric_covariate(right_name, numeric_covariates)?;
    let (numeric_name, scalar_text, scalar_on_left, covariate) = match (
        left_covariate,
        right_covariate,
    ) {
        (Some(covariate), None) => (covariate.name, right, false, covariate),
        (None, Some(covariate)) => (covariate.name, left, true, covariate),
        (Some(left_covariate), Some(right_covariate)) => {
            return formula_numeric_binary_expression_term(
                expression,
                left_covariate.name,
                op,
                right_covariate.name,
                left_covariate,
                right_covariate,
            );
        }
        (None, None) => {
            return Err(DeseqError::InvalidOptions {
                reason: format!(
                    "formula transform 'I({expression})' must contain one supplied numeric covariate"
                ),
            });
        }
    };
    let scalar = parse_formula_scalar(scalar_text, expression)?;
    if op == '/' && !scalar_on_left && scalar == 0.0 {
        return Err(DeseqError::InvalidOptions {
            reason: format!(
                "formula transform 'I({expression})' scalar '{scalar_text}' must be non-zero for division"
            ),
        });
    }
    let input = formula_numeric_value(covariate.values);
    let scalar_value = super::expression::Value::scalar(scalar);
    let arguments = if scalar_on_left {
        vec![(None, scalar_value), (None, input)]
    } else {
        vec![(None, input), (None, scalar_value)]
    };
    let values = formula_call_values(&op.to_string(), arguments, &format!("I({expression})"))?;
    Ok((
        format!(
            "{}_{}_{}",
            numeric_name,
            formula_scalar_operator_label(op, scalar_on_left),
            formula_scalar_value_label(scalar_text)
        ),
        values,
    ))
}

fn formula_numeric_binary_expression_term(
    expression: &str,
    left_name: &str,
    op: char,
    right_name: &str,
    left: &ExpandedNumericSpec<'_>,
    right: &ExpandedNumericSpec<'_>,
) -> Result<(String, Vec<f64>), DeseqError> {
    if left.values.len() != right.values.len() {
        return Err(DeseqError::InvalidOptions {
            reason: format!(
                "formula transform 'I({expression})' received numeric covariates with different sample counts"
            ),
        });
    }
    let values = formula_call_values(
        &op.to_string(),
        vec![
            (None, formula_numeric_value(left.values)),
            (None, formula_numeric_value(right.values)),
        ],
        &format!("I({expression})"),
    )?;
    Ok((
        format!(
            "{}_{}_{}",
            left_name,
            formula_binary_operator_label(op),
            right_name
        ),
        values,
    ))
}

fn parse_formula_scalar(value: &str, expression: &str) -> Result<f64, DeseqError> {
    let value = strip_formula_outer_parentheses(value.trim())?;
    let scalar = value
        .parse::<f64>()
        .map_err(|_| DeseqError::InvalidOptions {
            reason: format!("formula transform 'I({expression})' scalar '{value}' must be finite"),
        })?;
    if !scalar.is_finite() {
        return Err(DeseqError::InvalidOptions {
            reason: format!("formula transform 'I({expression})' scalar '{value}' must be finite"),
        });
    }
    Ok(scalar)
}

fn formula_binary_operator_label(op: char) -> &'static str {
    match op {
        '+' => "plus",
        '-' => "minus",
        '*' => "times",
        '/' => "div",
        _ => unreachable!("validated formula binary operator"),
    }
}

fn formula_scalar_operator_label(op: char, scalar_on_left: bool) -> &'static str {
    match (scalar_on_left, op) {
        (_, '+') => "plus",
        (false, '-') => "minus",
        (true, '-') => "rminus",
        (_, '*') => "times",
        (false, '/') => "div",
        (true, '/') => "rdiv",
        _ => unreachable!("validated formula scalar operator"),
    }
}

fn formula_scalar_value_label(value: &str) -> String {
    let mut label = String::new();
    for ch in value.trim().chars() {
        match ch {
            '-' => label.push_str("neg"),
            '+' => label.push_str("pos"),
            '.' => label.push('p'),
            ch if ch.is_ascii_alphanumeric() => label.push(ch),
            _ => label.push('_'),
        }
    }
    label
}

fn formula_numeric_identity_term(
    expression: &str,
    numeric_covariates: &[ExpandedNumericSpec<'_>],
) -> Result<(String, Vec<f64>), DeseqError> {
    let (numeric_term, sign, label) = match expression.trim().as_bytes().first().copied() {
        Some(b'+') => (expression.trim()[1..].trim(), 1.0, "identity"),
        Some(b'-') => (expression.trim()[1..].trim(), -1.0, "neg"),
        _ => (expression.trim(), 1.0, "identity"),
    };
    let numeric_name = formula_variable_name(numeric_term)?;
    let covariate = find_formula_numeric_covariate(numeric_name, numeric_covariates)?;
    let numeric_name = covariate.name;
    let values = formula_call_values(
        if sign < 0.0 { "-" } else { "+" },
        vec![(None, formula_numeric_value(covariate.values))],
        &format!("I({expression})"),
    )?;
    Ok((format!("{numeric_name}_{label}"), values))
}

fn formula_numeric_power_term(
    expression: &str,
    numeric_name: &str,
    exponent_text: &str,
    numeric_covariates: &[ExpandedNumericSpec<'_>],
) -> Result<(String, Vec<f64>), DeseqError> {
    let numeric_name = formula_variable_name(numeric_name.trim())?;
    let exponent = exponent_text
        .trim()
        .parse::<i32>()
        .map_err(|_| DeseqError::InvalidOptions {
            reason: format!("formula transform exponent '{exponent_text}' must be an integer"),
        })?;
    if !(2..=16).contains(&exponent) {
        return Err(DeseqError::InvalidOptions {
            reason: "formula transform powers must be integers from 2 through 16".to_string(),
        });
    }
    let covariate = find_formula_numeric_covariate(numeric_name, numeric_covariates)?;
    let numeric_name = covariate.name;
    let values = formula_call_values(
        "^",
        vec![
            (None, formula_numeric_value(covariate.values)),
            (None, super::expression::Value::scalar(exponent as f64)),
        ],
        &format!("I({expression})"),
    )?;
    Ok((format!("{numeric_name}_pow_{exponent}"), values))
}

fn formula_variable_term(name: &str) -> String {
    if validate_formula_variable(name).is_ok() {
        name.to_string()
    } else {
        format!("`{name}`")
    }
}
