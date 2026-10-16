// Raw and orthogonal polynomial column construction.

fn formula_numeric_raw_poly_columns(
    expression: &str,
    numeric_name: &str,
    values: &[f64],
    degree: i32,
) -> Result<FormulaNumericTransformExpansion, DeseqError> {
    let mut replacement_terms = Vec::with_capacity(degree as usize);
    let mut derived = Vec::with_capacity(degree as usize);
    for exponent in 1..=degree {
        let name = format!("{numeric_name}_poly_{exponent}");
        let mut column_values = Vec::with_capacity(values.len());
        for (idx, value) in values.iter().copied().enumerate() {
            let transformed = value.powi(exponent);
            if !transformed.is_finite() {
                return Err(DeseqError::InvalidOptions {
                    reason: format!(
                        "formula transform 'poly({expression})' produced non-finite value at sample {idx}"
                    ),
                });
            }
            column_values.push(transformed);
        }
        replacement_terms.push(formula_variable_term(&name));
        derived.push((name, column_values));
    }
    Ok((format!("({})", replacement_terms.join(" + ")), derived))
}

fn formula_numeric_orthogonal_poly_columns(
    expression: &str,
    numeric_name: &str,
    values: &[f64],
    degree: i32,
) -> Result<FormulaNumericTransformExpansion, DeseqError> {
    validate_formula_poly_values(expression, values, degree)?;
    let mut basis = Vec::with_capacity(degree as usize + 1);
    let constant_norm = (values.len() as f64).sqrt();
    basis.push(vec![1.0 / constant_norm; values.len()]);
    let mut replacement_terms = Vec::with_capacity(degree as usize);
    let mut derived = Vec::with_capacity(degree as usize);
    for exponent in 1..=degree {
        let mut column_values = Vec::with_capacity(values.len());
        for (idx, value) in values.iter().copied().enumerate() {
            let transformed = value.powi(exponent);
            if !transformed.is_finite() {
                return Err(DeseqError::InvalidOptions {
                    reason: format!(
                        "formula transform 'poly({expression})' produced non-finite value at sample {idx}"
                    ),
                });
            }
            column_values.push(transformed);
        }
        for previous in &basis {
            let projection = formula_dot_product(&column_values, previous);
            for (value, previous_value) in column_values.iter_mut().zip(previous) {
                *value -= projection * previous_value;
            }
        }
        let norm = formula_dot_product(&column_values, &column_values).sqrt();
        if !norm.is_finite() || norm <= 0.0 {
            return Err(DeseqError::InvalidOptions {
                reason: format!(
                    "formula transform 'poly({expression})' produced a dependent polynomial column"
                ),
            });
        }
        for value in &mut column_values {
            *value /= norm;
        }
        let name = format!("poly({numeric_name}, {degree}){exponent}");
        replacement_terms.push(formula_variable_term(&name));
        derived.push((name, column_values.clone()));
        basis.push(column_values);
    }
    Ok((format!("({})", replacement_terms.join(" + ")), derived))
}

fn validate_formula_poly_values(
    expression: &str,
    values: &[f64],
    degree: i32,
) -> Result<(), DeseqError> {
    let mut unique = Vec::new();
    for (idx, value) in values.iter().copied().enumerate() {
        if !value.is_finite() {
            return Err(DeseqError::InvalidOptions {
                reason: format!(
                    "formula transform 'poly({expression})' received non-finite value at sample {idx}"
                ),
            });
        }
        if !unique.contains(&value) {
            unique.push(value);
        }
    }
    if unique.len() <= degree as usize {
        return Err(DeseqError::InvalidOptions {
            reason: format!(
                "formula transform 'poly({expression})' degree must be less than the number of unique values"
            ),
        });
    }
    Ok(())
}

fn formula_dot_product(left: &[f64], right: &[f64]) -> f64 {
    left.iter()
        .zip(right)
        .map(|(left, right)| left * right)
        .sum()
}
