// Formula labels around the standalone polynomial operation.

fn formula_numeric_raw_poly_columns(
    expression: &str,
    numeric_name: &str,
    values: &[f64],
    degree: i32,
) -> Result<FormulaNumericTransformExpansion, DeseqError> {
    formula_polynomial_columns(expression, numeric_name, values, degree, true)
}

fn formula_numeric_orthogonal_poly_columns(
    expression: &str,
    numeric_name: &str,
    values: &[f64],
    degree: i32,
) -> Result<FormulaNumericTransformExpansion, DeseqError> {
    formula_polynomial_columns(expression, numeric_name, values, degree, false)
}

fn formula_polynomial_columns(
    expression: &str,
    numeric_name: &str,
    values: &[f64],
    degree: i32,
    raw: bool,
) -> Result<FormulaNumericTransformExpansion, DeseqError> {
    let basis = formula_call_values(
        "poly",
        vec![
            (Some("x"), formula_numeric_value(values)),
            (
                Some("degree"),
                super::expression::Value::scalar(degree as f64),
            ),
            (
                Some("raw"),
                super::expression::Value::new(super::expression::Data::Logical(vec![Some(raw)])),
            ),
        ],
        &format!("poly({expression})"),
    )?;
    let mut replacements = Vec::with_capacity(degree as usize);
    let mut derived = Vec::with_capacity(degree as usize);
    for column in 0..degree as usize {
        let name = if raw {
            format!("{numeric_name}_poly_{}", column + 1)
        } else {
            format!("poly({numeric_name}, {degree}){}", column + 1)
        };
        replacements.push(formula_variable_term(&name));
        derived.push((
            name,
            (0..values.len())
                .map(|row| basis[row * degree as usize + column])
                .collect(),
        ));
    }
    Ok((format!("({})", replacements.join(" + ")), derived))
}
