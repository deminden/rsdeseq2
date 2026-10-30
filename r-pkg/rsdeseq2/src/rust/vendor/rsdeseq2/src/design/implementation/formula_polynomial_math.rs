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
    let basis = super::numeric_transform::polynomial_basis(values, degree as usize, raw, None)
        .map_err(|error| DeseqError::InvalidOptions {
            reason: format!("formula transform 'poly({expression})': {error}"),
        })?;
    let mut replacements = Vec::with_capacity(degree as usize);
    let mut derived = Vec::with_capacity(degree as usize);
    for column in 0..degree as usize {
        let name = if raw {
            format!("{numeric_name}_poly_{}", column + 1)
        } else {
            format!("poly({numeric_name}, {degree}){}", column + 1)
        };
        replacements.push(formula_variable_term(&name));
        derived.push((name, (0..values.len()).map(|row| basis.values.as_slice()[row * degree as usize + column]).collect()));
    }
    Ok((format!("({})", replacements.join(" + ")), derived))
}
