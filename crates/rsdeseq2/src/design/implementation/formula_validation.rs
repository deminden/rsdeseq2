fn validate_formula_variable(variable: &str) -> Result<(), DeseqError> {
    if variable.is_empty()
        || variable
            .chars()
            .any(|character| character.is_whitespace() || "+-*/:^()".contains(character))
    {
        return Err(DeseqError::InvalidOptions {
            reason: format!("formula variable '{variable}' is not a supported bare variable name"),
        });
    }
    Ok(())
}

fn formula_variable_name(variable: &str) -> Result<&str, DeseqError> {
    let variable = variable.trim();
    if let Some(stripped) = variable
        .strip_prefix('`')
        .and_then(|value| value.strip_suffix('`'))
    {
        if stripped.is_empty() || stripped.contains('`') {
            return Err(DeseqError::InvalidOptions {
                reason: format!(
                    "formula variable '{variable}' is not a supported quoted variable name"
                ),
            });
        }
        return Ok(stripped);
    }
    validate_formula_variable(variable)?;
    Ok(variable)
}

fn validate_formula_model_frame_column_name(name: &str) -> Result<(), DeseqError> {
    if name.is_empty() || name.contains('`') {
        return Err(DeseqError::InvalidOptions {
            reason: format!("formula model-frame column name '{name}' is not supported"),
        });
    }
    Ok(())
}
