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
