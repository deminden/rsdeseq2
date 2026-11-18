use super::*;
use crate::errors::invalid_formula;
pub(super) fn validate_formula_variable(variable: &str) -> Result<(), Error> {
    if variable.is_empty()
        || variable
            .chars()
            .any(|character| character.is_whitespace() || "+-*/:^()".contains(character))
    {
        return Err(invalid_formula(format!(
            "formula variable '{variable}' is not a supported bare variable name"
        )));
    }
    Ok(())
}

pub(super) fn formula_variable_name(variable: &str) -> Result<&str, Error> {
    let variable = variable.trim();
    if let Some(stripped) = variable
        .strip_prefix('`')
        .and_then(|value| value.strip_suffix('`'))
    {
        if stripped.is_empty() || stripped.contains('`') {
            return Err(invalid_formula(format!(
                "formula variable '{variable}' is not a supported quoted variable name"
            )));
        }
        return Ok(stripped);
    }
    validate_formula_variable(variable)?;
    Ok(variable)
}

pub(super) fn validate_formula_model_frame_column_name(name: &str) -> Result<(), Error> {
    if name.is_empty() || name.contains('`') {
        return Err(invalid_formula(format!(
            "formula model-frame column name '{name}' is not supported"
        )));
    }
    Ok(())
}
