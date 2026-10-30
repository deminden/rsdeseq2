// Factor operations and recoding execute in the shared expression engine.
// This adapter translates labels and the compatibility API's reference metadata.
fn formula_text_value(values: &[String]) -> super::expression::Value {
    super::expression::Value::new(super::expression::Data::Text(
        values.iter().cloned().map(Some).collect(),
    ))
}

fn formula_factor_result(
    name: String,
    primitive: super::expression::Value,
    preferred_reference: Option<&str>,
) -> Result<FormulaDerivedFactorCovariate, DeseqError> {
    let levels = primitive
        .attribute("levels")
        .ok_or_else(|| DeseqError::InvalidOptions {
            reason: format!("formula factor '{name}' has no levels"),
        })?
        .text()?
        .into_iter()
        .map(|v| {
            v.ok_or_else(|| DeseqError::InvalidOptions {
                reason: format!("formula factor '{name}' contains a missing level"),
            })
        })
        .collect::<Result<Vec<_>, _>>()?;
    let sample_levels = primitive
        .numeric()?
        .into_iter()
        .enumerate()
        .map(|(i, code)| {
            if !code.is_finite() || code.fract() != 0.0 || code < 1.0 || code > levels.len() as f64
            {
                return Err(DeseqError::InvalidOptions {
                    reason: format!(
                        "formula factor '{name}' has a missing or invalid level at sample {i}"
                    ),
                });
            }
            Ok(levels[code as usize - 1].clone())
        })
        .collect::<Result<Vec<_>, _>>()?;
    let reference = preferred_reference
        .filter(|r| levels.iter().any(|v| v == *r))
        .or_else(|| levels.first().map(String::as_str))
        .ok_or_else(|| DeseqError::InvalidOptions {
            reason: format!("formula factor '{name}' requires at least one level"),
        })?
        .to_string();
    validate_factor_design_inputs_with_levels(&name, &sample_levels, &reference, Some(&levels))?;
    Ok(FormulaDerivedFactorCovariate {
        name,
        sample_levels,
        reference,
        levels: Some(levels),
        primitive,
    })
}
