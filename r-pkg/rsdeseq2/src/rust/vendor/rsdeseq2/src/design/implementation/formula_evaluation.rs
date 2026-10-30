// The string parser translates its validated arguments and retains its public
// labels. Numerical execution uses the same evaluator as the primitive bridge.
fn formula_numeric_value(values: &[f64]) -> super::expression::Value {
    super::expression::Value::new(super::expression::Data::Number {
        values: values.to_vec(),
        integer: false,
    })
}

fn formula_call_values(
    function: &str,
    arguments: Vec<(Option<&str>, super::expression::Value)>,
    context: &str,
) -> Result<Vec<f64>, DeseqError> {
    let result = formula_call_value(function, arguments, context)?;
    formula_value_numbers(result, context)
}

fn formula_call_value(
    function: &str,
    arguments: Vec<(Option<&str>, super::expression::Value)>,
    context: &str,
) -> Result<super::expression::Value, DeseqError> {
    use super::expression::{Argument, Expression};
    let expression = Expression::Call {
        function: function.to_string(),
        arguments: arguments
            .into_iter()
            .map(|(name, value)| Argument {
                name: name.map(str::to_string),
                expression: Expression::Literal(value),
            })
            .collect(),
    };
    formula_evaluate_expression(&expression, context)
}

fn formula_expression_values(
    expression: &super::expression::Expression,
    context: &str,
) -> Result<Vec<f64>, DeseqError> {
    formula_value_numbers(formula_evaluate_expression(expression, context)?, context)
}

fn formula_evaluate_expression(
    expression: &super::expression::Expression,
    context: &str,
) -> Result<super::expression::Value, DeseqError> {
    use super::expression::EvaluationOptions;
    super::expression::evaluate(
        expression,
        &[],
        &EvaluationOptions {
            collation: "C".into(),
        },
    )
    .map_err(|error| DeseqError::InvalidOptions {
        reason: format!("formula transform '{context}': {error}"),
    })
}

fn formula_value_numbers(
    result: super::expression::Value,
    context: &str,
) -> Result<Vec<f64>, DeseqError> {
    use super::expression::Data;
    let Data::Number { values, .. } = result.data else {
        return Err(DeseqError::InvalidOptions {
            reason: format!("formula transform '{context}' must produce numeric values"),
        });
    };
    formula_finite_values(values, context)
}

fn formula_finite_values(values: Vec<f64>, context: &str) -> Result<Vec<f64>, DeseqError> {
    if let Some(sample) = values.iter().position(|value| !value.is_finite()) {
        return Err(DeseqError::InvalidOptions {
            reason: format!(
                "formula transform '{context}' produced non-finite value at sample {sample}"
            ),
        });
    }
    Ok(values)
}
