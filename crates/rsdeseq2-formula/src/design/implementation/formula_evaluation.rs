use super::*;
// The string parser translates its validated arguments and retains its public
// labels. Numerical execution uses the same evaluator as the primitive bridge.
pub(super) fn formula_numeric_value(values: &[f64]) -> super::super::expression::Value {
    super::super::expression::Value::new(super::super::expression::Data::Number {
        values: values.to_vec(),
        integer: false,
    })
}

pub(super) fn formula_call_values(
    function: &str,
    arguments: Vec<(Option<&str>, super::super::expression::Value)>,
    context: &str,
) -> Result<Vec<f64>, DeseqError> {
    let result = formula_call_value(function, arguments, context)?;
    formula_value_numbers(result, context)
}

pub(super) fn formula_call_value(
    function: &str,
    arguments: Vec<(Option<&str>, super::super::expression::Value)>,
    context: &str,
) -> Result<super::super::expression::Value, DeseqError> {
    use super::super::expression::{Argument, Expression};
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

pub(super) fn formula_expression_values(
    expression: &super::super::expression::Expression,
    context: &str,
) -> Result<Vec<f64>, DeseqError> {
    formula_value_numbers(formula_evaluate_expression(expression, context)?, context)
}

pub(super) fn formula_evaluate_expression(
    expression: &super::super::expression::Expression,
    context: &str,
) -> Result<super::super::expression::Value, DeseqError> {
    use super::super::expression::EvaluationOptions;
    super::super::expression::evaluate(
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

pub(super) fn formula_value_numbers(
    result: super::super::expression::Value,
    context: &str,
) -> Result<Vec<f64>, DeseqError> {
    use super::super::expression::Data;
    let Data::Number { values, .. } = result.data else {
        return Err(DeseqError::InvalidOptions {
            reason: format!("formula transform '{context}' must produce numeric values"),
        });
    };
    formula_finite_values(values, context)
}

pub(super) fn formula_finite_values(
    values: Vec<f64>,
    context: &str,
) -> Result<Vec<f64>, DeseqError> {
    if let Some(sample) = values.iter().position(|value| !value.is_finite()) {
        return Err(DeseqError::InvalidOptions {
            reason: format!(
                "formula transform '{context}' produced non-finite value at sample {sample}"
            ),
        });
    }
    Ok(values)
}
