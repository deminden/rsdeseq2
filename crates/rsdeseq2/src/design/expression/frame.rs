//! Shared model-frame row policy, categorical preparation and offset arithmetic.
use super::{
    Argument, Data, EvaluationOptions, Expression, Value, evaluate_many, invalid,
    validation::validate_value,
};
use crate::errors::DeseqError;

/// Missing-value action supplied as metadata by a frontend.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MissingPolicy {
    /// Retain all rows, including missing values.
    Pass,
    /// Reject missing covariates.
    Fail,
    /// Remove rows having any missing entry in any model variable.
    Omit,
}
/// Evaluated columns and row selection, before frontend object construction.
#[derive(Clone, Debug, PartialEq)]
pub struct Frame {
    /// Columns in the original variable order.
    pub values: Vec<Value>,
    /// Number of rows before missing-value selection.
    pub original_rows: usize,
    /// Retained zero-based row indices.
    pub kept_rows: Vec<usize>,
    /// Omitted zero-based row indices.
    pub omitted_rows: Vec<usize>,
}
/// Evaluate a frame in one shared numerical transaction, before frontend object
/// construction. Prediction metadata remains attached to the selected values.
pub fn evaluate(
    expressions: &[Expression],
    variables: &[(String, Value)],
    options: &EvaluationOptions,
    data_rows: usize,
    policy: MissingPolicy,
) -> Result<Frame, DeseqError> {
    assemble(
        evaluate_many(expressions, variables, options)?,
        data_rows,
        policy,
    )
}
fn rows(value: &Value) -> usize {
    value.dimensions.map_or(value.len(), |d| d[0])
}
fn missing(value: &Value, row: usize) -> bool {
    let width = value.dimensions.map_or(1, |d| d[1]);
    (0..width).any(|column| {
        let i = row * width + column;
        match &value.data {
            Data::Number { values, .. } => values[i].is_nan(),
            Data::Logical(v) => v[i].is_none(),
            Data::Text(v) => v[i].is_none(),
            _ => false,
        }
    })
}
fn select(value: &Value, kept: &[usize]) -> Value {
    let mut out = value.clone();
    let width = value.dimensions.map_or(1, |d| d[1]);
    let indices = kept
        .iter()
        .flat_map(|r| (0..width).map(move |c| r * width + c))
        .collect::<Vec<_>>();
    out.data = match &value.data {
        Data::Number { values, integer } => Data::Number {
            values: indices.iter().map(|i| values[*i]).collect(),
            integer: *integer,
        },
        Data::Logical(v) => Data::Logical(indices.iter().map(|i| v[*i]).collect()),
        Data::Text(v) => Data::Text(indices.iter().map(|i| v[*i].clone()).collect()),
        _ => unreachable!("frame columns were checked before selection"),
    };
    if let Some([_, columns]) = value.dimensions {
        out.names = None;
        out.dimensions = Some([kept.len(), columns]);
        out.dimnames[0] = value.dimnames[0]
            .as_ref()
            .map(|n| kept.iter().map(|i| n[*i].clone()).collect());
        // Matrix row subsetting constructs fresh label vectors in R, dropping
        // names on those label vectors even for an all-row subset.
        out.dimname_names = [None, None];
    } else {
        out.names = value
            .names
            .as_ref()
            .map(|n| kept.iter().map(|i| n[*i].clone()).collect());
    }
    out
}
/// Validate column row counts and apply the same missing-value policy in every
/// frontend. With no variables, the data's row count determines an intercept.
pub fn assemble(
    values: Vec<Value>,
    data_rows: usize,
    policy: MissingPolicy,
) -> Result<Frame, DeseqError> {
    for value in &values {
        validate_value(value, 0)?;
        if matches!(value.data, Data::Null | Data::List(_)) {
            return Err(invalid(
                "model variables must be atomic vectors or matrices",
            ));
        }
    }
    let original_rows = values.first().map_or(data_rows, rows);
    if values.iter().any(|v| rows(v) != original_rows) {
        return Err(invalid("model variable lengths differ"));
    }
    let mut kept_rows = Vec::new();
    let mut omitted_rows = Vec::new();
    for row in 0..original_rows {
        let absent = values.iter().any(|v| missing(v, row));
        if absent && policy == MissingPolicy::Fail {
            return Err(invalid("missing values in model variables"));
        }
        if absent && policy == MissingPolicy::Omit {
            omitted_rows.push(row);
        } else {
            kept_rows.push(row);
        }
    }
    let update = values
        .iter()
        .any(|v| v.names.is_some() || v.dimensions.is_some() || !v.attributes.is_empty());
    let values = if !omitted_rows.is_empty() || policy == MissingPolicy::Omit && update {
        values.iter().map(|v| select(v, &kept_rows)).collect()
    } else {
        values
    };
    Ok(Frame {
        values,
        original_rows,
        kept_rows,
        omitted_rows,
    })
}
/// Convert characters/logicals and the response for numeric model construction,
/// preserving declared factor levels and performing no foreign calculation.
pub fn prepare_model_values(
    values: &[Value],
    response: Option<usize>,
    options: &EvaluationOptions,
) -> Result<Vec<Value>, DeseqError> {
    if response.is_some_and(|i| i >= values.len()) {
        return Err(invalid("model response index is outside columns"));
    }
    let variables = values
        .iter()
        .enumerate()
        .map(|(i, v)| (format!("v{i}"), v.clone()))
        .collect::<Vec<_>>();
    let wrap = |function: &str, expression: Expression| Expression::Call {
        function: function.into(),
        arguments: vec![Argument {
            name: None,
            expression,
        }],
    };
    let expressions = values
        .iter()
        .enumerate()
        .map(|(i, v)| {
            let mut x = Expression::Symbol(format!("v{i}"));
            if matches!(v.data, Data::Text(_)) {
                x = wrap("factor", x);
            }
            if response == Some(i) {
                x = wrap("as.double", x);
            } else if matches!(v.data, Data::Logical(_)) {
                x = Expression::Call {
                    function: "factor".into(),
                    arguments: vec![
                        Argument {
                            name: None,
                            expression: x,
                        },
                        Argument {
                            name: Some("levels".into()),
                            expression: Expression::Literal(Value::new(Data::Logical(vec![
                                Some(false),
                                Some(true),
                            ]))),
                        },
                    ],
                };
            }
            x
        })
        .collect::<Vec<_>>();
    evaluate_many(&expressions, &variables, options)
}
/// Sum explicit model offsets in their declared order. The scalar zero seed is
/// intentional, including its effect on signed zero. Indices are zero-based.
pub fn offsets(
    values: &[Value],
    indices: &[usize],
    base: Option<&Value>,
) -> Result<Option<Value>, DeseqError> {
    if indices.is_empty() {
        return Ok(base.cloned());
    }
    let mut result = base.cloned().unwrap_or_else(|| Value::scalar(0.0));
    validate_value(&result, 0)?;
    for i in indices {
        let value = values
            .get(*i)
            .ok_or_else(|| invalid("model offset index is outside columns"))?;
        validate_value(value, 0)?;
        result = super::arithmetic::binary("+", &result, value)?;
    }
    Ok(Some(result))
}
