use super::{Data, Value, invalid};
use crate::errors::DeseqError;

pub(super) fn parallel(name: &str, args: &[(Option<String>, Value)]) -> Result<Value, DeseqError> {
    let mut remove_missing = false;
    let mut flag_seen = false;
    let mut inputs = Vec::new();
    for (label, value) in args {
        if label.as_deref() == Some("na.rm") {
            if flag_seen {
                return Err(invalid("duplicate na.rm argument"));
            }
            flag_seen = true;
            remove_missing = value.bool_scalar()?;
        } else {
            if value.has_class("factor")
                || !matches!(value.data, Data::Number { .. } | Data::Logical(_))
            {
                return Err(invalid(
                    "parallel extrema require numeric or logical vectors",
                ));
            }
            inputs.push(value);
        }
    }
    if inputs.is_empty() {
        return Err(invalid("parallel extrema require at least one input"));
    }
    let integer = inputs.iter().all(|v| {
        matches!(
            v.data,
            Data::Logical(_) | Data::Number { integer: true, .. }
        )
    });
    let n = if inputs.iter().any(|v| v.is_empty()) {
        0
    } else {
        inputs.iter().map(|v| v.len()).max().unwrap()
    };
    let mut result = if n > 0 && inputs[0].len() == n {
        inputs[0].clone()
    } else {
        Value::numbers(vec![])
    };
    let arrays = inputs
        .iter()
        .map(|v| v.numeric())
        .collect::<Result<Vec<_>, _>>()?;
    let mut values = Vec::with_capacity(n);
    for i in 0..n {
        let mut value = arrays[0][i % arrays[0].len()];
        for array in &arrays[1..] {
            let next = array[i % array.len()];
            if value.is_nan() && !next.is_nan() {
                if remove_missing {
                    value = next;
                }
            } else if next.is_nan() {
                if !remove_missing || value.is_nan() {
                    value = next;
                }
            } else if if name == "pmax" {
                next > value
            } else {
                next < value
            } {
                value = next;
            }
        }
        values.push(value);
    }
    result.data = Data::Number { values, integer };
    Ok(result)
}
