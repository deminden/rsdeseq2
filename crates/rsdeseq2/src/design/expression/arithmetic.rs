use super::{Data, Value, invalid};
use crate::errors::DeseqError;

pub(super) fn power(x: f64, y: f64) -> f64 {
    if y == 2.0 {
        return x * x;
    }
    if x == 1.0 || y == 0.0 {
        return 1.0;
    }
    if x == 0.0 {
        return if y > 0.0 {
            0.0
        } else if y < 0.0 {
            f64::INFINITY
        } else {
            y
        };
    }
    if (-11.0..=11.0).contains(&x) {
        if y == 3.0 {
            return x * x * x;
        }
        if y == 4.0 {
            return x * x * x * x;
        }
    }
    if x.is_finite() && y.is_finite() {
        return x.powf(y);
    }
    if x.is_nan() || y.is_nan() {
        return x + y;
    }
    if x.is_infinite() {
        if x > 0.0 {
            return if y < 0.0 { 0.0 } else { f64::INFINITY };
        }
        if y.is_finite() && y.fract() == 0.0 {
            return if y < 0.0 {
                0.0
            } else if y % 2.0 != 0.0 {
                f64::NEG_INFINITY
            } else {
                f64::INFINITY
            };
        }
    }
    if y.is_infinite() && x >= 0.0 {
        return if y > 0.0 {
            if x >= 1.0 { f64::INFINITY } else { 0.0 }
        } else if x < 1.0 {
            f64::INFINITY
        } else {
            0.0
        };
    }
    f64::NAN
}
pub(super) fn unary(name: &str, value: &Value) -> Result<Value, DeseqError> {
    if value.has_class("factor") {
        return Err(invalid(
            "numeric formula operations require numeric covariates",
        ));
    }
    let input = value.numeric()?;
    if name == "!" {
        let mut result = value.clone();
        result.data = Data::Logical(
            input
                .iter()
                .map(|v| if v.is_nan() { None } else { Some(*v == 0.0) })
                .collect(),
        );
        return Ok(result);
    }
    let integer = matches!(
        value.data,
        Data::Number { integer: true, .. } | Data::Logical(_)
    );
    Ok(value.with_numbers(
        input
            .iter()
            .map(|v| if name == "-" { -v } else { *v })
            .collect(),
        integer,
    ))
}
pub(super) fn binary(name: &str, left: &Value, right: &Value) -> Result<Value, DeseqError> {
    let a = left.numeric()?;
    let b = right.numeric()?;
    let n = if a.is_empty() || b.is_empty() {
        0
    } else {
        a.len().max(b.len())
    };
    if left.dimensions.is_some()
        && right.dimensions.is_some()
        && left.dimensions != right.dimensions
    {
        return Err(invalid("non-conformable formula arrays"));
    }
    let template = if left.dimensions.is_some() {
        left
    } else if right.dimensions.is_some() || b.len() > a.len() {
        right
    } else {
        left
    };
    if template.dimensions.is_some() && template.len() != n {
        return Err(invalid("formula vector is longer than its matrix"));
    }
    let logical = matches!(name, "<" | "<=" | ">" | ">=" | "==" | "!=" | "&" | "|");
    let mut result = template.clone();
    if a.len() == b.len() {
        result.attributes = right.attributes.clone();
        for (name, value) in &left.attributes {
            result.set_attribute(name, value.clone());
        }
        result.names = left.names.clone().or_else(|| right.names.clone());
        if left.dimensions.is_some() && right.dimensions.is_some() {
            for i in 0..2 {
                result.dimnames[i] = left.dimnames[i]
                    .clone()
                    .or_else(|| right.dimnames[i].clone());
                result.dimname_names[i] = if left.dimnames[i].is_some() {
                    left.dimname_names[i].clone()
                } else {
                    right.dimname_names[i].clone()
                };
            }
        }
    }
    if logical {
        let bools = (0..n)
            .map(|i| {
                let (x, y) = (a[i % a.len()], b[i % b.len()]);
                if name == "&" {
                    if x == 0.0 || y == 0.0 {
                        return Some(false);
                    }
                    if x.is_nan() || y.is_nan() {
                        return None;
                    }
                    return Some(true);
                }
                if name == "|" {
                    if (!x.is_nan() && x != 0.0) || (!y.is_nan() && y != 0.0) {
                        return Some(true);
                    }
                    if x.is_nan() || y.is_nan() {
                        return None;
                    }
                    return Some(false);
                }
                if x.is_nan() || y.is_nan() {
                    None
                } else {
                    Some(match name {
                        "<" => x < y,
                        "<=" => x <= y,
                        ">" => x > y,
                        ">=" => x >= y,
                        "==" => x == y,
                        "!=" => x != y,
                        _ => unreachable!(),
                    })
                }
            })
            .collect();
        result.data = Data::Logical(bools);
    } else {
        let mut integer = matches!(
            &left.data,
            Data::Number { integer: true, .. } | Data::Logical(_)
        ) && matches!(
            &right.data,
            Data::Number { integer: true, .. } | Data::Logical(_)
        ) && !matches!(name, "/" | "^");
        let mut out = Vec::with_capacity(n);
        for i in 0..n {
            let (x, y) = (a[i % a.len()], b[i % b.len()]);
            let v = match name {
                "+" => x + y,
                "-" => x - y,
                "*" => x * y,
                "/" => x / y,
                "^" => power(x, y),
                "%%" => x - (x / y).floor() * y,
                "%/%" => {
                    let quotient = (x / y).floor();
                    if quotient == 0.0 { 0.0 } else { quotient }
                }
                _ => return Err(invalid("unsupported formula operator")),
            };
            out.push(
                if integer && (!v.is_finite() || v <= i32::MIN as f64 || v > i32::MAX as f64) {
                    super::NA_REAL
                } else {
                    if integer {
                        super::integer_storage(v)
                    } else {
                        v
                    }
                },
            );
        }
        if n == 0 {
            integer = false
        }
        result.data = Data::Number {
            values: out,
            integer,
        };
    }
    if result.dimensions.is_some() {
        result.names = None;
    }
    if n == 0 {
        result.dimensions = None;
        result.names = None;
        result.dimnames = [None, None];
        result.dimname_names = [None, None];
    }
    Ok(result)
}
