use super::{Data, Expression, Value, functions, invalid};
use crate::errors::DeseqError;

pub(super) fn function_name(function: &str) -> Result<&str, DeseqError> {
    let name = function
        .strip_prefix("base::")
        .or_else(|| function.strip_prefix("stats::"))
        .unwrap_or(function);
    if !functions::supported(name) {
        return Err(invalid(format!(
            "custom or unsupported formula function '{function}'; provide a supported built-in transform"
        )));
    }
    let namespace = if matches!(name, "poly" | "offset" | "relevel" | "sd") {
        "stats::"
    } else {
        "base::"
    };
    if function.contains("::") && !function.starts_with(namespace) {
        return Err(invalid(format!(
            "unsupported formula namespace for '{function}'"
        )));
    }
    Ok(name)
}

pub(super) fn validate_expression(expression: &Expression, depth: usize) -> Result<(), DeseqError> {
    if depth > 128 {
        return Err(invalid("formula expression nesting exceeds 128"));
    }
    match expression {
        Expression::Literal(value) => validate_value(value, depth),
        Expression::Symbol(_) => Ok(()),
        Expression::Call {
            function,
            arguments,
        } => {
            function_name(function)?;
            for argument in arguments {
                validate_expression(&argument.expression, depth + 1)?;
            }
            Ok(())
        }
    }
}

// Public primitive inputs must be checked before any matrix indexing or factor
// lookup. This also bounds recursive metadata independently of AST recursion.
pub(super) fn validate_value(value: &Value, depth: usize) -> Result<(), DeseqError> {
    validate_value_with_limit(value, depth, 128)
}
pub(super) fn validate_value_with_limit(
    value: &Value,
    depth: usize,
    limit: usize,
) -> Result<(), DeseqError> {
    if depth > limit {
        return Err(invalid(format!("formula value nesting exceeds {limit}")));
    }
    if let Some([rows, columns]) = value.dimensions {
        if rows.checked_mul(columns) != Some(value.len()) {
            return Err(invalid("formula dimensions do not match stored values"));
        }
        for (names, length) in value.dimnames.iter().zip([rows, columns]) {
            if names.as_ref().is_some_and(|n| n.len() != length) {
                return Err(invalid("formula dimension names do not match dimensions"));
            }
        }
        for i in 0..2 {
            if let Some(names) = &value.dimname_names[i]
                && value.dimnames[i]
                    .as_ref()
                    .is_none_or(|v| v.len() != names.len())
            {
                return Err(invalid(
                    "formula label names require matching dimension labels",
                ));
            }
        }
    } else if value
        .dimnames
        .iter()
        .chain(&value.dimname_names)
        .any(Option::is_some)
    {
        return Err(invalid("formula dimension names require dimensions"));
    }
    if value.names.as_ref().is_some_and(|n| n.len() != value.len()) {
        return Err(invalid("formula names do not match stored values"));
    }
    if let Data::Number {
        values,
        integer: true,
    } = &value.data
        && values.iter().any(|v| {
            !v.is_nan()
                && (!v.is_finite()
                    || v.fract() != 0.0
                    || *v <= i32::MIN as f64
                    || *v > i32::MAX as f64)
        })
    {
        return Err(invalid("formula integer values are outside their range"));
    }
    for (index, (name, attribute)) in value.attributes.iter().enumerate() {
        if name.is_empty() || value.attributes[..index].iter().any(|(n, _)| n == name) {
            return Err(invalid("formula attributes require unique nonempty names"));
        }
        validate_value_with_limit(attribute, depth + 1, limit)?;
    }
    if let Data::List(values) = &value.data {
        for v in values {
            validate_value_with_limit(v, depth + 1, limit)?;
        }
    }
    if value.has_class("factor") {
        let Some(Value {
            data: Data::Text(levels),
            ..
        }) = value.attribute("levels")
        else {
            return Err(invalid("formula factor requires text levels"));
        };
        let Data::Number {
            values,
            integer: true,
        } = &value.data
        else {
            return Err(invalid("formula factor requires integer codes"));
        };
        if levels
            .iter()
            .enumerate()
            .any(|(i, v)| levels[..i].contains(v))
            || values
                .iter()
                .any(|v| !v.is_nan() && (*v < 1.0 || *v > levels.len() as f64))
        {
            return Err(invalid("formula factor codes or levels are invalid"));
        }
    }
    Ok(())
}
