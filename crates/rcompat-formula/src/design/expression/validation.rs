use super::{Data, Expression, Value, functions};
use crate::errors::invalid_formula;
use crate::errors::{Error, FormulaResource};

const MAX_EXPRESSION_DEPTH: usize = 128;

pub(super) fn function_name(function: &str) -> Result<&str, Error> {
    let name = function
        .strip_prefix("base::")
        .or_else(|| function.strip_prefix("stats::"))
        .unwrap_or(function);
    if !functions::supported(name) {
        return Err(Error::UnsupportedFormulaFunction {
            function: function.into(),
            expected_namespace: None,
        });
    }
    let namespace = if matches!(name, "poly" | "offset" | "relevel" | "sd") {
        "stats::"
    } else {
        "base::"
    };
    if function.contains("::") && !function.starts_with(namespace) {
        return Err(Error::UnsupportedFormulaFunction {
            function: function.into(),
            expected_namespace: Some(namespace),
        });
    }
    Ok(name)
}

pub(super) fn validate_expression_depth(depth: usize) -> Result<(), Error> {
    if depth > MAX_EXPRESSION_DEPTH {
        return Err(Error::FormulaLimitExceeded {
            resource: FormulaResource::ExpressionDepth,
            limit: MAX_EXPRESSION_DEPTH,
            observed: Some(depth),
        });
    }
    Ok(())
}

pub(super) fn validate_expression(expression: &Expression, depth: usize) -> Result<(), Error> {
    validate_expression_depth(depth)?;
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
pub(super) fn validate_value(value: &Value, depth: usize) -> Result<(), Error> {
    validate_value_with_limit(value, depth, MAX_EXPRESSION_DEPTH)
}
pub(super) fn validate_value_with_limit(
    value: &Value,
    depth: usize,
    limit: usize,
) -> Result<(), Error> {
    if depth > limit {
        return Err(Error::FormulaLimitExceeded {
            resource: FormulaResource::ValueDepth,
            limit,
            observed: Some(depth),
        });
    }
    if let Some([rows, columns]) = value.dimensions {
        if rows.checked_mul(columns) != Some(value.len()) {
            return Err(invalid_formula(
                "formula dimensions do not match stored values",
            ));
        }
        for (names, length) in value.dimnames.iter().zip([rows, columns]) {
            if names.as_ref().is_some_and(|n| n.len() != length) {
                return Err(invalid_formula(
                    "formula dimension names do not match dimensions",
                ));
            }
        }
        for i in 0..2 {
            if let Some(names) = &value.dimname_names[i]
                && value.dimnames[i]
                    .as_ref()
                    .is_none_or(|v| v.len() != names.len())
            {
                return Err(invalid_formula(
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
        return Err(invalid_formula(
            "formula dimension names require dimensions",
        ));
    }
    if value.names.as_ref().is_some_and(|n| n.len() != value.len()) {
        return Err(invalid_formula("formula names do not match stored values"));
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
        return Err(invalid_formula(
            "formula integer values are outside their range",
        ));
    }
    for (index, (name, attribute)) in value.attributes.iter().enumerate() {
        if name.is_empty() || value.attributes[..index].iter().any(|(n, _)| n == name) {
            return Err(invalid_formula(
                "formula attributes require unique nonempty names",
            ));
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
            return Err(invalid_formula("formula factor requires text levels"));
        };
        let Data::Number {
            values,
            integer: true,
        } = &value.data
        else {
            return Err(invalid_formula("formula factor requires integer codes"));
        };
        if levels
            .iter()
            .enumerate()
            .any(|(i, v)| levels[..i].contains(v))
            || values
                .iter()
                .any(|v| !v.is_nan() && (*v < 1.0 || *v > levels.len() as f64))
        {
            return Err(invalid_formula(
                "formula factor codes or levels are invalid",
            ));
        }
    }
    Ok(())
}
