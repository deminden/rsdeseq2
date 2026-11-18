use super::*;
use crate::errors::invalid_formula;
// Shared transform arguments, options and string parsing.

pub(super) fn reject_duplicate_formula_poly_argument(
    seen: bool,
    argument: &str,
    expression: &str,
) -> Result<(), Error> {
    reject_duplicate_formula_transform_argument(seen, "poly", argument, expression)
}

pub(super) fn reject_duplicate_formula_transform_argument(
    seen: bool,
    transform: &str,
    argument: &str,
    expression: &str,
) -> Result<(), Error> {
    if seen {
        return Err(invalid_formula(format!(
            "formula transform '{transform}({expression})' argument '{argument}' matched multiple times"
        )));
    }
    Ok(())
}

pub(super) fn split_formula_named_argument(argument: &str) -> Option<(String, &str)> {
    let equals = formula_top_level_equals_index(argument)?;
    let (key, value_with_equals) = argument.split_at(equals);
    let value = &value_with_equals[1..];
    let key = key
        .split_whitespace()
        .collect::<String>()
        .to_ascii_lowercase();
    Some((key, value.trim()))
}

pub(super) fn formula_top_level_equals_index(argument: &str) -> Option<usize> {
    let mut depth = 0_i32;
    let mut in_backticks = false;
    let mut quote = None;
    let mut escaped_quote = false;
    for (idx, character) in argument.char_indices() {
        if let Some(active_quote) = quote {
            if escaped_quote {
                escaped_quote = false;
            } else if character == '\\' {
                escaped_quote = true;
            } else if character == active_quote {
                quote = None;
            }
            continue;
        }
        if matches!(character, '"' | '\'') && !in_backticks {
            quote = Some(character);
            continue;
        }
        if character == '`' {
            in_backticks = !in_backticks;
            continue;
        }
        if in_backticks {
            continue;
        }
        match character {
            '(' => depth += 1,
            ')' => depth -= 1,
            '=' if depth == 0 => return Some(idx),
            _ => {}
        }
    }
    None
}

pub(super) fn parse_formula_bool_option(
    value: &str,
    argument: &str,
    expression: &str,
) -> Result<bool, Error> {
    match value
        .split_whitespace()
        .collect::<String>()
        .to_ascii_lowercase()
        .as_str()
    {
        "true" | "t" => Ok(true),
        "false" | "f" => Ok(false),
        _ => Err(invalid_formula(format!(
            "formula transform 'poly({expression})' argument '{argument}' must be TRUE or FALSE"
        ))),
    }
}

pub(super) fn parse_formula_string_argument(
    value: &str,
    transform: &str,
    argument: &str,
    expression: &str,
) -> Result<String, Error> {
    let value = strip_formula_outer_parentheses(value.trim())?.trim();
    let parsed = if value.starts_with('"') || value.starts_with('\'') {
        parse_formula_quoted_string(value, transform, argument, expression)?
    } else {
        value.to_string()
    };
    if parsed.is_empty() || parsed.contains('`') {
        return Err(invalid_formula(format!(
            "formula transform '{transform}({expression})' argument '{argument}' must be a simple level string"
        )));
    }
    Ok(parsed)
}

pub(super) fn parse_formula_quoted_string(
    value: &str,
    transform: &str,
    argument: &str,
    expression: &str,
) -> Result<String, Error> {
    let quote = value
        .chars()
        .next()
        .ok_or_else(|| crate::errors::invariant("quoted formula string is not empty"))?;
    let Some(inner) = value
        .strip_prefix(quote)
        .and_then(|candidate| candidate.strip_suffix(quote))
    else {
        return Err(invalid_formula(format!(
            "formula transform '{transform}({expression})' argument '{argument}' must be a balanced quoted string"
        )));
    };
    let mut parsed = String::with_capacity(inner.len());
    let mut chars = inner.chars();
    while let Some(character) = chars.next() {
        if character == '\\' {
            let Some(escaped) = chars.next() else {
                return Err(invalid_formula(format!(
                    "formula transform '{transform}({expression})' argument '{argument}' has a trailing string escape"
                )));
            };
            parsed.push(match escaped {
                'n' => '\n',
                'r' => '\r',
                't' => '\t',
                other => other,
            });
        } else {
            parsed.push(character);
        }
    }
    Ok(parsed)
}

pub(super) fn parse_formula_string_vector_argument(
    value: &str,
    transform: &str,
    argument: &str,
    expression: &str,
) -> Result<Vec<String>, Error> {
    let value = strip_formula_outer_parentheses(value.trim())?.trim();
    if value.starts_with('\'') || value.starts_with('"') {
        return Ok(vec![parse_formula_string_argument(
            value, transform, argument, expression,
        )?]);
    }
    let Some(inner) = value
        .strip_prefix("c(")
        .and_then(|candidate| candidate.strip_suffix(')'))
    else {
        return Err(invalid_formula(format!(
            "formula transform '{transform}({expression})' argument '{argument}' must be c(\"level\", ...)"
        )));
    };
    let mut levels = Vec::new();
    for level in split_formula_transform_arguments(inner)? {
        let parsed = parse_formula_string_argument(&level, transform, argument, expression)?;
        if argument != "labels" && levels.iter().any(|candidate| candidate == &parsed) {
            return Err(invalid_formula(format!(
                "formula transform '{transform}({expression})' argument '{argument}' contains duplicate level '{parsed}'"
            )));
        }
        levels.push(parsed);
    }
    if levels.is_empty() {
        return Err(invalid_formula(format!(
            "formula transform '{transform}({expression})' argument '{argument}' must contain at least one level"
        )));
    }
    Ok(levels)
}

pub(super) fn formula_string_vector_label(values: &[String]) -> String {
    let quoted = values
        .iter()
        .map(|value| format!("\"{}\"", formula_quote_string(value)))
        .collect::<Vec<_>>();
    format!("c({})", quoted.join(", "))
}

pub(super) fn formula_quote_string(value: &str) -> String {
    let mut quoted = String::with_capacity(value.len());
    for character in value.chars() {
        match character {
            '\\' => quoted.push_str("\\\\"),
            '"' => quoted.push_str("\\\""),
            '\n' => quoted.push_str("\\n"),
            '\r' => quoted.push_str("\\r"),
            '\t' => quoted.push_str("\\t"),
            other => quoted.push(other),
        }
    }
    quoted
}
