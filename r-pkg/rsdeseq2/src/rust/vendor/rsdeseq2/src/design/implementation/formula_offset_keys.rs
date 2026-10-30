// R terms deduplicates language expressions before model.offset sums them.
// Compare lexical primitives, not evaluated vectors: different variables can
// contain identical values and must still contribute independently.
#[derive(PartialEq)]
enum FormulaOffsetToken {
    Symbol(String),
    Number(u64, bool),
    Text(String),
    Punctuation(char),
}

fn formula_offset_key(expression: &str) -> Result<Vec<FormulaOffsetToken>, DeseqError> {
    let mut tokens = Vec::new();
    let mut at = 0;
    while at < expression.len() {
        let c = expression[at..]
            .chars()
            .next()
            .expect("valid character boundary");
        if c.is_whitespace() {
            at += c.len_utf8();
            continue;
        }
        if matches!(c, '\'' | '"' | '`') {
            let start = at;
            at += 1;
            let mut escaped = false;
            let mut closed = false;
            while at < expression.len() {
                let next = expression[at..]
                    .chars()
                    .next()
                    .expect("valid character boundary");
                at += next.len_utf8();
                if escaped {
                    escaped = false;
                } else if next == '\\' {
                    escaped = true;
                } else if next == c {
                    closed = true;
                    break;
                }
            }
            if !closed {
                return Err(DeseqError::InvalidOptions {
                    reason: "formula offset has an unbalanced quoted token".into(),
                });
            }
            let text = &expression[start..at];
            tokens.push(if c == '`' {
                FormulaOffsetToken::Symbol(formula_variable_name(text)?.into())
            } else {
                FormulaOffsetToken::Text(parse_formula_quoted_string(
                    text,
                    "offset",
                    "expression",
                    expression,
                )?)
            });
            continue;
        }
        if c.is_ascii_digit()
            || c == '.' && expression[at + 1..].starts_with(|v: char| v.is_ascii_digit())
        {
            let start = at;
            at += c.len_utf8();
            while at < expression.len() {
                let next = expression.as_bytes()[at];
                if next.is_ascii_digit() || next == b'.' {
                    at += 1;
                } else if matches!(next, b'e' | b'E') {
                    at += 1;
                    if at < expression.len() && matches!(expression.as_bytes()[at], b'+' | b'-') {
                        at += 1;
                    }
                } else {
                    break;
                }
            }
            let number =
                expression[start..at]
                    .parse::<f64>()
                    .map_err(|_| DeseqError::InvalidOptions {
                        reason: "formula offset contains an invalid numeric token".into(),
                    })?;
            let integer = expression[at..].starts_with('L');
            if integer {
                at += 1;
            }
            tokens.push(FormulaOffsetToken::Number(number.to_bits(), integer));
            continue;
        }
        if c.is_alphabetic() || matches!(c, '.' | '_') {
            let start = at;
            at += c.len_utf8();
            while at < expression.len() {
                let next = expression[at..]
                    .chars()
                    .next()
                    .expect("valid character boundary");
                if !next.is_alphanumeric() && !matches!(next, '.' | '_') {
                    break;
                }
                at += next.len_utf8();
            }
            tokens.push(FormulaOffsetToken::Symbol(expression[start..at].into()));
            continue;
        }
        tokens.push(FormulaOffsetToken::Punctuation(c));
        at += c.len_utf8();
    }
    Ok(tokens)
}
