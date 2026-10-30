//! Identifier extraction for the supported native formula grammar.
use crate::errors::DeseqError;

// Run only after the numerical formula parser has accepted the expression.
// Quoted strings/argument names/functions/constants are not variable uses;
// removed terms and offset arguments remain uses of their original variables.
pub(super) fn formula_variables(formula: &str) -> Result<Vec<String>, DeseqError> {
    let chars: Vec<char> = formula.chars().collect();
    let mut variables = Vec::new();
    let mut i = 0;
    while i < chars.len() {
        let c = chars[i];
        if c == '\'' || c == '"' {
            let quote = c;
            i += 1;
            while i < chars.len() && chars[i] != quote {
                if chars[i] == '\\' {
                    i += 1;
                }
                i += 1;
            }
            if i == chars.len() {
                return Err(DeseqError::InvalidOptions {
                    reason: "unterminated formula origin string".into(),
                });
            }
            i += 1;
            continue;
        }
        if c == '%' {
            i += 1;
            while i < chars.len() && chars[i] != '%' {
                i += 1;
            }
            i += 1;
            continue;
        }
        if c.is_ascii_digit() || (c == '.' && chars.get(i + 1).is_some_and(char::is_ascii_digit)) {
            i += 1;
            while i < chars.len() {
                let c = chars[i];
                if c.is_ascii_alphanumeric()
                    || c == '.'
                    || ((c == '+' || c == '-') && matches!(chars[i - 1], 'e' | 'E' | 'p' | 'P'))
                {
                    i += 1;
                } else {
                    break;
                }
            }
            continue;
        }
        let quoted = c == '`';
        if quoted || c.is_alphabetic() || c == '_' || c == '.' {
            let start = i + usize::from(quoted);
            i = start;
            if quoted {
                while i < chars.len() && chars[i] != '`' {
                    i += 1;
                }
                if i == chars.len() {
                    return Err(DeseqError::InvalidOptions {
                        reason: "unterminated formula origin identifier".into(),
                    });
                }
            } else {
                while i < chars.len()
                    && (chars[i].is_alphanumeric() || chars[i] == '_' || chars[i] == '.')
                {
                    i += 1;
                }
            }
            let name: String = chars[start..i].iter().collect();
            i += usize::from(quoted);
            let mut next = i;
            while next < chars.len() && chars[next].is_whitespace() {
                next += 1;
            }
            let function_or_argument = matches!(chars.get(next), Some('(' | '='));
            let constant = !quoted
                && matches!(
                    name.as_str(),
                    "TRUE"
                        | "FALSE"
                        | "NULL"
                        | "NA"
                        | "NA_real_"
                        | "NA_integer_"
                        | "NA_character_"
                        | "NA_complex_"
                        | "Inf"
                        | "NaN"
                );
            if !function_or_argument && !constant && !variables.contains(&name) {
                variables.push(name);
            }
            continue;
        }
        i += 1;
    }
    Ok(variables)
}

#[cfg(test)]
mod tests {
    use super::formula_variables;

    #[test]
    fn original_identifiers_include_removed_terms_and_offsets_but_exclude_call_syntax() {
        // Expected identifier sets checked against base R 4.6.2 all.vars().
        let cases: &[(&str, &[&str])] = &[
            (
                "~ condition + nuisance - nuisance",
                &["condition", "nuisance"],
            ),
            (
                "~ condition + offset(log(exposure))",
                &["condition", "exposure"],
            ),
            (
                "~ relevel(x=factor(condition, levels=c('A','B')), ref='B')",
                &["condition"],
            ),
            (
                "~ I(1e-3 * dose + 2L) + poly(age, degree=2, raw=TRUE)",
                &["dose", "age"],
            ),
            (
                "~ `cell type` %in% batch + `cell type`",
                &["cell type", "batch"],
            ),
            ("~ 0 + condition", &["condition"]),
            ("~ . - condition", &[".", "condition"]),
            ("~ 1", &[]),
            ("~ factor(numeric)", &["numeric"]),
            (
                r#"~ factor(condition, labels=c("a\"b", "c\\d"))"#,
                &["condition"],
            ),
        ];
        for (formula, expected) in cases {
            assert_eq!(formula_variables(formula).unwrap(), *expected, "{formula}");
        }
        assert!(formula_variables("~ `unterminated").is_err());
        assert!(formula_variables("~ factor(x, 'unterminated").is_err());
    }
}
