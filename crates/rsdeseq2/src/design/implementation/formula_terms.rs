fn expand_parenthesized_formula_terms(rhs: &str) -> Result<String, DeseqError> {
    let signed_terms = split_formula_signed_terms(rhs)?;
    let mut expanded = Vec::new();
    for (sign, term) in signed_terms {
        for expanded_term in expand_parenthesized_formula_term(&term)? {
            expanded.push(apply_formula_term_sign(sign, expanded_term));
        }
    }
    Ok(join_formula_terms(&expanded))
}

fn apply_formula_term_sign(sign: i8, term: String) -> String {
    if sign >= 0 {
        return term;
    }
    term.strip_prefix("- ")
        .map_or_else(|| format!("- {term}"), ToOwned::to_owned)
}

fn join_formula_terms(terms: &[String]) -> String {
    let mut joined = String::new();
    for term in terms {
        if joined.is_empty() {
            joined.push_str(term);
        } else if term.trim_start().starts_with('-') {
            joined.push(' ');
            joined.push_str(term);
        } else {
            joined.push_str(" + ");
            joined.push_str(term);
        }
    }
    joined
}

fn split_formula_signed_terms(rhs: &str) -> Result<Vec<(i8, String)>, DeseqError> {
    let mut terms = Vec::new();
    let mut depth = 0_i32;
    let mut in_backticks = false;
    let mut quote = None;
    let mut escaped_quote = false;
    let mut sign = 1_i8;
    let mut start = 0_usize;
    let mut preserve_negative_one_sentinel = false;
    for (idx, character) in rhs.char_indices() {
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
            ')' => {
                depth -= 1;
                if depth < 0 {
                    return Err(DeseqError::InvalidOptions {
                        reason: "formula parentheses are unbalanced".to_string(),
                    });
                }
            }
            '+' | '-' if depth == 0 => {
                let term = rhs[start..idx].trim();
                if term.is_empty() {
                    if character == '-' {
                        preserve_negative_one_sentinel = sign < 0;
                        sign = -1;
                    } else {
                        preserve_negative_one_sentinel = false;
                    }
                    start = idx + character.len_utf8();
                    continue;
                }
                let term = if sign < 0 && preserve_negative_one_sentinel && term == "1" {
                    "-1"
                } else {
                    term
                };
                terms.push((sign, term.to_string()));
                sign = if character == '-' { -1 } else { 1 };
                preserve_negative_one_sentinel = false;
                start = idx + character.len_utf8();
            }
            _ => {}
        }
    }
    if depth != 0 {
        return Err(DeseqError::InvalidOptions {
            reason: "formula parentheses are unbalanced".to_string(),
        });
    }
    if in_backticks {
        return Err(DeseqError::InvalidOptions {
            reason: "formula backtick-quoted variable name is unbalanced".to_string(),
        });
    }
    if quote.is_some() {
        return Err(DeseqError::InvalidOptions {
            reason: "formula quoted string is unbalanced".to_string(),
        });
    }
    let term = rhs[start..].trim();
    if term.is_empty() {
        return Err(DeseqError::InvalidOptions {
            reason: "formula contains an empty term".to_string(),
        });
    }
    let term = if sign < 0 && preserve_negative_one_sentinel && term == "1" {
        "-1"
    } else {
        term
    };
    terms.push((sign, term.to_string()));
    Ok(terms)
}

fn expand_parenthesized_formula_term(term: &str) -> Result<Vec<String>, DeseqError> {
    let term = strip_formula_outer_parentheses(term.trim())?;
    if term == "-1" {
        return Ok(vec![term.to_string()]);
    }
    let power_pieces = split_formula_top_level(term, '^')?;
    if power_pieces.len() > 1 {
        return expand_parenthesized_formula_power(term, &power_pieces);
    }
    let nested_in_pieces = split_formula_top_level_operator(term, "%in%")?;
    if nested_in_pieces.len() > 1 {
        return expand_parenthesized_formula_in_operator(term, &nested_in_pieces);
    }
    let additive_pieces = split_formula_top_level(term, '+')?;
    if additive_pieces.len() > 1 {
        return split_formula_additive_group(term);
    }
    for delimiter in ['*', ':', '/'] {
        let pieces = split_formula_top_level(term, delimiter)?;
        if pieces.len() > 1 {
            return expand_parenthesized_formula_operator(&pieces, delimiter);
        }
    }
    if formula_contains_unquoted_char(term, '-') {
        return expand_signed_additive_formula_group(term);
    }
    if formula_contains_unquoted_parentheses(term) {
        return split_formula_additive_group(term);
    }
    Ok(vec![term.to_string()])
}

fn expand_parenthesized_formula_in_operator(
    term: &str,
    pieces: &[String],
) -> Result<Vec<String>, DeseqError> {
    if pieces.len() != 2 {
        return Err(DeseqError::InvalidOptions {
            reason: format!("formula nesting term '{term}' must use one '%in%' operator"),
        });
    }
    let inner = split_formula_additive_group(&pieces[0])?;
    let outer = split_formula_additive_group(&pieces[1])?;
    let alternatives = [outer, inner];
    Ok(formula_alternative_products(&alternatives)
        .into_iter()
        .map(|product| product.join(":"))
        .collect())
}

fn expand_parenthesized_formula_operator(
    pieces: &[String],
    delimiter: char,
) -> Result<Vec<String>, DeseqError> {
    let mut alternatives = Vec::with_capacity(pieces.len());
    for piece in pieces {
        alternatives.push(split_formula_additive_group(piece)?);
    }
    match delimiter {
        '*' => expand_parenthesized_star(&alternatives),
        ':' => Ok(formula_alternative_products(&alternatives)
            .into_iter()
            .map(|product| product.join(":"))
            .collect()),
        '/' => {
            let mut terms = Vec::new();
            for product in formula_alternative_products(&alternatives) {
                for prefix_len in 1..=product.len() {
                    push_unique_formula_term(&mut terms, product[..prefix_len].join(":"));
                }
            }
            Ok(terms)
        }
        _ => unreachable!("unsupported formula delimiter"),
    }
}

fn expand_parenthesized_star(alternatives: &[Vec<String>]) -> Result<Vec<String>, DeseqError> {
    let mut terms = Vec::new();
    for group in alternatives {
        for term in group {
            push_unique_formula_term(&mut terms, term.clone());
        }
    }
    for order in 2..=alternatives.len() {
        for group_subset in formula_group_combinations(alternatives, order) {
            for product in formula_alternative_products(&group_subset) {
                push_unique_formula_term(&mut terms, product.join(":"));
            }
        }
    }
    Ok(terms)
}

fn expand_parenthesized_formula_power(
    term: &str,
    pieces: &[String],
) -> Result<Vec<String>, DeseqError> {
    if pieces.len() != 2 {
        return Err(DeseqError::InvalidOptions {
            reason: format!("formula power term '{term}' must have one exponent"),
        });
    }
    let terms = split_formula_power_base_terms(&pieces[0])?;
    if terms.iter().any(|term| term == ".") {
        return Ok(vec![term.to_string()]);
    }
    let order = parse_formula_interaction_power(term, &pieces[1])?;
    let max_order = order.min(terms.len());
    let mut expanded = Vec::new();
    for interaction_order in 1..=max_order {
        for combination in formula_term_combinations(&terms, interaction_order) {
            push_unique_formula_term(&mut expanded, combination.join(":"));
        }
    }
    Ok(expanded)
}

fn split_formula_power_base_terms(term: &str) -> Result<Vec<String>, DeseqError> {
    let stripped = strip_formula_outer_parentheses(term.trim())?;
    if !formula_contains_unquoted_char(stripped, '-') {
        return split_formula_additive_group(stripped);
    }
    simplify_signed_additive_formula_group(stripped)
}

fn parse_formula_interaction_power(term: &str, exponent: &str) -> Result<usize, DeseqError> {
    let exponent = exponent.trim();
    if exponent.is_empty() {
        return Err(DeseqError::InvalidOptions {
            reason: format!("formula power term '{term}' has an empty exponent"),
        });
    }
    let order = exponent
        .parse::<usize>()
        .map_err(|_| DeseqError::InvalidOptions {
            reason: format!(
                "formula power term '{term}' requires a positive integer exponent"
            ),
        })?;
    if order == 0 {
        return Err(DeseqError::InvalidOptions {
            reason: format!(
                "formula power term '{term}' requires a positive integer exponent"
            ),
        });
    }
    Ok(order)
}

fn split_formula_additive_group(term: &str) -> Result<Vec<String>, DeseqError> {
    let stripped = strip_formula_outer_parentheses(term.trim())?;
    if formula_contains_unquoted_char(stripped, '-') {
        return simplify_signed_additive_formula_group(stripped);
    }
    let pieces = split_formula_top_level(stripped, '+')?;
    if pieces.len() == 1 {
        if formula_contains_unquoted_parentheses(stripped) {
            return expand_parenthesized_formula_term(stripped);
        }
        return Ok(vec![stripped.to_string()]);
    }
    let mut terms = Vec::new();
    for piece in pieces {
        if formula_contains_unquoted_char(&piece, '-') {
            for expanded in simplify_signed_additive_formula_group(&piece)? {
                push_unique_formula_term(&mut terms, expanded);
            }
            continue;
        }
        for expanded in expand_parenthesized_formula_term(&piece)? {
            push_unique_formula_term(&mut terms, expanded);
        }
    }
    Ok(terms)
}

fn simplify_signed_additive_formula_group(term: &str) -> Result<Vec<String>, DeseqError> {
    let stripped = strip_formula_outer_parentheses(term.trim())?;
    let signed_terms = split_formula_signed_terms(stripped)?;
    let mut terms = Vec::new();
    for (sign, signed_term) in signed_terms {
        if sign < 0 {
            if let Ok(expanded_terms) = expand_parenthesized_formula_term(&signed_term) {
                for expanded in expanded_terms {
                    remove_signed_formula_group_term(&mut terms, &expanded);
                }
            }
            continue;
        }
        for expanded in expand_parenthesized_formula_term(&signed_term)? {
            push_unique_formula_term(&mut terms, expanded);
        }
    }
    Ok(terms)
}

fn remove_signed_formula_group_term(terms: &mut Vec<String>, term: &str) {
    terms.retain(|candidate| candidate != term);
}

fn expand_signed_additive_formula_group(term: &str) -> Result<Vec<String>, DeseqError> {
    let stripped = strip_formula_outer_parentheses(term.trim())?;
    let signed_terms = split_formula_signed_terms(stripped)?;
    if signed_terms.len() == 1 && signed_terms[0].0 >= 0 && signed_terms[0].1.trim() == stripped {
        return Err(DeseqError::InvalidOptions {
            reason: format!("formula term '{term}' contains unsupported nested subtraction"),
        });
    }
    let mut terms = Vec::new();
    for (sign, signed_term) in signed_terms {
        for expanded in expand_parenthesized_formula_term(&signed_term)? {
            push_unique_formula_term(&mut terms, apply_formula_term_sign(sign, expanded));
        }
    }
    Ok(terms)
}

fn strip_formula_outer_parentheses(term: &str) -> Result<&str, DeseqError> {
    let mut stripped = term.trim();
    loop {
        if !(stripped.starts_with('(') && stripped.ends_with(')')) {
            return Ok(stripped);
        }
        let mut depth = 0_i32;
        let mut quote = None;
        let mut escaped_quote = false;
        let mut encloses_whole_term = true;
        for (idx, character) in stripped.char_indices() {
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
            if matches!(character, '"' | '\'') {
                quote = Some(character);
                continue;
            }
            match character {
                '(' => depth += 1,
                ')' => {
                    depth -= 1;
                    if depth < 0 {
                        return Err(DeseqError::InvalidOptions {
                            reason: "formula parentheses are unbalanced".to_string(),
                        });
                    }
                    if depth == 0 && idx + character.len_utf8() < stripped.len() {
                        encloses_whole_term = false;
                        break;
                    }
                }
                _ => {}
            }
        }
        if depth != 0 {
            return Err(DeseqError::InvalidOptions {
                reason: "formula parentheses are unbalanced".to_string(),
            });
        }
        if quote.is_some() {
            return Err(DeseqError::InvalidOptions {
                reason: "formula quoted string is unbalanced".to_string(),
            });
        }
        if !encloses_whole_term {
            return Ok(stripped);
        }
        stripped = stripped[1..stripped.len() - 1].trim();
    }
}

fn split_formula_top_level(term: &str, delimiter: char) -> Result<Vec<String>, DeseqError> {
    let mut pieces = Vec::new();
    let mut depth = 0_i32;
    let mut in_backticks = false;
    let mut quote = None;
    let mut escaped_quote = false;
    let mut start = 0_usize;
    for (idx, character) in term.char_indices() {
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
            ')' => {
                depth -= 1;
                if depth < 0 {
                    return Err(DeseqError::InvalidOptions {
                        reason: "formula parentheses are unbalanced".to_string(),
                    });
                }
            }
            _ if character == delimiter && depth == 0 => {
                let piece = term[start..idx].trim();
                if piece.is_empty() {
                    return Err(DeseqError::InvalidOptions {
                        reason: format!("formula term '{term}' contains an empty component"),
                    });
                }
                pieces.push(piece.to_string());
                start = idx + character.len_utf8();
            }
            _ => {}
        }
    }
    if depth != 0 {
        return Err(DeseqError::InvalidOptions {
            reason: "formula parentheses are unbalanced".to_string(),
        });
    }
    if in_backticks {
        return Err(DeseqError::InvalidOptions {
            reason: "formula backtick-quoted variable name is unbalanced".to_string(),
        });
    }
    if quote.is_some() {
        return Err(DeseqError::InvalidOptions {
            reason: "formula quoted string is unbalanced".to_string(),
        });
    }
    let piece = term[start..].trim();
    if piece.is_empty() {
        return Err(DeseqError::InvalidOptions {
            reason: format!("formula term '{term}' contains an empty component"),
        });
    }
    pieces.push(piece.to_string());
    Ok(pieces)
}

fn formula_contains_unquoted_parentheses(term: &str) -> bool {
    let mut in_backticks = false;
    let mut quote = None;
    let mut escaped_quote = false;
    for character in term.chars() {
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
        if !in_backticks && matches!(character, '(' | ')') {
            return true;
        }
    }
    false
}

fn formula_contains_unquoted_char(term: &str, needle: char) -> bool {
    let mut in_backticks = false;
    let mut quote = None;
    let mut escaped_quote = false;
    for character in term.chars() {
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
        if !in_backticks && character == needle {
            return true;
        }
    }
    false
}

fn formula_contains_top_level(term: &str, delimiter: char) -> Result<bool, DeseqError> {
    Ok(split_formula_top_level(term, delimiter)?.len() > 1)
}

fn split_formula_top_level_operator(
    term: &str,
    operator: &str,
) -> Result<Vec<String>, DeseqError> {
    let mut pieces = Vec::new();
    let mut depth = 0_i32;
    let mut in_backticks = false;
    let mut quote = None;
    let mut escaped_quote = false;
    let mut start = 0_usize;
    let mut idx = 0_usize;
    while idx < term.len() {
        let Some(character) = term[idx..].chars().next() else {
            break;
        };
        if let Some(active_quote) = quote {
            if escaped_quote {
                escaped_quote = false;
            } else if character == '\\' {
                escaped_quote = true;
            } else if character == active_quote {
                quote = None;
            }
            idx += character.len_utf8();
            continue;
        }
        if matches!(character, '"' | '\'') && !in_backticks {
            quote = Some(character);
            idx += character.len_utf8();
            continue;
        }
        if character == '`' {
            in_backticks = !in_backticks;
            idx += character.len_utf8();
            continue;
        }
        if in_backticks {
            idx += character.len_utf8();
            continue;
        }
        match character {
            '(' => depth += 1,
            ')' => {
                depth -= 1;
                if depth < 0 {
                    return Err(DeseqError::InvalidOptions {
                        reason: "formula parentheses are unbalanced".to_string(),
                    });
                }
            }
            _ if depth == 0 && term[idx..].starts_with(operator) => {
                let piece = term[start..idx].trim();
                if piece.is_empty() {
                    return Err(DeseqError::InvalidOptions {
                        reason: format!("formula term '{term}' contains an empty component"),
                    });
                }
                pieces.push(piece.to_string());
                idx += operator.len();
                start = idx;
                continue;
            }
            _ => {}
        }
        idx += character.len_utf8();
    }
    if depth != 0 {
        return Err(DeseqError::InvalidOptions {
            reason: "formula parentheses are unbalanced".to_string(),
        });
    }
    if in_backticks {
        return Err(DeseqError::InvalidOptions {
            reason: "formula backtick-quoted variable name is unbalanced".to_string(),
        });
    }
    if quote.is_some() {
        return Err(DeseqError::InvalidOptions {
            reason: "formula quoted string is unbalanced".to_string(),
        });
    }
    let piece = term[start..].trim();
    if piece.is_empty() {
        return Err(DeseqError::InvalidOptions {
            reason: format!("formula term '{term}' contains an empty component"),
        });
    }
    pieces.push(piece.to_string());
    Ok(pieces)
}

fn formula_alternative_products(alternatives: &[Vec<String>]) -> Vec<Vec<String>> {
    let mut products: Vec<Vec<String>> = vec![Vec::new()];
    for group in alternatives {
        let mut next = Vec::new();
        for prefix in &products {
            for term in group {
                let mut product = prefix.clone();
                product.push(term.clone());
                next.push(product);
            }
        }
        products = next;
    }
    products
}

fn formula_group_combinations(groups: &[Vec<String>], order: usize) -> Vec<Vec<Vec<String>>> {
    fn push_group_combinations(
        groups: &[Vec<String>],
        order: usize,
        start: usize,
        current: &mut Vec<Vec<String>>,
        output: &mut Vec<Vec<Vec<String>>>,
    ) {
        if current.len() == order {
            output.push(current.clone());
            return;
        }
        let remaining = order - current.len();
        for idx in start..=groups.len() - remaining {
            current.push(groups[idx].clone());
            push_group_combinations(groups, order, idx + 1, current, output);
            current.pop();
        }
    }
    if order == 0 || order > groups.len() {
        return Vec::new();
    }
    let mut output = Vec::new();
    push_group_combinations(groups, order, 0, &mut Vec::new(), &mut output);
    output
}

fn formula_term_combinations(terms: &[String], order: usize) -> Vec<Vec<String>> {
    fn push_term_combinations(
        terms: &[String],
        order: usize,
        start: usize,
        current: &mut Vec<String>,
        output: &mut Vec<Vec<String>>,
    ) {
        if current.len() == order {
            output.push(current.clone());
            return;
        }
        let remaining = order - current.len();
        for idx in start..=terms.len() - remaining {
            current.push(terms[idx].clone());
            push_term_combinations(terms, order, idx + 1, current, output);
            current.pop();
        }
    }
    if order == 0 || order > terms.len() {
        return Vec::new();
    }
    let mut output = Vec::new();
    push_term_combinations(terms, order, 0, &mut Vec::new(), &mut output);
    output
}

fn push_unique_formula_term(terms: &mut Vec<String>, term: String) {
    if !terms.iter().any(|candidate| candidate == &term) {
        terms.push(term);
    }
}
