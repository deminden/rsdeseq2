fn remove_formula_term<'a>(
    term: &str,
    factors: &'a [ExpandedFactorSpec<'a>],
    numeric_covariates: &'a [ExpandedNumericSpec<'a>],
    state: &mut ExpandedFormulaDesignState<'a>,
) -> Result<(), DeseqError> {
    if term.is_empty() {
        return Err(DeseqError::InvalidOptions {
            reason: "formula subtraction must be followed by a term".to_string(),
        });
    }
    if term == "1" {
        state.has_intercept = false;
        return Ok(());
    }
    if term == "0" {
        state.has_intercept = true;
        return Ok(());
    }
    if term == "-1" {
        state.has_intercept = true;
        return Ok(());
    }
    if formula_contains_top_level(term, '^')? {
        remove_power_formula_term(term, factors, numeric_covariates, state)?;
        return Ok(());
    }
    if formula_contains_top_level(term, '(')? || formula_contains_top_level(term, ')')? {
        return Err(DeseqError::InvalidOptions {
            reason: format!("formula term '{term}' is not supported by the primitive parser"),
        });
    }
    if formula_contains_top_level(term, '/')? {
        remove_nested_formula_term(term, factors, numeric_covariates, state)?;
        return Ok(());
    }
    if formula_contains_top_level(term, '*')? {
        remove_star_formula_term(term, factors, numeric_covariates, state)?;
        return Ok(());
    }
    if formula_contains_top_level(term, ':')? {
        remove_interaction_formula_term(term, factors, numeric_covariates, state)?;
        return Ok(());
    }
    remove_main_formula_term(term, factors, numeric_covariates, state)
}

fn add_star_formula_term<'a>(
    term: &str,
    factors: &'a [ExpandedFactorSpec<'a>],
    numeric_covariates: &'a [ExpandedNumericSpec<'a>],
    state: &mut ExpandedFormulaDesignState<'a>,
) -> Result<(), DeseqError> {
    let pieces = split_formula_pieces(term, '*')?;
    let alternatives =
        formula_piece_variable_alternatives(&pieces, factors, numeric_covariates)?;
    for alternative in &alternatives {
        for variable in alternative {
            add_main_formula_term(
                variable,
                factors,
                numeric_covariates,
                &mut state.selected_factors,
                &mut state.selected_numeric_covariates,
            )?;
        }
    }
    for left in 0..alternatives.len() {
        for right in (left + 1)..alternatives.len() {
            for left_variable in &alternatives[left] {
                for right_variable in &alternatives[right] {
                    if left_variable == right_variable {
                        continue;
                    }
                    add_pairwise_formula_interaction(
                        left_variable,
                        right_variable,
                        term,
                        factors,
                        numeric_covariates,
                        state,
                    )?;
                }
            }
        }
    }
    for order in 3..=alternatives.len() {
        for subset in formula_alternative_combinations(&alternatives, order) {
            for product in formula_variable_alternative_products(&subset) {
                if formula_variables_are_unique(&product) {
                    add_higher_order_formula_interaction(
                        &product,
                        term,
                        factors,
                        numeric_covariates,
                        state,
                    )?;
                }
            }
        }
    }
    Ok(())
}

fn add_nested_formula_term<'a>(
    term: &str,
    factors: &'a [ExpandedFactorSpec<'a>],
    numeric_covariates: &'a [ExpandedNumericSpec<'a>],
    state: &mut ExpandedFormulaDesignState<'a>,
) -> Result<(), DeseqError> {
    let pieces = split_formula_pieces(term, '/')?;
    let alternatives =
        formula_piece_variable_alternatives(&pieces, factors, numeric_covariates)?;
    for variable in &alternatives[0] {
        add_main_formula_term(
            variable,
            factors,
            numeric_covariates,
            &mut state.selected_factors,
            &mut state.selected_numeric_covariates,
        )?;
    }
    for prefix_len in 2..=alternatives.len() {
        for product in formula_variable_alternative_products(&alternatives[..prefix_len]) {
            if !formula_variables_are_unique(&product) {
                continue;
            }
            if product.len() == 2 {
                add_pairwise_formula_interaction(
                    product[0],
                    product[1],
                    term,
                    factors,
                    numeric_covariates,
                    state,
                )?;
            } else {
                add_higher_order_formula_interaction(
                    &product,
                    term,
                    factors,
                    numeric_covariates,
                    state,
                )?;
            }
        }
    }
    Ok(())
}

fn remove_star_formula_term<'a>(
    term: &str,
    factors: &'a [ExpandedFactorSpec<'a>],
    numeric_covariates: &'a [ExpandedNumericSpec<'a>],
    state: &mut ExpandedFormulaDesignState<'a>,
) -> Result<(), DeseqError> {
    let pieces = split_formula_pieces(term, '*')?;
    let alternatives =
        formula_piece_variable_alternatives(&pieces, factors, numeric_covariates)?;
    for alternative in &alternatives {
        for variable in alternative {
            remove_main_formula_term(variable, factors, numeric_covariates, state)?;
        }
    }
    for left in 0..alternatives.len() {
        for right in (left + 1)..alternatives.len() {
            for left_variable in &alternatives[left] {
                for right_variable in &alternatives[right] {
                    if left_variable == right_variable {
                        continue;
                    }
                    remove_pairwise_formula_interaction(
                        left_variable,
                        right_variable,
                        factors,
                        numeric_covariates,
                        state,
                    )?;
                }
            }
        }
    }
    for order in 3..=alternatives.len() {
        for subset in formula_alternative_combinations(&alternatives, order) {
            for product in formula_variable_alternative_products(&subset) {
                if formula_variables_are_unique(&product) {
                    remove_higher_order_formula_interaction(
                        &product,
                        factors,
                        numeric_covariates,
                        state,
                    )?;
                }
            }
        }
    }
    Ok(())
}

fn remove_nested_formula_term<'a>(
    term: &str,
    factors: &'a [ExpandedFactorSpec<'a>],
    numeric_covariates: &'a [ExpandedNumericSpec<'a>],
    state: &mut ExpandedFormulaDesignState<'a>,
) -> Result<(), DeseqError> {
    let pieces = split_formula_pieces(term, '/')?;
    let alternatives =
        formula_piece_variable_alternatives(&pieces, factors, numeric_covariates)?;
    for variable in &alternatives[0] {
        remove_main_formula_term(variable, factors, numeric_covariates, state)?;
    }
    for prefix_len in 2..=alternatives.len() {
        for product in formula_variable_alternative_products(&alternatives[..prefix_len]) {
            if !formula_variables_are_unique(&product) {
                continue;
            }
            if product.len() == 2 {
                remove_pairwise_formula_interaction(
                    product[0],
                    product[1],
                    factors,
                    numeric_covariates,
                    state,
                )?;
            } else {
                remove_higher_order_formula_interaction(
                    &product,
                    factors,
                    numeric_covariates,
                    state,
                )?;
            }
        }
    }
    Ok(())
}

fn add_power_formula_term<'a>(
    term: &str,
    factors: &'a [ExpandedFactorSpec<'a>],
    numeric_covariates: &'a [ExpandedNumericSpec<'a>],
    state: &mut ExpandedFormulaDesignState<'a>,
) -> Result<(), DeseqError> {
    apply_power_formula_term(term, factors, numeric_covariates, state, true)
}

fn remove_power_formula_term<'a>(
    term: &str,
    factors: &'a [ExpandedFactorSpec<'a>],
    numeric_covariates: &'a [ExpandedNumericSpec<'a>],
    state: &mut ExpandedFormulaDesignState<'a>,
) -> Result<(), DeseqError> {
    apply_power_formula_term(term, factors, numeric_covariates, state, false)
}

fn apply_power_formula_term<'a>(
    term: &str,
    factors: &'a [ExpandedFactorSpec<'a>],
    numeric_covariates: &'a [ExpandedNumericSpec<'a>],
    state: &mut ExpandedFormulaDesignState<'a>,
    add: bool,
) -> Result<(), DeseqError> {
    let pieces = split_formula_pieces(term, '^')?;
    if pieces.len() != 2 {
        return Err(DeseqError::InvalidOptions {
            reason: format!("formula power term '{term}' must have one exponent"),
        });
    }
    let base = strip_formula_outer_parentheses(&pieces[0])?;
    let variables = formula_piece_variables(base, factors, numeric_covariates)?;
    let order = parse_formula_state_interaction_power(term, &pieces[1])?;
    let max_order = order.min(variables.len());
    for interaction_order in 1..=max_order {
        for product in formula_variable_combinations(&variables, interaction_order) {
            if product.len() == 1 {
                if add {
                    add_main_formula_term(
                        product[0],
                        factors,
                        numeric_covariates,
                        &mut state.selected_factors,
                        &mut state.selected_numeric_covariates,
                    )?;
                } else {
                    remove_main_formula_term(product[0], factors, numeric_covariates, state)?;
                }
            } else if product.len() == 2 {
                if add {
                    add_pairwise_formula_interaction(
                        product[0],
                        product[1],
                        term,
                        factors,
                        numeric_covariates,
                        state,
                    )?;
                } else {
                    remove_pairwise_formula_interaction(
                        product[0],
                        product[1],
                        factors,
                        numeric_covariates,
                        state,
                    )?;
                }
            } else if add {
                add_higher_order_formula_interaction(
                    &product,
                    term,
                    factors,
                    numeric_covariates,
                    state,
                )?;
            } else {
                remove_higher_order_formula_interaction(&product, factors, numeric_covariates, state)?;
            }
        }
    }
    Ok(())
}

fn parse_formula_state_interaction_power(term: &str, exponent: &str) -> Result<usize, DeseqError> {
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

fn add_main_formula_term<'a>(
    term: &str,
    factors: &'a [ExpandedFactorSpec<'a>],
    numeric_covariates: &'a [ExpandedNumericSpec<'a>],
    selected_factors: &mut Vec<ExpandedFactorSpec<'a>>,
    selected_numeric_covariates: &mut Vec<ExpandedNumericSpec<'a>>,
) -> Result<(), DeseqError> {
    if term == "." {
        add_formula_dot_main_terms(
            factors,
            numeric_covariates,
            selected_factors,
            selected_numeric_covariates,
        );
        return Ok(());
    }
    let variable = resolved_formula_variable_name(term, factors, numeric_covariates)?;
    if let Some(factor) = factors
        .iter()
        .find(|candidate| candidate.factor == variable)
    {
        if selected_factors
            .iter()
            .any(|candidate| candidate.factor == factor.factor)
        {
            return Ok(());
        }
        selected_factors.push(factor.clone());
        return Ok(());
    }
    if let Some(covariate) = numeric_covariates
        .iter()
        .find(|candidate| candidate.name == variable)
    {
        if selected_numeric_covariates
            .iter()
            .any(|candidate| candidate.name == covariate.name)
        {
            return Ok(());
        }
        selected_numeric_covariates.push(covariate.clone());
        return Ok(());
    }
    Err(DeseqError::InvalidOptions {
        reason: format!("formula variable '{variable}' is not present in supplied design metadata"),
    })
}

fn remove_main_formula_term<'a>(
    term: &str,
    factors: &'a [ExpandedFactorSpec<'a>],
    numeric_covariates: &'a [ExpandedNumericSpec<'a>],
    state: &mut ExpandedFormulaDesignState<'a>,
) -> Result<(), DeseqError> {
    if term == "." {
        state.selected_factors.clear();
        state.selected_numeric_covariates.clear();
        return Ok(());
    }
    let variable = resolved_formula_variable_name(term, factors, numeric_covariates)?;
    if factors.iter().any(|candidate| candidate.factor == variable) {
        state
            .selected_factors
            .retain(|candidate| candidate.factor != variable);
        return Ok(());
    }
    if numeric_covariates
        .iter()
        .any(|candidate| candidate.name == variable)
    {
        state
            .selected_numeric_covariates
            .retain(|candidate| candidate.name != variable);
        return Ok(());
    }
    Err(DeseqError::InvalidOptions {
        reason: format!("formula variable '{variable}' is not present in supplied design metadata"),
    })
}

fn add_formula_dot_main_terms<'a>(
    factors: &'a [ExpandedFactorSpec<'a>],
    numeric_covariates: &'a [ExpandedNumericSpec<'a>],
    selected_factors: &mut Vec<ExpandedFactorSpec<'a>>,
    selected_numeric_covariates: &mut Vec<ExpandedNumericSpec<'a>>,
) {
    for factor in factors {
        if !selected_factors
            .iter()
            .any(|selected| selected.factor == factor.factor)
        {
            selected_factors.push(factor.clone());
        }
    }
    for covariate in numeric_covariates {
        if !selected_numeric_covariates
            .iter()
            .any(|selected| selected.name == covariate.name)
        {
            selected_numeric_covariates.push(covariate.clone());
        }
    }
}

fn formula_piece_variable_alternatives<'a>(
    pieces: &[String],
    factors: &'a [ExpandedFactorSpec<'a>],
    numeric_covariates: &'a [ExpandedNumericSpec<'a>],
) -> Result<Vec<Vec<&'a str>>, DeseqError> {
    pieces
        .iter()
        .map(|piece| formula_piece_variables(piece, factors, numeric_covariates))
        .collect()
}

fn formula_piece_variables<'a>(
    piece: &str,
    factors: &'a [ExpandedFactorSpec<'a>],
    numeric_covariates: &'a [ExpandedNumericSpec<'a>],
) -> Result<Vec<&'a str>, DeseqError> {
    if piece == "." {
        let mut variables = Vec::with_capacity(factors.len() + numeric_covariates.len());
        variables.extend(factors.iter().map(|factor| factor.factor));
        variables.extend(numeric_covariates.iter().map(|covariate| covariate.name));
        if variables.is_empty() {
            return Err(DeseqError::InvalidOptions {
                reason: "formula '.' requires at least one supplied design column".to_string(),
            });
        }
        return Ok(variables);
    }
    let variable = resolve_formula_variable(
        resolved_formula_variable_name(piece, factors, numeric_covariates)?,
        factors,
        numeric_covariates,
    )?;
    match variable {
        FormulaVariableRef::Factor(name) | FormulaVariableRef::Numeric(name) => Ok(vec![name]),
    }
}

fn formula_alternative_combinations<'a>(
    alternatives: &'a [Vec<&'a str>],
    order: usize,
) -> Vec<Vec<Vec<&'a str>>> {
    fn push_combinations<'a>(
        alternatives: &'a [Vec<&'a str>],
        order: usize,
        start: usize,
        current: &mut Vec<Vec<&'a str>>,
        output: &mut Vec<Vec<Vec<&'a str>>>,
    ) {
        if current.len() == order {
            output.push(current.clone());
            return;
        }
        let remaining = order - current.len();
        for idx in start..=alternatives.len() - remaining {
            current.push(alternatives[idx].clone());
            push_combinations(alternatives, order, idx + 1, current, output);
            current.pop();
        }
    }

    if order == 0 || order > alternatives.len() {
        return Vec::new();
    }
    let mut output = Vec::new();
    push_combinations(alternatives, order, 0, &mut Vec::new(), &mut output);
    output
}

fn formula_variable_alternative_products<'a>(alternatives: &[Vec<&'a str>]) -> Vec<Vec<&'a str>> {
    let mut products: Vec<Vec<&'a str>> = vec![Vec::new()];
    for alternative in alternatives {
        let mut next = Vec::new();
        for product in &products {
            for variable in alternative {
                let mut extended = product.clone();
                extended.push(*variable);
                next.push(extended);
            }
        }
        products = next;
    }
    products
}

fn formula_variable_combinations<'a>(variables: &[&'a str], order: usize) -> Vec<Vec<&'a str>> {
    fn push_combinations<'a>(
        variables: &[&'a str],
        order: usize,
        start: usize,
        current: &mut Vec<&'a str>,
        output: &mut Vec<Vec<&'a str>>,
    ) {
        if current.len() == order {
            output.push(current.clone());
            return;
        }
        let remaining = order - current.len();
        for idx in start..=variables.len() - remaining {
            current.push(variables[idx]);
            push_combinations(variables, order, idx + 1, current, output);
            current.pop();
        }
    }
    if order == 0 || order > variables.len() {
        return Vec::new();
    }
    let mut output = Vec::new();
    push_combinations(variables, order, 0, &mut Vec::new(), &mut output);
    output
}

fn formula_variables_are_unique(variables: &[&str]) -> bool {
    variables
        .iter()
        .enumerate()
        .all(|(idx, variable)| !variables[..idx].contains(variable))
}

fn add_interaction_formula_term<'a>(
    term: &str,
    factors: &'a [ExpandedFactorSpec<'a>],
    numeric_covariates: &'a [ExpandedNumericSpec<'a>],
    state: &mut ExpandedFormulaDesignState<'a>,
) -> Result<(), DeseqError> {
    let pieces = split_formula_pieces(term, ':')?;
    let alternatives =
        formula_piece_variable_alternatives(&pieces, factors, numeric_covariates)?;
    let has_dot = pieces.iter().any(|piece| piece == ".");
    let mut added_any = false;
    for product in formula_variable_alternative_products(&alternatives) {
        if !formula_variables_are_unique(&product) {
            continue;
        }
        if product.len() == 2 {
            add_pairwise_formula_interaction(
                product[0],
                product[1],
                term,
                factors,
                numeric_covariates,
                state,
            )?;
        } else {
            add_higher_order_formula_interaction(
                &product,
                term,
                factors,
                numeric_covariates,
                state,
            )?;
        }
        added_any = true;
    }
    if !added_any && !has_dot {
        return Err(DeseqError::InvalidOptions {
            reason: format!("formula interaction '{term}' cannot use one variable twice"),
        });
    }
    Ok(())
}

fn add_pairwise_formula_interaction<'a>(
    left: &str,
    right: &str,
    display_term: &str,
    factors: &'a [ExpandedFactorSpec<'a>],
    numeric_covariates: &'a [ExpandedNumericSpec<'a>],
    state: &mut ExpandedFormulaDesignState<'a>,
) -> Result<(), DeseqError> {
    let left_variable_name = resolved_formula_variable_name(left, factors, numeric_covariates)?;
    let right_variable_name = resolved_formula_variable_name(right, factors, numeric_covariates)?;
    if left_variable_name == right_variable_name {
        return Err(DeseqError::InvalidOptions {
            reason: format!("formula interaction '{display_term}' cannot use one variable twice"),
        });
    }

    let left_numeric = numeric_covariates
        .iter()
        .find(|candidate| candidate.name == left_variable_name);
    let right_numeric = numeric_covariates
        .iter()
        .find(|candidate| candidate.name == right_variable_name);
    let left_factor = factors
        .iter()
        .find(|candidate| candidate.factor == left_variable_name);
    let right_factor = factors
        .iter()
        .find(|candidate| candidate.factor == right_variable_name);

    match (left_factor, right_factor, left_numeric, right_numeric) {
        (Some(left_factor), Some(right_factor), None, None) => {
            if state.factor_interactions.iter().any(|interaction| {
                (interaction.left_factor == left_factor.factor
                    && interaction.right_factor == right_factor.factor)
                    || (interaction.left_factor == right_factor.factor
                        && interaction.right_factor == left_factor.factor)
            }) {
                return Ok(());
            }
            state
                .factor_interactions
                .push(ExpandedFactorInteractionSpec {
                    left_factor: left_factor.factor,
                    right_factor: right_factor.factor,
                });
        }
        (Some(factor), None, None, Some(numeric)) => {
            add_factor_numeric_formula_interaction(
                factor.factor,
                numeric.name,
                &mut state.factor_numeric_interactions,
                display_term,
            )?;
        }
        (None, Some(factor), Some(numeric), None) => {
            add_factor_numeric_formula_interaction(
                factor.factor,
                numeric.name,
                &mut state.factor_numeric_interactions,
                display_term,
            )?;
        }
        (None, None, Some(left_numeric), Some(right_numeric)) => {
            if state.numeric_interactions.iter().any(|interaction| {
                (interaction.left_numeric == left_numeric.name
                    && interaction.right_numeric == right_numeric.name)
                    || (interaction.left_numeric == right_numeric.name
                        && interaction.right_numeric == left_numeric.name)
            }) {
                return Ok(());
            }
            state
                .numeric_interactions
                .push(ExpandedNumericInteractionSpec {
                    left_numeric: left_numeric.name,
                    right_numeric: right_numeric.name,
                });
        }
        _ => {
            return Err(DeseqError::InvalidOptions {
                reason: format!(
                    "formula interaction '{display_term}' references variables missing from supplied design metadata"
                ),
            });
        }
    }
    Ok(())
}

fn remove_interaction_formula_term<'a>(
    term: &str,
    factors: &'a [ExpandedFactorSpec<'a>],
    numeric_covariates: &'a [ExpandedNumericSpec<'a>],
    state: &mut ExpandedFormulaDesignState<'a>,
) -> Result<(), DeseqError> {
    let pieces = split_formula_pieces(term, ':')?;
    let alternatives =
        formula_piece_variable_alternatives(&pieces, factors, numeric_covariates)?;
    let has_dot = pieces.iter().any(|piece| piece == ".");
    let mut removed_any = false;
    for product in formula_variable_alternative_products(&alternatives) {
        if !formula_variables_are_unique(&product) {
            continue;
        }
        if product.len() == 2 {
            remove_pairwise_formula_interaction(
                product[0],
                product[1],
                factors,
                numeric_covariates,
                state,
            )?;
        } else {
            remove_higher_order_formula_interaction(
                &product,
                factors,
                numeric_covariates,
                state,
            )?;
        }
        removed_any = true;
    }
    if !removed_any && !has_dot {
        return Err(DeseqError::InvalidOptions {
            reason: format!("formula interaction '{term}' cannot use one variable twice"),
        });
    }
    Ok(())
}

fn remove_pairwise_formula_interaction<'a>(
    left: &str,
    right: &str,
    factors: &'a [ExpandedFactorSpec<'a>],
    numeric_covariates: &'a [ExpandedNumericSpec<'a>],
    state: &mut ExpandedFormulaDesignState<'a>,
) -> Result<(), DeseqError> {
    let left_variable_name = resolved_formula_variable_name(left, factors, numeric_covariates)?;
    let right_variable_name = resolved_formula_variable_name(right, factors, numeric_covariates)?;
    if left_variable_name == right_variable_name {
        return Err(DeseqError::InvalidOptions {
            reason: format!("formula interaction '{left}:{right}' cannot use one variable twice"),
        });
    }
    let left_variable = resolve_formula_variable(left_variable_name, factors, numeric_covariates)?;
    let right_variable = resolve_formula_variable(right_variable_name, factors, numeric_covariates)?;
    match (&left_variable, &right_variable) {
        (FormulaVariableRef::Factor(left_factor), FormulaVariableRef::Factor(right_factor)) => {
            state.factor_interactions.retain(|interaction| {
                !((interaction.left_factor == *left_factor
                    && interaction.right_factor == *right_factor)
                    || (interaction.left_factor == *right_factor
                        && interaction.right_factor == *left_factor))
            });
        }
        (FormulaVariableRef::Factor(factor), FormulaVariableRef::Numeric(numeric))
        | (FormulaVariableRef::Numeric(numeric), FormulaVariableRef::Factor(factor)) => {
            state.factor_numeric_interactions.retain(|interaction| {
                !(interaction.factor == *factor && interaction.numeric == *numeric)
            });
        }
        (FormulaVariableRef::Numeric(left_numeric), FormulaVariableRef::Numeric(right_numeric)) => {
            state.numeric_interactions.retain(|interaction| {
                !((interaction.left_numeric == *left_numeric
                    && interaction.right_numeric == *right_numeric)
                    || (interaction.left_numeric == *right_numeric
                        && interaction.right_numeric == *left_numeric))
            });
        }
    }
    Ok(())
}

fn add_higher_order_formula_interaction<'a>(
    pieces: &[&str],
    display_term: &str,
    factors: &'a [ExpandedFactorSpec<'a>],
    numeric_covariates: &'a [ExpandedNumericSpec<'a>],
    state: &mut ExpandedFormulaDesignState<'a>,
) -> Result<(), DeseqError> {
    if pieces.len() < 3 {
        return Err(DeseqError::InvalidOptions {
            reason: format!(
                "formula interaction '{display_term}' must use at least three variables"
            ),
        });
    }
    let mut variables = Vec::with_capacity(pieces.len());
    for piece in pieces {
        let variable_name = resolved_formula_variable_name(piece, factors, numeric_covariates)?;
        if variables
            .iter()
            .any(|variable: &FormulaVariableRef<'_>| variable.name() == variable_name)
        {
            return Err(DeseqError::InvalidOptions {
                reason: format!(
                    "formula interaction '{display_term}' cannot use one variable twice"
                ),
            });
        }
        variables.push(resolve_formula_variable(variable_name, factors, numeric_covariates)?);
    }
    if state
        .higher_order_interactions
        .iter()
        .any(|interaction| same_formula_variable_set(&interaction.variables, &variables))
    {
        return Ok(());
    }
    state
        .higher_order_interactions
        .push(FormulaHigherOrderInteractionSpec { variables });
    Ok(())
}

fn remove_higher_order_formula_interaction<'a>(
    pieces: &[&str],
    factors: &'a [ExpandedFactorSpec<'a>],
    numeric_covariates: &'a [ExpandedNumericSpec<'a>],
    state: &mut ExpandedFormulaDesignState<'a>,
) -> Result<(), DeseqError> {
    if pieces.len() < 3 {
        return Err(DeseqError::InvalidOptions {
            reason: "formula interaction removal must use at least three variables".to_string(),
        });
    }
    let mut variables = Vec::with_capacity(pieces.len());
    for piece in pieces {
        let variable_name = resolved_formula_variable_name(piece, factors, numeric_covariates)?;
        if variables
            .iter()
            .any(|variable: &FormulaVariableRef<'_>| variable.name() == variable_name)
        {
            return Err(DeseqError::InvalidOptions {
                reason: format!(
                    "formula interaction '{}' cannot use one variable twice",
                    pieces.join(":")
                ),
            });
        }
        variables.push(resolve_formula_variable(variable_name, factors, numeric_covariates)?);
    }
    state
        .higher_order_interactions
        .retain(|interaction| !same_formula_variable_set(&interaction.variables, &variables));
    Ok(())
}

fn resolve_formula_variable<'a>(
    variable: &str,
    factors: &'a [ExpandedFactorSpec<'a>],
    numeric_covariates: &'a [ExpandedNumericSpec<'a>],
) -> Result<FormulaVariableRef<'a>, DeseqError> {
    let exact_factors = factors
        .iter()
        .filter(|candidate| candidate.factor == variable)
        .collect::<Vec<_>>();
    match exact_factors.as_slice() {
        [_] => {}
        [] => {}
        _ => {
            return Err(DeseqError::InvalidOptions {
                reason: format!("formula variable '{variable}' appears more than once"),
            });
        }
    }
    let exact_numeric = numeric_covariates
        .iter()
        .filter(|candidate| candidate.name == variable)
        .collect::<Vec<_>>();
    match exact_numeric.as_slice() {
        [_] => {}
        [] => {}
        _ => {
            return Err(DeseqError::InvalidOptions {
                reason: format!("formula variable '{variable}' appears more than once"),
            });
        }
    }
    if exact_factors.len() + exact_numeric.len() > 1 {
        return Err(DeseqError::InvalidOptions {
            reason: format!("formula variable '{variable}' appears more than once"),
        });
    }
    if let [factor] = exact_factors.as_slice() {
        return Ok(FormulaVariableRef::Factor(factor.factor));
    }
    if let [covariate] = exact_numeric.as_slice() {
        return Ok(FormulaVariableRef::Numeric(covariate.name));
    }

    let factor_aliases = factors
        .iter()
        .filter(|candidate| {
            r_like_name_candidates(candidate.factor)
                .into_iter()
                .any(|candidate| candidate == variable)
        })
        .collect::<Vec<_>>();
    let numeric_aliases = numeric_covariates
        .iter()
        .filter(|candidate| {
            r_like_name_candidates(candidate.name)
                .into_iter()
                .any(|candidate| candidate == variable)
        })
        .collect::<Vec<_>>();
    if factor_aliases.len() + numeric_aliases.len() > 1 {
        return Err(DeseqError::InvalidOptions {
            reason: format!("formula variable '{variable}' resolves ambiguously after R-style cleanup"),
        });
    }
    if let [factor] = factor_aliases.as_slice() {
        return Ok(FormulaVariableRef::Factor(factor.factor));
    }
    if let [covariate] = numeric_aliases.as_slice() {
        return Ok(FormulaVariableRef::Numeric(covariate.name));
    }
    Err(DeseqError::InvalidOptions {
        reason: format!("formula variable '{variable}' is not present in supplied design metadata"),
    })
}

fn resolved_formula_variable_name<'a>(
    variable: &str,
    factors: &'a [ExpandedFactorSpec<'a>],
    numeric_covariates: &'a [ExpandedNumericSpec<'a>],
) -> Result<&'a str, DeseqError> {
    let variable_name = formula_variable_name(variable).unwrap_or(variable.trim());
    let resolved = resolve_formula_variable(variable_name, factors, numeric_covariates)?;
    match resolved {
        FormulaVariableRef::Factor(name) | FormulaVariableRef::Numeric(name) => Ok(name),
    }
}

fn same_formula_variable_set(
    left: &[FormulaVariableRef<'_>],
    right: &[FormulaVariableRef<'_>],
) -> bool {
    left.len() == right.len()
        && left.iter().all(|left_variable| {
            right
                .iter()
                .any(|right_variable| right_variable == left_variable)
        })
}
