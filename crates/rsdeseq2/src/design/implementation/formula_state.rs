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
    if term == "0" || term == "-1" {
        return Ok(());
    }
    if term.contains('^') || term.contains('(') || term.contains(')') {
        return Err(DeseqError::InvalidOptions {
            reason: format!("formula term '{term}' is not supported by the primitive parser"),
        });
    }
    if term.contains('/') {
        remove_nested_formula_term(term, factors, numeric_covariates, state)?;
        return Ok(());
    }
    if term.contains('*') {
        remove_star_formula_term(term, factors, numeric_covariates, state)?;
        return Ok(());
    }
    if term.contains(':') {
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
    for piece in &pieces {
        add_main_formula_term(
            piece,
            factors,
            numeric_covariates,
            &mut state.selected_factors,
            &mut state.selected_numeric_covariates,
        )?;
    }
    for left in 0..pieces.len() {
        for right in (left + 1)..pieces.len() {
            add_pairwise_formula_interaction(
                pieces[left],
                pieces[right],
                term,
                factors,
                numeric_covariates,
                state,
            )?;
        }
    }
    for order in 3..=pieces.len() {
        for subset in formula_piece_combinations(&pieces, order) {
            add_higher_order_formula_interaction(
                &subset,
                term,
                factors,
                numeric_covariates,
                state,
            )?;
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
    add_main_formula_term(
        pieces[0],
        factors,
        numeric_covariates,
        &mut state.selected_factors,
        &mut state.selected_numeric_covariates,
    )?;
    add_pairwise_formula_interaction(
        pieces[0],
        pieces[1],
        term,
        factors,
        numeric_covariates,
        state,
    )?;
    for prefix_len in 3..=pieces.len() {
        add_higher_order_formula_interaction(
            &pieces[..prefix_len],
            term,
            factors,
            numeric_covariates,
            state,
        )?;
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
    for piece in &pieces {
        remove_main_formula_term(piece, factors, numeric_covariates, state)?;
    }
    for left in 0..pieces.len() {
        for right in (left + 1)..pieces.len() {
            remove_pairwise_formula_interaction(
                pieces[left],
                pieces[right],
                factors,
                numeric_covariates,
                state,
            )?;
        }
    }
    for order in 3..=pieces.len() {
        for subset in formula_piece_combinations(&pieces, order) {
            remove_higher_order_formula_interaction(&subset, factors, numeric_covariates, state)?;
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
    remove_main_formula_term(pieces[0], factors, numeric_covariates, state)?;
    remove_pairwise_formula_interaction(pieces[0], pieces[1], factors, numeric_covariates, state)?;
    for prefix_len in 3..=pieces.len() {
        remove_higher_order_formula_interaction(
            &pieces[..prefix_len],
            factors,
            numeric_covariates,
            state,
        )?;
    }
    Ok(())
}

fn add_main_formula_term<'a>(
    term: &str,
    factors: &'a [ExpandedFactorSpec<'a>],
    numeric_covariates: &'a [ExpandedNumericSpec<'a>],
    selected_factors: &mut Vec<ExpandedFactorSpec<'a>>,
    selected_numeric_covariates: &mut Vec<ExpandedNumericSpec<'a>>,
) -> Result<(), DeseqError> {
    validate_formula_variable(term)?;
    if let Some(factor) = factors.iter().find(|candidate| candidate.factor == term) {
        if selected_factors
            .iter()
            .any(|candidate| candidate.factor == term)
        {
            return Err(DeseqError::InvalidOptions {
                reason: format!("formula main effect '{term}' appears more than once"),
            });
        }
        selected_factors.push(factor.clone());
        return Ok(());
    }
    if let Some(covariate) = numeric_covariates
        .iter()
        .find(|candidate| candidate.name == term)
    {
        if selected_numeric_covariates
            .iter()
            .any(|candidate| candidate.name == term)
        {
            return Err(DeseqError::InvalidOptions {
                reason: format!("formula main effect '{term}' appears more than once"),
            });
        }
        selected_numeric_covariates.push(covariate.clone());
        return Ok(());
    }
    Err(DeseqError::InvalidOptions {
        reason: format!("formula variable '{term}' is not present in supplied design metadata"),
    })
}

fn remove_main_formula_term<'a>(
    term: &str,
    factors: &'a [ExpandedFactorSpec<'a>],
    numeric_covariates: &'a [ExpandedNumericSpec<'a>],
    state: &mut ExpandedFormulaDesignState<'a>,
) -> Result<(), DeseqError> {
    validate_formula_variable(term)?;
    if factors.iter().any(|candidate| candidate.factor == term) {
        state
            .selected_factors
            .retain(|candidate| candidate.factor != term);
        return Ok(());
    }
    if numeric_covariates
        .iter()
        .any(|candidate| candidate.name == term)
    {
        state
            .selected_numeric_covariates
            .retain(|candidate| candidate.name != term);
        return Ok(());
    }
    Err(DeseqError::InvalidOptions {
        reason: format!("formula variable '{term}' is not present in supplied design metadata"),
    })
}

fn add_interaction_formula_term<'a>(
    term: &str,
    factors: &'a [ExpandedFactorSpec<'a>],
    numeric_covariates: &'a [ExpandedNumericSpec<'a>],
    state: &mut ExpandedFormulaDesignState<'a>,
) -> Result<(), DeseqError> {
    let pieces = split_formula_pieces(term, ':')?;
    if pieces.len() == 2 {
        add_pairwise_formula_interaction(
            pieces[0],
            pieces[1],
            term,
            factors,
            numeric_covariates,
            state,
        )
    } else {
        add_higher_order_formula_interaction(&pieces, term, factors, numeric_covariates, state)
    }
}

fn add_pairwise_formula_interaction<'a>(
    left: &str,
    right: &str,
    display_term: &str,
    factors: &'a [ExpandedFactorSpec<'a>],
    numeric_covariates: &'a [ExpandedNumericSpec<'a>],
    state: &mut ExpandedFormulaDesignState<'a>,
) -> Result<(), DeseqError> {
    validate_formula_variable(left)?;
    validate_formula_variable(right)?;
    if left == right {
        return Err(DeseqError::InvalidOptions {
            reason: format!("formula interaction '{display_term}' cannot use one variable twice"),
        });
    }

    let left_factor = factors.iter().find(|candidate| candidate.factor == left);
    let right_factor = factors.iter().find(|candidate| candidate.factor == right);
    let left_numeric = numeric_covariates
        .iter()
        .find(|candidate| candidate.name == left);
    let right_numeric = numeric_covariates
        .iter()
        .find(|candidate| candidate.name == right);

    match (left_factor, right_factor, left_numeric, right_numeric) {
        (Some(left_factor), Some(right_factor), None, None) => {
            if state.factor_interactions.iter().any(|interaction| {
                (interaction.left_factor == left_factor.factor
                    && interaction.right_factor == right_factor.factor)
                    || (interaction.left_factor == right_factor.factor
                        && interaction.right_factor == left_factor.factor)
            }) {
                return Err(DeseqError::InvalidOptions {
                    reason: format!("formula interaction '{display_term}' appears more than once"),
                });
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
                return Err(DeseqError::InvalidOptions {
                    reason: format!("formula interaction '{display_term}' appears more than once"),
                });
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
    if pieces.len() == 2 {
        remove_pairwise_formula_interaction(
            pieces[0],
            pieces[1],
            factors,
            numeric_covariates,
            state,
        )
    } else {
        remove_higher_order_formula_interaction(&pieces, factors, numeric_covariates, state)
    }
}

fn remove_pairwise_formula_interaction<'a>(
    left: &str,
    right: &str,
    factors: &'a [ExpandedFactorSpec<'a>],
    numeric_covariates: &'a [ExpandedNumericSpec<'a>],
    state: &mut ExpandedFormulaDesignState<'a>,
) -> Result<(), DeseqError> {
    validate_formula_variable(left)?;
    validate_formula_variable(right)?;
    if left == right {
        return Err(DeseqError::InvalidOptions {
            reason: format!("formula interaction '{left}:{right}' cannot use one variable twice"),
        });
    }
    let left_variable = resolve_formula_variable(left, factors, numeric_covariates)?;
    let right_variable = resolve_formula_variable(right, factors, numeric_covariates)?;
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
        validate_formula_variable(piece)?;
        if variables
            .iter()
            .any(|variable: &FormulaVariableRef<'_>| variable.name() == *piece)
        {
            return Err(DeseqError::InvalidOptions {
                reason: format!(
                    "formula interaction '{display_term}' cannot use one variable twice"
                ),
            });
        }
        variables.push(resolve_formula_variable(
            piece,
            factors,
            numeric_covariates,
        )?);
    }
    if state
        .higher_order_interactions
        .iter()
        .any(|interaction| same_formula_variable_set(&interaction.variables, &variables))
    {
        return Err(DeseqError::InvalidOptions {
            reason: format!("formula interaction '{display_term}' appears more than once"),
        });
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
        validate_formula_variable(piece)?;
        if variables
            .iter()
            .any(|variable: &FormulaVariableRef<'_>| variable.name() == *piece)
        {
            return Err(DeseqError::InvalidOptions {
                reason: format!(
                    "formula interaction '{}' cannot use one variable twice",
                    pieces.join(":")
                ),
            });
        }
        variables.push(resolve_formula_variable(
            piece,
            factors,
            numeric_covariates,
        )?);
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
    if let Some(factor) = factors
        .iter()
        .find(|candidate| candidate.factor == variable)
    {
        return Ok(FormulaVariableRef::Factor(factor.factor));
    }
    if let Some(covariate) = numeric_covariates
        .iter()
        .find(|candidate| candidate.name == variable)
    {
        return Ok(FormulaVariableRef::Numeric(covariate.name));
    }
    Err(DeseqError::InvalidOptions {
        reason: format!("formula variable '{variable}' is not present in supplied design metadata"),
    })
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
