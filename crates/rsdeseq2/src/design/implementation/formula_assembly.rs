fn expanded_formula_design_from_state<'a>(
    state: &ExpandedFormulaDesignState<'a>,
    factors: &'a [ExpandedFactorSpec<'a>],
    numeric_covariates: &'a [ExpandedNumericSpec<'a>],
) -> Result<ExpandedAdditiveFactorDesign, DeseqError> {
    let used_factors = formula_used_factors(state, factors)?;
    let used_numeric_covariates = formula_used_numeric_covariates(state, numeric_covariates)?;
    let n_samples = formula_design_sample_count(
        state,
        &used_factors,
        &used_numeric_covariates,
        factors,
        numeric_covariates,
    )?;
    validate_formula_used_factor_specs(&used_factors, n_samples)?;
    validate_numeric_covariate_specs(&used_numeric_covariates, n_samples)?;

    // Formula state is first normalized to the actually referenced variables;
    // coefficient names and matrix columns are then emitted in the same order
    // so expanded and reported designs stay aligned.
    let factor_levels = used_factors
        .iter()
        .map(|factor| {
            ordered_levels_with_declared(factor.sample_levels, factor.reference, factor.levels)
        })
        .collect::<Vec<_>>();

    let mut expanded_names = Vec::new();
    if state.has_intercept {
        expanded_names.push("Intercept".to_string());
    }
    for factor in &state.selected_factors {
        let levels = formula_factor_levels(factor.factor, &used_factors, &factor_levels)?;
        expanded_names.extend(
            levels
                .iter()
                .map(|level| expanded_factor_coefficient_name(factor.factor, level)),
        );
    }
    expanded_names.extend(
        state
            .selected_numeric_covariates
            .iter()
            .map(|covariate| covariate.name.to_string()),
    );
    for interaction in &state.factor_interactions {
        let left_levels =
            formula_factor_levels(interaction.left_factor, &used_factors, &factor_levels)?;
        let right_levels =
            formula_factor_levels(interaction.right_factor, &used_factors, &factor_levels)?;
        for left_level in left_levels {
            for right_level in right_levels {
                expanded_names.push(interaction_coefficient_name(
                    interaction.left_factor,
                    left_level,
                    interaction.right_factor,
                    right_level,
                ));
            }
        }
    }
    for interaction in &state.factor_numeric_interactions {
        let levels = formula_factor_levels(interaction.factor, &used_factors, &factor_levels)?;
        for level in levels {
            expanded_names.push(factor_numeric_interaction_coefficient_name(
                interaction.factor,
                level,
                interaction.numeric,
            ));
        }
    }
    for interaction in &state.numeric_interactions {
        expanded_names.push(numeric_interaction_coefficient_name(
            interaction.left_numeric,
            interaction.right_numeric,
        ));
    }
    for interaction in &state.higher_order_interactions {
        expanded_names.extend(formula_higher_order_expanded_names(
            interaction,
            &used_factors,
            &factor_levels,
        )?);
    }
    validate_unique_coefficient_names(&expanded_names, "expanded formula design")?;

    let mut expanded_values = Vec::with_capacity(n_samples * expanded_names.len());
    for sample in 0..n_samples {
        if state.has_intercept {
            expanded_values.push(1.0);
        }
        for factor in &state.selected_factors {
            let levels = formula_factor_levels(factor.factor, &used_factors, &factor_levels)?;
            let sample_level = factor.sample_levels[sample].as_str();
            for level in levels {
                expanded_values.push((sample_level == level) as u8 as f64);
            }
        }
        for covariate in &state.selected_numeric_covariates {
            expanded_values.push(covariate.values[sample]);
        }
        for interaction in &state.factor_interactions {
            let left = formula_factor_spec(interaction.left_factor, &used_factors)?;
            let right = formula_factor_spec(interaction.right_factor, &used_factors)?;
            let left_levels =
                formula_factor_levels(interaction.left_factor, &used_factors, &factor_levels)?;
            let right_levels =
                formula_factor_levels(interaction.right_factor, &used_factors, &factor_levels)?;
            let left_sample_level = left.sample_levels[sample].as_str();
            let right_sample_level = right.sample_levels[sample].as_str();
            for left_level in left_levels {
                for right_level in right_levels {
                    expanded_values.push(
                        (left_sample_level == left_level && right_sample_level == right_level) as u8
                            as f64,
                    );
                }
            }
        }
        for interaction in &state.factor_numeric_interactions {
            let factor = formula_factor_spec(interaction.factor, &used_factors)?;
            let numeric = formula_numeric_spec(interaction.numeric, &used_numeric_covariates)?;
            let levels = formula_factor_levels(interaction.factor, &used_factors, &factor_levels)?;
            let sample_level = factor.sample_levels[sample].as_str();
            let numeric_value = numeric.values[sample];
            for level in levels {
                expanded_values.push(((sample_level == level) as u8 as f64) * numeric_value);
            }
        }
        for interaction in &state.numeric_interactions {
            expanded_values.push(
                formula_numeric_spec(interaction.left_numeric, &used_numeric_covariates)?.values
                    [sample]
                    * formula_numeric_spec(interaction.right_numeric, &used_numeric_covariates)?
                        .values[sample],
            );
        }
        for interaction in &state.higher_order_interactions {
            expanded_values.extend(formula_higher_order_expanded_values(
                interaction,
                &used_factors,
                &used_numeric_covariates,
                &factor_levels,
                sample,
            )?);
        }
    }
    let expanded_design = DesignMatrix::from_row_major(
        n_samples,
        expanded_names.len(),
        expanded_values,
        Some(expanded_names.clone()),
    )?;

    let mut standard_names = Vec::new();
    let mut standard_expanded_names = Vec::new();
    if state.has_intercept {
        standard_names.push("Intercept".to_string());
        standard_expanded_names.push("Intercept".to_string());
    }
    for (factor_idx, factor) in state.selected_factors.iter().enumerate() {
        let levels = formula_factor_levels(factor.factor, &used_factors, &factor_levels)?;
        let use_all_levels = !state.has_intercept && factor_idx == 0;
        for level in levels
            .iter()
            .filter(|level| use_all_levels || level.as_str() != factor.reference)
        {
            if use_all_levels {
                standard_names.push(expanded_factor_coefficient_name(factor.factor, level));
            } else {
                standard_names.push(standard_factor_coefficient_name(
                    factor.factor,
                    level,
                    factor.reference,
                ));
            }
            standard_expanded_names.push(expanded_factor_coefficient_name(factor.factor, level));
        }
    }
    for covariate in &state.selected_numeric_covariates {
        standard_names.push(covariate.name.to_string());
        standard_expanded_names.push(covariate.name.to_string());
    }
    for interaction in &state.factor_interactions {
        push_formula_factor_interaction_standard_names(
            interaction,
            &used_factors,
            &factor_levels,
            state,
            &mut standard_names,
            &mut standard_expanded_names,
        )?;
    }
    for interaction in &state.factor_numeric_interactions {
        push_formula_factor_numeric_standard_names(
            interaction,
            &used_factors,
            &factor_levels,
            state,
            &mut standard_names,
            &mut standard_expanded_names,
        )?;
    }
    for interaction in &state.numeric_interactions {
        let name = numeric_interaction_coefficient_name(
            interaction.left_numeric,
            interaction.right_numeric,
        );
        standard_names.push(name.clone());
        standard_expanded_names.push(name);
    }
    for interaction in &state.higher_order_interactions {
        push_formula_higher_order_standard_names(
            interaction,
            &used_factors,
            &factor_levels,
            state,
            &mut standard_names,
            &mut standard_expanded_names,
        )?;
    }
    validate_unique_coefficient_names(&standard_names, "standard formula design")?;

    let mut standard_values = Vec::with_capacity(n_samples * standard_expanded_names.len());
    for sample in 0..n_samples {
        for expanded_name in &standard_expanded_names {
            let col = expanded_names
                .iter()
                .position(|candidate| candidate == expanded_name)
                .ok_or_else(|| DeseqError::InvalidOptions {
                    reason: format!(
                        "formula standard column '{expanded_name}' is not present in expanded design"
                    ),
                })?;
            standard_values.push(
                expanded_design
                    .matrix()
                    .get(sample, col)
                    .copied()
                    .unwrap_or(0.0),
            );
        }
    }
    let standard_design = DesignMatrix::from_row_major(
        n_samples,
        standard_names.len(),
        standard_values,
        Some(standard_names),
    )?;

    let mut coefficient_groups = Vec::with_capacity(standard_design.n_coefficients());
    for expanded_name in &standard_expanded_names {
        let col = expanded_names
            .iter()
            .position(|candidate| candidate == expanded_name)
            .ok_or_else(|| DeseqError::InvalidOptions {
                reason: format!(
                    "formula standard column '{expanded_name}' is not present in expanded design"
                ),
            })?;
        coefficient_groups.push(vec![col]);
    }

    Ok(ExpandedAdditiveFactorDesign {
        expanded_design,
        standard_design,
        coefficient_groups,
        factor_levels: state
            .selected_factors
            .iter()
            .map(|factor| {
                ordered_levels_with_declared(factor.sample_levels, factor.reference, factor.levels)
            })
            .collect(),
        numeric_covariates: state
            .selected_numeric_covariates
            .iter()
            .map(|covariate| covariate.name.to_string())
            .collect(),
        interactions: state
            .factor_interactions
            .iter()
            .map(|interaction| format!("{}:{}", interaction.left_factor, interaction.right_factor))
            .collect(),
        factor_numeric_interactions: state
            .factor_numeric_interactions
            .iter()
            .map(|interaction| format!("{}:{}", interaction.factor, interaction.numeric))
            .collect(),
        numeric_interactions: state
            .numeric_interactions
            .iter()
            .map(|interaction| {
                format!("{}:{}", interaction.left_numeric, interaction.right_numeric)
            })
            .collect(),
        higher_order_interactions: state
            .higher_order_interactions
            .iter()
            .map(formula_higher_order_interaction_name)
            .collect(),
    })
}

fn formula_used_factors<'a>(
    state: &ExpandedFormulaDesignState<'a>,
    factors: &'a [ExpandedFactorSpec<'a>],
) -> Result<Vec<ExpandedFactorSpec<'a>>, DeseqError> {
    let mut used = Vec::new();
    for factor in &state.selected_factors {
        push_unique_factor_spec(&mut used, factor.clone())?;
    }
    for interaction in &state.factor_interactions {
        push_unique_factor_spec(
            &mut used,
            formula_factor_spec(interaction.left_factor, factors)?,
        )?;
        push_unique_factor_spec(
            &mut used,
            formula_factor_spec(interaction.right_factor, factors)?,
        )?;
    }
    for interaction in &state.factor_numeric_interactions {
        push_unique_factor_spec(&mut used, formula_factor_spec(interaction.factor, factors)?)?;
    }
    for interaction in &state.higher_order_interactions {
        for variable in &interaction.variables {
            if let FormulaVariableRef::Factor(factor) = variable {
                push_unique_factor_spec(&mut used, formula_factor_spec(factor, factors)?)?;
            }
        }
    }
    Ok(used)
}

fn formula_used_numeric_covariates<'a>(
    state: &ExpandedFormulaDesignState<'a>,
    numeric_covariates: &'a [ExpandedNumericSpec<'a>],
) -> Result<Vec<ExpandedNumericSpec<'a>>, DeseqError> {
    let mut used = Vec::new();
    for covariate in &state.selected_numeric_covariates {
        push_unique_numeric_spec(&mut used, covariate.clone())?;
    }
    for interaction in &state.factor_numeric_interactions {
        push_unique_numeric_spec(
            &mut used,
            formula_numeric_spec(interaction.numeric, numeric_covariates)?,
        )?;
    }
    for interaction in &state.numeric_interactions {
        push_unique_numeric_spec(
            &mut used,
            formula_numeric_spec(interaction.left_numeric, numeric_covariates)?,
        )?;
        push_unique_numeric_spec(
            &mut used,
            formula_numeric_spec(interaction.right_numeric, numeric_covariates)?,
        )?;
    }
    for interaction in &state.higher_order_interactions {
        for variable in &interaction.variables {
            if let FormulaVariableRef::Numeric(numeric) = variable {
                push_unique_numeric_spec(
                    &mut used,
                    formula_numeric_spec(numeric, numeric_covariates)?,
                )?;
            }
        }
    }
    Ok(used)
}

fn push_unique_factor_spec<'a>(
    used: &mut Vec<ExpandedFactorSpec<'a>>,
    factor: ExpandedFactorSpec<'a>,
) -> Result<(), DeseqError> {
    if used
        .iter()
        .any(|candidate| candidate.factor == factor.factor)
    {
        return Ok(());
    }
    used.push(factor);
    Ok(())
}

fn push_unique_numeric_spec<'a>(
    used: &mut Vec<ExpandedNumericSpec<'a>>,
    covariate: ExpandedNumericSpec<'a>,
) -> Result<(), DeseqError> {
    if used
        .iter()
        .any(|candidate| candidate.name == covariate.name)
    {
        return Ok(());
    }
    used.push(covariate);
    Ok(())
}

fn formula_design_sample_count(
    state: &ExpandedFormulaDesignState<'_>,
    used_factors: &[ExpandedFactorSpec<'_>],
    used_numeric_covariates: &[ExpandedNumericSpec<'_>],
    supplied_factors: &[ExpandedFactorSpec<'_>],
    supplied_numeric_covariates: &[ExpandedNumericSpec<'_>],
) -> Result<usize, DeseqError> {
    if let Some(factor) = used_factors.first() {
        return Ok(factor.sample_levels.len());
    }
    if let Some(covariate) = used_numeric_covariates.first() {
        if covariate.values.is_empty() {
            return Err(DeseqError::InvalidOptions {
                reason: "numeric covariate values must be non-empty".to_string(),
            });
        }
        return Ok(covariate.values.len());
    }
    if !state.has_intercept {
        return Err(DeseqError::InvalidOptions {
            reason: "formula without an intercept must include at least one design term"
                .to_string(),
        });
    }
    if let Some(factor) = supplied_factors.first() {
        if factor.sample_levels.is_empty() {
            return Err(DeseqError::InvalidOptions {
                reason: "intercept-only formula requires non-empty sample metadata".to_string(),
            });
        }
        return Ok(factor.sample_levels.len());
    }
    if let Some(covariate) = supplied_numeric_covariates.first() {
        if covariate.values.is_empty() {
            return Err(DeseqError::InvalidOptions {
                reason: "intercept-only formula requires non-empty sample metadata".to_string(),
            });
        }
        return Ok(covariate.values.len());
    }
    Err(DeseqError::InvalidOptions {
        reason: "intercept-only formula requires supplied sample metadata".to_string(),
    })
}

fn validate_formula_used_factor_specs(
    factors: &[ExpandedFactorSpec<'_>],
    n_samples: usize,
) -> Result<(), DeseqError> {
    for (idx, factor) in factors.iter().enumerate() {
        validate_factor_design_inputs_with_levels(
            factor.factor,
            factor.sample_levels,
            factor.reference,
            factor.levels,
        )?;
        if factor.sample_levels.len() != n_samples {
            return Err(invalid_dimensions(
                "formula factor sample levels",
                n_samples,
                factor.sample_levels.len(),
            ));
        }
        if factors[..idx]
            .iter()
            .any(|previous| previous.factor == factor.factor)
        {
            return Err(DeseqError::InvalidOptions {
                reason: format!("factor '{}' appears more than once", factor.factor),
            });
        }
    }
    Ok(())
}

fn formula_factor_spec<'a>(
    factor: &str,
    factors: &'a [ExpandedFactorSpec<'a>],
) -> Result<ExpandedFactorSpec<'a>, DeseqError> {
    factors
        .iter()
        .find(|candidate| candidate.factor == factor)
        .cloned()
        .ok_or_else(|| DeseqError::InvalidOptions {
            reason: format!("formula factor '{factor}' is not present in supplied design metadata"),
        })
}

fn formula_numeric_spec<'a>(
    numeric: &str,
    numeric_covariates: &'a [ExpandedNumericSpec<'a>],
) -> Result<ExpandedNumericSpec<'a>, DeseqError> {
    numeric_covariates
        .iter()
        .find(|candidate| candidate.name == numeric)
        .cloned()
        .ok_or_else(|| DeseqError::InvalidOptions {
            reason: format!(
                "formula numeric covariate '{numeric}' is not present in supplied design metadata"
            ),
        })
}

fn formula_factor_levels<'a>(
    factor: &str,
    factors: &[ExpandedFactorSpec<'_>],
    factor_levels: &'a [Vec<String>],
) -> Result<&'a [String], DeseqError> {
    let idx = factors
        .iter()
        .position(|candidate| candidate.factor == factor)
        .ok_or_else(|| DeseqError::InvalidOptions {
            reason: format!("formula factor '{factor}' is not present in supplied design metadata"),
        })?;
    Ok(&factor_levels[idx])
}

fn formula_has_factor_main_effect(state: &ExpandedFormulaDesignState<'_>, factor: &str) -> bool {
    state
        .selected_factors
        .iter()
        .any(|candidate| candidate.factor == factor)
}

fn formula_has_numeric_main_effect(state: &ExpandedFormulaDesignState<'_>, numeric: &str) -> bool {
    state
        .selected_numeric_covariates
        .iter()
        .any(|candidate| candidate.name == numeric)
}

fn push_formula_factor_interaction_standard_names(
    interaction: &ExpandedFactorInteractionSpec<'_>,
    used_factors: &[ExpandedFactorSpec<'_>],
    factor_levels: &[Vec<String>],
    state: &ExpandedFormulaDesignState<'_>,
    standard_names: &mut Vec<String>,
    standard_expanded_names: &mut Vec<String>,
) -> Result<(), DeseqError> {
    let left = formula_factor_spec(interaction.left_factor, used_factors)?;
    let right = formula_factor_spec(interaction.right_factor, used_factors)?;
    let left_levels = formula_factor_levels(interaction.left_factor, used_factors, factor_levels)?;
    let right_levels =
        formula_factor_levels(interaction.right_factor, used_factors, factor_levels)?;
    let left_main = formula_has_factor_main_effect(state, interaction.left_factor);
    let right_main = formula_has_factor_main_effect(state, interaction.right_factor);
    let left_treatment = right_main;
    let right_treatment = left_main;
    for left_level in left_levels
        .iter()
        .filter(|level| !left_treatment || level.as_str() != left.reference)
    {
        for right_level in right_levels
            .iter()
            .filter(|level| !right_treatment || level.as_str() != right.reference)
        {
            standard_names.push(formula_factor_interaction_standard_name(
                FormulaFactorInteractionPiece {
                    factor: interaction.left_factor,
                    level: left_level,
                    reference: left.reference,
                    use_treatment_name: left_treatment,
                },
                FormulaFactorInteractionPiece {
                    factor: interaction.right_factor,
                    level: right_level,
                    reference: right.reference,
                    use_treatment_name: right_treatment,
                },
            ));
            standard_expanded_names.push(interaction_coefficient_name(
                interaction.left_factor,
                left_level,
                interaction.right_factor,
                right_level,
            ));
        }
    }
    Ok(())
}

fn push_formula_factor_numeric_standard_names(
    interaction: &ExpandedFactorNumericInteractionSpec<'_>,
    used_factors: &[ExpandedFactorSpec<'_>],
    factor_levels: &[Vec<String>],
    state: &ExpandedFormulaDesignState<'_>,
    standard_names: &mut Vec<String>,
    standard_expanded_names: &mut Vec<String>,
) -> Result<(), DeseqError> {
    let factor = formula_factor_spec(interaction.factor, used_factors)?;
    let levels = formula_factor_levels(interaction.factor, used_factors, factor_levels)?;
    let numeric_main = formula_has_numeric_main_effect(state, interaction.numeric);
    let use_treatment_levels = numeric_main;
    for level in levels
        .iter()
        .filter(|level| !use_treatment_levels || level.as_str() != factor.reference)
    {
        if use_treatment_levels {
            standard_names.push(standard_factor_numeric_interaction_coefficient_name(
                interaction.factor,
                level,
                factor.reference,
                interaction.numeric,
            ));
        } else {
            standard_names.push(factor_numeric_interaction_coefficient_name(
                interaction.factor,
                level,
                interaction.numeric,
            ));
        }
        standard_expanded_names.push(factor_numeric_interaction_coefficient_name(
            interaction.factor,
            level,
            interaction.numeric,
        ));
    }
    Ok(())
}

fn formula_higher_order_expanded_names(
    interaction: &FormulaHigherOrderInteractionSpec<'_>,
    used_factors: &[ExpandedFactorSpec<'_>],
    factor_levels: &[Vec<String>],
) -> Result<Vec<String>, DeseqError> {
    formula_higher_order_names(interaction, used_factors, factor_levels, &|_| false)
}

fn formula_higher_order_expanded_values(
    interaction: &FormulaHigherOrderInteractionSpec<'_>,
    used_factors: &[ExpandedFactorSpec<'_>],
    used_numeric_covariates: &[ExpandedNumericSpec<'_>],
    factor_levels: &[Vec<String>],
    sample: usize,
) -> Result<Vec<f64>, DeseqError> {
    let mut values = Vec::new();
    let level_grid = formula_higher_order_level_grid(interaction, used_factors, factor_levels)?;
    for levels in level_grid {
        let mut value = 1.0;
        for (variable, level) in interaction.variables.iter().zip(levels) {
            match variable {
                FormulaVariableRef::Factor(factor) => {
                    let factor = formula_factor_spec(factor, used_factors)?;
                    let sample_level = factor.sample_levels[sample].as_str();
                    let Some(level) = level else {
                        return Err(DeseqError::InvalidOptions {
                            reason: "formula factor interaction level is missing".to_string(),
                        });
                    };
                    value *= (sample_level == level) as u8 as f64;
                }
                FormulaVariableRef::Numeric(numeric) => {
                    value *= formula_numeric_spec(numeric, used_numeric_covariates)?.values[sample];
                }
            }
        }
        values.push(value);
    }
    Ok(values)
}

fn push_formula_higher_order_standard_names(
    interaction: &FormulaHigherOrderInteractionSpec<'_>,
    used_factors: &[ExpandedFactorSpec<'_>],
    factor_levels: &[Vec<String>],
    state: &ExpandedFormulaDesignState<'_>,
    standard_names: &mut Vec<String>,
    standard_expanded_names: &mut Vec<String>,
) -> Result<(), DeseqError> {
    let treatment_filter = |factor: &str| formula_has_factor_main_effect(state, factor);
    let standard =
        formula_higher_order_names(interaction, used_factors, factor_levels, &treatment_filter)?;
    let expanded_by_name = formula_higher_order_names_with_choice_filter(
        interaction,
        used_factors,
        factor_levels,
        &|_| false,
        &treatment_filter,
    )?;

    for (standard_name, expanded_name) in standard.into_iter().zip(expanded_by_name) {
        standard_names.push(standard_name);
        standard_expanded_names.push(expanded_name);
    }
    Ok(())
}

fn formula_higher_order_names<F>(
    interaction: &FormulaHigherOrderInteractionSpec<'_>,
    used_factors: &[ExpandedFactorSpec<'_>],
    factor_levels: &[Vec<String>],
    use_treatment_for_factor: &F,
) -> Result<Vec<String>, DeseqError>
where
    F: Fn(&str) -> bool,
{
    formula_higher_order_names_with_choice_filter(
        interaction,
        used_factors,
        factor_levels,
        use_treatment_for_factor,
        use_treatment_for_factor,
    )
}

fn formula_higher_order_names_with_choice_filter<F, G>(
    interaction: &FormulaHigherOrderInteractionSpec<'_>,
    used_factors: &[ExpandedFactorSpec<'_>],
    factor_levels: &[Vec<String>],
    use_treatment_for_factor: &F,
    use_treatment_choices_for_factor: &G,
) -> Result<Vec<String>, DeseqError>
where
    F: Fn(&str) -> bool,
    G: Fn(&str) -> bool,
{
    let mut names = Vec::new();
    let level_grid = formula_higher_order_level_grid_with_treatment(
        interaction,
        used_factors,
        factor_levels,
        use_treatment_choices_for_factor,
    )?;
    for levels in level_grid {
        let mut pieces = Vec::with_capacity(interaction.variables.len());
        for (variable, level) in interaction.variables.iter().zip(levels) {
            match variable {
                FormulaVariableRef::Factor(factor) => {
                    let factor_spec = formula_factor_spec(factor, used_factors)?;
                    let Some(level) = level else {
                        return Err(DeseqError::InvalidOptions {
                            reason: "formula factor interaction level is missing".to_string(),
                        });
                    };
                    if use_treatment_for_factor(factor) {
                        pieces.push(standard_factor_coefficient_name(
                            factor,
                            level,
                            factor_spec.reference,
                        ));
                    } else {
                        pieces.push(expanded_factor_coefficient_name(factor, level));
                    }
                }
                FormulaVariableRef::Numeric(numeric) => pieces.push((*numeric).to_string()),
            }
        }
        names.push(pieces.join(":"));
    }
    Ok(names)
}

fn formula_higher_order_level_grid<'a>(
    interaction: &FormulaHigherOrderInteractionSpec<'_>,
    used_factors: &[ExpandedFactorSpec<'_>],
    factor_levels: &'a [Vec<String>],
) -> Result<Vec<Vec<Option<&'a str>>>, DeseqError> {
    formula_higher_order_level_grid_with_treatment(
        interaction,
        used_factors,
        factor_levels,
        &|_| false,
    )
}

fn formula_higher_order_level_grid_with_treatment<'a, F>(
    interaction: &FormulaHigherOrderInteractionSpec<'_>,
    used_factors: &[ExpandedFactorSpec<'_>],
    factor_levels: &'a [Vec<String>],
    use_treatment_for_factor: &F,
) -> Result<Vec<Vec<Option<&'a str>>>, DeseqError>
where
    F: Fn(&str) -> bool,
{
    let mut grids: Vec<Vec<Option<&str>>> = vec![Vec::new()];
    for variable in &interaction.variables {
        let choices = match variable {
            FormulaVariableRef::Factor(factor) => {
                let factor_spec = formula_factor_spec(factor, used_factors)?;
                formula_factor_levels(factor, used_factors, factor_levels)?
                    .iter()
                    .filter(|level| {
                        !use_treatment_for_factor(factor) || level.as_str() != factor_spec.reference
                    })
                    .map(|level| Some(level.as_str()))
                    .collect::<Vec<_>>()
            }
            FormulaVariableRef::Numeric(_) => vec![None],
        };
        let mut next = Vec::with_capacity(grids.len() * choices.len());
        for prefix in &grids {
            for choice in &choices {
                let mut row = prefix.clone();
                row.push(*choice);
                next.push(row);
            }
        }
        grids = next;
    }
    Ok(grids)
}

fn formula_higher_order_interaction_name(
    interaction: &FormulaHigherOrderInteractionSpec<'_>,
) -> String {
    interaction
        .variables
        .iter()
        .map(FormulaVariableRef::name)
        .collect::<Vec<_>>()
        .join(":")
}

struct FormulaFactorInteractionPiece<'a> {
    factor: &'a str,
    level: &'a str,
    reference: &'a str,
    use_treatment_name: bool,
}

fn formula_factor_interaction_standard_name(
    left: FormulaFactorInteractionPiece<'_>,
    right: FormulaFactorInteractionPiece<'_>,
) -> String {
    let left = if left.use_treatment_name {
        standard_factor_coefficient_name(left.factor, left.level, left.reference)
    } else {
        expanded_factor_coefficient_name(left.factor, left.level)
    };
    let right = if right.use_treatment_name {
        standard_factor_coefficient_name(right.factor, right.level, right.reference)
    } else {
        expanded_factor_coefficient_name(right.factor, right.level)
    };
    format!("{left}:{right}")
}

fn add_factor_numeric_formula_interaction<'a>(
    factor: &'a str,
    numeric: &'a str,
    factor_numeric_interactions: &mut Vec<ExpandedFactorNumericInteractionSpec<'a>>,
    _term: &str,
) -> Result<(), DeseqError> {
    if factor_numeric_interactions
        .iter()
        .any(|interaction| interaction.factor == factor && interaction.numeric == numeric)
    {
        return Ok(());
    }
    factor_numeric_interactions.push(ExpandedFactorNumericInteractionSpec { factor, numeric });
    Ok(())
}

fn split_formula_pieces(term: &str, delimiter: char) -> Result<Vec<String>, DeseqError> {
    let pieces = split_formula_top_level(term, delimiter)?;
    if pieces.len() < 2 || pieces.iter().any(|piece| piece.is_empty()) {
        return Err(DeseqError::InvalidOptions {
            reason: format!("formula term '{term}' must contain non-empty variables"),
        });
    }
    Ok(pieces)
}
