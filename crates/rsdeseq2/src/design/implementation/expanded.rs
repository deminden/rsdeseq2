/// Build a primitive expanded one-factor design and matching standard design.
///
/// The expanded matrix has columns `Intercept`, then one indicator per factor
/// level. The standard matrix has `Intercept`, then one treatment-style
/// indicator for every non-reference level. This helper owns only the matrix
/// construction step; callers still choose how to use the returned
/// `coefficient_groups` in beta-prior fitting and result assembly.
pub fn expanded_factor_design<S: AsRef<str>>(
    factor: &str,
    sample_levels: &[S],
    reference: &str,
) -> Result<ExpandedFactorDesign, DeseqError> {
    validate_factor_design_inputs(factor, sample_levels, reference)?;
    let levels = ordered_levels(sample_levels, reference);
    let n_samples = sample_levels.len();

    let mut expanded_names = Vec::with_capacity(levels.len() + 1);
    expanded_names.push("Intercept".to_string());
    expanded_names.extend(
        levels
            .iter()
            .map(|level| expanded_factor_coefficient_name(factor, level)),
    );
    let mut expanded_values = Vec::with_capacity(n_samples * expanded_names.len());
    for level in sample_levels {
        let level = level.as_ref();
        expanded_values.push(1.0);
        for candidate in &levels {
            expanded_values.push((level == candidate) as u8 as f64);
        }
    }
    let expanded_design = DesignMatrix::from_row_major(
        n_samples,
        expanded_names.len(),
        expanded_values,
        Some(expanded_names),
    )?;

    let non_reference_levels = levels
        .iter()
        .filter(|level| level.as_str() != reference)
        .collect::<Vec<_>>();
    let mut standard_names = Vec::with_capacity(non_reference_levels.len() + 1);
    standard_names.push("Intercept".to_string());
    standard_names.extend(
        non_reference_levels
            .iter()
            .map(|level| standard_factor_coefficient_name(factor, level, reference)),
    );
    let mut standard_values = Vec::with_capacity(n_samples * standard_names.len());
    for level in sample_levels {
        let level = level.as_ref();
        standard_values.push(1.0);
        for candidate in &non_reference_levels {
            standard_values.push((level == candidate.as_str()) as u8 as f64);
        }
    }
    let standard_design = DesignMatrix::from_row_major(
        n_samples,
        standard_names.len(),
        standard_values,
        Some(standard_names),
    )?;

    let mut coefficient_groups = Vec::with_capacity(standard_design.n_coefficients());
    coefficient_groups.push(vec![0]);
    for level in non_reference_levels {
        let expanded_column = levels
            .iter()
            .position(|candidate| candidate == level)
            .map(|idx| idx + 1)
            .ok_or_else(|| DeseqError::InvalidOptions {
                reason: format!("factor level '{}' is not present in expanded levels", level),
            })?;
        coefficient_groups.push(vec![expanded_column]);
    }

    Ok(ExpandedFactorDesign {
        expanded_design,
        standard_design,
        coefficient_groups,
        levels,
    })
}

/// Build an expanded additive-factor design and matching standard design.
///
/// This covers primitive `~ factor1 + factor2 + ...` model-matrix
/// construction for categorical terms. The expanded matrix has one intercept
/// plus one indicator per level of each factor. The standard matrix has one
/// intercept plus one treatment-style indicator per non-reference level of
/// each factor. Interactions, nested terms, continuous covariates, and formula
/// parsing remain caller or wrapper responsibilities.
pub fn expanded_additive_factor_design(
    factors: &[ExpandedFactorSpec<'_>],
) -> Result<ExpandedAdditiveFactorDesign, DeseqError> {
    expanded_additive_design(factors, &[])
}

/// Build an expanded additive design with categorical factors and numeric covariates.
///
/// Categorical terms use expanded indicator columns and treatment-style
/// reported columns as in [`expanded_additive_factor_design`]. Numeric
/// covariates are included unchanged in both expanded and standard designs and
/// are mapped one-to-one in the coefficient groups. Interactions, splines,
/// nested terms, transformed variables, and formula parsing remain caller or
/// wrapper responsibilities.
pub fn expanded_additive_design(
    factors: &[ExpandedFactorSpec<'_>],
    numeric_covariates: &[ExpandedNumericSpec<'_>],
) -> Result<ExpandedAdditiveFactorDesign, DeseqError> {
    expanded_additive_design_with_interactions(factors, numeric_covariates, &[])
}

/// Build an expanded additive design with categorical factors, numerics, and factor interactions.
///
/// Factor-by-factor interactions add all level-pair products to the expanded
/// design and non-reference-by-non-reference treatment-style products to the
/// standard design. The interaction terms reference factor names already
/// supplied in `factors`.
pub fn expanded_additive_design_with_interactions(
    factors: &[ExpandedFactorSpec<'_>],
    numeric_covariates: &[ExpandedNumericSpec<'_>],
    interactions: &[ExpandedFactorInteractionSpec<'_>],
) -> Result<ExpandedAdditiveFactorDesign, DeseqError> {
    expanded_additive_design_with_all_interactions(
        factors,
        numeric_covariates,
        interactions,
        &[],
        &[],
    )
}

/// Build an expanded additive design with all primitive pairwise interaction types.
pub fn expanded_additive_design_with_all_interactions(
    factors: &[ExpandedFactorSpec<'_>],
    numeric_covariates: &[ExpandedNumericSpec<'_>],
    factor_interactions: &[ExpandedFactorInteractionSpec<'_>],
    factor_numeric_interactions: &[ExpandedFactorNumericInteractionSpec<'_>],
    numeric_interactions: &[ExpandedNumericInteractionSpec<'_>],
) -> Result<ExpandedAdditiveFactorDesign, DeseqError> {
    let n_samples = additive_design_sample_count(factors, numeric_covariates)?;
    if !factors.is_empty() {
        validate_additive_factor_specs(factors)?;
    }
    validate_numeric_covariate_specs(numeric_covariates, n_samples)?;
    let interaction_indices = resolve_interaction_indices(factors, factor_interactions)?;
    let factor_numeric_indices = resolve_factor_numeric_interaction_indices(
        factors,
        numeric_covariates,
        factor_numeric_interactions,
    )?;
    let numeric_interaction_indices =
        resolve_numeric_interaction_indices(numeric_covariates, numeric_interactions)?;
    let factor_levels = factors
        .iter()
        .map(|factor| {
            ordered_levels_with_declared(factor.sample_levels, factor.reference, factor.levels)
        })
        .collect::<Vec<_>>();

    let mut expanded_names = vec!["Intercept".to_string()];
    for (factor, levels) in factors.iter().zip(&factor_levels) {
        expanded_names.extend(
            levels
                .iter()
                .map(|level| expanded_factor_coefficient_name(factor.factor, level)),
        );
    }
    expanded_names.extend(
        numeric_covariates
            .iter()
            .map(|covariate| covariate.name.to_string()),
    );
    for (left_idx, right_idx) in &interaction_indices {
        let left = &factors[*left_idx];
        let right = &factors[*right_idx];
        for left_level in &factor_levels[*left_idx] {
            for right_level in &factor_levels[*right_idx] {
                expanded_names.push(interaction_coefficient_name(
                    left.factor,
                    left_level,
                    right.factor,
                    right_level,
                ));
            }
        }
    }
    for (factor_idx, numeric_idx) in &factor_numeric_indices {
        let factor = &factors[*factor_idx];
        let numeric = &numeric_covariates[*numeric_idx];
        for level in &factor_levels[*factor_idx] {
            expanded_names.push(factor_numeric_interaction_coefficient_name(
                factor.factor,
                level,
                numeric.name,
            ));
        }
    }
    for (left_idx, right_idx) in &numeric_interaction_indices {
        expanded_names.push(numeric_interaction_coefficient_name(
            numeric_covariates[*left_idx].name,
            numeric_covariates[*right_idx].name,
        ));
    }
    validate_unique_coefficient_names(&expanded_names, "expanded additive design")?;
    let mut expanded_values = Vec::with_capacity(n_samples * expanded_names.len());
    for sample in 0..n_samples {
        expanded_values.push(1.0);
        for (factor, levels) in factors.iter().zip(&factor_levels) {
            let sample_level = factor.sample_levels[sample].as_str();
            for level in levels {
                expanded_values.push((sample_level == level) as u8 as f64);
            }
        }
        for covariate in numeric_covariates {
            expanded_values.push(covariate.values[sample]);
        }
        for (left_idx, right_idx) in &interaction_indices {
            let left = &factors[*left_idx];
            let right = &factors[*right_idx];
            let left_sample_level = left.sample_levels[sample].as_str();
            let right_sample_level = right.sample_levels[sample].as_str();
            for left_level in &factor_levels[*left_idx] {
                for right_level in &factor_levels[*right_idx] {
                    expanded_values.push(
                        (left_sample_level == left_level && right_sample_level == right_level) as u8
                            as f64,
                    );
                }
            }
        }
        for (factor_idx, numeric_idx) in &factor_numeric_indices {
            let factor = &factors[*factor_idx];
            let numeric = &numeric_covariates[*numeric_idx];
            let sample_level = factor.sample_levels[sample].as_str();
            let numeric_value = numeric.values[sample];
            for level in &factor_levels[*factor_idx] {
                expanded_values.push(((sample_level == level) as u8 as f64) * numeric_value);
            }
        }
        for (left_idx, right_idx) in &numeric_interaction_indices {
            expanded_values.push(
                numeric_covariates[*left_idx].values[sample]
                    * numeric_covariates[*right_idx].values[sample],
            );
        }
    }
    let expanded_design = DesignMatrix::from_row_major(
        n_samples,
        expanded_names.len(),
        expanded_values,
        Some(expanded_names),
    )?;

    let mut standard_names = vec!["Intercept".to_string()];
    for (factor, levels) in factors.iter().zip(&factor_levels) {
        standard_names.extend(
            levels
                .iter()
                .filter(|level| level.as_str() != factor.reference)
                .map(|level| {
                    standard_factor_coefficient_name(factor.factor, level, factor.reference)
                }),
        );
    }
    standard_names.extend(
        numeric_covariates
            .iter()
            .map(|covariate| covariate.name.to_string()),
    );
    for (left_idx, right_idx) in &interaction_indices {
        let left = &factors[*left_idx];
        let right = &factors[*right_idx];
        for left_level in factor_levels[*left_idx]
            .iter()
            .filter(|level| level.as_str() != left.reference)
        {
            for right_level in factor_levels[*right_idx]
                .iter()
                .filter(|level| level.as_str() != right.reference)
            {
                standard_names.push(standard_interaction_coefficient_name(
                    left.factor,
                    left_level,
                    left.reference,
                    right.factor,
                    right_level,
                    right.reference,
                ));
            }
        }
    }
    for (factor_idx, numeric_idx) in &factor_numeric_indices {
        let factor = &factors[*factor_idx];
        let numeric = &numeric_covariates[*numeric_idx];
        for level in factor_levels[*factor_idx]
            .iter()
            .filter(|level| level.as_str() != factor.reference)
        {
            standard_names.push(standard_factor_numeric_interaction_coefficient_name(
                factor.factor,
                level,
                factor.reference,
                numeric.name,
            ));
        }
    }
    for (left_idx, right_idx) in &numeric_interaction_indices {
        standard_names.push(numeric_interaction_coefficient_name(
            numeric_covariates[*left_idx].name,
            numeric_covariates[*right_idx].name,
        ));
    }
    validate_unique_coefficient_names(&standard_names, "standard additive design")?;
    let mut standard_values = Vec::with_capacity(n_samples * standard_names.len());
    for sample in 0..n_samples {
        standard_values.push(1.0);
        for (factor, levels) in factors.iter().zip(&factor_levels) {
            let sample_level = factor.sample_levels[sample].as_str();
            for level in levels
                .iter()
                .filter(|level| level.as_str() != factor.reference)
            {
                standard_values.push((sample_level == level) as u8 as f64);
            }
        }
        for covariate in numeric_covariates {
            standard_values.push(covariate.values[sample]);
        }
        for (left_idx, right_idx) in &interaction_indices {
            let left = &factors[*left_idx];
            let right = &factors[*right_idx];
            let left_sample_level = left.sample_levels[sample].as_str();
            let right_sample_level = right.sample_levels[sample].as_str();
            for left_level in factor_levels[*left_idx]
                .iter()
                .filter(|level| level.as_str() != left.reference)
            {
                for right_level in factor_levels[*right_idx]
                    .iter()
                    .filter(|level| level.as_str() != right.reference)
                {
                    standard_values.push(
                        (left_sample_level == left_level && right_sample_level == right_level) as u8
                            as f64,
                    );
                }
            }
        }
        for (factor_idx, numeric_idx) in &factor_numeric_indices {
            let factor = &factors[*factor_idx];
            let numeric = &numeric_covariates[*numeric_idx];
            let sample_level = factor.sample_levels[sample].as_str();
            let numeric_value = numeric.values[sample];
            for level in factor_levels[*factor_idx]
                .iter()
                .filter(|level| level.as_str() != factor.reference)
            {
                standard_values.push(((sample_level == level) as u8 as f64) * numeric_value);
            }
        }
        for (left_idx, right_idx) in &numeric_interaction_indices {
            standard_values.push(
                numeric_covariates[*left_idx].values[sample]
                    * numeric_covariates[*right_idx].values[sample],
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
    coefficient_groups.push(vec![0]);
    let mut expanded_offset = 1;
    for (factor, levels) in factors.iter().zip(&factor_levels) {
        for (level_idx, level) in levels.iter().enumerate() {
            if level.as_str() != factor.reference {
                coefficient_groups.push(vec![expanded_offset + level_idx]);
            }
        }
        expanded_offset += levels.len();
    }
    for idx in 0..numeric_covariates.len() {
        coefficient_groups.push(vec![expanded_offset + idx]);
    }
    expanded_offset += numeric_covariates.len();
    for (left_idx, right_idx) in &interaction_indices {
        let left_levels = &factor_levels[*left_idx];
        let right_levels = &factor_levels[*right_idx];
        for (left_level_idx, left_level) in left_levels.iter().enumerate() {
            if left_level.as_str() == factors[*left_idx].reference {
                continue;
            }
            for (right_level_idx, right_level) in right_levels.iter().enumerate() {
                if right_level.as_str() == factors[*right_idx].reference {
                    continue;
                }
                let expanded_column =
                    expanded_offset + left_level_idx * right_levels.len() + right_level_idx;
                coefficient_groups.push(vec![expanded_column]);
            }
        }
        expanded_offset += left_levels.len() * right_levels.len();
    }
    for (factor_idx, numeric_idx) in &factor_numeric_indices {
        let levels = &factor_levels[*factor_idx];
        let factor = &factors[*factor_idx];
        for (level_idx, level) in levels.iter().enumerate() {
            if level.as_str() != factor.reference {
                coefficient_groups.push(vec![expanded_offset + level_idx]);
            }
        }
        let _ = numeric_idx;
        expanded_offset += levels.len();
    }
    for idx in 0..numeric_interaction_indices.len() {
        coefficient_groups.push(vec![expanded_offset + idx]);
    }

    Ok(ExpandedAdditiveFactorDesign {
        expanded_design,
        standard_design,
        coefficient_groups,
        factor_levels,
        numeric_covariates: numeric_covariates
            .iter()
            .map(|covariate| covariate.name.to_string())
            .collect(),
        interactions: factor_interactions
            .iter()
            .map(|interaction| format!("{}:{}", interaction.left_factor, interaction.right_factor))
            .collect(),
        factor_numeric_interactions: factor_numeric_interactions
            .iter()
            .map(|interaction| format!("{}:{}", interaction.factor, interaction.numeric))
            .collect(),
        numeric_interactions: numeric_interactions
            .iter()
            .map(|interaction| {
                format!("{}:{}", interaction.left_numeric, interaction.right_numeric)
            })
            .collect(),
        higher_order_interactions: Vec::new(),
    })
}
