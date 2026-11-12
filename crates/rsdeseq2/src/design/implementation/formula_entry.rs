// Public formula design entry points.

/// Build an expanded additive design from a primitive DESeq2-style formula subset.
///
/// Supported right-hand-side terms are intercept-preserving main effects
/// (`condition`, `dose`), intercept-only `1`, pairwise interactions
/// (`condition:dose`), nested shorthand (`condition/batch`), R nesting
/// operator terms (`batch %in% condition`), and `*` shorthand for main effects
/// plus interactions (`condition*dose`).
/// Interaction variables can appear without corresponding main effects. The
/// reported standard-design interaction columns then follow R model-matrix
/// treatment-coding shape for the supported primitive terms. Intercept removal
/// with `0` or `-1` and R-style intercept restoration with `+1` or `-0` are
/// supported for these primitive terms. Formula interactions can contain two
/// or more variables. Primitive `- term` subtraction, signed additive
/// parenthesized groups such as `(condition + batch - batch)`, and R-like
/// unary `+` / `-` sign sequences are supported for the same term subset.
/// Signed additive groups are reduced inside supported `*`, `:`, `/`, `%in%`,
/// and power expressions before operator products are formed.
/// Plain `I(numeric)` and signed `I(-numeric)` terms, simple
/// `I(numeric op scalar)` arithmetic, and integer numeric power transforms,
/// raw polynomial transforms, and default orthogonal polynomial transforms are
/// materialized as derived numeric covariates. Supported
/// parenthesized formula powers such as `(condition + batch + dose)^2` expand
/// into main effects plus interactions up to the requested order, with signed
/// additive bases normalized before expansion. Dot formula terms such as `.`,
/// `condition:.`, `condition*.`, and `.^k` use the supplied factors and
/// numeric covariates as the base variable set.
/// Supported
/// `offset(numeric)` and single-vector `offset(transform(numeric))` terms are
/// parsed by [`expanded_formula_design_with_offsets`]; this compatibility
/// helper returns only the design matrix.
pub fn expanded_formula_design<'a>(
    formula: &str,
    factors: &'a [ExpandedFactorSpec<'a>],
    numeric_covariates: &'a [ExpandedNumericSpec<'a>],
) -> Result<ExpandedAdditiveFactorDesign, DeseqError> {
    Ok(expanded_formula_design_with_offsets(formula, factors, numeric_covariates)?.design)
}

/// Build an expanded formula design and return supported per-sample offsets.
///
/// Supported offsets are `offset(numeric)` terms where `numeric` names a
/// supplied numeric covariate, plus single-vector offsets from the supported
/// numeric transform subset such as `offset(log2(numeric))` or
/// `offset(I(numeric + other_numeric))`. Multiple offset terms are summed
/// sample-wise.
pub fn expanded_formula_design_with_offsets<'a>(
    formula: &str,
    factors: &'a [ExpandedFactorSpec<'a>],
    numeric_covariates: &'a [ExpandedNumericSpec<'a>],
) -> Result<ExpandedFormulaDesignWithOffsets, DeseqError> {
    formula_design_with_column_order(formula, factors, numeric_covariates, false)
}

fn formula_design_with_column_order<'a>(
    formula: &str,
    factors: &'a [ExpandedFactorSpec<'a>],
    numeric_covariates: &'a [ExpandedNumericSpec<'a>],
    fitting_order: bool,
) -> Result<ExpandedFormulaDesignWithOffsets, DeseqError> {
    validate_expansion_input(formula)?;
    let rhs = formula_rhs(formula)?;
    let (rhs, offsets) = extract_formula_offsets(rhs, factors, numeric_covariates)?;
    let (rhs, derived_factor_covariates) = expand_formula_factor_transform_terms(&rhs, factors)?;
    let mut all_factors = Vec::with_capacity(factors.len() + derived_factor_covariates.len());
    all_factors.extend(factors.iter().cloned());
    all_factors.extend(
        derived_factor_covariates
            .iter()
            .map(|factor| ExpandedFactorSpec {
                factor: factor.name.as_str(),
                sample_levels: factor.sample_levels.as_slice(),
                reference: factor.reference.as_str(),
                levels: factor.levels.as_deref(),
            }),
    );
    validate_formula_factor_name_collisions(&all_factors)?;
    let (rhs, derived_numeric_covariates) = expand_formula_numeric_transform_terms(
        &rhs,
        &all_factors,
        numeric_covariates,
        fitting_order,
    )?;
    let mut all_numeric_covariates =
        Vec::with_capacity(numeric_covariates.len() + derived_numeric_covariates.len());
    all_numeric_covariates.extend(numeric_covariates.iter().cloned());
    all_numeric_covariates.extend(derived_numeric_covariates.iter().map(|covariate| {
        ExpandedNumericSpec {
            name: covariate.name.as_str(),
            values: covariate.values.as_slice(),
        }
    }));
    validate_formula_numeric_covariate_name_collisions(&all_numeric_covariates)?;
    validate_formula_factor_numeric_name_collisions(&all_factors, &all_numeric_covariates)?;
    let mut state = ExpandedFormulaDesignState {
        term_order: fitting_order.then(Vec::new),
        ..Default::default()
    };

    let rhs = expand_parenthesized_formula_terms(&rhs)?;
    state.frame_order = super::formula_origin::formula_variables(&format!("~{rhs}"))?;
    for (sign, raw_term) in split_formula_signed_terms(&rhs)? {
        (|| -> Result<(), DeseqError> {
            let term = raw_term.trim();
            if term.is_empty() {
                return Err(DeseqError::InvalidOptions {
                    reason: format!("formula '{formula}' contains an empty term"),
                });
            }
            if sign < 0 {
                remove_formula_term(term, &all_factors, &all_numeric_covariates, &mut state)?;
                return Ok(());
            }
            if term == "1" {
                state.has_intercept = true;
                return Ok(());
            }
            if term == "0" || term == "-1" {
                state.has_intercept = false;
                return Ok(());
            }
            if formula_contains_top_level(term, '^')? {
                add_power_formula_term(term, &all_factors, &all_numeric_covariates, &mut state)?;
                return Ok(());
            }
            if formula_contains_top_level(term, '(')? || formula_contains_top_level(term, ')')? {
                return Err(DeseqError::InvalidOptions {
                    reason: format!(
                        "formula term '{term}' is not supported by the primitive parser"
                    ),
                });
            }
            if formula_contains_top_level(term, '/')? {
                add_nested_formula_term(term, &all_factors, &all_numeric_covariates, &mut state)?;
                return Ok(());
            }
            if formula_contains_top_level(term, '*')? {
                add_star_formula_term(term, &all_factors, &all_numeric_covariates, &mut state)?;
                return Ok(());
            }
            if formula_contains_top_level(term, ':')? {
                add_interaction_formula_term(
                    term,
                    &all_factors,
                    &all_numeric_covariates,
                    &mut state,
                )?;
                return Ok(());
            }
            add_main_formula_term(
                term,
                &all_factors,
                &all_numeric_covariates,
                &mut state.selected_factors,
                &mut state.selected_numeric_covariates,
            )?;
            Ok(())
        })()?;
        refresh_formula_term_order(&mut state);
    }

    let model_frame = formula_model_frame_from_specs(&all_factors, &all_numeric_covariates);
    let mut design = if fitting_order {
        formula_fitting_design(
            &state,
            &all_factors,
            &all_numeric_covariates,
            &derived_factor_covariates,
            &derived_numeric_covariates,
        )?
    } else {
        expanded_formula_design_from_state(
            &state,
            &all_factors,
            &all_numeric_covariates,
            &derived_factor_covariates,
        )?
    };
    let origin = DesignOrigin::formula(super::formula_origin::formula_variables(formula)?)?
        .with_ordered_factors(
            derived_factor_covariates
                .iter()
                .filter(|f| f.primitive.has_class("ordered"))
                .map(|f| f.name.clone())
                .collect(),
        );
    design.standard_design = design.standard_design.with_origin(origin.clone());
    design.expanded_design = design.expanded_design.with_origin(origin);
    Ok(ExpandedFormulaDesignWithOffsets {
        design,
        offsets,
        model_frame,
    })
}
