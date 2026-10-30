// Ordered factors use the same polynomial coding and product plan as the R
// primitive interface. Unordered expanded columns retain full indicators.
fn formula_ordered_design(
    state: &ExpandedFormulaDesignState<'_>,
    factors: &[ExpandedFactorSpec<'_>],
    numeric: &[ExpandedNumericSpec<'_>],
    derived: &[FormulaDerivedFactorCovariate],
    samples: usize,
) -> Result<ExpandedAdditiveFactorDesign, DeseqError> {
    use super::{FactorCoding, ModelMatrixPlan, ModelMatrixVariable};
    use crate::matrix::{MatrixOrder, MatrixView};
    let ordered = factors
        .iter()
        .map(|factor| {
            derived
                .iter()
                .any(|d| d.name == factor.factor && d.primitive.has_class("ordered"))
        })
        .collect::<Vec<_>>();
    let levels = factors
        .iter()
        .map(|f| ordered_levels_with_declared(f.sample_levels, f.reference, f.levels))
        .collect::<Vec<_>>();
    let codes = factors
        .iter()
        .zip(&levels)
        .map(|(factor, levels)| {
            factor
                .sample_levels
                .iter()
                .map(|sample| {
                    levels
                        .iter()
                        .position(|level| level == sample)
                        .and_then(|i| i32::try_from(i + 1).ok())
                        .ok_or_else(|| DeseqError::InvalidOptions {
                            reason: "formula factor code is outside levels".into(),
                        })
                })
                .collect::<Result<Vec<_>, _>>()
        })
        .collect::<Result<Vec<_>, DeseqError>>()?;
    let mut variables = factors
        .iter()
        .enumerate()
        .map(|(i, factor)| ModelMatrixVariable::Factor {
            name: factor.factor,
            codes: &codes[i],
            levels: &levels[i],
            coding: if ordered[i] {
                FactorCoding::Polynomial
            } else {
                FactorCoding::Treatment { base: 0 }
            },
        })
        .collect::<Vec<_>>();
    for value in numeric {
        variables.push(ModelMatrixVariable::Numeric {
            name: value.name,
            values: MatrixView::new(samples, 1, value.values, MatrixOrder::RowMajor)?,
            names: None,
        });
    }
    // Multiplication follows frame-variable order, including numeric variables
    // preceding ordered factors. Reassociation can otherwise change bits.
    let mut variable_order = (0..variables.len()).collect::<Vec<_>>();
    variable_order.sort_by_key(|&i| {
        let name = if i < factors.len() {
            factors[i].factor
        } else {
            numeric[i - factors.len()].name
        };
        state
            .frame_order
            .iter()
            .position(|v| v == name)
            .unwrap_or(usize::MAX)
    });
    let mut positions = vec![0; variables.len()];
    for (position, &source) in variable_order.iter().enumerate() {
        positions[source] = position;
    }
    let variables = variable_order
        .iter()
        .map(|&i| variables[i])
        .collect::<Vec<_>>();
    let factor_index = |name: &str| {
        factors
            .iter()
            .position(|f| f.factor == name)
            .map(|i| positions[i])
            .ok_or_else(|| DeseqError::InvalidOptions {
                reason: format!("formula factor '{name}' is absent from the plan"),
            })
    };
    let numeric_index = |name: &str| {
        numeric
            .iter()
            .position(|v| v.name == name)
            .map(|i| positions[factors.len() + i])
            .ok_or_else(|| DeseqError::InvalidOptions {
                reason: format!("formula covariate '{name}' is absent from the plan"),
            })
    };
    let empty = || vec![0; variables.len()];
    let mut terms = Vec::new();
    for factor in &state.selected_factors {
        let mut term = empty();
        term[factor_index(factor.factor)?] = 1;
        terms.push(term);
    }
    for value in &state.selected_numeric_covariates {
        let mut term = empty();
        term[numeric_index(value.name)?] = 1;
        terms.push(term);
    }
    for interaction in &state.factor_interactions {
        let mut term = empty();
        term[factor_index(interaction.left_factor)?] =
            if formula_has_factor_main_effect(state, interaction.right_factor) {
                1
            } else {
                2
            };
        term[factor_index(interaction.right_factor)?] =
            if formula_has_factor_main_effect(state, interaction.left_factor) {
                1
            } else {
                2
            };
        terms.push(term);
    }
    for interaction in &state.factor_numeric_interactions {
        let mut term = empty();
        term[factor_index(interaction.factor)?] =
            if formula_has_numeric_main_effect(state, interaction.numeric) {
                1
            } else {
                2
            };
        term[numeric_index(interaction.numeric)?] = 1;
        terms.push(term);
    }
    for interaction in &state.numeric_interactions {
        let mut term = empty();
        term[numeric_index(interaction.left_numeric)?] = 1;
        term[numeric_index(interaction.right_numeric)?] = 1;
        terms.push(term);
    }
    for interaction in &state.higher_order_interactions {
        let mut term = empty();
        for variable in &interaction.variables {
            match variable {
                FormulaVariableRef::Factor(name) => {
                    term[factor_index(name)?] = if formula_has_factor_main_effect(state, name) {
                        1
                    } else {
                        2
                    }
                }
                FormulaVariableRef::Numeric(name) => term[numeric_index(name)?] = 1,
            }
        }
        terms.push(term);
    }
    // A factor in an interaction uses contrasts when the complementary term
    // is already in the model. Its own main effect alone is insufficient for
    // higher-order products (e.g. a+b+x+a:b:x still uses full a/b indicators).
    let memberships = formula_term_variables(state);
    for (term_index, term) in terms.iter_mut().enumerate() {
        if term.iter().filter(|&&code| code != 0).count() < 2 {
            continue;
        }
        for factor in 0..factors.len() {
            let slot = positions[factor];
            if term[slot] == 0 {
                continue;
            }
            let has_complement = formula_factor_has_complement(
                state,
                factors[factor].factor,
                &memberships[term_index],
            );
            term[slot] = if has_complement { 1 } else { 2 };
        }
    }
    let standard = ModelMatrixPlan::new(samples, &variables, &terms, state.has_intercept)?;
    let expanded_terms = terms
        .iter()
        .map(|term| {
            term.iter()
                .enumerate()
                .map(|(i, &code)| {
                    if code != 0 && variable_order[i] < factors.len() && !ordered[variable_order[i]]
                    {
                        2
                    } else {
                        code
                    }
                })
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let expanded = ModelMatrixPlan::new(samples, &variables, &expanded_terms, state.has_intercept)?;
    let standard_layout = standard.column_layout();
    let expanded_layout = expanded.column_layout();
    let standard_names = standard_layout
        .iter()
        .map(|pieces| {
            if pieces.is_empty() {
                return "Intercept".to_string();
            }
            pieces
                .iter()
                .map(|piece| {
                    let source = variable_order[piece.variable];
                    if source < factors.len() && !ordered[source] && !piece.full {
                        let factor = &factors[source];
                        standard_factor_coefficient_name(
                            factor.factor,
                            &levels[source][piece.column + 1],
                            factor.reference,
                        )
                    } else {
                        standard.piece_name(piece)
                    }
                })
                .collect::<Vec<_>>()
                .join(":")
        })
        .collect::<Vec<_>>();
    let expanded_names = expanded
        .column_names()
        .iter()
        .map(|name| {
            if name == "(Intercept)" {
                "Intercept".to_string()
            } else {
                name.clone()
            }
        })
        .collect::<Vec<_>>();
    let groups = standard_layout
        .into_iter()
        .map(|mut pieces| {
            for piece in &mut pieces {
                let source = variable_order[piece.variable];
                if source < factors.len() && !ordered[source] && !piece.full {
                    piece.column += 1;
                    piece.full = true;
                }
            }
            expanded_layout
                .iter()
                .position(|column| *column == pieces)
                .map(|i| vec![i])
                .ok_or_else(|| DeseqError::InvalidOptions {
                    reason: "standard formula column is absent from the expanded plan".into(),
                })
        })
        .collect::<Result<Vec<_>, _>>()?;
    let standard_design = DesignMatrix::from_row_major(
        samples,
        standard.n_cols(),
        standard.values(),
        Some(standard_names),
    )?;
    let expanded_design = DesignMatrix::from_row_major(
        samples,
        expanded.n_cols(),
        expanded.values(),
        Some(expanded_names),
    )?;
    Ok(formula_design_output(
        state,
        expanded_design,
        standard_design,
        groups,
    ))
}
