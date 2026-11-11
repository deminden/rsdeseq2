// Complete workflows preserve term and matrix-valued variable order while the
// expanded compatibility APIs retain their established explicit column layout.

/// Build a formula design for a complete fitting workflow in R term order.
/// Products and factor coding use the shared model-matrix implementation.
pub fn formula_design_for_fitting_from_model_frame(
    formula: &str,
    model_frame: &FormulaModelFrame,
) -> Result<ExpandedFormulaDesignWithOffsets, DeseqError> {
    let resolved = ResolvedFormulaModelFrame::new(model_frame)?;
    formula_design_with_column_order(
        formula,
        &resolved.factor_specs,
        &resolved.numeric_specs,
        true,
    )
}

fn refresh_formula_term_order(state: &mut ExpandedFormulaDesignState<'_>) {
    if state.term_order.is_none() {
        return;
    }
    let current = formula_term_variables(state)
        .into_iter()
        .map(|names| {
            let mut names = names.into_iter().map(str::to_owned).collect::<Vec<_>>();
            names.sort();
            names
        })
        .collect::<Vec<_>>();
    let order = state.term_order.as_mut().expect("term order enabled");
    order.retain(|term| current.contains(term));
    for term in current {
        if !order.contains(&term) {
            order.push(term);
        }
    }
}

struct FittingNumericGroup<'a> {
    key: String,
    group: Option<usize>,
    columns: Vec<&'a ExpandedNumericSpec<'a>>,
}

fn formula_fitting_design<'a>(
    state: &ExpandedFormulaDesignState<'a>,
    factors: &'a [ExpandedFactorSpec<'a>],
    numeric: &'a [ExpandedNumericSpec<'a>],
    derived_factors: &[FormulaDerivedFactorCovariate],
    derived_numeric: &[FormulaDerivedNumericCovariate],
) -> Result<ExpandedAdditiveFactorDesign, DeseqError> {
    use super::{FactorCoding, ModelMatrixPlan, ModelMatrixVariable};
    use crate::matrix::{MatrixOrder, MatrixView};
    let used_factors = formula_used_factors(state, factors)?;
    let used_numeric = formula_used_numeric_covariates(state, numeric)?;
    let samples =
        formula_design_sample_count(state, &used_factors, &used_numeric, factors, numeric)?;
    validate_formula_used_factor_specs(&used_factors, samples)?;
    validate_numeric_covariate_specs(&used_numeric, samples)?;
    let levels = used_factors
        .iter()
        .map(|f| ordered_levels_with_declared(f.sample_levels, f.reference, f.levels))
        .collect::<Vec<_>>();
    let ordered = used_factors
        .iter()
        .map(|f| {
            derived_factors
                .iter()
                .any(|d| d.name == f.factor && d.primitive.has_class("ordered"))
        })
        .collect::<Vec<_>>();
    let codes = used_factors
        .iter()
        .zip(&levels)
        .map(|(f, levels)| {
            f.sample_levels
                .iter()
                .map(|v| {
                    levels
                        .iter()
                        .position(|level| level == v)
                        .and_then(|i| i32::try_from(i + 1).ok())
                        .ok_or_else(|| DeseqError::InvalidOptions {
                            reason: "formula factor code is outside levels".into(),
                        })
                })
                .collect::<Result<Vec<_>, _>>()
        })
        .collect::<Result<Vec<_>, _>>()?;
    // poly() is one matrix-valued variable even though the compatibility
    // parser materializes its columns separately. Preserve that grouping so
    // interaction columns vary in the same order as the primitive R plan.
    let mut groups: Vec<FittingNumericGroup<'_>> = Vec::new();
    for v in &used_numeric {
        let group = derived_numeric
            .iter()
            .find(|d| d.name == v.name)
            .and_then(|d| d.group);
        let index = groups
            .iter()
            .position(|g| g.group == group && (group.is_some() || g.key == v.name));
        if let Some(index) = index {
            groups[index].columns.push(v);
        } else {
            groups.push(FittingNumericGroup {
                key: v.name.to_owned(),
                group,
                columns: vec![v],
            });
        }
    }
    let values = groups
        .iter()
        .map(|g| {
            (0..samples)
                .flat_map(|row| g.columns.iter().map(move |c| c.values[row]))
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let mut variables = Vec::new();
    for (i, f) in used_factors.iter().enumerate() {
        variables.push(ModelMatrixVariable::Factor {
            name: f.factor,
            codes: &codes[i],
            levels: &levels[i],
            coding: if ordered[i] {
                FactorCoding::Polynomial
            } else {
                FactorCoding::Treatment { base: 0 }
            },
        });
    }
    for (i, g) in groups.iter().enumerate() {
        variables.push(ModelMatrixVariable::Numeric {
            name: &g.key,
            values: MatrixView::new(samples, g.columns.len(), &values[i], MatrixOrder::RowMajor)?,
            names: None,
        });
    }
    let frame_position = |name: &str| {
        state
            .frame_order
            .iter()
            .position(|v| v == name)
            .unwrap_or(usize::MAX)
    };
    let mut sources = (0..variables.len()).collect::<Vec<_>>();
    sources.sort_by_key(|&i| {
        if i < used_factors.len() {
            frame_position(used_factors[i].factor)
        } else {
            groups[i - used_factors.len()]
                .columns
                .iter()
                .map(|c| frame_position(c.name))
                .min()
                .unwrap_or(usize::MAX)
        }
    });
    let mut positions = vec![0; variables.len()];
    for (position, &source) in sources.iter().enumerate() {
        positions[source] = position;
    }
    let variables = sources.iter().map(|&i| variables[i]).collect::<Vec<_>>();
    let variable_index = |name: &str| -> Result<usize, DeseqError> {
        if let Some(i) = used_factors.iter().position(|f| f.factor == name) {
            return Ok(positions[i]);
        }
        groups
            .iter()
            .position(|g| g.columns.iter().any(|c| c.name == name))
            .map(|i| positions[used_factors.len() + i])
            .ok_or_else(|| DeseqError::InvalidOptions {
                reason: format!("formula variable '{name}' is absent from fitting plan"),
            })
    };
    let resolve = |names: &[&str]| -> Result<Vec<usize>, DeseqError> {
        let mut indices = names
            .iter()
            .map(|name| variable_index(name))
            .collect::<Result<Vec<_>, _>>()?;
        indices.sort();
        indices.dedup();
        Ok(indices)
    };
    let mut memberships = Vec::new();
    for names in formula_term_variables(state) {
        let indices = resolve(&names)?;
        if !memberships.contains(&indices) {
            memberships.push(indices);
        }
    }
    let order = state
        .term_order
        .as_ref()
        .ok_or_else(|| DeseqError::InvalidOptions {
            reason: "fitting formula term order is absent".into(),
        })?
        .iter()
        .map(|names| resolve(&names.iter().map(String::as_str).collect::<Vec<_>>()))
        .collect::<Result<Vec<_>, _>>()?;
    memberships.sort_by_key(|term| {
        (
            term.len(),
            order.iter().position(|v| v == term).unwrap_or(usize::MAX),
        )
    });
    let terms = memberships
        .iter()
        .enumerate()
        .map(|(term_index, members)| {
            let mut term = vec![0; variables.len()];
            for &index in members {
                term[index] = if matches!(variables[index], ModelMatrixVariable::Factor { .. })
                    && members.len() > 1
                {
                    let complement = members
                        .iter()
                        .copied()
                        .filter(|&i| i != index)
                        .collect::<Vec<_>>();
                    if memberships[..term_index]
                        .iter()
                        .any(|preceding| complement.iter().all(|v| preceding.contains(v)))
                    {
                        1
                    } else {
                        2
                    }
                } else {
                    1
                };
            }
            term
        })
        .collect::<Vec<_>>();
    let expanded_terms = terms
        .iter()
        .map(|term| {
            term.iter()
                .enumerate()
                .map(|(i, &code)| {
                    if code != 0 && sources[i] < used_factors.len() && !ordered[sources[i]] {
                        2
                    } else {
                        code
                    }
                })
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let standard = ModelMatrixPlan::new(samples, &variables, &terms, state.has_intercept)?;
    let expanded = ModelMatrixPlan::new(samples, &variables, &expanded_terms, state.has_intercept)?;
    let standard_layout = standard.column_layout();
    let expanded_layout = expanded.column_layout();
    let names = |plan: &ModelMatrixPlan<'_>, layout: &[Vec<super::model_matrix::ColumnPiece>]| {
        layout
            .iter()
            .map(|pieces| {
                if pieces.is_empty() {
                    return "Intercept".into();
                }
                pieces
                    .iter()
                    .map(|p| {
                        let source = sources[p.variable];
                        if source >= used_factors.len() {
                            groups[source - used_factors.len()].columns[p.column]
                                .name
                                .to_owned()
                        } else if !ordered[source] && !p.full {
                            standard_factor_coefficient_name(
                                used_factors[source].factor,
                                &levels[source][p.column + 1],
                                used_factors[source].reference,
                            )
                        } else {
                            plan.piece_name(p)
                        }
                    })
                    .collect::<Vec<_>>()
                    .join(":")
            })
            .collect::<Vec<String>>()
    };
    let standard_names = names(&standard, &standard_layout);
    let expanded_names = names(&expanded, &expanded_layout);
    validate_unique_coefficient_names(&standard_names, "fitting formula design")?;
    validate_unique_coefficient_names(&expanded_names, "expanded fitting formula design")?;
    let groups = standard_layout
        .into_iter()
        .map(|mut pieces| {
            for p in &mut pieces {
                let source = sources[p.variable];
                if source < used_factors.len() && !ordered[source] && !p.full {
                    p.column += 1;
                    p.full = true;
                }
            }
            expanded_layout
                .iter()
                .position(|p| *p == pieces)
                .map(|i| vec![i])
                .ok_or_else(|| DeseqError::InvalidOptions {
                    reason: "standard formula column is absent from expanded fitting plan".into(),
                })
        })
        .collect::<Result<Vec<_>, _>>()?;
    Ok(formula_design_output(
        state,
        DesignMatrix::from_row_major(
            samples,
            expanded.n_cols(),
            expanded.values(),
            Some(expanded_names),
        )?,
        DesignMatrix::from_row_major(
            samples,
            standard.n_cols(),
            standard.values(),
            Some(standard_names),
        )?,
        groups,
    ))
}
