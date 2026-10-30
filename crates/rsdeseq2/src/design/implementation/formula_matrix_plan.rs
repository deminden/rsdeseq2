// Translate the compatibility parser's column order and labels into the shared
// model-matrix product plan. No factor indicators/products are calculated here.
fn formula_expanded_matrix_values(
    state: &ExpandedFormulaDesignState<'_>,
    factors: &[ExpandedFactorSpec<'_>],
    numeric: &[ExpandedNumericSpec<'_>],
    levels: &[Vec<String>],
    names: &[String],
    samples: usize,
) -> Result<Vec<f64>, DeseqError> {
    use super::model_matrix::{ColumnPiece, PlannedColumn};
    use super::{FactorCoding, ModelMatrixPlan, ModelMatrixVariable};
    use crate::matrix::{MatrixOrder, MatrixView};
    let codes = factors
        .iter()
        .zip(levels)
        .map(|(factor, levels)| {
            factor
                .sample_levels
                .iter()
                .map(|sample| {
                    levels
                        .iter()
                        .position(|level| sample == level)
                        .and_then(|i| i32::try_from(i + 1).ok())
                        .ok_or_else(|| DeseqError::InvalidOptions {
                            reason: "formula factor code is outside declared levels".into(),
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
            coding: FactorCoding::Treatment { base: 0 },
        })
        .collect::<Vec<_>>();
    for value in numeric {
        variables.push(ModelMatrixVariable::Numeric {
            name: value.name,
            values: MatrixView::new(samples, 1, value.values, MatrixOrder::RowMajor)?,
            names: None,
        });
    }
    let factor_piece = |name: &str, level: &str| -> Result<ColumnPiece, DeseqError> {
        let i = factors
            .iter()
            .position(|f| f.factor == name)
            .ok_or_else(|| DeseqError::InvalidOptions {
                reason: format!("formula factor '{name}' is absent from the matrix plan"),
            })?;
        let column = levels[i].iter().position(|v| v == level).ok_or_else(|| {
            DeseqError::InvalidOptions {
                reason: format!("formula factor level '{level}' is absent from the matrix plan"),
            }
        })?;
        Ok(ColumnPiece {
            variable: i,
            column,
            full: true,
        })
    };
    let numeric_piece = |name: &str| -> Result<ColumnPiece, DeseqError> {
        let i = numeric.iter().position(|v| v.name == name).ok_or_else(|| {
            DeseqError::InvalidOptions {
                reason: format!("formula numeric '{name}' is absent from the matrix plan"),
            }
        })?;
        Ok(ColumnPiece {
            variable: factors.len() + i,
            column: 0,
            full: false,
        })
    };
    let mut pieces = Vec::new();
    if state.has_intercept {
        pieces.push(Vec::new());
    }
    for factor in &state.selected_factors {
        for level in formula_factor_levels(factor.factor, factors, levels)? {
            pieces.push(vec![factor_piece(factor.factor, level)?]);
        }
    }
    for value in &state.selected_numeric_covariates {
        pieces.push(vec![numeric_piece(value.name)?]);
    }
    for interaction in &state.factor_interactions {
        for left in formula_factor_levels(interaction.left_factor, factors, levels)? {
            for right in formula_factor_levels(interaction.right_factor, factors, levels)? {
                pieces.push(vec![
                    factor_piece(interaction.left_factor, left)?,
                    factor_piece(interaction.right_factor, right)?,
                ]);
            }
        }
    }
    for interaction in &state.factor_numeric_interactions {
        for level in formula_factor_levels(interaction.factor, factors, levels)? {
            pieces.push(vec![
                factor_piece(interaction.factor, level)?,
                numeric_piece(interaction.numeric)?,
            ]);
        }
    }
    for interaction in &state.numeric_interactions {
        pieces.push(vec![
            numeric_piece(interaction.left_numeric)?,
            numeric_piece(interaction.right_numeric)?,
        ]);
    }
    for interaction in &state.higher_order_interactions {
        for grid in formula_higher_order_level_grid(interaction, factors, levels)? {
            let column = interaction
                .variables
                .iter()
                .zip(grid)
                .map(|(variable, level)| match variable {
                    FormulaVariableRef::Factor(name) => {
                        factor_piece(name, level.expect("factor grid contains a level"))
                    }
                    FormulaVariableRef::Numeric(name) => numeric_piece(name),
                })
                .collect::<Result<Vec<_>, DeseqError>>()?;
            pieces.push(column);
        }
    }
    if pieces.len() != names.len() {
        return Err(DeseqError::InvalidOptions {
            reason: "formula names and matrix plan columns differ".into(),
        });
    }
    let columns = pieces
        .into_iter()
        .zip(names)
        .enumerate()
        .map(|(index, (mut pieces, name))| {
            pieces.sort_by_key(|piece| {
                let name = if piece.variable < factors.len() {
                    factors[piece.variable].factor
                } else {
                    numeric[piece.variable - factors.len()].name
                };
                state
                    .frame_order
                    .iter()
                    .position(|v| v == name)
                    .unwrap_or(usize::MAX)
            });
            PlannedColumn {
                name: name.clone(),
                assignment: if pieces.is_empty() { 0 } else { index + 1 },
                pieces,
            }
        })
        .collect();
    Ok(ModelMatrixPlan::from_columns(samples, &variables, columns)?.values())
}
