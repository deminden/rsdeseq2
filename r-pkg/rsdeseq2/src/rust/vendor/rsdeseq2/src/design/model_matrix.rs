//! Shared matrix assembly from formula-term metadata and evaluated variables.
//!
//! Formula parsing and expression evaluation are separate from this operation.
//! The same validated plan can be used by native and foreign-language callers.
use crate::errors::DeseqError;
use crate::matrix::MatrixView;

mod contrasts;
mod polynomial;
use contrasts::CodingMatrix;

fn invalid(reason: impl Into<String>) -> DeseqError {
    DeseqError::InvalidOptions {
        reason: reason.into(),
    }
}

/// Factor contrast policy. Explicit matrices are caller-supplied primitives.
#[derive(Clone, Copy, Debug)]
pub enum FactorCoding<'a> {
    /// Omit a zero-based reference level.
    Treatment { base: usize },
    /// Sum-to-zero contrasts, with the last level represented by minus ones.
    Sum,
    /// Helmert contrasts in declared level order.
    Helmert,
    /// Orthogonal polynomial contrasts on scores 1 through the level count.
    Polynomial,
    /// A level-by-column matrix with optional suffixes for its columns.
    Matrix {
        values: MatrixView<'a, f64>,
        names: Option<&'a [String]>,
    },
}

/// A variable already aligned to the sample order.
#[derive(Clone, Copy, Debug)]
pub enum ModelMatrixVariable<'a> {
    /// A numeric vector or matrix. One-column variables retain the term name;
    /// wider variables append column names, or one-based column numbers.
    Numeric {
        name: &'a str,
        values: MatrixView<'a, f64>,
        names: Option<&'a [String]>,
    },
    /// One-based factor codes and their declared level order, including unused
    /// levels. Logical variables can be transported as two-level factors.
    Factor {
        name: &'a str,
        codes: &'a [i32],
        levels: &'a [String],
        coding: FactorCoding<'a>,
    },
}

struct Variable<'a> {
    name: &'a str,
    data: VariableData<'a>,
    suffixes: Vec<String>,
}
enum VariableData<'a> {
    Numeric(MatrixView<'a, f64>),
    Factor {
        codes: &'a [i32],
        levels: &'a [String],
        contrasted: CodingMatrix,
    },
}
#[derive(Clone)]
struct Piece {
    variable: usize,
    column: usize,
    full: bool,
}

/// An explicit column layout supplied by a native term parser. It uses the
/// same variable preparation and product evaluation as the R-style term plan.
pub(crate) struct PlannedColumn {
    pub name: String,
    pub assignment: usize,
    pub pieces: Vec<ColumnPiece>,
}
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct ColumnPiece {
    pub variable: usize,
    pub column: usize,
    pub full: bool,
}

/// Validated model-matrix schema. Terms contain one code per variable: zero
/// omits it, one uses contrasts, and two uses full factor indicators. Within an
/// interaction the first variable's columns vary fastest, as in R model.matrix.
/// Matrix construction does not perform a rank test or change coefficient units.
pub struct ModelMatrixPlan<'a> {
    samples: usize,
    variables: Vec<Variable<'a>>,
    columns: Vec<Vec<Piece>>,
    names: Vec<String>,
    assignments: Vec<usize>,
    response_dropped: bool,
}
impl<'a> ModelMatrixPlan<'a> {
    pub(crate) fn column_layout(&self) -> Vec<Vec<ColumnPiece>> {
        self.columns
            .iter()
            .map(|column| {
                column
                    .iter()
                    .map(|piece| ColumnPiece {
                        variable: piece.variable,
                        column: piece.column,
                        full: piece.full,
                    })
                    .collect()
            })
            .collect()
    }
    pub(crate) fn piece_name(&self, piece: &ColumnPiece) -> String {
        let variable = &self.variables[piece.variable];
        let suffix = match &variable.data {
            VariableData::Factor { levels, .. } if piece.full => &levels[piece.column],
            _ => &variable.suffixes[piece.column],
        };
        format!("{}{suffix}", variable.name)
    }
    pub(crate) fn from_columns(
        samples: usize,
        variables: &[ModelMatrixVariable<'a>],
        columns: Vec<PlannedColumn>,
    ) -> Result<Self, DeseqError> {
        let mut plan = Self::new(samples, variables, &[], false)?;
        samples
            .checked_mul(columns.len())
            .ok_or_else(|| invalid("model matrix shape overflow"))?;
        for column in columns {
            if column.name.contains('\0') {
                return Err(invalid("model matrix column names must contain no NUL"));
            }
            let pieces = column
                .pieces
                .into_iter()
                .map(|piece| {
                    let variable = plan.variables.get(piece.variable).ok_or_else(|| {
                        invalid("model matrix piece variable is outside the plan")
                    })?;
                    let width = match &variable.data {
                        VariableData::Numeric(matrix) => matrix.n_cols(),
                        VariableData::Factor {
                            levels, contrasted, ..
                        } => {
                            if piece.full {
                                levels.len()
                            } else {
                                contrasted.names.len()
                            }
                        }
                    };
                    if piece.column >= width {
                        return Err(invalid("model matrix piece column is outside its variable"));
                    }
                    Ok(Piece {
                        variable: piece.variable,
                        column: piece.column,
                        full: piece.full,
                    })
                })
                .collect::<Result<Vec<_>, DeseqError>>()?;
            plan.assignments.push(column.assignment);
            plan.columns.push(pieces);
            plan.names.push(column.name);
        }
        Ok(plan)
    }
    /// Prepare a matrix without copying numeric input matrices. An absent
    /// intercept promotes the first participating factor to full indicators.
    /// Terms may be empty, yielding an intercept-only or zero-column matrix.
    pub fn new(
        samples: usize,
        variables: &[ModelMatrixVariable<'a>],
        terms: &[Vec<i32>],
        intercept: bool,
    ) -> Result<Self, DeseqError> {
        Self::with_response(samples, variables, terms, intercept, None)
    }
    /// Assemble term metadata with an optional zero-based response variable.
    /// The response is transported numerically, contributes to interactions,
    /// and is omitted when it appears as a standalone right-hand-side term.
    pub fn with_response(
        samples: usize,
        variables: &[ModelMatrixVariable<'a>],
        terms: &[Vec<i32>],
        intercept: bool,
        response: Option<usize>,
    ) -> Result<Self, DeseqError> {
        if response.is_some_and(|r| {
            r >= variables.len() || matches!(variables[r], ModelMatrixVariable::Factor { .. })
        }) {
            return Err(invalid(
                "model matrix response must identify a numeric variable",
            ));
        }
        if samples == 0 {
            return Err(invalid("model matrix requires samples"));
        }
        for term in terms {
            if term.len() != variables.len() || term.iter().any(|c| !(0..=2).contains(c)) {
                return Err(invalid(
                    "model matrix terms must contain one code (0, 1 or 2) per variable",
                ));
            }
        }
        let promoted = if intercept {
            None
        } else {
            terms.iter().enumerate().find_map(|(t, term)| {
                variables.iter().enumerate().find_map(|(i, v)| {
                    (term[i] != 0 && matches!(v, ModelMatrixVariable::Factor { .. }))
                        .then_some((t, i))
                })
            })
        };
        let mut prepared = Vec::with_capacity(variables.len());
        for variable in variables {
            let (name, data, suffixes) = match *variable {
                ModelMatrixVariable::Numeric {
                    name,
                    values,
                    names,
                } => {
                    if values.n_rows() != samples
                        || values.as_slice().iter().any(|v| !v.is_finite())
                    {
                        return Err(invalid(
                            "model matrix numeric variables must have finite values and one row per sample",
                        ));
                    }
                    let suffixes = numeric_suffixes(values.n_cols(), names)?;
                    (name, VariableData::Numeric(values), suffixes)
                }
                ModelMatrixVariable::Factor {
                    name,
                    codes,
                    levels,
                    coding,
                } => {
                    if codes.len() != samples
                        || levels.len() < 2
                        || codes.iter().any(|&v| v < 1 || v as usize > levels.len())
                    {
                        return Err(invalid(
                            "model matrix factors require valid one-based codes and at least two declared levels",
                        ));
                    }
                    if levels
                        .iter()
                        .enumerate()
                        .any(|(i, v)| levels[..i].contains(v))
                    {
                        return Err(invalid("model matrix factor levels must be unique"));
                    }
                    let contrasted = CodingMatrix::new(levels, coding)?;
                    let suffixes = contrasted.names.clone();
                    (
                        name,
                        VariableData::Factor {
                            codes,
                            levels,
                            contrasted,
                        },
                        suffixes,
                    )
                }
            };
            if name.is_empty() || name.contains('\0') {
                return Err(invalid(
                    "model matrix variable names must be nonempty and contain no NUL",
                ));
            }
            prepared.push(Variable {
                name,
                data,
                suffixes,
            });
        }
        let mut plan = Self {
            samples,
            variables: prepared,
            columns: Vec::new(),
            names: Vec::new(),
            assignments: Vec::new(),
            response_dropped: false,
        };
        if intercept {
            plan.columns.push(Vec::new());
            plan.names.push("(Intercept)".into());
            plan.assignments.push(0);
        }
        for (term_index, term) in terms.iter().enumerate() {
            if let Some(response) = response
                && term[response] != 0
                && term
                    .iter()
                    .enumerate()
                    .all(|(i, &code)| i == response || code == 0)
            {
                plan.response_dropped = true;
                continue;
            }
            let mut combinations = vec![Vec::new()];
            for (index, &code) in term.iter().enumerate().filter(|(_, code)| **code != 0) {
                let variable = &plan.variables[index];
                let full = code == 2 || promoted == Some((term_index, index));
                let width = match &variable.data {
                    VariableData::Factor { levels, .. } if full => levels.len(),
                    _ => variable.suffixes.len(),
                };
                let count = combinations
                    .len()
                    .checked_mul(width)
                    .ok_or_else(|| invalid("model matrix interaction shape overflow"))?;
                samples
                    .checked_mul(count)
                    .ok_or_else(|| invalid("model matrix shape overflow"))?;
                let mut next = Vec::with_capacity(count);
                for column in 0..width {
                    for previous in &combinations {
                        let mut pieces = previous.clone();
                        pieces.push(Piece {
                            variable: index,
                            column,
                            full,
                        });
                        next.push(pieces);
                    }
                }
                combinations = next;
            }
            // A zero-code column is not a formula term; it must not create a
            // second intercept with an unrelated assignment.
            if combinations[0].is_empty() {
                return Err(invalid("model matrix term contains no variables"));
            }
            for pieces in combinations {
                let name = pieces
                    .iter()
                    .map(|piece| {
                        plan.piece_name(&ColumnPiece {
                            variable: piece.variable,
                            column: piece.column,
                            full: piece.full,
                        })
                    })
                    .collect::<Vec<_>>()
                    .join(":");
                plan.names.push(name);
                plan.assignments.push(term_index + 1);
                plan.columns.push(pieces);
            }
        }
        samples
            .checked_mul(plan.names.len())
            .ok_or_else(|| invalid("model matrix shape overflow"))?;
        if plan.names.iter().any(|name| name.contains('\0')) {
            return Err(invalid("model matrix column names must contain no NUL"));
        }
        Ok(plan)
    }
    /// Sample count.
    pub fn n_rows(&self) -> usize {
        self.samples
    }
    /// Coefficient count; zero is allowed before fitting validation.
    pub fn n_cols(&self) -> usize {
        self.names.len()
    }
    /// Exact unsanitized column labels in matrix order.
    pub fn column_names(&self) -> &[String] {
        &self.names
    }
    /// Zero for intercept, otherwise the one-based formula-term index.
    pub fn assignments(&self) -> &[usize] {
        &self.assignments
    }
    /// Whether response-only right-hand-side terms were omitted.
    pub fn response_term_dropped(&self) -> bool {
        self.response_dropped
    }
    /// Materialize in row-major order, retaining multiplication order and
    /// signed-zero behavior within each interaction.
    pub fn values(&self) -> Vec<f64> {
        let mut values = Vec::with_capacity(self.samples * self.n_cols());
        for sample in 0..self.samples {
            for pieces in &self.columns {
                let mut value = None;
                for piece in pieces {
                    let next = match &self.variables[piece.variable].data {
                        VariableData::Numeric(matrix) => *matrix.at(sample, piece.column),
                        VariableData::Factor {
                            codes, contrasted, ..
                        } => {
                            let level = codes[sample] as usize - 1;
                            if piece.full {
                                (level == piece.column) as u8 as f64
                            } else {
                                contrasted.value(level, piece.column)
                            }
                        }
                    };
                    value = Some(value.map_or(next, |previous: f64| previous * next));
                }
                values.push(value.unwrap_or(1.0));
            }
        }
        values
    }
}
fn numeric_suffixes(columns: usize, names: Option<&[String]>) -> Result<Vec<String>, DeseqError> {
    if let Some(names) = names
        && names.len() != columns
    {
        return Err(invalid(
            "model matrix numeric column names have incompatible length",
        ));
    }
    Ok(if columns == 1 {
        vec![String::new()]
    } else {
        names
            .map(<[String]>::to_vec)
            .unwrap_or_else(|| (1..=columns).map(|n| n.to_string()).collect())
    })
}

#[cfg(test)]
mod explicit_columns_tests {
    use super::*;
    use crate::matrix::MatrixOrder;

    #[test]
    fn explicit_column_pieces_are_validated_before_materialization() {
        let values = [1.0, -0.0];
        let variable = ModelMatrixVariable::Numeric {
            name: "x",
            values: MatrixView::new(2, 1, &values, MatrixOrder::RowMajor).unwrap(),
            names: None,
        };
        for (name, variable_index, column) in [("x", 1, 0), ("x", 0, 1), ("x\0", 0, 0)] {
            assert!(
                ModelMatrixPlan::from_columns(
                    2,
                    &[variable],
                    vec![PlannedColumn {
                        name: name.into(),
                        assignment: 1,
                        pieces: vec![ColumnPiece {
                            variable: variable_index,
                            column,
                            full: false
                        }],
                    }]
                )
                .is_err()
            );
        }
        let result = ModelMatrixPlan::from_columns(
            2,
            &[variable],
            vec![
                PlannedColumn {
                    name: "intercept".into(),
                    assignment: 0,
                    pieces: vec![],
                },
                PlannedColumn {
                    name: "x".into(),
                    assignment: 1,
                    pieces: vec![ColumnPiece {
                        variable: 0,
                        column: 0,
                        full: false,
                    }],
                },
            ],
        )
        .unwrap();
        assert_eq!(
            result
                .values()
                .iter()
                .map(|v| v.to_bits())
                .collect::<Vec<_>>(),
            [1.0f64, 1.0, 1.0, -0.0].map(f64::to_bits)
        );
    }
}
