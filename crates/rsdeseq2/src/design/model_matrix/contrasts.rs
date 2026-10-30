use super::{FactorCoding, invalid};
use crate::errors::DeseqError;

pub(super) struct CodingMatrix {
    pub names: Vec<String>,
    pub values: Vec<f64>,
    treatment_base: Option<usize>,
}
impl CodingMatrix {
    pub fn new(levels: &[String], coding: FactorCoding<'_>) -> Result<Self, DeseqError> {
        let n = levels.len();
        // Treatment columns are indicators, including in expanded native plans.
        // Retaining n*(n-1) zeros/ones would add quadratic storage to a parser
        // migration that only needs level lookup and shared product evaluation.
        if let FactorCoding::Treatment { base } = coding {
            if base >= n {
                return Err(invalid(
                    "factor treatment base is outside its declared levels",
                ));
            }
            let names = levels
                .iter()
                .enumerate()
                .filter(|(i, _)| *i != base)
                .map(|(_, name)| name.clone())
                .collect::<Vec<_>>();
            if names.iter().any(|name| name.contains('\0')) {
                return Err(invalid("factor contrast names/values are invalid"));
            }
            return Ok(Self {
                names,
                values: Vec::new(),
                treatment_base: Some(base),
            });
        }
        let width = match coding {
            FactorCoding::Matrix { values, .. } => values.n_cols(),
            _ => n - 1,
        };
        let mut names = (1..=width).map(|n| n.to_string()).collect::<Vec<_>>();
        let mut values = vec![
            0.0;
            n.checked_mul(width)
                .ok_or_else(|| invalid("factor coding shape overflow"))?
        ];
        match coding {
            FactorCoding::Treatment { .. } => unreachable!("treatment indicators returned above"),
            FactorCoding::Sum => {
                for column in 0..width {
                    values[column * width + column] = 1.0;
                    values[(n - 1) * width + column] = -1.0;
                }
            }
            FactorCoding::Helmert => {
                for column in 0..width {
                    for level in 0..=column {
                        values[level * width + column] = -1.0;
                    }
                    values[(column + 1) * width + column] = (column + 1) as f64;
                }
            }
            FactorCoding::Polynomial => {
                if n > 95 {
                    return Err(invalid(
                        "orthogonal polynomial contrasts require at most 95 levels",
                    ));
                }
                values = super::polynomial::contrasts(n);
                for (column, name) in names.iter_mut().enumerate() {
                    *name = match column {
                        0 => ".L".into(),
                        1 => ".Q".into(),
                        2 => ".C".into(),
                        _ => format!("^{}", column + 1),
                    };
                }
            }
            FactorCoding::Matrix {
                values: matrix,
                names: supplied,
            } => {
                if matrix.n_rows() != n || matrix.as_slice().iter().any(|v| !v.is_finite()) {
                    return Err(invalid(
                        "explicit factor contrasts must be a finite level-by-column matrix",
                    ));
                }
                if let Some(supplied) = supplied {
                    if supplied.len() != width {
                        return Err(invalid(
                            "factor contrast column names have incompatible length",
                        ));
                    }
                    names = supplied.to_vec();
                }
                for level in 0..n {
                    for column in 0..width {
                        values[level * width + column] = *matrix.at(level, column);
                    }
                }
            }
        }
        if names.iter().any(|s| s.contains('\0')) || values.iter().any(|v| !v.is_finite()) {
            return Err(invalid("factor contrast names/values are invalid"));
        }
        Ok(Self {
            names,
            values,
            treatment_base: None,
        })
    }
    pub fn value(&self, level: usize, column: usize) -> f64 {
        if let Some(base) = self.treatment_base {
            let selected = if column < base { column } else { column + 1 };
            (level == selected) as u8 as f64
        } else {
            self.values[level * self.names.len() + column]
        }
    }
}
