//! Numerical formula transforms shared independently of formula parsing.
//!
//! Inputs are already evaluated finite covariates in sample order. These
//! routines do not alter fitting rank thresholds or dispatch to a frontend.
use crate::errors::DeseqError;
use crate::matrix::{MatrixView, RowMajorMatrix};

mod polynomial;
pub(super) mod polynomial_qr;
mod reduction;
pub use polynomial::{PolynomialBasis, PolynomialCoefficients, polynomial_basis};

fn invalid(reason: impl Into<String>) -> DeseqError {
    DeseqError::InvalidOptions {
        reason: reason.into(),
    }
}

/// Centering or scaling policy, with one supplied value per matrix column.
#[derive(Clone, Copy, Debug)]
pub enum ScalePolicy<'a> {
    /// Column mean for centering; sqrt(sum of squares / max(1, n-1)) for scaling.
    Automatic,
    /// Omit this operation, retaining its absence in output metadata.
    Disabled,
    /// Explicit centers or divisors. Divisors may be negative, but not zero.
    Supplied(&'a [f64]),
}

/// Scaled values and numerical metadata for interface adapters.
#[derive(Clone, Debug, PartialEq)]
pub struct ScaledCovariates {
    /// Values in row-major sample-by-column order.
    pub values: RowMajorMatrix<f64>,
    /// Applied centers; absent when centering was disabled.
    pub center: Option<Vec<f64>>,
    /// Applied divisors; absent when scaling was disabled.
    pub scale: Option<Vec<f64>>,
}

/// Center/scale finite model variables using the same operation order for
/// borrowed row-major and column-major inputs. Non-finite outputs are rejected
/// because they cannot be passed to fitting as valid design covariates.
pub fn scale_covariates(
    input: MatrixView<'_, f64>,
    center: ScalePolicy<'_>,
    scale: ScalePolicy<'_>,
) -> Result<ScaledCovariates, DeseqError> {
    let rows = input.n_rows();
    let columns = input.n_cols();
    if rows == 0 || columns == 0 || input.as_slice().iter().any(|v| !v.is_finite()) {
        return Err(invalid("scaling requires nonempty finite covariates"));
    }
    for policy in [center, scale] {
        if let ScalePolicy::Supplied(values) = policy
            && (values.len() != columns || values.iter().any(|v| !v.is_finite()))
        {
            return Err(invalid(
                "scale controls require one finite value per column",
            ));
        }
    }
    let mut values = vec![0.0; rows * columns];
    let mut centers =
        (!matches!(center, ScalePolicy::Disabled)).then(|| Vec::with_capacity(columns));
    let mut divisors =
        (!matches!(scale, ScalePolicy::Disabled)).then(|| Vec::with_capacity(columns));
    for column in 0..columns {
        // R's column mean is a single reduction; poly's scalar mean has a
        // correction pass. They are different operations, not frontend modes.
        let center_value = match center {
            ScalePolicy::Automatic => {
                (0..rows).fold(0.0, |sum, row| sum + input.at(row, column)) / rows as f64
            }
            ScalePolicy::Disabled => 0.0,
            ScalePolicy::Supplied(v) => v[column],
        };
        for row in 0..rows {
            values[row * columns + column] = if matches!(center, ScalePolicy::Disabled) {
                *input.at(row, column)
            } else {
                *input.at(row, column) - center_value
            };
        }
        let divisor = match scale {
            ScalePolicy::Automatic => {
                let squares = (0..rows)
                    .map(|r| {
                        let x = values[r * columns + column];
                        x * x
                    })
                    .sum::<f64>();
                (squares / rows.saturating_sub(1).max(1) as f64).sqrt()
            }
            ScalePolicy::Disabled => 1.0,
            ScalePolicy::Supplied(v) => v[column],
        };
        if !center_value.is_finite() || !divisor.is_finite() || divisor == 0.0 {
            return Err(invalid(
                "scaling produced a non-finite center or zero/non-finite divisor",
            ));
        }
        if !matches!(scale, ScalePolicy::Disabled) {
            for row in 0..rows {
                values[row * columns + column] /= divisor;
            }
        }
        if let Some(v) = &mut centers {
            v.push(center_value);
        }
        if let Some(v) = &mut divisors {
            v.push(divisor);
        }
    }
    if values.iter().any(|v| !v.is_finite()) {
        return Err(invalid("scaling produced non-finite model covariates"));
    }
    Ok(ScaledCovariates {
        values: RowMajorMatrix::from_row_major(rows, columns, values)?,
        center: centers,
        scale: divisors,
    })
}

/// R-compatible finite integer coercion represented as a double model column.
/// Truncation is toward zero. The reserved i32 minimum and out-of-range values
/// would become missing in R and are rejected for a fitting covariate.
pub fn integer_covariate(value: f64) -> Result<f64, DeseqError> {
    let integer = value.trunc();
    if !integer.is_finite() || integer <= i32::MIN as f64 || integer > i32::MAX as f64 {
        return Err(invalid(
            "integer formula coercion produced a missing/out-of-range value",
        ));
    }
    Ok((integer as i32) as f64)
}

// Scalar mean used by polynomial construction. Binary64 reduction and the
// correction/overflow fallback are explicit and shared across all frontends.
fn scalar_mean(values: &[f64]) -> f64 {
    let n = values.len() as f64;
    let total = values.iter().fold(0.0, |sum, x| sum + x);
    let ordinary = total.is_finite();
    let mut mean = if ordinary {
        total / n
    } else {
        values.iter().map(|x| x / n).sum::<f64>()
    };
    if mean.is_finite() {
        mean += if ordinary {
            values.iter().map(|x| x - mean).sum::<f64>() / n
        } else {
            values.iter().map(|x| (x - mean) / n).sum::<f64>()
        };
    }
    mean
}
