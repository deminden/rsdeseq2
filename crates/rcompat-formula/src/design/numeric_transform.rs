//! Numerical formula transforms shared independently of formula parsing.
//!
//! Inputs are already evaluated finite covariates in sample order. These
//! routines do not alter fitting rank thresholds or dispatch to a frontend.
use crate::errors::Error;
use crate::errors::{ErrorKind, missing_prerequisite};
use crate::errors::{invalid_dimensions, invalid_options};
use crate::matrix::{MatrixView, RowMajorMatrix, try_alloc_matrix};
use crate::numeric::profile::{NumericProfile, numeric_profile, r_linux_arm64, r_linux_x86};

pub(in crate::design) mod binary128;

// The x87 extended-precision kernels are the only unsafe code in the core.
// Elsewhere a binary64 stand-in keeps the R-profile branches compiling; that
// profile cannot be selected on those hosts.
#[cfg(target_arch = "x86_64")]
#[allow(unsafe_code)]
mod extended;
#[cfg(not(target_arch = "x86_64"))]
#[path = "numeric_transform/extended_unavailable.rs"]
mod extended;
mod polynomial;
pub(super) mod polynomial_qr;
mod reduction;
pub(super) use polynomial::polynomial_expression;
pub use polynomial::{PolynomialBasis, PolynomialCoefficients, polynomial_basis};

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
) -> Result<ScaledCovariates, Error> {
    let rows = input.n_rows();
    let columns = input.n_cols();
    if rows == 0 || columns == 0 || input.as_slice().iter().any(|v| !v.is_finite()) {
        return Err(missing_prerequisite(
            ErrorKind::Options,
            "scaling",
            "nonempty finite covariates",
        ));
    }
    for policy in [center, scale] {
        if let ScalePolicy::Supplied(values) = policy
            && (values.len() != columns || values.iter().any(|v| !v.is_finite()))
        {
            return Err(missing_prerequisite(
                ErrorKind::Options,
                "scale control",
                "one finite value per column",
            ));
        }
    }
    let output = scale_expression(input, center, scale)?;
    if output.center.iter().flatten().any(|v| !v.is_finite())
        || output
            .scale
            .iter()
            .flatten()
            .any(|v| !v.is_finite() || *v == 0.0)
    {
        return Err(invalid_options(
            "scaling produced a non-finite center or zero/non-finite divisor",
        ));
    }
    if output.values.as_slice().iter().any(|v| !v.is_finite()) {
        return Err(invalid_options(
            "scaling produced non-finite model covariates",
        ));
    }
    Ok(output)
}

// Expression evaluation may produce missing/non-finite intermediates that a
// later built-in replaces. The finite fitting API above validates the final
// result; both APIs use these same column reductions and element operations.
pub(super) fn scale_expression(
    input: MatrixView<'_, f64>,
    center: ScalePolicy<'_>,
    scale: ScalePolicy<'_>,
) -> Result<ScaledCovariates, Error> {
    let rows = input.n_rows();
    let columns = input.n_cols();
    for policy in [center, scale] {
        if let ScalePolicy::Supplied(values) = policy
            && values.len() != columns
        {
            return Err(missing_prerequisite(
                ErrorKind::Options,
                "scale control",
                "one value per column",
            ));
        }
    }
    if rows == 0 {
        return Err(invalid_dimensions("row-major matrix rows", 1, 0));
    }
    let mut values = try_alloc_matrix(rows, columns, "scaled covariates")?;
    values.resize(rows * columns, 0.0);
    let mut centers = (!matches!(center, ScalePolicy::Disabled))
        .then(|| try_alloc_matrix(columns, 1, "scale centers"))
        .transpose()?;
    let mut divisors = (!matches!(scale, ScalePolicy::Disabled))
        .then(|| try_alloc_matrix(columns, 1, "scale divisors"))
        .transpose()?;
    let portable = numeric_profile() == NumericProfile::Portable;
    for column in 0..columns {
        // R's column mean is a single reduction; poly's scalar mean has a
        // correction pass. They are different operations, not frontend modes.
        let center_value = match center {
            ScalePolicy::Automatic => {
                let present = (0..rows)
                    .map(|row| *input.at(row, column))
                    .filter(|v| !v.is_nan());
                if r_linux_arm64() {
                    let count = present.clone().count();
                    binary128::column_mean(present, count)
                } else if r_linux_x86() {
                    let count = present.clone().count();
                    let mut total = extended::sum(present);
                    total.divide(count as f64);
                    total.value()
                } else {
                    let (total, count) = present.fold((0.0, 0usize), |(sum, count), value| {
                        (sum + value, count + 1)
                    });
                    total / count as f64
                }
            }
            ScalePolicy::Disabled => 0.0,
            ScalePolicy::Supplied(v) => v[column],
        };
        let center_value =
            if portable && matches!(center, ScalePolicy::Automatic) && center_value.is_nan() {
                f64::NAN
            } else {
                center_value
            };
        for row in 0..rows {
            let original = *input.at(row, column);
            values[row * columns + column] = if matches!(center, ScalePolicy::Disabled) {
                original
            } else {
                let value = original - center_value;
                if portable && value.is_nan() {
                    first_nan([original, center_value])
                } else {
                    value
                }
            };
        }
        let divisor = match scale {
            ScalePolicy::Automatic => {
                let present = (0..rows)
                    .map(|r| values[r * columns + column])
                    .filter(|v| !v.is_nan());
                if r_linux_arm64() || r_linux_x86() {
                    (column_sum(present.clone().map(|x| x * x))
                        / present.count().saturating_sub(1).max(1) as f64)
                        .sqrt()
                } else {
                    let (squares, count) =
                        present.fold((0.0, 0usize), |(sum, count), x| (sum + x * x, count + 1));
                    (squares / count.saturating_sub(1).max(1) as f64).sqrt()
                }
            }
            ScalePolicy::Disabled => 1.0,
            ScalePolicy::Supplied(v) => v[column],
        };
        if !matches!(scale, ScalePolicy::Disabled) {
            for row in 0..rows {
                let entry = &mut values[row * columns + column];
                let value = *entry / divisor;
                *entry = if portable && value.is_nan() {
                    first_nan([*entry, divisor])
                } else {
                    value
                };
            }
        }
        if let Some(v) = &mut centers {
            v.push(center_value);
        }
        if let Some(v) = &mut divisors {
            v.push(divisor);
        }
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
pub fn integer_covariate(value: f64) -> Result<f64, Error> {
    let integer = value.trunc();
    if !integer.is_finite() || integer <= i32::MIN as f64 || integer > i32::MAX as f64 {
        return Err(invalid_options(
            "integer formula coercion produced a missing/out-of-range value",
        ));
    }
    Ok((integer as i32) as f64)
}

/// R's `mean()` for doubles under the active numeric profile (two-pass,
/// extended precision under the Linux R profiles).
#[doc(hidden)]
pub fn scalar_mean(values: &[f64]) -> f64 {
    if r_linux_arm64() {
        binary128::mean(values)
    } else if r_linux_x86() {
        extended::mean(values)
    } else {
        let mean = scalar_mean_binary64(values);
        if mean.is_nan() && numeric_profile() == NumericProfile::Portable {
            first_nan(values.iter().copied())
        } else {
            mean
        }
    }
}

fn first_nan(values: impl IntoIterator<Item = f64>) -> f64 {
    values
        .into_iter()
        .find(|value| value.is_nan())
        .map_or(f64::NAN, |value| {
            f64::from_bits(value.to_bits() | (1 << 51))
        })
}

fn scalar_mean_binary64(values: &[f64]) -> f64 {
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

pub(super) fn column_sum(values: impl Iterator<Item = f64>) -> f64 {
    if r_linux_arm64() {
        binary128::sum_as_double(values)
    } else if r_linux_x86() {
        extended::sum(values).value()
    } else {
        values.sum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::matrix::MatrixOrder;

    #[test]
    fn scaling_rejects_empty_views_before_allocating_column_metadata() {
        let empty = MatrixView::new(0, usize::MAX, &[], MatrixOrder::ColumnMajor).unwrap();
        let error =
            scale_expression(empty, ScalePolicy::Automatic, ScalePolicy::Automatic).unwrap_err();
        std::assert_matches!(
            error,
            Error::InvalidDimensions { context, expected: 1, actual: 0 }
                if context == "row-major matrix rows"
        );
    }
}
