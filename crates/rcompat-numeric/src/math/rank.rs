use crate::errors::Error;
use crate::errors::{ErrorKind, Requirement, invalid_parameter};
use crate::matrix::try_alloc_matrix;

/// Default absolute tolerance for deterministic matrix-rank checks.
pub const DEFAULT_RANK_TOLERANCE: f64 = 1e-10;

/// Compute a deterministic row-major matrix rank with partial-pivot elimination.
///
/// This is intentionally small and inspectable. It is used for DESeq2-style
/// full-rank checks before GLM fitting and for observation-weight rank checks.
pub fn matrix_rank(
    values: &[f64],
    n_rows: usize,
    n_cols: usize,
    tolerance: f64,
) -> Result<usize, Error> {
    validate_rank_input(values, n_rows, n_cols, tolerance)?;
    let mut matrix = try_alloc_matrix(n_rows, n_cols, "rank matrix")?;
    matrix.extend_from_slice(values);
    eliminate_rank(&mut matrix, n_rows, n_cols, tolerance)
}

/// Apply the same rank calculation to a caller-owned work matrix.
#[doc(hidden)]
pub fn matrix_rank_in_place(
    matrix: &mut [f64],
    n_rows: usize,
    n_cols: usize,
    tolerance: f64,
) -> Result<usize, Error> {
    validate_rank_input(matrix, n_rows, n_cols, tolerance)?;
    eliminate_rank(matrix, n_rows, n_cols, tolerance)
}

fn validate_rank_input(
    values: &[f64],
    n_rows: usize,
    n_cols: usize,
    tolerance: f64,
) -> Result<(), Error> {
    if !tolerance.is_finite() || tolerance < 0.0 {
        return Err(invalid_parameter(
            ErrorKind::Options,
            "rank tolerance",
            Requirement::FiniteAndNonNegative,
        ));
    }
    let expected = n_rows
        .checked_mul(n_cols)
        .ok_or_else(|| Error::InvalidDimensions {
            context: "rank matrix shape overflow".to_string(),
            expected: usize::MAX,
            actual: values.len(),
        })?;
    if values.len() != expected {
        return Err(Error::InvalidDimensions {
            context: "rank matrix values".to_string(),
            expected,
            actual: values.len(),
        });
    }
    for (index, value) in values.iter().copied().enumerate() {
        if !value.is_finite() {
            return Err(Error::NonFiniteValue {
                context: "rank matrix".to_string(),
                index: Some(index),
                value,
            });
        }
    }
    Ok(())
}

fn eliminate_rank(
    matrix: &mut [f64],
    n_rows: usize,
    n_cols: usize,
    tolerance: f64,
) -> Result<usize, Error> {
    if n_rows == 0 || n_cols == 0 {
        return Ok(0);
    }
    let mut rank = 0_usize;
    let mut row = 0_usize;

    for col in 0..n_cols {
        let pivot = (row..n_rows)
            .max_by(|left, right| {
                matrix[*left * n_cols + col]
                    .abs()
                    .total_cmp(&matrix[*right * n_cols + col].abs())
            })
            .unwrap_or(row);
        let pivot_value = matrix[pivot * n_cols + col];
        if pivot_value.abs() <= tolerance {
            continue;
        }
        if pivot != row {
            for c in col..n_cols {
                matrix.swap(row * n_cols + c, pivot * n_cols + c);
            }
        }
        let pivot_value = matrix[row * n_cols + col];
        for r in row + 1..n_rows {
            let factor = matrix[r * n_cols + col] / pivot_value;
            if !factor.is_finite() {
                return Err(Error::NonFiniteValue {
                    context: "rank elimination factor".to_string(),
                    index: Some(r * n_cols + col),
                    value: factor,
                });
            }
            if factor == 0.0 {
                continue;
            }
            matrix[r * n_cols + col] = 0.0;
            for c in col + 1..n_cols {
                let source = matrix[row * n_cols + c];
                let target = matrix[r * n_cols + c];
                let product = factor * source;
                let updated = target - product;
                if !product.is_finite() || !updated.is_finite() {
                    return Err(Error::NonFiniteValue {
                        context: "rank elimination update".to_string(),
                        index: Some(r * n_cols + c),
                        value: updated,
                    });
                }
                matrix[r * n_cols + c] = updated;
            }
        }
        rank += 1;
        row += 1;
        if row == n_rows {
            break;
        }
    }
    Ok(rank)
}

/// Rank after scaling each nonzero column by its maximum absolute entry.
///
/// Scaling is confined to validation; fitting retains the caller's units.
/// The tolerance is relative to each column's original scale.
pub fn matrix_rank_relative(
    values: &[f64],
    n_rows: usize,
    n_cols: usize,
    tolerance: f64,
) -> Result<usize, Error> {
    let expected = n_rows
        .checked_mul(n_cols)
        .ok_or_else(|| Error::InvalidDimensions {
            context: "relative rank matrix shape overflow".into(),
            expected: usize::MAX,
            actual: values.len(),
        })?;
    if values.len() != expected {
        return Err(Error::InvalidDimensions {
            context: "relative rank matrix values".into(),
            expected,
            actual: values.len(),
        });
    }
    let mut scaled = try_alloc_matrix(n_rows, n_cols, "relative rank matrix")?;
    for (index, value) in values.iter().copied().enumerate() {
        if !value.is_finite() {
            return Err(Error::NonFiniteValue {
                context: "relative rank matrix".into(),
                index: Some(index),
                value,
            });
        }
    }
    scaled.extend_from_slice(values);
    for col in 0..n_cols {
        let scale = (0..n_rows)
            .map(|row| values[row * n_cols + col].abs())
            .fold(0.0_f64, f64::max);
        if scale > 0.0 {
            for row in 0..n_rows {
                scaled[row * n_cols + col] /= scale;
            }
        }
    }
    matrix_rank_in_place(&mut scaled, n_rows, n_cols, tolerance)
}

/// R's default `qr()` tolerance.
pub const R_QR_TOLERANCE: f64 = 1e-7;

/// Rank of a row-major matrix as R's default `qr(x)$rank` reports it.
///
/// Not available in this revision: the previous implementation was removed
/// and is being replaced by one written from a specification.
pub fn matrix_rank_r_qr(
    values: &[f64],
    n_rows: usize,
    n_cols: usize,
    tolerance: f64,
) -> Result<usize, Error> {
    let _ = (values, n_rows, n_cols, tolerance);
    Err(Error::UnsupportedFeature {
        feature: "R qr() rank (being rewritten)".to_string(),
    })
}

#[cfg(test)]
mod tests {
    use super::matrix_rank;
    use crate::errors::Error;

    #[test]
    fn matrix_rank_rejects_nonfinite_elimination_update() {
        let values = [
            1.0,
            f64::MAX,
            1.0,
            -f64::MAX, //
        ];

        let err = matrix_rank(&values, 2, 2, 0.0).unwrap_err();

        std::assert_matches!(
            err,
            Error::NonFiniteValue { context, index: Some(3), value: f64::INFINITY }
                if context == "rank elimination update"
        );
    }
}
