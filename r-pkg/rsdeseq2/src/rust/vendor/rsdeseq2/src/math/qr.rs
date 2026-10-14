use crate::errors::DeseqError;

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
) -> Result<usize, DeseqError> {
    if !tolerance.is_finite() || tolerance < 0.0 {
        return Err(DeseqError::InvalidOptions {
            reason: "rank tolerance must be finite and non-negative".to_string(),
        });
    }
    let expected = n_rows
        .checked_mul(n_cols)
        .ok_or_else(|| DeseqError::InvalidDimensions {
            context: "rank matrix shape overflow".to_string(),
            expected: usize::MAX,
            actual: values.len(),
        })?;
    if values.len() != expected {
        return Err(DeseqError::InvalidDimensions {
            context: "rank matrix values".to_string(),
            expected,
            actual: values.len(),
        });
    }
    for (index, value) in values.iter().copied().enumerate() {
        if !value.is_finite() {
            return Err(DeseqError::NonFiniteValue {
                context: "rank matrix".to_string(),
                index: Some(index),
                value,
            });
        }
    }
    if n_rows == 0 || n_cols == 0 {
        return Ok(0);
    }

    let mut matrix = values.to_vec();
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
                return Err(DeseqError::NonFiniteValue {
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
                    return Err(DeseqError::NonFiniteValue {
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

#[cfg(test)]
mod tests {
    use super::matrix_rank;

    #[test]
    fn matrix_rank_rejects_nonfinite_elimination_update() {
        let values = [
            1.0,
            f64::MAX,
            1.0,
            -f64::MAX, //
        ];

        let err = matrix_rank(&values, 2, 2, 0.0).unwrap_err();

        assert!(err.to_string().contains("rank elimination update"));
    }
}
