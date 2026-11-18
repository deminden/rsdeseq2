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

/// Rank after scaling each nonzero column by its maximum absolute entry.
///
/// Scaling is confined to validation; fitting retains the caller's units.
/// The tolerance is relative to each column's original scale.
pub fn matrix_rank_relative(
    values: &[f64],
    n_rows: usize,
    n_cols: usize,
    tolerance: f64,
) -> Result<usize, DeseqError> {
    let expected = n_rows
        .checked_mul(n_cols)
        .ok_or_else(|| DeseqError::InvalidDimensions {
            context: "relative rank matrix shape overflow".into(),
            expected: usize::MAX,
            actual: values.len(),
        })?;
    if values.len() != expected {
        return Err(DeseqError::InvalidDimensions {
            context: "relative rank matrix values".into(),
            expected,
            actual: values.len(),
        });
    }
    let mut scaled = values.to_vec();
    for (index, value) in values.iter().copied().enumerate() {
        if !value.is_finite() {
            return Err(DeseqError::NonFiniteValue {
                context: "relative rank matrix".into(),
                index: Some(index),
                value,
            });
        }
    }
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
    matrix_rank(&scaled, n_rows, n_cols, tolerance)
}

/// R's default `qr()` tolerance (LINPACK `dqrdc2`).
pub const R_QR_TOLERANCE: f64 = 1e-7;

/// Rank of a row-major matrix by R's default `qr(x)$rank`: LINPACK `dqrdc2`
/// Householder QR with limited column pivoting.
///
/// A column is moved to the end, and does not count towards the rank, when
/// the norm of its part not yet explained by earlier columns falls below
/// `tolerance` times its original norm. DESeq2's `checkFullRank` compares
/// this rank, with `tolerance = 1e-7`, to the number of columns. Follows the
/// Fortran in R's `src/appl/dqrdc2.f`, including R's 1999 change that
/// recomputes a column norm once the downdate ratio drops below 1e-6.
pub fn matrix_rank_r_qr(
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
    if let Some((index, value)) = values
        .iter()
        .copied()
        .enumerate()
        .find(|(_, value)| !value.is_finite())
    {
        return Err(DeseqError::NonFiniteValue {
            context: "rank matrix".to_string(),
            index: Some(index),
            value,
        });
    }
    if n_rows == 0 || n_cols == 0 {
        return Ok(0);
    }
    let n = n_rows;
    let p = n_cols;
    // Column-major working copy, as LINPACK stores it.
    let mut x: Vec<Vec<f64>> = (0..p)
        .map(|col| (0..n).map(|row| values[row * p + col]).collect())
        .collect();
    let mut qraux: Vec<f64> = x.iter().map(|col| dnrm2(col)).collect();
    let mut work1 = qraux.clone();
    let mut work2: Vec<f64> = qraux
        .iter()
        .map(|&norm| if norm == 0.0 { 1.0 } else { norm })
        .collect();
    let lup = n.min(p);
    // `k` is one past the last column still in play (LINPACK's k = p + 1).
    let mut k = p + 1;
    for l in 0..lup {
        // Cycle negligible columns to the end until a usable one sits at l.
        while l + 1 < k && qraux[l] < work2[l] * tolerance {
            let column = x.remove(l);
            x.push(column);
            let moved = qraux.remove(l);
            qraux.push(moved);
            let moved = work1.remove(l);
            work1.push(moved);
            let moved = work2.remove(l);
            work2.push(moved);
            k -= 1;
        }
        if l + 1 == n {
            break;
        }
        let nrmxl = dnrm2(&x[l][l..]);
        if nrmxl == 0.0 {
            continue;
        }
        let nrmxl = if x[l][l] != 0.0 {
            nrmxl.copysign(x[l][l])
        } else {
            nrmxl
        };
        for value in &mut x[l][l..] {
            *value *= 1.0 / nrmxl;
        }
        x[l][l] += 1.0;
        let (pivot, rest) = x.split_at_mut(l + 1);
        let householder = &pivot[l][l..];
        for (offset, column) in rest.iter_mut().enumerate() {
            let j = l + 1 + offset;
            let dot = householder
                .iter()
                .zip(&column[l..])
                .map(|(a, b)| a * b)
                .sum::<f64>();
            let t = -dot / householder[0];
            for (value, h) in column[l..].iter_mut().zip(householder) {
                *value += t * h;
            }
            if qraux[j] != 0.0 {
                let ratio = column[l].abs() / qraux[j];
                let tt = (1.0 - ratio * ratio).max(0.0);
                if tt.abs() < 1e-6 {
                    qraux[j] = dnrm2(&column[l + 1..]);
                    work1[j] = qraux[j];
                } else {
                    qraux[j] *= tt.sqrt();
                }
            }
        }
        qraux[l] = x[l][l];
        x[l][l] = -nrmxl;
    }
    Ok((k - 1).min(n))
}

/// Euclidean norm with BLAS `dnrm2`'s scaling, so no intermediate overflows.
fn dnrm2(values: &[f64]) -> f64 {
    let mut scale = 0.0_f64;
    let mut ssq = 1.0_f64;
    for value in values.iter().map(|v| v.abs()).filter(|&v| v != 0.0) {
        if scale < value {
            ssq = 1.0 + ssq * (scale / value) * (scale / value);
            scale = value;
        } else {
            ssq += (value / scale) * (value / scale);
        }
    }
    scale * ssq.sqrt()
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
