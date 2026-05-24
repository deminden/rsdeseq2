//! Helpers for expanding compact non-all-zero rows back to full gene order.

use crate::errors::{invalid_dimensions, DeseqError};
use crate::matrix::RowMajorMatrix;

/// Return indices for rows that should be fit or copied from compact outputs.
pub(crate) fn nonzero_gene_indices(all_zero: &[bool]) -> Vec<usize> {
    all_zero
        .iter()
        .copied()
        .enumerate()
        .filter_map(|(gene, is_zero)| (!is_zero).then_some(gene))
        .collect()
}

/// Expand a compact row matrix into full gene order, filling all-zero rows with `NaN`.
pub(crate) fn expand_matrix_with_nan_rows(
    matrix: &RowMajorMatrix<f64>,
    all_zero: &[bool],
) -> Result<RowMajorMatrix<f64>, DeseqError> {
    let n_cols = matrix.n_cols();
    let mut values = vec![f64::NAN; all_zero.len() * n_cols];
    let mut compact_row = 0_usize;
    for (row, is_zero) in all_zero.iter().copied().enumerate() {
        if is_zero {
            continue;
        }
        let src = matrix.row(compact_row)?;
        let start = row * n_cols;
        values[start..start + n_cols].copy_from_slice(src);
        compact_row += 1;
    }
    if compact_row != matrix.n_rows() {
        return Err(invalid_dimensions(
            "expanded matrix non-zero rows",
            compact_row,
            matrix.n_rows(),
        ));
    }
    RowMajorMatrix::from_row_major(all_zero.len(), n_cols, values)
}

/// Expand compact per-gene floating-point values, filling all-zero rows with `NaN`.
pub(crate) fn expand_gene_values_with_nan_rows(
    values: &[f64],
    all_zero: &[bool],
) -> Result<Vec<f64>, DeseqError> {
    expand_gene_values_with_fill_rows(values, all_zero, f64::NAN)
}

/// Replace full-length all-zero rows with `NaN`, preserving non-all-zero values.
pub(crate) fn mask_all_zero_values_with_nan_rows(
    values: &[f64],
    all_zero: &[bool],
) -> Result<Vec<f64>, DeseqError> {
    if values.len() != all_zero.len() {
        return Err(invalid_dimensions(
            "all-zero masked vector rows",
            all_zero.len(),
            values.len(),
        ));
    }
    Ok(values
        .iter()
        .copied()
        .zip(all_zero.iter().copied())
        .map(|(value, is_zero)| if is_zero { f64::NAN } else { value })
        .collect())
}

/// Expand compact per-gene values, filling all-zero rows with a caller-supplied value.
pub(crate) fn expand_gene_values_with_fill_rows<T: Clone>(
    values: &[T],
    all_zero: &[bool],
    fill: T,
) -> Result<Vec<T>, DeseqError> {
    let mut expanded = vec![fill; all_zero.len()];
    let mut compact_row = 0_usize;
    for (row, is_zero) in all_zero.iter().copied().enumerate() {
        if is_zero {
            continue;
        }
        let Some(value) = values.get(compact_row) else {
            return Err(invalid_dimensions(
                "expanded vector non-zero rows",
                compact_row + 1,
                values.len(),
            ));
        };
        expanded[row] = value.clone();
        compact_row += 1;
    }
    if compact_row != values.len() {
        return Err(invalid_dimensions(
            "expanded vector non-zero rows",
            compact_row,
            values.len(),
        ));
    }
    Ok(expanded)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn nonzero_indices_skip_all_zero_rows() {
        assert_eq!(
            nonzero_gene_indices(&[true, false, true, false]),
            vec![1, 3]
        );
    }

    #[test]
    fn expand_values_restores_full_gene_order() {
        assert_eq!(
            expand_gene_values_with_fill_rows(&[10, 20], &[true, false, true, false], 0).unwrap(),
            vec![0, 10, 0, 20]
        );
    }

    #[test]
    fn expand_values_validates_compact_row_count() {
        assert!(expand_gene_values_with_fill_rows(&[10], &[false, false], 0).is_err());
        assert!(expand_gene_values_with_fill_rows(&[10, 20], &[false], 0).is_err());
    }

    #[test]
    fn expand_matrix_restores_full_gene_order_with_nan_rows() {
        let compact = RowMajorMatrix::from_row_major(2, 2, vec![1.0, 2.0, 3.0, 4.0]).unwrap();
        let expanded = expand_matrix_with_nan_rows(&compact, &[true, false, true, false]).unwrap();

        assert!(expanded.row(0).unwrap()[0].is_nan());
        assert_eq!(expanded.row(1).unwrap(), &[1.0, 2.0]);
        assert!(expanded.row(2).unwrap()[1].is_nan());
        assert_eq!(expanded.row(3).unwrap(), &[3.0, 4.0]);
    }

    #[test]
    fn mask_full_length_values_marks_all_zero_rows() {
        let masked =
            mask_all_zero_values_with_nan_rows(&[1.0, 2.0, 3.0], &[false, true, false]).unwrap();

        assert_eq!(masked[0], 1.0);
        assert!(masked[1].is_nan());
        assert_eq!(masked[2], 3.0);
        assert!(mask_all_zero_values_with_nan_rows(&[1.0], &[false, true]).is_err());
    }
}
