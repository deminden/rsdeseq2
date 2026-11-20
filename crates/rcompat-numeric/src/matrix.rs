//! Dense row-major matrices and borrowed matrix views.

pub mod borrowed;
pub use borrowed::{MatrixOrder, MatrixView};

use core::ops::{Bound, RangeBounds};

use crate::errors::{Error, invalid_dimensions};

/// A dense row-major matrix.
///
/// The value at row `r`, column `c` is stored at `values[r * n_cols + c]`.
#[derive(Clone, Debug, PartialEq)]
pub struct RowMajorMatrix<T> {
    n_rows: usize,
    n_cols: usize,
    values: Vec<T>,
}

impl<T> RowMajorMatrix<T> {
    /// Create a row-major matrix after validating dimensions.
    pub fn from_row_major(n_rows: usize, n_cols: usize, values: Vec<T>) -> Result<Self, Error> {
        if n_rows == 0 {
            return Err(Error::InvalidDimensions {
                context: "row-major matrix rows".to_string(),
                expected: 1,
                actual: 0,
            });
        }
        if n_cols == 0 {
            return Err(Error::InvalidDimensions {
                context: "row-major matrix columns".to_string(),
                expected: 1,
                actual: 0,
            });
        }
        let expected = n_rows
            .checked_mul(n_cols)
            .ok_or_else(|| Error::InvalidDimensions {
                context: "row-major matrix shape overflow".to_string(),
                expected: usize::MAX,
                actual: values.len(),
            })?;
        if values.len() != expected {
            return Err(invalid_dimensions(
                "row-major matrix values",
                expected,
                values.len(),
            ));
        }
        Ok(Self {
            n_rows,
            n_cols,
            values,
        })
    }

    /// Number of rows.
    pub fn n_rows(&self) -> usize {
        self.n_rows
    }

    /// Number of columns.
    pub fn n_cols(&self) -> usize {
        self.n_cols
    }

    /// Number of stored values.
    pub fn len(&self) -> usize {
        self.values.len()
    }

    /// Reusable row-index span.
    pub fn row_indices(&self) -> core::range::Range<usize> {
        core::range::Range {
            start: 0,
            end: self.n_rows,
        }
    }

    /// Reusable column-index span.
    pub fn col_indices(&self) -> core::range::Range<usize> {
        core::range::Range {
            start: 0,
            end: self.n_cols,
        }
    }

    /// Whether the matrix has no stored values.
    pub fn is_empty(&self) -> bool {
        self.values.is_empty()
    }

    /// Matrix values in row-major order.
    pub fn as_slice(&self) -> &[T] {
        &self.values
    }

    #[doc(hidden)]
    pub fn as_mut_slice(&mut self) -> &mut [T] {
        &mut self.values
    }

    /// Consume the matrix and return row-major values.
    pub fn into_values(self) -> Vec<T> {
        self.values
    }

    /// Return a row by index.
    pub fn row(&self, row: usize) -> Result<&[T], Error> {
        if row >= self.n_rows {
            return Err(Error::InvalidDimensions {
                context: "row index".to_string(),
                expected: self.n_rows.saturating_sub(1),
                actual: row,
            });
        }
        let start = row * self.n_cols;
        Ok(&self.values[start..start + self.n_cols])
    }

    /// Return a contiguous row block.
    ///
    /// The range accepts both legacy range syntax (`1..3`) and the newer
    /// `core::range` types. The returned slice is still row-major and contains
    /// `n_rows_in_range * n_cols` values.
    pub fn rows<R: RangeBounds<usize>>(&self, rows: R) -> Result<&[T], Error> {
        let (start_row, end_row) = normalize_index_range(rows, self.n_rows, "row range")?;
        let start = start_row
            .checked_mul(self.n_cols)
            .ok_or_else(|| Error::InvalidDimensions {
                context: "row range start overflow".to_string(),
                expected: self.len(),
                actual: usize::MAX,
            })?;
        let end = end_row
            .checked_mul(self.n_cols)
            .ok_or_else(|| Error::InvalidDimensions {
                context: "row range end overflow".to_string(),
                expected: self.len(),
                actual: usize::MAX,
            })?;
        Ok(&self.values[start..end])
    }

    /// Return a value by row and column.
    pub fn get(&self, row: usize, col: usize) -> Option<&T> {
        if row >= self.n_rows || col >= self.n_cols {
            return None;
        }
        self.values.get(row * self.n_cols + col)
    }
}

impl<T: Clone> RowMajorMatrix<T> {
    /// Create a matrix filled with one repeated value.
    pub fn from_elem(n_rows: usize, n_cols: usize, value: T) -> Result<Self, Error> {
        let len = n_rows
            .checked_mul(n_cols)
            .ok_or_else(|| Error::InvalidDimensions {
                context: "row-major matrix shape overflow".to_string(),
                expected: usize::MAX,
                actual: 0,
            })?;
        let mut values = try_alloc_matrix(len, 1, "row-major matrix")?;
        values.resize(len, value);
        Self::from_row_major(n_rows, n_cols, values)
    }
}

/// Reserve an empty buffer for `rows × columns` values.
///
/// Overflowing sizes and failed allocations become a
/// [`Error::InvalidDimensions`] naming `context` instead of a panic or an
/// abort, so very large inputs fail cleanly through the R bridge.
#[doc(hidden)]
pub fn try_alloc_matrix<T>(rows: usize, columns: usize, context: &str) -> Result<Vec<T>, Error> {
    let max = isize::MAX as usize / std::mem::size_of::<T>().max(1);
    let len = rows
        .checked_mul(columns)
        .filter(|&len| len <= max)
        .ok_or_else(|| Error::InvalidDimensions {
            context: format!("{context} allocation overflow"),
            expected: max,
            actual: rows.saturating_mul(columns),
        })?;
    let mut values = Vec::new();
    values
        .try_reserve_exact(len)
        .map_err(|_| Error::InvalidDimensions {
            context: format!(
                "{context} allocation failed ({} bytes requested)",
                len * std::mem::size_of::<T>()
            ),
            expected: len,
            actual: 0,
        })?;
    Ok(values)
}

impl RowMajorMatrix<f64> {
    /// Validate that all matrix values are finite.
    pub fn validate_finite(&self, context: &str) -> Result<(), Error> {
        for (idx, value) in self.values.iter().copied().enumerate() {
            if !value.is_finite() {
                return Err(Error::NonFiniteValue {
                    context: context.to_string(),
                    index: Some(idx),
                    value,
                });
            }
        }
        Ok(())
    }
}

#[doc(hidden)]
pub fn normalize_index_range<R: RangeBounds<usize>>(
    range: R,
    len: usize,
    context: &str,
) -> Result<(usize, usize), Error> {
    let start = match range.start_bound() {
        Bound::Included(value) => *value,
        Bound::Excluded(value) => value
            .checked_add(1)
            .ok_or_else(|| invalid_range_bound(context, len, usize::MAX))?,
        Bound::Unbounded => 0,
    };
    let end = match range.end_bound() {
        Bound::Included(value) => value
            .checked_add(1)
            .ok_or_else(|| invalid_range_bound(context, len, usize::MAX))?,
        Bound::Excluded(value) => *value,
        Bound::Unbounded => len,
    };
    if start > end {
        return Err(invalid_range_bound(context, end, start));
    }
    if end > len {
        return Err(invalid_range_bound(context, len, end));
    }
    Ok((start, end))
}

fn invalid_range_bound(context: &str, expected: usize, actual: usize) -> Error {
    Error::InvalidDimensions {
        context: context.to_string(),
        expected,
        actual,
    }
}

#[cfg(test)]
mod tests {
    use super::RowMajorMatrix;
    use crate::errors::Error;
    use std::assert_matches;

    #[test]
    fn matrix_index_spans_are_copy_and_reusable() {
        let matrix = RowMajorMatrix::from_row_major(2, 3, vec![1, 2, 3, 4, 5, 6]).unwrap();
        let rows = matrix.row_indices();
        let cols = matrix.col_indices();

        let first_rows = rows.into_iter().collect::<Vec<_>>();
        let second_rows = rows.into_iter().collect::<Vec<_>>();
        let first_cols = cols.into_iter().collect::<Vec<_>>();
        let second_cols = cols.into_iter().collect::<Vec<_>>();

        assert_eq!(first_rows, vec![0, 1]);
        assert_eq!(second_rows, first_rows);
        assert_eq!(first_cols, vec![0, 1, 2]);
        assert_eq!(second_cols, first_cols);
    }

    #[test]
    fn matrix_rows_accept_legacy_and_new_ranges() {
        let matrix = RowMajorMatrix::from_row_major(3, 2, vec![1, 2, 3, 4, 5, 6]).unwrap();

        assert_eq!(matrix.rows(1..3).unwrap(), &[3, 4, 5, 6]);
        assert_eq!(
            matrix
                .rows(core::range::Range { start: 0, end: 2 })
                .unwrap(),
            &[1, 2, 3, 4]
        );
        assert_eq!(
            matrix.rows(core::range::RangeFrom { start: 2 }).unwrap(),
            &[5, 6]
        );
    }

    #[test]
    fn matrix_rows_reports_range_errors_with_debuggable_match() {
        let matrix = RowMajorMatrix::from_row_major(2, 2, vec![1, 2, 3, 4]).unwrap();
        let start = 2;
        let end = 1;

        assert_matches!(
            matrix.rows(start..end).unwrap_err(),
            Error::InvalidDimensions { .. }
        );
        assert_matches!(
            matrix
                .rows(core::range::Range { start: 0, end: 3 })
                .unwrap_err(),
            Error::InvalidDimensions { .. }
        );
    }
}
