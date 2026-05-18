use crate::errors::{invalid_dimensions, DeseqError};

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
    pub fn from_row_major(
        n_rows: usize,
        n_cols: usize,
        values: Vec<T>,
    ) -> Result<Self, DeseqError> {
        if n_rows == 0 {
            return Err(DeseqError::InvalidDimensions {
                context: "row-major matrix rows".to_string(),
                expected: 1,
                actual: 0,
            });
        }
        if n_cols == 0 {
            return Err(DeseqError::InvalidDimensions {
                context: "row-major matrix columns".to_string(),
                expected: 1,
                actual: 0,
            });
        }
        let expected = n_rows
            .checked_mul(n_cols)
            .ok_or_else(|| DeseqError::InvalidDimensions {
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

    /// Whether the matrix has no stored values.
    pub fn is_empty(&self) -> bool {
        self.values.is_empty()
    }

    /// Matrix values in row-major order.
    pub fn as_slice(&self) -> &[T] {
        &self.values
    }

    /// Consume the matrix and return row-major values.
    pub fn into_values(self) -> Vec<T> {
        self.values
    }

    /// Return a row by index.
    pub fn row(&self, row: usize) -> Result<&[T], DeseqError> {
        if row >= self.n_rows {
            return Err(DeseqError::InvalidDimensions {
                context: "row index".to_string(),
                expected: self.n_rows.saturating_sub(1),
                actual: row,
            });
        }
        let start = row * self.n_cols;
        Ok(&self.values[start..start + self.n_cols])
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
    pub fn from_elem(n_rows: usize, n_cols: usize, value: T) -> Result<Self, DeseqError> {
        let len = n_rows
            .checked_mul(n_cols)
            .ok_or_else(|| DeseqError::InvalidDimensions {
                context: "row-major matrix shape overflow".to_string(),
                expected: usize::MAX,
                actual: 0,
            })?;
        Self::from_row_major(n_rows, n_cols, vec![value; len])
    }
}

impl RowMajorMatrix<f64> {
    /// Validate that all matrix values are finite.
    pub fn validate_finite(&self, context: &str) -> Result<(), DeseqError> {
        for (idx, value) in self.values.iter().copied().enumerate() {
            if !value.is_finite() {
                return Err(DeseqError::NonFiniteValue {
                    context: context.to_string(),
                    index: Some(idx),
                    value,
                });
            }
        }
        Ok(())
    }
}
