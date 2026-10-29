//! Validated immutable matrix storage; layout never chooses arithmetic.
use super::RowMajorMatrix;
use crate::errors::{DeseqError, invalid_dimensions};

/// Storage order of a borrowed dense matrix.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum MatrixOrder {
    /// Consecutive columns within each row.
    RowMajor,
    /// Consecutive rows within each column.
    ColumnMajor,
}

/// Shape-checked dense storage borrowed for a synchronous calculation.
///
/// Zero rows are allowed; columns must be nonzero. No ownership is transferred
/// and no storage is transposed or copied.
#[derive(Clone, Copy, Debug)]
pub struct MatrixView<'a, T> {
    values: &'a [T],
    rows: usize,
    columns: usize,
    order: MatrixOrder,
}
impl<'a, T> MatrixView<'a, T> {
    /// Validate the complete allocated slice against its logical shape.
    pub fn new(
        rows: usize,
        columns: usize,
        values: &'a [T],
        order: MatrixOrder,
    ) -> Result<Self, DeseqError> {
        if columns == 0 {
            return Err(invalid_dimensions("borrowed matrix columns", 1, 0));
        }
        let length = rows.checked_mul(columns).ok_or_else(|| {
            invalid_dimensions("borrowed matrix shape overflow", usize::MAX, values.len())
        })?;
        if length != values.len() {
            return Err(invalid_dimensions(
                "borrowed matrix values",
                length,
                values.len(),
            ));
        }
        Ok(Self {
            values,
            rows,
            columns,
            order,
        })
    }
    /// Logical row count.
    pub fn n_rows(&self) -> usize {
        self.rows
    }
    /// Logical column count.
    pub fn n_cols(&self) -> usize {
        self.columns
    }
    /// Physical storage order.
    pub fn order(&self) -> MatrixOrder {
        self.order
    }
    /// Read-only underlying allocation.
    pub fn as_slice(&self) -> &'a [T] {
        self.values
    }
    /// Read one cell after checking both logical indices.
    pub fn get(&self, row: usize, column: usize) -> Result<&'a T, DeseqError> {
        if row >= self.rows {
            return Err(invalid_dimensions("borrowed matrix row", self.rows, row));
        }
        if column >= self.columns {
            return Err(invalid_dimensions(
                "borrowed matrix column",
                self.columns,
                column,
            ));
        }
        Ok(self.at(row, column))
    }
    #[inline]
    pub(crate) fn at(&self, row: usize, column: usize) -> &'a T {
        debug_assert!(row < self.rows && column < self.columns);
        &self.values[match self.order {
            MatrixOrder::RowMajor => row * self.columns + column,
            MatrixOrder::ColumnMajor => row + column * self.rows,
        }]
    }
}
impl<'a, T> From<&'a RowMajorMatrix<T>> for MatrixView<'a, T> {
    fn from(matrix: &'a RowMajorMatrix<T>) -> Self {
        Self {
            values: matrix.as_slice(),
            rows: matrix.n_rows(),
            columns: matrix.n_cols(),
            order: MatrixOrder::RowMajor,
        }
    }
}
