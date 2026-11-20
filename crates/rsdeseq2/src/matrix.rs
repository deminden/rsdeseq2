//! Count and numeric matrices, borrowed views and row flags.

use nalgebra::{DMatrix, DVector};

use crate::errors::DeseqError;

mod count_matrix;
mod counts;
mod row_flags;
pub use count_matrix::{CountMatrix, CountsSummary};
pub use counts::{CountRows, CountView};
pub use row_flags::RowFlags;

pub use rcompat_numeric::matrix::{MatrixOrder, MatrixView, RowMajorMatrix, borrowed};
pub(crate) use rcompat_numeric::matrix::{normalize_index_range, try_alloc_matrix};

/// Allocate an initialized column-major workspace through the shared matrix allocator.
pub(crate) fn try_zero_dmatrix(
    rows: usize,
    columns: usize,
    context: &str,
) -> Result<DMatrix<f64>, DeseqError> {
    let mut values = try_alloc_matrix(rows, columns, context)?;
    values.resize(rows * columns, 0.0);
    Ok(DMatrix::from_vec(rows, columns, values))
}

/// Allocate and fill a vector through the shared matrix allocator.
pub(crate) fn try_filled_vec<T: Clone>(
    len: usize,
    value: T,
    context: &str,
) -> Result<Vec<T>, DeseqError> {
    let mut values = try_alloc_matrix(len, 1, context)?;
    values.resize(len, value);
    Ok(values)
}

pub(crate) fn try_zero_dvector(len: usize, context: &str) -> Result<DVector<f64>, DeseqError> {
    Ok(DVector::from_vec(try_filled_vec(len, 0.0, context)?))
}

/// Copy validated row-major storage into a fallibly allocated column-major matrix.
pub(crate) fn try_dmatrix_from_row_major(
    source: &RowMajorMatrix<f64>,
    context: &str,
) -> Result<DMatrix<f64>, DeseqError> {
    let rows = source.n_rows();
    let columns = source.n_cols();
    let mut values = try_alloc_matrix(rows, columns, context)?;
    for column in 0..columns {
        values.extend(source.as_slice()[column..].iter().step_by(columns).copied());
    }
    Ok(DMatrix::from_vec(rows, columns, values))
}
