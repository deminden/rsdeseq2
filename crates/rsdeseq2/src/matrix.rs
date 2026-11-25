//! Count and numeric matrices, borrowed views and row flags.

use crate::errors::DeseqError;

mod count_matrix;
mod counts;
mod row_flags;
pub use count_matrix::{CountMatrix, CountsSummary};
pub use counts::{CountRows, CountView};
pub use row_flags::RowFlags;

pub(crate) use rcompat_numeric::math::dense::DenseMatrix;
pub use rcompat_numeric::matrix::{MatrixOrder, MatrixView, RowMajorMatrix, borrowed};
pub(crate) use rcompat_numeric::matrix::{normalize_index_range, try_alloc_matrix};

/// Allocate an initialized column-major workspace through the shared matrix allocator.
pub(crate) fn try_zero_matrix(
    rows: usize,
    columns: usize,
    context: &str,
) -> Result<DenseMatrix, DeseqError> {
    DenseMatrix::try_zeros(rows, columns, context)
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

pub(crate) fn try_zero_vector(len: usize, context: &str) -> Result<Vec<f64>, DeseqError> {
    try_filled_vec(len, 0.0, context)
}

/// Copy validated row-major storage into a fallibly allocated column-major matrix.
pub(crate) fn try_dense_from_row_major(
    source: &RowMajorMatrix<f64>,
    context: &str,
) -> Result<DenseMatrix, DeseqError> {
    let columns = source.n_cols();
    DenseMatrix::try_from_fn(source.n_rows(), columns, context, |row, column| {
        source.as_slice()[row * columns + column]
    })
}
