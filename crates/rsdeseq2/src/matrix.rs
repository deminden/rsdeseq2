//! Count and numeric matrices, borrowed views and row flags.

mod count_matrix;
mod counts;
mod row_flags;
pub use count_matrix::{CountMatrix, CountsSummary};
pub use counts::{CountRows, CountView};
pub use row_flags::RowFlags;

pub use rcompat_numeric::matrix::{MatrixOrder, MatrixView, RowMajorMatrix, borrowed};
pub(crate) use rcompat_numeric::matrix::{normalize_index_range, try_alloc_matrix};
