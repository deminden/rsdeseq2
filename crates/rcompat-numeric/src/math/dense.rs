//! Small dense linear algebra on column-major `f64` matrices.
//!
//! [`DenseMatrix`] stores its entries column by column, so a column is a
//! contiguous slice. The kernels follow Golub and Van Loan, *Matrix
//! Computations* (4th ed., 2013), cited as GVL with the section of each
//! algorithm:
//!
//! - [`matmul_to`], [`gemv_to`], [`tr_gemv_to`]: column (gaxpy) and dot
//!   product forms of matrix products (GVL 1.1);
//! - [`triangular`]: forward and back substitution (GVL 3.1);
//! - [`lu`]: Gaussian elimination with partial pivoting, determinants and
//!   general inverses (GVL 3.2 and 3.4);
//! - [`inverse`]: closed-form inverses of matrices up to 4 by 4 (adjugate
//!   divided by the determinant) and the LU inverse beyond;
//! - [`cholesky`]: the gaxpy Cholesky factorization (GVL 4.2);
//! - [`householder`]: Householder reflections and QR with column pivoting
//!   (GVL 5.1, 5.2 and 5.4);
//! - [`eigen`]: the cyclic Jacobi method for symmetric matrices (GVL 8.5);
//! - [`svd`]: the one-sided Jacobi singular value decomposition (GVL 8.6).
//!
//! The kernels work in place on caller-owned storage wherever the callers run
//! them in a loop, so repeated solves do not allocate.

use core::ops::{Index, IndexMut};

use crate::errors::Error;
use crate::matrix::try_alloc_matrix;

/// Cholesky factorization of symmetric positive definite matrices.
pub mod cholesky;
/// Symmetric eigen-decomposition by the cyclic Jacobi method.
pub mod eigen;
/// Householder reflections and column-pivoted Householder QR.
pub mod householder;
/// Closed-form and LU inverses and determinants.
pub mod inverse;
/// LU factorization with partial pivoting.
pub mod lu;
/// Singular value decomposition by one-sided Jacobi rotations.
pub mod svd;
/// Forward and back substitution.
pub mod triangular;

pub use cholesky::Cholesky;
pub use eigen::SymmetricEigen;
pub use householder::ColPivQr;
pub use inverse::{determinant, try_inverse_in_place};
pub use lu::{Lu, lu_solve_in_place};
pub use svd::Svd;

/// A dense matrix of `f64` stored in column-major order.
///
/// The default value is the empty 0 by 0 matrix.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct DenseMatrix {
    nrows: usize,
    ncols: usize,
    data: Vec<f64>,
}

fn shape_error(context: &str, expected: (usize, usize), actual: (usize, usize)) -> Error {
    Error::MatrixShapeMismatch {
        context: context.to_string(),
        expected,
        actual,
    }
}

impl DenseMatrix {
    /// A `nrows` by `ncols` matrix of zeros, allocated fallibly; `context`
    /// names the matrix in an allocation error.
    pub fn try_zeros(nrows: usize, ncols: usize, context: &str) -> Result<Self, Error> {
        Self::try_from_element(nrows, ncols, 0.0, context)
    }

    /// A `nrows` by `ncols` matrix of zeros.
    pub fn zeros(nrows: usize, ncols: usize) -> Result<Self, Error> {
        Self::try_zeros(nrows, ncols, "dense matrix")
    }

    /// A `nrows` by `ncols` matrix with every entry equal to `value`.
    pub fn from_element(nrows: usize, ncols: usize, value: f64) -> Result<Self, Error> {
        Self::try_from_element(nrows, ncols, value, "dense matrix")
    }

    fn try_from_element(
        nrows: usize,
        ncols: usize,
        value: f64,
        context: &str,
    ) -> Result<Self, Error> {
        let mut data = try_alloc_matrix(nrows, ncols, context)?;
        data.resize(nrows * ncols, value);
        Ok(Self { nrows, ncols, data })
    }

    /// The `n` by `n` identity matrix.
    pub fn identity(n: usize) -> Result<Self, Error> {
        let mut matrix = Self::zeros(n, n)?;
        matrix.fill_with_identity();
        Ok(matrix)
    }

    /// A matrix whose entry `(i, j)` is `entry(i, j)`, evaluated column by
    /// column; `context` names the matrix in an allocation error.
    pub fn try_from_fn(
        nrows: usize,
        ncols: usize,
        context: &str,
        mut entry: impl FnMut(usize, usize) -> f64,
    ) -> Result<Self, Error> {
        let mut data = try_alloc_matrix(nrows, ncols, context)?;
        for j in 0..ncols {
            data.extend((0..nrows).map(|i| entry(i, j)));
        }
        Ok(Self { nrows, ncols, data })
    }

    /// Wrap `values`, stored column by column, as a `nrows` by `ncols` matrix.
    pub fn from_column_major(nrows: usize, ncols: usize, values: Vec<f64>) -> Result<Self, Error> {
        if nrows.checked_mul(ncols) != Some(values.len()) {
            return Err(Error::InvalidDimensions {
                context: "dense matrix values".to_string(),
                expected: nrows.saturating_mul(ncols),
                actual: values.len(),
            });
        }
        Ok(Self {
            nrows,
            ncols,
            data: values,
        })
    }

    /// Copy `values`, stored row by row, into a `nrows` by `ncols` matrix.
    pub fn from_row_major(nrows: usize, ncols: usize, values: &[f64]) -> Result<Self, Error> {
        if nrows.checked_mul(ncols) != Some(values.len()) {
            return Err(Error::InvalidDimensions {
                context: "dense matrix values".to_string(),
                expected: nrows.saturating_mul(ncols),
                actual: values.len(),
            });
        }
        Self::try_from_fn(nrows, ncols, "dense matrix", |i, j| values[i * ncols + j])
    }

    /// Number of rows.
    pub fn nrows(&self) -> usize {
        self.nrows
    }

    /// Number of columns.
    pub fn ncols(&self) -> usize {
        self.ncols
    }

    /// `(rows, columns)`.
    pub fn shape(&self) -> (usize, usize) {
        (self.nrows, self.ncols)
    }

    /// Whether the matrix has no entries.
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    /// The entries in column-major order.
    pub fn as_slice(&self) -> &[f64] {
        &self.data
    }

    /// The entries in column-major order, mutably.
    pub fn as_mut_slice(&mut self) -> &mut [f64] {
        &mut self.data
    }

    /// The column-major storage.
    pub fn into_vec(self) -> Vec<f64> {
        self.data
    }

    /// Iterate over the entries in column-major order.
    pub fn iter(&self) -> core::slice::Iter<'_, f64> {
        self.data.iter()
    }

    /// Column `j` as a slice. Panics if `j` is out of range.
    pub fn column(&self, j: usize) -> &[f64] {
        &self.data[j * self.nrows..(j + 1) * self.nrows]
    }

    /// Column `j` as a mutable slice. Panics if `j` is out of range.
    pub fn column_mut(&mut self, j: usize) -> &mut [f64] {
        &mut self.data[j * self.nrows..(j + 1) * self.nrows]
    }

    /// Entry `(i, j)`, or `None` outside the matrix.
    pub fn get(&self, i: usize, j: usize) -> Option<f64> {
        (i < self.nrows && j < self.ncols).then(|| self.data[j * self.nrows + i])
    }

    /// Set every entry to `value`.
    pub fn fill(&mut self, value: f64) {
        self.data.fill(value);
    }

    /// Set the diagonal to one and every other entry to zero.
    pub fn fill_with_identity(&mut self) {
        self.data.fill(0.0);
        for i in 0..self.nrows.min(self.ncols) {
            self.data[i * self.nrows + i] = 1.0;
        }
    }

    /// Make `self` an entry-wise copy of `source`, taking its shape.
    ///
    /// Storage is reused when it is large enough.
    pub fn copy_from(&mut self, source: &Self) {
        self.nrows = source.nrows;
        self.ncols = source.ncols;
        self.data.clear();
        self.data.extend_from_slice(&source.data);
    }

    /// Swap rows `a` and `b` in every column.
    pub fn swap_rows(&mut self, a: usize, b: usize) {
        if a == b {
            return;
        }
        for column in self.data.chunks_exact_mut(self.nrows.max(1)) {
            column.swap(a, b);
        }
    }

    /// Sum of the diagonal entries, added in order.
    pub fn trace(&self) -> f64 {
        (0..self.nrows.min(self.ncols)).map(|i| self[(i, i)]).sum()
    }

    /// The transpose as a new matrix.
    pub fn transpose(&self) -> Result<Self, Error> {
        Self::try_from_fn(self.ncols, self.nrows, "dense matrix transpose", |i, j| {
            self[(j, i)]
        })
    }

    /// Write the transpose into `out`, which must be `ncols` by `nrows`.
    pub fn transpose_to(&self, out: &mut Self) -> Result<(), Error> {
        if out.shape() != (self.ncols, self.nrows) {
            return Err(shape_error(
                "dense matrix transpose",
                (self.ncols, self.nrows),
                out.shape(),
            ));
        }
        for j in 0..self.ncols {
            for (i, &value) in self.column(j).iter().enumerate() {
                out[(j, i)] = value;
            }
        }
        Ok(())
    }
}

impl Index<(usize, usize)> for DenseMatrix {
    type Output = f64;

    fn index(&self, (i, j): (usize, usize)) -> &f64 {
        debug_assert!(i < self.nrows && j < self.ncols);
        &self.data[j * self.nrows + i]
    }
}

impl IndexMut<(usize, usize)> for DenseMatrix {
    fn index_mut(&mut self, (i, j): (usize, usize)) -> &mut f64 {
        debug_assert!(i < self.nrows && j < self.ncols);
        &mut self.data[j * self.nrows + i]
    }
}

/// Smallest dimension from which [`matmul_to`] accumulates with fused
/// multiply-adds.
const BLOCKED_PRODUCT_DIMENSION: usize = 6;
/// Number of products [`matmul_to`] fuses into one partial sum.
const PRODUCT_BLOCK: usize = 256;

/// `out = a * b` for `a` of shape `m × k` and `b` of shape `k × n`
/// (GVL 1.1.11, the gaxpy form).
///
/// When any of `m`, `k` and `n` is at most 5, column `j` of the product
/// starts as `b[(0, j)] * a[.., 0]` and accumulates `b[(l, j)] * a[.., l]`
/// for `l = 1, 2, ...` in order, so every entry is the left-to-right sum of
/// its products, each rounded separately. When all three exceed 5, every entry is accumulated in
/// blocks of 256 consecutive products: each block is a left-to-right chain
/// of fused multiply-adds starting from zero, and the block sums are added
/// to the entry in order. `out` must be `m × n`.
pub fn matmul_to(a: &DenseMatrix, b: &DenseMatrix, out: &mut DenseMatrix) -> Result<(), Error> {
    if a.ncols != b.nrows {
        return Err(shape_error(
            "matrix product inner dimension",
            (a.nrows, a.ncols),
            b.shape(),
        ));
    }
    if out.shape() != (a.nrows, b.ncols) {
        return Err(shape_error(
            "matrix product output",
            (a.nrows, b.ncols),
            out.shape(),
        ));
    }
    let blocked = a.nrows >= BLOCKED_PRODUCT_DIMENSION
        && a.ncols >= BLOCKED_PRODUCT_DIMENSION
        && b.ncols >= BLOCKED_PRODUCT_DIMENSION;
    for j in 0..b.ncols {
        let column = &mut out.data[j * a.nrows..(j + 1) * a.nrows];
        let factors = &b.data[j * b.nrows..(j + 1) * b.nrows];
        if blocked {
            for (i, target) in column.iter_mut().enumerate() {
                let mut entry = 0.0;
                for (block, chunk) in factors.chunks(PRODUCT_BLOCK).enumerate() {
                    let start = block * PRODUCT_BLOCK;
                    let partial = chunk.iter().enumerate().fold(0.0_f64, |sum, (l, &factor)| {
                        a.data[(start + l) * a.nrows + i].mul_add(factor, sum)
                    });
                    entry = if block == 0 { partial } else { entry + partial };
                }
                *target = entry;
            }
        } else {
            gaxpy(a, factors, column);
        }
    }
    Ok(())
}

/// The product `a * b` as a new matrix (see [`matmul_to`]).
pub fn matmul(a: &DenseMatrix, b: &DenseMatrix) -> Result<DenseMatrix, Error> {
    let mut out = DenseMatrix::try_zeros(a.nrows, b.ncols, "matrix product")?;
    matmul_to(a, b, &mut out)?;
    Ok(out)
}

/// `y = a * x` (GVL 1.1.4, column-oriented gaxpy).
///
/// `y` starts as `x[0] * a[.., 0]` and then receives `x[j] * a[.., j]` for
/// `j = 1, 2, ...` in order. `x` must have `a.ncols()` entries and `y` `a.nrows()`.
pub fn gemv_to(a: &DenseMatrix, x: &[f64], y: &mut [f64]) -> Result<(), Error> {
    if x.len() != a.ncols || y.len() != a.nrows {
        return Err(shape_error(
            "matrix-vector product",
            (a.nrows, a.ncols),
            (y.len(), x.len()),
        ));
    }
    gaxpy(a, x, y);
    Ok(())
}

/// `y = a x` column by column: `y` starts as `x[0] * a[.., 0]` (zero when
/// `a` has no columns) and receives `x[j] * a[.., j]` for `j = 1, 2, ...`.
fn gaxpy(a: &DenseMatrix, x: &[f64], y: &mut [f64]) {
    let Some((&first, rest)) = x.split_first() else {
        y.fill(0.0);
        return;
    };
    for (target, &value) in y.iter_mut().zip(a.column(0)) {
        *target = value * first;
    }
    for (j, &factor) in rest.iter().enumerate() {
        for (target, &value) in y.iter_mut().zip(a.column(j + 1)) {
            *target += value * factor;
        }
    }
}

/// The product `a * x` as a new vector (see [`gemv_to`]).
pub fn gemv(a: &DenseMatrix, x: &[f64]) -> Result<Vec<f64>, Error> {
    let mut y = try_alloc_matrix(a.nrows, 1, "matrix-vector product")?;
    y.resize(a.nrows, 0.0);
    gemv_to(a, x, &mut y)?;
    Ok(y)
}

/// `y = aᵀ * x` (GVL 1.1.3, dot product form): `y[j]` is the left-to-right
/// sum of the products of column `j` with `x`, each rounded separately (the
/// entries of the column gaxpy [`gemv_to`] applied to the transpose).
pub fn tr_gemv_to(a: &DenseMatrix, x: &[f64], y: &mut [f64]) -> Result<(), Error> {
    if x.len() != a.nrows || y.len() != a.ncols {
        return Err(shape_error(
            "transposed matrix-vector product",
            (a.ncols, a.nrows),
            (y.len(), x.len()),
        ));
    }
    for (j, target) in y.iter_mut().enumerate() {
        let column = a.column(j);
        *target = match column.split_first() {
            Some((&head, tail)) => tail
                .iter()
                .zip(&x[1..])
                .fold(head * x[0], |sum, (a, b)| sum + a * b),
            None => 0.0,
        };
    }
    Ok(())
}

/// Number of interleaved accumulators in the blocked part of [`dot`].
const DOT_LANES: usize = 8;

/// Dot product of two slices (only the common prefix is used if the lengths
/// differ), every product and sum rounded separately.
///
/// The longest prefix whose length is a multiple of 8 is accumulated into
/// eight interleaved accumulators (position `p` goes to accumulator
/// `p mod 8`). They are combined as `(s0 + s4) + (s1 + s5) + (s2 + s6) +
/// (s3 + s7)`, left to right, and the remaining tail is added in order.
/// Fewer than 8 entries give the left-to-right sum.
pub fn dot(x: &[f64], y: &[f64]) -> f64 {
    let n = x.len().min(y.len());
    let blocked = DOT_LANES * (n / DOT_LANES);
    let mut lanes = [0.0_f64; DOT_LANES];
    for (p, (&a, &b)) in x[..blocked].iter().zip(&y[..blocked]).enumerate() {
        lanes[p % DOT_LANES] += a * b;
    }
    let s = lanes;
    let mut sum = (s[0] + s[4]) + (s[1] + s[5]) + (s[2] + s[6]) + (s[3] + s[7]);
    for (&a, &b) in x[blocked..n].iter().zip(&y[blocked..n]) {
        sum += a * b;
    }
    sum
}

#[cfg(test)]
pub(crate) mod test_support {
    use super::DenseMatrix;

    /// A small deterministic generator of uniform values in `[-1, 1)`.
    pub(crate) struct Lcg(pub(crate) u64);

    impl Lcg {
        pub(crate) fn next(&mut self) -> f64 {
            self.0 = self
                .0
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            ((self.0 >> 11) as f64 / (1_u64 << 53) as f64) * 2.0 - 1.0
        }

        pub(crate) fn matrix(&mut self, rows: usize, cols: usize) -> DenseMatrix {
            DenseMatrix::try_from_fn(rows, cols, "test", |_, _| self.next()).unwrap()
        }

        /// A diagonally dominant, hence well-conditioned, square matrix.
        pub(crate) fn well_conditioned(&mut self, n: usize) -> DenseMatrix {
            let mut a = self.matrix(n, n);
            for i in 0..n {
                a[(i, i)] += n as f64 + 1.0;
            }
            a
        }

        /// A symmetric positive definite matrix `bᵀb + n I`.
        pub(crate) fn spd(&mut self, n: usize) -> DenseMatrix {
            let b = self.matrix(n, n);
            let mut a = super::matmul(&b.transpose().unwrap(), &b).unwrap();
            for i in 0..n {
                a[(i, i)] += n as f64;
            }
            a
        }
    }

    pub(crate) fn max_abs_diff(a: &[f64], b: &[f64]) -> f64 {
        assert_eq!(a.len(), b.len());
        a.iter()
            .zip(b)
            .map(|(x, y)| (x - y).abs())
            .fold(0.0, f64::max)
    }
}

#[cfg(test)]
mod tests {
    use super::test_support::{Lcg, max_abs_diff};
    use super::*;

    #[test]
    fn construction_and_layout() {
        let a = DenseMatrix::from_row_major(2, 3, &[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]).unwrap();
        assert_eq!(a.shape(), (2, 3));
        assert_eq!(a.as_slice(), &[1.0, 4.0, 2.0, 5.0, 3.0, 6.0]);
        assert_eq!(a.column(1), &[2.0, 5.0]);
        assert_eq!(a[(1, 2)], 6.0);
        assert_eq!(a.get(2, 0), None);
        let t = a.transpose().unwrap();
        assert_eq!(t.shape(), (3, 2));
        assert_eq!(t[(2, 1)], 6.0);
        let mut out = DenseMatrix::zeros(3, 2).unwrap();
        a.transpose_to(&mut out).unwrap();
        assert_eq!(out, t);
        assert!(
            a.transpose_to(&mut DenseMatrix::zeros(2, 3).unwrap())
                .is_err()
        );
        assert!(DenseMatrix::from_column_major(2, 2, vec![1.0; 3]).is_err());
        assert!(DenseMatrix::zeros(usize::MAX, 2).is_err());
        let mut i = DenseMatrix::from_element(3, 3, 7.0).unwrap();
        i.fill_with_identity();
        assert_eq!(i, DenseMatrix::identity(3).unwrap());
        assert_eq!(i.trace(), 3.0);
        let mut s = a.clone();
        s.swap_rows(0, 1);
        assert_eq!(s.as_slice(), &[4.0, 1.0, 5.0, 2.0, 6.0, 3.0]);
        let mut c = DenseMatrix::default();
        c.copy_from(&a);
        assert_eq!(c, a);
    }

    #[test]
    fn dot_uses_eight_lanes() {
        let x: Vec<f64> = (0..19).map(|i| 1.0 + f64::from(i) * 1e-3).collect();
        let y: Vec<f64> = (0..19).map(|i| 1.0 - f64::from(i) * 7e-4).collect();
        let mut lanes = [0.0; 8];
        for p in 0..16 {
            lanes[p % 8] += x[p] * y[p];
        }
        let s = lanes;
        let expected = (s[0] + s[4])
            + (s[1] + s[5])
            + (s[2] + s[6])
            + (s[3] + s[7])
            + x[16] * y[16]
            + x[17] * y[17]
            + x[18] * y[18];
        assert_eq!(dot(&x, &y).to_bits(), expected.to_bits());
        assert_eq!(
            dot(&x[..3], &y[..5]),
            x[0] * y[0] + x[1] * y[1] + x[2] * y[2]
        );
    }

    #[test]
    fn products_match_their_definitions() {
        let mut rng = Lcg(7);
        let a = rng.matrix(4, 3);
        let b = rng.matrix(3, 5);
        let c = matmul(&a, &b).unwrap();
        for i in 0..4 {
            for j in 0..5 {
                let expected = (0..3).fold(0.0, |s, l| s + a[(i, l)] * b[(l, j)]);
                assert_eq!(c[(i, j)].to_bits(), expected.to_bits());
            }
        }
        assert!(matmul(&b, &b).is_err());
        let x = [0.5, -1.0, 2.0];
        let y = gemv(&a, &x).unwrap();
        for i in 0..4 {
            let expected = (0..3).fold(0.0, |s, l| s + a[(i, l)] * x[l]);
            assert_eq!(y[i].to_bits(), expected.to_bits());
        }
        // All dimensions above 5: fused blocks of 256 products.
        let a = rng.matrix(7, 300);
        let b = rng.matrix(300, 6);
        let c = matmul(&a, &b).unwrap();
        for i in 0..7 {
            for j in 0..6 {
                let block = |range: core::ops::Range<usize>| {
                    range.fold(0.0_f64, |s, l| a[(i, l)].mul_add(b[(l, j)], s))
                };
                let expected = block(0..256) + block(256..300);
                assert_eq!(c[(i, j)].to_bits(), expected.to_bits());
            }
        }
        let a = rng.matrix(4, 3);
        let y = gemv(&a, &x).unwrap();
        let mut z = [0.0; 3];
        tr_gemv_to(&a, &y, &mut z).unwrap();
        let at = a.transpose().unwrap();
        let expected = gemv(&at, &y).unwrap();
        assert!(max_abs_diff(&z, &expected) < 1e-14);
        assert!(gemv_to(&a, &[1.0], &mut [0.0; 4]).is_err());
    }
}
