//! Gaussian elimination with partial pivoting (GVL 3.2 and 3.4.2).
//!
//! Step `k` picks the first row `p >= k` whose entry in column `k` has the
//! largest magnitude, swaps rows `k` and `p` across the whole matrix, scales
//! the entries below the pivot by the pivot's reciprocal (the multipliers of
//! `L`) and subtracts the outer product of the multipliers and the pivot row
//! from the trailing submatrix, column by column (GVL Algorithm 3.4.1, outer
//! product form). The factors overwrite the matrix: `U` on and above the
//! diagonal, the unit lower factor `L` below it.

use super::DenseMatrix;
use super::triangular::{solve_lower_unit_in_place, solve_upper_in_place};
use crate::errors::Error;

/// Position of the first entry of largest magnitude in `values`.
fn first_largest(values: &[f64]) -> usize {
    let mut best = 0;
    let mut largest = values.first().map_or(0.0, |value| value.abs());
    for (i, value) in values.iter().enumerate().skip(1) {
        if value.abs() > largest {
            largest = value.abs();
            best = i;
        }
    }
    best
}

/// One elimination step at column `k` with the (already swapped) nonzero
/// `pivot` at `(k, k)`.
fn eliminate(a: &mut DenseMatrix, k: usize, pivot: f64) {
    let n = a.nrows();
    let reciprocal = 1.0 / pivot;
    for value in &mut a.column_mut(k)[k + 1..] {
        *value *= reciprocal;
    }
    let ncols = a.ncols();
    let (left, right) = a.as_mut_slice().split_at_mut((k + 1) * n);
    let multipliers = &left[k * n + k + 1..k * n + n];
    for column in right.chunks_exact_mut(n).take(ncols - k - 1) {
        let factor = column[k];
        for (target, &multiplier) in column[k + 1..].iter_mut().zip(multipliers) {
            *target -= factor * multiplier;
        }
    }
}

/// Factorize the square matrix `a` in place and solve `a x = b` for the
/// right-hand side columns stored column by column in `rhs`.
///
/// Every row interchange is applied to `rhs` as soon as it is chosen; the
/// unit lower and the upper triangular systems are then solved in turn.
/// Returns `false` at the first zero pivot (the matrix is singular), leaving
/// `a` and `rhs` partly updated. `rhs.len()` must be a multiple of the order.
pub fn lu_solve_in_place(a: &mut DenseMatrix, rhs: &mut [f64]) -> bool {
    let n = a.ncols();
    if a.nrows() != n {
        return false;
    }
    for k in 0..n {
        let pivot_row = k + first_largest(&a.column(k)[k..]);
        let pivot = a[(pivot_row, k)];
        if pivot == 0.0 {
            return false;
        }
        if pivot_row != k {
            a.swap_rows(k, pivot_row);
            for column in rhs.chunks_exact_mut(n) {
                column.swap(k, pivot_row);
            }
        }
        eliminate(a, k, pivot);
    }
    solve_lower_unit_in_place(a, rhs);
    solve_upper_in_place(a, rhs)
}

/// An LU factorization `P A = L U` with partial pivoting.
///
/// A column without a nonzero pivot is left unreduced; the factorization of
/// a singular matrix then has a zero on the diagonal of `U`, which
/// [`Lu::solve_in_place`] reports and [`Lu::determinant`] turns into zero.
#[derive(Clone, Debug, Default)]
pub struct Lu {
    factors: DenseMatrix,
    /// `swaps[k]` is the row exchanged with row `k` at step `k`.
    swaps: Vec<usize>,
}

impl Lu {
    /// Factorize the square matrix `a`, consuming it.
    pub fn new(mut a: DenseMatrix) -> Result<Self, Error> {
        let n = a.ncols();
        if a.nrows() != n {
            return Err(Error::MatrixShapeMismatch {
                context: "LU factorization of a non-square matrix".to_string(),
                expected: (n, n),
                actual: a.shape(),
            });
        }
        let mut swaps = crate::matrix::try_alloc_matrix(n, 1, "LU row interchanges")?;
        for k in 0..n {
            let pivot_row = k + first_largest(&a.column(k)[k..]);
            swaps.push(pivot_row);
            let pivot = a[(pivot_row, k)];
            if pivot == 0.0 {
                continue;
            }
            a.swap_rows(k, pivot_row);
            eliminate(&mut a, k, pivot);
        }
        Ok(Self { factors: a, swaps })
    }

    /// The matrix order.
    pub fn order(&self) -> usize {
        self.swaps.len()
    }

    /// The combined factors: `U` on and above the diagonal, the multipliers
    /// of the unit lower factor below it.
    pub fn factors(&self) -> &DenseMatrix {
        &self.factors
    }

    /// Solve `A x = b` in place for the right-hand side columns stored column
    /// by column in `rhs`: the interchanges are applied in order, then the
    /// unit lower and the upper triangular systems are solved. Returns
    /// `false` if `U` has a zero on its diagonal.
    pub fn solve_in_place(&self, rhs: &mut [f64]) -> bool {
        let n = self.order();
        for column in rhs.chunks_exact_mut(n.max(1)) {
            for (k, &row) in self.swaps.iter().enumerate() {
                column.swap(k, row);
            }
        }
        solve_lower_unit_in_place(&self.factors, rhs);
        solve_upper_in_place(&self.factors, rhs)
    }

    /// The determinant: the product of the diagonal of `U`, taken in order,
    /// with the sign of the row permutation (GVL 3.2.1).
    pub fn determinant(&self) -> f64 {
        let mut determinant = 1.0;
        for (k, &row) in self.swaps.iter().enumerate() {
            determinant *= self.factors[(k, k)];
            if row != k {
                determinant = -determinant;
            }
        }
        determinant
    }

    /// Write `A⁻¹` into `out` (resized to the order) by solving against the
    /// identity. Returns `false` if `A` is singular.
    pub fn try_inverse_to(&self, out: &mut DenseMatrix) -> Result<bool, Error> {
        let n = self.order();
        if out.shape() != (n, n) {
            *out = DenseMatrix::try_zeros(n, n, "LU inverse")?;
        }
        out.fill_with_identity();
        Ok(self.solve_in_place(out.as_mut_slice()))
    }
}

#[cfg(test)]
mod tests {
    use super::super::test_support::{Lcg, max_abs_diff};
    use super::super::{gemv, matmul};
    use super::*;

    #[test]
    fn random_systems_have_small_residuals() {
        let mut rng = Lcg(3);
        for n in [1, 2, 3, 6, 12] {
            let a = rng.well_conditioned(n);
            let x: Vec<f64> = (0..n).map(|_| rng.next()).collect();
            let b = gemv(&a, &x).unwrap();

            let mut work = a.clone();
            let mut solution = b.clone();
            assert!(lu_solve_in_place(&mut work, &mut solution));
            assert!(max_abs_diff(&solution, &x) < 1e-13);

            let lu = Lu::new(a.clone()).unwrap();
            let mut solution = b.clone();
            assert!(lu.solve_in_place(&mut solution));
            assert!(max_abs_diff(&solution, &x) < 1e-13);
            // Both routes perform the same operations.
            assert_eq!(lu.factors(), &work);

            let mut inverse = DenseMatrix::default();
            assert!(lu.try_inverse_to(&mut inverse).unwrap());
            let identity = matmul(&a, &inverse).unwrap();
            let expected = DenseMatrix::identity(n).unwrap();
            assert!(max_abs_diff(identity.as_slice(), expected.as_slice()) < 1e-13);
        }
    }

    #[test]
    fn known_factorization_and_determinant() {
        // Partial pivoting swaps the rows of [[1, 2], [3, 4]]:
        // P A = [[1, 0], [1/3, 1]] [[3, 4], [0, 2/3]].
        let a = DenseMatrix::from_row_major(2, 2, &[1.0, 2.0, 3.0, 4.0]).unwrap();
        let lu = Lu::new(a).unwrap();
        let f = lu.factors();
        assert_eq!((f[(0, 0)], f[(0, 1)], f[(1, 0)]), (3.0, 4.0, 1.0 / 3.0));
        assert!((f[(1, 1)] - 2.0 / 3.0).abs() < 1e-15);
        assert!((lu.determinant() + 2.0).abs() < 1e-15);
        let a = DenseMatrix::from_row_major(3, 3, &[2.0, 0.0, 1.0, 1.0, 3.0, 2.0, 1.0, 1.0, 2.0])
            .unwrap();
        assert!((Lu::new(a).unwrap().determinant() - 6.0).abs() < 1e-14);
    }

    #[test]
    fn singular_matrices_are_detected() {
        let a = DenseMatrix::from_row_major(3, 3, &[1.0, 2.0, 3.0, 2.0, 4.0, 6.0, 1.0, 0.0, 1.0])
            .unwrap();
        let mut work = a.clone();
        assert!(!lu_solve_in_place(&mut work, &mut [1.0, 2.0, 3.0]));
        let lu = Lu::new(a).unwrap();
        assert_eq!(lu.determinant(), 0.0);
        assert!(!lu.solve_in_place(&mut [1.0, 2.0, 3.0]));
        let mut out = DenseMatrix::default();
        assert!(!lu.try_inverse_to(&mut out).unwrap());
        let zero = DenseMatrix::zeros(2, 2).unwrap();
        assert!(!lu_solve_in_place(&mut zero.clone(), &mut [1.0, 1.0]));
        assert!(Lu::new(DenseMatrix::zeros(2, 3).unwrap()).is_err());
        assert!(lu_solve_in_place(&mut DenseMatrix::default(), &mut []));
    }
}
