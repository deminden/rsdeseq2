//! Cholesky factorization `A = L Lᵀ` of a symmetric positive definite matrix
//! (GVL 4.2, gaxpy Cholesky, Algorithm 4.2.1).
//!
//! Column `j` of `L` is formed from column `j` of the lower triangle of `A`
//! by subtracting `L[j, k] * L[j.., k]` for every earlier column `k`, in
//! order, and then dividing by the square root of the resulting diagonal
//! entry. Only the lower triangle of `A` is read.

use super::DenseMatrix;
use super::triangular::{solve_lower_in_place, solve_lower_transpose_in_place};
use crate::errors::Error;

/// The lower-triangular Cholesky factor of a symmetric positive definite
/// matrix.
#[derive(Clone, Debug)]
pub struct Cholesky {
    factor: DenseMatrix,
}

impl Cholesky {
    /// Factorize `a`, consuming it.
    ///
    /// Returns `None` if `a` is not square or a diagonal entry of the
    /// reduction is not strictly positive (including `NaN`), that is, if `a`
    /// is not numerically positive definite.
    pub fn new(mut a: DenseMatrix) -> Option<Self> {
        let n = a.nrows();
        if a.ncols() != n {
            return None;
        }
        let values = a.as_mut_slice();
        for j in 0..n {
            let (done, rest) = values.split_at_mut(j * n);
            let column = &mut rest[..n];
            for k in 0..j {
                let previous = &done[k * n..(k + 1) * n];
                let factor = previous[j];
                for (target, &value) in column[j..].iter_mut().zip(&previous[j..]) {
                    *target -= factor * value;
                }
            }
            let diagonal = column[j];
            if diagonal.is_nan() || diagonal <= 0.0 {
                return None;
            }
            let root = diagonal.sqrt();
            column[j] = root;
            for value in &mut column[j + 1..] {
                *value /= root;
            }
            column[..j].fill(0.0);
        }
        Some(Self { factor: a })
    }

    /// The lower-triangular factor `L` (zero above the diagonal).
    pub fn l(&self) -> &DenseMatrix {
        &self.factor
    }

    /// Solve `A x = b` in place: `L y = b` by forward substitution, then
    /// `Lᵀ x = y` by back substitution. `rhs` holds one or more right-hand
    /// sides column by column.
    pub fn solve_in_place(&self, rhs: &mut [f64]) {
        // The diagonal of L is positive, so neither solve can fail.
        solve_lower_in_place(&self.factor, rhs);
        solve_lower_transpose_in_place(&self.factor, rhs);
    }

    /// The solution of `A x = b` as a new vector.
    pub fn solve(&self, rhs: &[f64]) -> Result<Vec<f64>, Error> {
        let mut x = crate::matrix::try_alloc_matrix(rhs.len(), 1, "Cholesky solution")?;
        x.extend_from_slice(rhs);
        self.solve_in_place(&mut x);
        Ok(x)
    }
}

#[cfg(test)]
mod tests {
    use super::super::test_support::{Lcg, max_abs_diff};
    use super::super::{gemv, matmul};
    use super::*;

    #[test]
    fn random_positive_definite_systems_have_small_residuals() {
        let mut rng = Lcg(17);
        for n in [1, 2, 4, 7, 15] {
            let a = rng.spd(n);
            let x: Vec<f64> = (0..n).map(|_| rng.next()).collect();
            let b = gemv(&a, &x).unwrap();
            let cholesky = Cholesky::new(a.clone()).unwrap();
            assert!(max_abs_diff(&cholesky.solve(&b).unwrap(), &x) < 1e-12);
            let l = cholesky.l();
            let product = matmul(l, &l.transpose().unwrap()).unwrap();
            assert!(max_abs_diff(product.as_slice(), a.as_slice()) < 1e-12);
            for j in 0..n {
                assert!(l[(j, j)] > 0.0);
                for i in 0..j {
                    assert_eq!(l[(i, j)], 0.0);
                }
            }
        }
    }

    #[test]
    fn known_factorization() {
        // [[4, 2], [2, 10]] = [[2, 0], [1, 3]] [[2, 1], [0, 3]].
        let a = DenseMatrix::from_row_major(2, 2, &[4.0, 2.0, 2.0, 10.0]).unwrap();
        let cholesky = Cholesky::new(a).unwrap();
        assert_eq!(cholesky.l().as_slice(), &[2.0, 1.0, 0.0, 3.0]);
        assert_eq!(cholesky.solve(&[6.0, 12.0]).unwrap(), vec![1.0, 1.0]);
    }

    #[test]
    fn indefinite_and_singular_matrices_are_rejected() {
        let indefinite = DenseMatrix::from_row_major(2, 2, &[1.0, 2.0, 2.0, 1.0]).unwrap();
        assert!(Cholesky::new(indefinite).is_none());
        let singular = DenseMatrix::from_row_major(2, 2, &[1.0, 1.0, 1.0, 1.0]).unwrap();
        assert!(Cholesky::new(singular).is_none());
        assert!(Cholesky::new(DenseMatrix::from_element(1, 1, f64::NAN).unwrap()).is_none());
        assert!(Cholesky::new(DenseMatrix::zeros(2, 3).unwrap()).is_none());
    }
}
