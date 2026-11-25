//! Eigen-decomposition of a symmetric matrix by the cyclic Jacobi method
//! (GVL 8.5.2 to 8.5.3, Algorithms 8.5.1 and 8.5.3).
//!
//! Sweeps visit the off-diagonal positions `(p, q)`, `p < q`, row by row.
//! Each nonzero `a[p, q]` is annihilated by the rotation `J(p, q, θ)` with
//! `τ = (a[q, q] - a[p, p]) / (2 a[p, q])`, `t = sign(τ) / (|τ| + √(1 + τ²))`,
//! `c = 1 / √(1 + t²)` and `s = t c` (the smaller rotation angle,
//! Algorithm 8.5.1), applied as `A ← Jᵀ A J` and accumulated into
//! `V ← V J`. The iteration stops when the off-diagonal Frobenius norm is at
//! most machine precision times the Frobenius norm of the matrix.

use super::DenseMatrix;
use crate::errors::Error;

/// Maximum number of sweeps; cyclic Jacobi converges quadratically, so
/// well-scaled matrices need far fewer.
const MAX_SWEEPS: usize = 100;

/// Eigenvalues and orthonormal eigenvectors of a symmetric matrix.
#[derive(Clone, Debug)]
pub struct SymmetricEigen {
    /// The eigenvalues, in the order of the converged diagonal (unsorted).
    pub eigenvalues: Vec<f64>,
    /// The eigenvectors, column `j` belonging to `eigenvalues[j]`.
    pub eigenvectors: DenseMatrix,
}

/// The rotation `(c, s)` that annihilates `a_pq` in the symmetric 2 by 2
/// block `[[a_pp, a_pq], [a_pq, a_qq]]` (GVL Algorithm 8.5.1).
pub(crate) fn symmetric_rotation(a_pp: f64, a_qq: f64, a_pq: f64) -> (f64, f64) {
    if a_pq == 0.0 {
        return (1.0, 0.0);
    }
    let tau = (a_qq - a_pp) / (2.0 * a_pq);
    let t = if tau.abs() > 1e150 {
        // √(1 + τ²) = |τ| to working precision; avoid squaring τ.
        0.5 / tau
    } else {
        1.0_f64.copysign(tau) / (tau.abs() + (1.0 + tau * tau).sqrt())
    };
    let c = 1.0 / (1.0 + t * t).sqrt();
    (c, t * c)
}

impl SymmetricEigen {
    /// Decompose the symmetric matrix `a` (only symmetric input is
    /// meaningful; both triangles are read and updated).
    ///
    /// Returns `Ok(None)` if `a` has a nonfinite entry or the sweeps do not
    /// converge, and an error if `a` is not square.
    pub fn new(mut a: DenseMatrix) -> Result<Option<Self>, Error> {
        let n = a.nrows();
        if a.ncols() != n {
            return Err(Error::MatrixShapeMismatch {
                context: "eigen-decomposition of a non-square matrix".to_string(),
                expected: (n, n),
                actual: a.shape(),
            });
        }
        if a.iter().any(|v| !v.is_finite()) {
            return Ok(None);
        }
        let mut v = DenseMatrix::try_zeros(n, n, "eigenvectors")?;
        v.fill_with_identity();
        let total = a.iter().fold(0.0, |sum, x| sum + x * x);
        let tolerance = f64::EPSILON * f64::EPSILON * total;
        let mut converged = false;
        for _ in 0..=MAX_SWEEPS {
            let off: f64 = (0..n)
                .flat_map(|j| (0..n).filter(move |&i| i != j).map(move |i| (i, j)))
                .fold(0.0, |sum, (i, j)| sum + a[(i, j)] * a[(i, j)]);
            if off <= tolerance {
                converged = true;
                break;
            }
            for p in 0..n {
                for q in p + 1..n {
                    let (c, s) = symmetric_rotation(a[(p, p)], a[(q, q)], a[(p, q)]);
                    if s == 0.0 {
                        continue;
                    }
                    rotate_columns(&mut a, p, q, c, s);
                    rotate_rows(&mut a, p, q, c, s);
                    a[(p, q)] = 0.0;
                    a[(q, p)] = 0.0;
                    rotate_columns(&mut v, p, q, c, s);
                }
            }
            if a.iter().any(|x| !x.is_finite()) {
                return Ok(None);
            }
        }
        if !converged {
            return Ok(None);
        }
        let eigenvalues = (0..n).map(|i| a[(i, i)]).collect();
        Ok(Some(Self {
            eigenvalues,
            eigenvectors: v,
        }))
    }
}

/// `M ← M J(p, q, θ)`: columns `p` and `q` become `c m_p - s m_q` and
/// `s m_p + c m_q`.
pub(crate) fn rotate_columns(m: &mut DenseMatrix, p: usize, q: usize, c: f64, s: f64) {
    for k in 0..m.nrows() {
        let (x, y) = (m[(k, p)], m[(k, q)]);
        m[(k, p)] = c * x - s * y;
        m[(k, q)] = s * x + c * y;
    }
}

/// `M ← J(p, q, θ)ᵀ M`: rows `p` and `q` become `c m_p - s m_q` and
/// `s m_p + c m_q`.
fn rotate_rows(m: &mut DenseMatrix, p: usize, q: usize, c: f64, s: f64) {
    for k in 0..m.ncols() {
        let (x, y) = (m[(p, k)], m[(q, k)]);
        m[(p, k)] = c * x - s * y;
        m[(q, k)] = s * x + c * y;
    }
}

#[cfg(test)]
mod tests {
    use super::super::matmul;
    use super::super::test_support::{Lcg, max_abs_diff};
    use super::*;

    fn reconstruct(e: &SymmetricEigen) -> DenseMatrix {
        let n = e.eigenvalues.len();
        let scaled = DenseMatrix::try_from_fn(n, n, "test", |i, j| {
            e.eigenvectors[(i, j)] * e.eigenvalues[j]
        })
        .unwrap();
        matmul(&scaled, &e.eigenvectors.transpose().unwrap()).unwrap()
    }

    #[test]
    fn random_symmetric_matrices_are_diagonalized() {
        let mut rng = Lcg(31);
        for n in [1, 2, 3, 6, 10] {
            let b = rng.matrix(n, n);
            let a = DenseMatrix::try_from_fn(n, n, "test", |i, j| b[(i, j)] + b[(j, i)]).unwrap();
            let e = SymmetricEigen::new(a.clone()).unwrap().unwrap();
            assert!(max_abs_diff(reconstruct(&e).as_slice(), a.as_slice()) < 1e-12);
            let v = &e.eigenvectors;
            let gram = matmul(&v.transpose().unwrap(), v).unwrap();
            let identity = DenseMatrix::identity(n).unwrap();
            assert!(max_abs_diff(gram.as_slice(), identity.as_slice()) < 1e-13);
        }
    }

    #[test]
    fn known_spectrum_and_singular_matrices() {
        let a = DenseMatrix::from_row_major(2, 2, &[2.0, 1.0, 1.0, 2.0]).unwrap();
        let mut values = SymmetricEigen::new(a).unwrap().unwrap().eigenvalues;
        values.sort_by(f64::total_cmp);
        assert!((values[0] - 1.0).abs() < 1e-15 && (values[1] - 3.0).abs() < 1e-15);
        // Rank one: eigenvalues 0 and 2.
        let a = DenseMatrix::from_element(2, 2, 1.0).unwrap();
        let mut values = SymmetricEigen::new(a).unwrap().unwrap().eigenvalues;
        values.sort_by(f64::total_cmp);
        assert!(values[0].abs() < 1e-15 && (values[1] - 2.0).abs() < 1e-15);
        let zero = SymmetricEigen::new(DenseMatrix::zeros(3, 3).unwrap())
            .unwrap()
            .unwrap();
        assert_eq!(zero.eigenvalues, vec![0.0; 3]);
        let nan = DenseMatrix::from_element(2, 2, f64::NAN).unwrap();
        assert!(SymmetricEigen::new(nan).unwrap().is_none());
        assert!(SymmetricEigen::new(DenseMatrix::zeros(2, 3).unwrap()).is_err());
    }
}
