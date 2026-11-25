//! Thin singular value decomposition by one-sided Jacobi rotations
//! (GVL 8.6, one-sided Jacobi SVD; the rotations are those of the symmetric
//! Jacobi method of GVL 8.5 applied implicitly to `AᵀA`).
//!
//! The columns of a working copy `U` of `A` (`m ≥ n`) are orthogonalized in
//! cyclic sweeps over the pairs `(p, q)`, `p < q`: with `α = ‖u_p‖²`,
//! `β = ‖u_q‖²` and `γ = u_p · u_q`, the pair is rotated unless
//! `|γ| ≤ ε √(α β)`, using the rotation that annihilates the off-diagonal
//! entry of `[[α, γ], [γ, β]]`; the same rotation accumulates into `V`. When
//! a sweep rotates nothing, the singular values are the column norms of `U`
//! and the left singular vectors its normalized columns. A wide matrix is
//! decomposed through its transpose.

use super::DenseMatrix;
use super::eigen::{rotate_columns, symmetric_rotation};
use crate::errors::Error;

/// Maximum number of sweeps.
const MAX_SWEEPS: usize = 100;

/// Thin singular value decomposition `A = U diag(σ) Vᵀ` of an `m × n`
/// matrix, with `r = min(m, n)` singular values in non-increasing order.
#[derive(Clone, Debug)]
pub struct Svd {
    /// Left singular vectors, `m × r`. The column of a zero singular value
    /// is zero.
    pub u: DenseMatrix,
    /// Singular values, non-increasing.
    pub singular_values: Vec<f64>,
    /// Right singular vectors, transposed: `r × n`.
    pub v_t: DenseMatrix,
}

impl Svd {
    /// Decompose `a`. Returns `Ok(None)` if `a` has a nonfinite entry or the
    /// sweeps do not converge.
    pub fn new(a: DenseMatrix) -> Result<Option<Self>, Error> {
        if a.iter().any(|v| !v.is_finite()) {
            return Ok(None);
        }
        if a.nrows() < a.ncols() {
            let Some(svd) = Self::tall(a.transpose()?)? else {
                return Ok(None);
            };
            // Aᵀ = U S Vᵀ, so A = V S Uᵀ.
            return Ok(Some(Self {
                u: svd.v_t.transpose()?,
                singular_values: svd.singular_values,
                v_t: svd.u.transpose()?,
            }));
        }
        Self::tall(a)
    }

    fn tall(mut u: DenseMatrix) -> Result<Option<Self>, Error> {
        let (m, n) = u.shape();
        let mut v = DenseMatrix::try_zeros(n, n, "right singular vectors")?;
        v.fill_with_identity();
        let mut converged = false;
        for _ in 0..MAX_SWEEPS {
            let mut rotated = false;
            for p in 0..n {
                for q in p + 1..n {
                    let alpha = super::dot(u.column(p), u.column(p));
                    let beta = super::dot(u.column(q), u.column(q));
                    let gamma = super::dot(u.column(p), u.column(q));
                    if gamma.abs() <= f64::EPSILON * (alpha * beta).sqrt() {
                        continue;
                    }
                    let (c, s) = symmetric_rotation(alpha, beta, gamma);
                    if s == 0.0 {
                        continue;
                    }
                    rotated = true;
                    rotate_columns(&mut u, p, q, c, s);
                    rotate_columns(&mut v, p, q, c, s);
                }
            }
            if u.iter().any(|x| !x.is_finite()) {
                return Ok(None);
            }
            if !rotated {
                converged = true;
                break;
            }
        }
        if !converged {
            return Ok(None);
        }
        let norms: Vec<f64> = (0..n)
            .map(|j| super::dot(u.column(j), u.column(j)).sqrt())
            .collect();
        let mut order: Vec<usize> = (0..n).collect();
        order.sort_by(|&i, &j| norms[j].total_cmp(&norms[i]));
        let singular_values: Vec<f64> = order.iter().map(|&j| norms[j]).collect();
        let left = DenseMatrix::try_from_fn(m, n, "left singular vectors", |i, k| {
            let sigma = singular_values[k];
            if sigma > 0.0 {
                u[(i, order[k])] / sigma
            } else {
                0.0
            }
        })?;
        let v_t =
            DenseMatrix::try_from_fn(n, n, "right singular vectors", |k, j| v[(j, order[k])])?;
        Ok(Some(Self {
            u: left,
            singular_values,
            v_t,
        }))
    }
}

#[cfg(test)]
mod tests {
    use super::super::matmul;
    use super::super::test_support::{Lcg, max_abs_diff};
    use super::*;

    fn reconstruct(svd: &Svd) -> DenseMatrix {
        let (m, r) = svd.u.shape();
        let scaled =
            DenseMatrix::try_from_fn(m, r, "test", |i, k| svd.u[(i, k)] * svd.singular_values[k])
                .unwrap();
        matmul(&scaled, &svd.v_t).unwrap()
    }

    #[test]
    fn random_matrices_are_reconstructed() {
        let mut rng = Lcg(37);
        for (m, n) in [(1, 1), (5, 3), (3, 5), (8, 8), (12, 10)] {
            let a = rng.matrix(m, n);
            let svd = Svd::new(a.clone()).unwrap().unwrap();
            let r = m.min(n);
            assert_eq!(svd.u.shape(), (m, r));
            assert_eq!(svd.v_t.shape(), (r, n));
            assert!(max_abs_diff(reconstruct(&svd).as_slice(), a.as_slice()) < 1e-12);
            assert!(svd.singular_values.windows(2).all(|w| w[0] >= w[1]));
            let vvt = matmul(&svd.v_t, &svd.v_t.transpose().unwrap()).unwrap();
            let identity = DenseMatrix::identity(r).unwrap();
            assert!(max_abs_diff(vvt.as_slice(), identity.as_slice()) < 1e-13);
        }
    }

    #[test]
    fn known_singular_values_and_rank() {
        let a = DenseMatrix::from_row_major(2, 2, &[3.0, 0.0, 0.0, -4.0]).unwrap();
        let svd = Svd::new(a).unwrap().unwrap();
        assert_eq!(svd.singular_values, vec![4.0, 3.0]);
        // A_ij = (i + 1) + (j + 1) / 10 has rank two.
        let a =
            DenseMatrix::try_from_fn(6, 5, "test", |i, j| (i + 1) as f64 + (j + 1) as f64 / 10.0)
                .unwrap();
        let svd = Svd::new(a.clone()).unwrap().unwrap();
        assert!(svd.singular_values[1] > 1e-3);
        assert!(svd.singular_values[2..].iter().all(|&s| s < 1e-12));
        assert!(max_abs_diff(reconstruct(&svd).as_slice(), a.as_slice()) < 1e-12);
        let nan = DenseMatrix::from_element(2, 2, f64::NAN).unwrap();
        assert!(Svd::new(nan).unwrap().is_none());
    }
}
