//! Householder reflections (GVL 5.1.2 and 5.1.3), Householder QR
//! (GVL 5.2.2) and QR with column pivoting (GVL 5.4.2, with the column
//! chosen as in complete pivoting, GVL 3.4.8).
//!
//! A reflector is stored as a unit vector `v`: `H = I - 2 v vᵀ`. For a
//! column `x` with `σ = ‖x‖` and `s = σ` carrying the sign of `x₁`, the
//! vector `x + s e₁`, scaled to unit length, maps `x` to `-s e₁`; adding `s`
//! to `x₁` avoids cancellation (GVL 5.1.3). The squared length of
//! `x + s e₁` is `2 (σ² + |x₁| σ)`.

use super::DenseMatrix;
use super::triangular::solve_upper_in_place;
use crate::errors::Error;

/// Turn `column` into the unit Householder vector that reduces it to a
/// multiple of the first unit vector, and return that multiple (the diagonal
/// entry of `R`).
///
/// For `x = column`, `σ = ‖x‖` is the square root of `x · x` (computed by
/// [`super::dot`]) and `s` is `σ` with the sign of `x₁`. The first entry becomes
/// `x₁ + s` and the column is divided by `sqrt(2 (σ² + |x₁| σ))` and then
/// by its own computed norm; the return value is `-s`. A zero column is left as it is (apart from the
/// signed zero added to its first entry) and `s` is returned: there is
/// nothing to reflect.
pub fn make_reflector(column: &mut [f64]) -> f64 {
    let Some(&first) = column.first() else {
        return 0.0;
    };
    let squared_norm = super::dot(column, column);
    let norm = squared_norm.sqrt();
    let signed_norm = norm.copysign(first);
    let length_squared = (squared_norm + first.abs() * norm) * 2.0;
    column[0] = first + signed_norm;
    if length_squared == 0.0 {
        return signed_norm;
    }
    let length = length_squared.sqrt();
    for value in column.iter_mut() {
        *value /= length;
    }
    // Normalize once more by the computed length, removing most of the
    // rounding of the scaling above.
    let length = super::dot(column, column).sqrt();
    for value in column.iter_mut() {
        *value /= length;
    }
    -signed_norm
}

/// Replace `target` by `sign * H target` for the reflector with unit vector
/// `axis`, where `sign` is `1` or `-1`: with `f = 2 (axis · target)` (by
/// [`super::dot`]), entry `i` becomes
/// `sign * target[i] + (-f * sign) * axis[i]`.
pub fn reflect_with_sign(axis: &[f64], target: &mut [f64], sign: f64) {
    let factor = super::dot(axis, target) * 2.0;
    let scaled = -factor * sign;
    for (value, &a) in target.iter_mut().zip(axis) {
        *value = sign * *value + scaled * a;
    }
}

/// One step of Householder QR on column `col` of `a`, from row `col` down.
///
/// The column below and on the diagonal is replaced by its unit reflector
/// ([`make_reflector`]) and every later column, from row `col` down, by
/// `sign * H` applied to it ([`reflect_with_sign`]), where `sign` is the
/// sign of the returned diagonal entry `d`. Row `col` of the later columns
/// then holds row `col` of an `R` factor whose diagonal entry is `|d|`.
pub fn reduce_column(a: &mut DenseMatrix, col: usize) -> f64 {
    let n = a.nrows();
    let ncols = a.ncols();
    let (left, right) = a.as_mut_slice().split_at_mut((col + 1) * n);
    let axis = &mut left[col * n + col..];
    let diagonal = make_reflector(axis);
    let sign = diagonal.signum();
    for column in right.chunks_exact_mut(n).take(ncols - col - 1) {
        reflect_with_sign(axis, &mut column[col..], sign);
    }
    diagonal
}

/// Householder QR with column pivoting, `A P = Q R` (GVL 5.4.2).
///
/// At step `k` the column holding the entry of largest magnitude in the
/// trailing submatrix (rows and columns `k..`; the first such entry in
/// column-major order) is swapped into position `k`, the choice of column
/// made by complete pivoting (GVL 3.4.8), and reduced with
/// [`reduce_column`], so the diagonal of `R` is nonnegative.
#[derive(Clone, Debug)]
pub struct ColPivQr {
    /// Reflectors below and on the diagonal, `R` above it.
    factors: DenseMatrix,
    /// Signed reflector values returned by [`reduce_column`]; `R` has
    /// their magnitudes on its diagonal.
    diagonal: Vec<f64>,
    /// `permutation[k]` is the original column moved to position `k`.
    permutation: Vec<usize>,
}

impl ColPivQr {
    /// Factorize `a` (with at least as many rows as columns), consuming it.
    pub fn new(mut a: DenseMatrix) -> Result<Self, Error> {
        let (m, n) = a.shape();
        if m < n {
            return Err(Error::MatrixShapeMismatch {
                context: "column-pivoted QR of a wide matrix".to_string(),
                expected: (n, n),
                actual: (m, n),
            });
        }
        let mut permutation: Vec<usize> = (0..n).collect();
        let mut diagonal = crate::matrix::try_alloc_matrix(n, 1, "QR diagonal")?;
        for k in 0..n {
            let mut best = k;
            let mut largest = -1.0;
            for j in k..n {
                for value in &a.column(j)[k..] {
                    if value.abs() > largest {
                        largest = value.abs();
                        best = j;
                    }
                }
            }
            if best != k {
                let values = a.as_mut_slice();
                for i in 0..m {
                    values.swap(k * m + i, best * m + i);
                }
                permutation.swap(k, best);
            }
            let diagonal_entry = reduce_column(&mut a, k);
            diagonal.push(diagonal_entry);
        }
        Ok(Self {
            factors: a,
            diagonal,
            permutation,
        })
    }

    /// The column permutation: entry `k` is the original column at position
    /// `k` of `R`.
    pub fn permutation(&self) -> &[usize] {
        &self.permutation
    }

    /// The diagonal of `R` (nonnegative).
    pub fn diagonal(&self) -> Vec<f64> {
        self.diagonal.iter().map(|d| d.abs()).collect()
    }

    /// Replace `b` (with one entry per row of `A`) by `Qᵀ b`, applying the
    /// reflectors in order, each with the sign of its diagonal entry as in
    /// [`reduce_column`].
    pub fn q_transpose_mul(&self, b: &mut [f64]) {
        let m = self.factors.nrows();
        for (k, d) in self.diagonal.iter().enumerate() {
            reflect_with_sign(&self.factors.column(k)[k..m], &mut b[k..m], d.signum());
        }
    }

    /// The thin orthogonal factor `Q` (`m × n`).
    ///
    /// Column `j` is the unit vector `e_j` with the reflectors `j, j-1, ...,
    /// 0` applied in turn, each with the sign of its diagonal entry.
    pub fn q(&self) -> Result<DenseMatrix, Error> {
        let m = self.factors.nrows();
        let n = self.diagonal.len();
        let mut q = DenseMatrix::try_zeros(m, n, "QR orthogonal factor")?;
        for j in 0..n {
            let column = q.column_mut(j);
            column[j] = 1.0;
            for k in (0..=j).rev() {
                let sign = self.diagonal[k].signum();
                reflect_with_sign(&self.factors.column(k)[k..m], &mut column[k..], sign);
            }
        }
        Ok(q)
    }

    /// Undo the column permutation: entry `k` of `z` (in pivoted order)
    /// becomes entry `permutation()[k]` of the result.
    pub fn unpermute(&self, z: &[f64]) -> Vec<f64> {
        let mut x = vec![0.0; z.len()];
        for (&column, &value) in self.permutation.iter().zip(z) {
            x[column] = value;
        }
        x
    }

    /// The upper-triangular factor `R` (`n × n`).
    pub fn r(&self) -> Result<DenseMatrix, Error> {
        let n = self.diagonal.len();
        DenseMatrix::try_from_fn(n, n, "QR triangular factor", |i, j| {
            if i == j {
                self.diagonal[i].abs()
            } else if i < j {
                self.factors[(i, j)]
            } else {
                0.0
            }
        })
    }

    /// The least-squares solution of `A x = b` (the solution when `A` is
    /// square): `R z = (Qᵀ b)[..n]` by back substitution, then `x = P z`.
    /// Returns `None` if `R` has a zero on its diagonal.
    pub fn solve(&self, b: &[f64]) -> Result<Option<Vec<f64>>, Error> {
        let m = self.factors.nrows();
        let n = self.diagonal.len();
        if b.len() != m {
            return Err(Error::InvalidDimensions {
                context: "QR right-hand side".to_string(),
                expected: m,
                actual: b.len(),
            });
        }
        let mut work = crate::matrix::try_alloc_matrix(m, 1, "QR right-hand side")?;
        work.extend_from_slice(b);
        self.q_transpose_mul(&mut work);
        work.truncate(n);
        if !solve_upper_in_place(&self.r()?, &mut work) {
            return Ok(None);
        }
        Ok(Some(self.unpermute(&work)))
    }
}

#[cfg(test)]
mod tests {
    use super::super::gemv;
    use super::super::test_support::{Lcg, max_abs_diff};
    use super::*;

    #[test]
    fn reflector_maps_a_column_onto_the_first_axis() {
        let mut x = [3.0, 4.0];
        let original = x;
        let d = make_reflector(&mut x);
        assert_eq!(d, -5.0);
        assert!((x[0] * x[0] + x[1] * x[1] - 1.0).abs() < 1e-15);
        let mut y = original;
        reflect_with_sign(&x, &mut y, 1.0);
        assert!((y[0] + 5.0).abs() < 1e-14 && y[1].abs() < 1e-14);
        let mut y = original;
        reflect_with_sign(&x, &mut y, -1.0);
        assert!((y[0] - 5.0).abs() < 1e-14 && y[1].abs() < 1e-14);
        let mut zero = [0.0, 0.0];
        assert_eq!(make_reflector(&mut zero), 0.0);
        assert_eq!(make_reflector(&mut []), 0.0);
    }

    #[test]
    fn reduced_columns_form_r_with_a_positive_diagonal() {
        let mut rng = Lcg(23);
        let (m, n) = (7, 4);
        let a = rng.matrix(m, n);
        let mut work = a.clone();
        let mut diagonal = Vec::new();
        for col in 0..n {
            diagonal.push(reduce_column(&mut work, col));
        }
        // Rebuild Qᵀ A column by column from the reflectors and compare
        // with R (rows scaled by the signs of the diagonal).
        for j in 0..n {
            let mut column = a.column(j).to_vec();
            for (k, &d) in diagonal.iter().enumerate() {
                reflect_with_sign(&work.column(k)[k..], &mut column[k..], d.signum());
            }
            for i in 0..n {
                let expected = match i.cmp(&j) {
                    core::cmp::Ordering::Less => work[(i, j)],
                    core::cmp::Ordering::Equal => diagonal[i].abs(),
                    core::cmp::Ordering::Greater => 0.0,
                };
                assert!((column[i] - expected).abs() < 1e-13);
            }
            assert!(column[n..].iter().all(|v| v.abs() < 1e-13));
        }
    }

    #[test]
    fn pivoted_qr_solves_square_and_least_squares_problems() {
        let mut rng = Lcg(29);
        for (m, n) in [(1, 1), (3, 3), (6, 6), (9, 4)] {
            let a = rng.well_conditioned(m);
            let a = DenseMatrix::try_from_fn(m, n, "test", |i, j| a[(i, j)]).unwrap();
            let x: Vec<f64> = (0..n).map(|_| rng.next()).collect();
            let b = gemv(&a, &x).unwrap();
            let qr = ColPivQr::new(a.clone()).unwrap();
            let solution = qr.solve(&b).unwrap().unwrap();
            assert!(max_abs_diff(&solution, &x) < 1e-12);
            let mut p = qr.permutation().to_vec();
            p.sort_unstable();
            assert_eq!(p, (0..n).collect::<Vec<_>>());
            let d = qr.diagonal();
            assert!(d.iter().all(|&v| v > 0.0));
            // Q has orthonormal columns and Qᵀ A P = R.
            let q = qr.q().unwrap();
            let qtq = super::super::matmul(&q.transpose().unwrap(), &q).unwrap();
            let identity = DenseMatrix::identity(n).unwrap();
            assert!(max_abs_diff(qtq.as_slice(), identity.as_slice()) < 1e-13);
            let mut qtb = vec![0.0; n];
            super::super::tr_gemv_to(&q, &b, &mut qtb).unwrap();
            assert!(solve_upper_in_place(&qr.r().unwrap(), &mut qtb));
            assert!(max_abs_diff(&qr.unpermute(&qtb), &x) < 1e-12);
            // Rebuild A P from Q R: apply the reflectors in reverse to R.
            let r = qr.r().unwrap();
            for (position, &column) in qr.permutation().iter().enumerate() {
                let mut v = vec![0.0; m];
                v[..n].copy_from_slice(r.column(position));
                for k in (0..n).rev() {
                    let sign = qr.diagonal[k].signum();
                    reflect_with_sign(&qr.factors.column(k)[k..], &mut v[k..], sign);
                }
                assert!(max_abs_diff(&v, a.column(column)) < 1e-12);
            }
        }
        // An inconsistent overdetermined system: the residual is orthogonal
        // to the columns.
        let a = DenseMatrix::from_row_major(3, 2, &[1.0, 0.0, 1.0, 1.0, 1.0, 2.0]).unwrap();
        let b = [6.0, 0.0, 0.0];
        let x = ColPivQr::new(a.clone())
            .unwrap()
            .solve(&b)
            .unwrap()
            .unwrap();
        assert!((x[0] - 5.0).abs() < 1e-12 && (x[1] + 3.0).abs() < 1e-12);
        let fitted = gemv(&a, &x).unwrap();
        let residual: Vec<f64> = b.iter().zip(&fitted).map(|(b, f)| b - f).collect();
        let normal = gemv(&a.transpose().unwrap(), &residual).unwrap();
        assert!(normal.iter().all(|v| v.abs() < 1e-12));
    }

    #[test]
    fn exactly_singular_systems_are_detected() {
        let a = DenseMatrix::from_row_major(2, 2, &[1.0, 0.0, 0.0, 0.0]).unwrap();
        assert!(
            ColPivQr::new(a)
                .unwrap()
                .solve(&[1.0, 1.0])
                .unwrap()
                .is_none()
        );
        assert!(ColPivQr::new(DenseMatrix::zeros(2, 3).unwrap()).is_err());
    }
}
