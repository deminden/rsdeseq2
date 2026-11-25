//! Column-oriented forward and back substitution (GVL 3.1.2 and 3.1.3).
//!
//! Each function solves in place for every column of a right-hand side
//! stored column by column in a slice whose length is a multiple of the
//! matrix order. Only the referenced triangle of the matrix is read.

use super::DenseMatrix;

fn columns(n: usize, rhs: &mut [f64]) -> core::slice::ChunksExactMut<'_, f64> {
    rhs.chunks_exact_mut(n.max(1))
}

/// Solve `L x = b` for lower-triangular `L` with every diagonal entry taken
/// to be one (the unit lower factor of an LU factorization).
///
/// Column-oriented forward substitution (GVL Algorithm 3.1.3): once `x[i]`
/// is known, `x[i] * L[i+1.., i]` is subtracted from the remaining entries.
pub fn solve_lower_unit_in_place(l: &DenseMatrix, rhs: &mut [f64]) {
    let n = l.nrows().min(l.ncols());
    for b in columns(n, rhs) {
        for i in 0..n {
            let coefficient = b[i];
            for (target, &value) in b[i + 1..].iter_mut().zip(&l.column(i)[i + 1..n]) {
                *target -= coefficient * value;
            }
        }
    }
}

/// Solve `L x = b` for lower-triangular `L`.
///
/// Column-oriented forward substitution (GVL Algorithm 3.1.3). Returns
/// `false`, leaving `rhs` partly updated, if a diagonal entry is zero.
pub fn solve_lower_in_place(l: &DenseMatrix, rhs: &mut [f64]) -> bool {
    let n = l.nrows().min(l.ncols());
    for b in columns(n, rhs) {
        for i in 0..n {
            let diagonal = l[(i, i)];
            if diagonal == 0.0 {
                return false;
            }
            let coefficient = b[i] / diagonal;
            b[i] = coefficient;
            for (target, &value) in b[i + 1..].iter_mut().zip(&l.column(i)[i + 1..n]) {
                *target -= coefficient * value;
            }
        }
    }
    true
}

/// Solve `Lᵀ x = b` for lower-triangular `L` without forming the transpose.
///
/// Row-oriented back substitution on `Lᵀ` (GVL Algorithm 3.1.2): `x[i]` is
/// `b[i]` minus the dot product of `L[i+1.., i]` with the known entries,
/// divided by `L[i, i]`. Returns `false` if a diagonal entry is zero.
pub fn solve_lower_transpose_in_place(l: &DenseMatrix, rhs: &mut [f64]) -> bool {
    let n = l.nrows().min(l.ncols());
    for b in columns(n, rhs) {
        for i in (0..n).rev() {
            let diagonal = l[(i, i)];
            if diagonal == 0.0 {
                return false;
            }
            let known = super::dot(&l.column(i)[i + 1..n], &b[i + 1..]);
            b[i] = (b[i] - known) / diagonal;
        }
    }
    true
}

/// Solve `U x = b` for upper-triangular `U`.
///
/// Column-oriented back substitution (GVL Algorithm 3.1.4): once `x[i]` is
/// known, `x[i] * U[..i, i]` is subtracted from the entries above it.
/// Returns `false`, leaving `rhs` partly updated, if a diagonal entry is
/// zero.
pub fn solve_upper_in_place(u: &DenseMatrix, rhs: &mut [f64]) -> bool {
    let n = u.nrows().min(u.ncols());
    for b in columns(n, rhs) {
        for i in (0..n).rev() {
            let diagonal = u[(i, i)];
            if diagonal == 0.0 {
                return false;
            }
            let coefficient = b[i] / diagonal;
            b[i] = coefficient;
            for (target, &value) in b[..i].iter_mut().zip(&u.column(i)[..i]) {
                *target -= coefficient * value;
            }
        }
    }
    true
}

#[cfg(test)]
mod tests {
    use super::super::test_support::{Lcg, max_abs_diff};
    use super::super::{gemv, matmul};
    use super::*;

    fn triangle(rng: &mut Lcg, n: usize, lower: bool) -> DenseMatrix {
        let mut a = rng.well_conditioned(n);
        for j in 0..n {
            for i in 0..n {
                if (lower && i < j) || (!lower && i > j) {
                    a[(i, j)] = 0.0;
                }
            }
        }
        a
    }

    #[test]
    fn substitutions_solve_random_triangular_systems() {
        let mut rng = Lcg(11);
        for n in [1, 2, 5, 9] {
            let x: Vec<f64> = (0..2 * n).map(|_| rng.next()).collect();
            let as_matrix = DenseMatrix::from_column_major(n, 2, x.clone()).unwrap();

            let l = triangle(&mut rng, n, true);
            let mut b = matmul(&l, &as_matrix).unwrap().into_vec();
            assert!(solve_lower_in_place(&l, &mut b));
            assert!(max_abs_diff(&b, &x) < 1e-13);

            let mut b = matmul(&l.transpose().unwrap(), &as_matrix)
                .unwrap()
                .into_vec();
            assert!(solve_lower_transpose_in_place(&l, &mut b));
            assert!(max_abs_diff(&b, &x) < 1e-13);

            let u = triangle(&mut rng, n, false);
            let mut b = matmul(&u, &as_matrix).unwrap().into_vec();
            assert!(solve_upper_in_place(&u, &mut b));
            assert!(max_abs_diff(&b, &x) < 1e-13);

            let mut unit = l.clone();
            for i in 0..n {
                unit[(i, i)] = 1.0;
            }
            let mut b = gemv(&unit, &x[..n]).unwrap();
            // The diagonal is never read.
            for i in 0..n {
                unit[(i, i)] = f64::NAN;
            }
            solve_lower_unit_in_place(&unit, &mut b);
            assert!(max_abs_diff(&b, &x[..n]) < 1e-13);
        }
    }

    #[test]
    fn zero_diagonals_are_reported() {
        let u = DenseMatrix::from_row_major(2, 2, &[1.0, 2.0, 0.0, 0.0]).unwrap();
        assert!(!solve_upper_in_place(&u, &mut [1.0, 1.0]));
        let l = u.transpose().unwrap();
        assert!(!solve_lower_in_place(&l, &mut [1.0, 1.0]));
        assert!(!solve_lower_transpose_in_place(&l, &mut [1.0, 1.0]));
    }
}
