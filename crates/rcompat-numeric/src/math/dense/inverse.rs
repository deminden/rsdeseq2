//! Inverses and determinants of small square matrices.
//!
//! Up to order 4 the inverse is the adjugate (the transposed matrix of
//! cofactors) divided by the determinant, the determinant being the cofactor
//! expansion along the first row; each cofactor is written out as a sum of
//! signed products of the entries. Larger matrices use the LU factorization
//! with partial pivoting (GVL 3.2 and 3.4).

use super::DenseMatrix;
use super::lu::Lu;
use crate::errors::Error;

/// Invert the square matrix `m` in place.
///
/// Orders 0 to 4 use the closed form described in the module documentation;
/// larger orders solve `A X = I` through [`Lu`]. Returns `Ok(false)`, with
/// `m` unchanged for orders up to 4, if the determinant (or a pivot) is
/// exactly zero. Returns an error for a non-square matrix.
pub fn try_inverse_in_place(m: &mut DenseMatrix) -> Result<bool, Error> {
    let n = m.nrows();
    if m.ncols() != n {
        return Err(Error::MatrixShapeMismatch {
            context: "inverse of a non-square matrix".to_string(),
            expected: (n, n),
            actual: m.shape(),
        });
    }
    let values = m.as_mut_slice();
    Ok(match n {
        0 => true,
        1 => invert1(values),
        2 => invert2(values),
        3 => invert3(values),
        4 => invert4(values),
        _ => {
            let lu = Lu::new(m.clone())?;
            lu.try_inverse_to(m)?
        }
    })
}

fn invert1(m: &mut [f64]) -> bool {
    if m[0] == 0.0 {
        return false;
    }
    m[0] = 1.0 / m[0];
    true
}

fn invert2(m: &mut [f64]) -> bool {
    // Column-major: m = [a c; b d] with a = m[0], b = m[1], c = m[2], d = m[3].
    let determinant = m[0] * m[3] - m[1] * m[2];
    if determinant == 0.0 {
        return false;
    }
    let (a, b, c, d) = (m[0], m[1], m[2], m[3]);
    m[0] = d / determinant;
    m[1] = -b / determinant;
    m[2] = -c / determinant;
    m[3] = a / determinant;
    true
}

fn invert3(m: &mut [f64]) -> bool {
    let e = |i: usize, j: usize| m[i + 3 * j];
    let minor_00 = e(1, 1) * e(2, 2) - e(2, 1) * e(1, 2);
    let minor_01 = e(1, 0) * e(2, 2) - e(2, 0) * e(1, 2);
    let minor_02 = e(1, 0) * e(2, 1) - e(2, 0) * e(1, 1);
    let determinant = e(0, 0) * minor_00 - e(0, 1) * minor_01 + e(0, 2) * minor_02;
    if determinant == 0.0 {
        return false;
    }
    let inverse = [
        minor_00 / determinant,
        -minor_01 / determinant,
        minor_02 / determinant,
        (e(0, 2) * e(2, 1) - e(2, 2) * e(0, 1)) / determinant,
        (e(0, 0) * e(2, 2) - e(2, 0) * e(0, 2)) / determinant,
        (e(0, 1) * e(2, 0) - e(2, 1) * e(0, 0)) / determinant,
        (e(0, 1) * e(1, 2) - e(1, 1) * e(0, 2)) / determinant,
        (e(0, 2) * e(1, 0) - e(1, 2) * e(0, 0)) / determinant,
        (e(0, 0) * e(1, 1) - e(1, 0) * e(0, 1)) / determinant,
    ];
    m.copy_from_slice(&inverse);
    true
}

/// Cofactors of a 4 by 4 matrix stored column by column, each the expansion
/// of its 3 by 3 minor into six signed triple products.
fn cofactors4(m: &[f64]) -> [f64; 16] {
    let mut c = [0.0; 16];
    c[0] = m[5] * m[10] * m[15] - m[5] * m[11] * m[14] - m[9] * m[6] * m[15]
        + m[9] * m[7] * m[14]
        + m[13] * m[6] * m[11]
        - m[13] * m[7] * m[10];
    c[4] = -m[4] * m[10] * m[15] + m[4] * m[11] * m[14] + m[8] * m[6] * m[15]
        - m[8] * m[7] * m[14]
        - m[12] * m[6] * m[11]
        + m[12] * m[7] * m[10];
    c[8] = m[4] * m[9] * m[15] - m[4] * m[11] * m[13] - m[8] * m[5] * m[15]
        + m[8] * m[7] * m[13]
        + m[12] * m[5] * m[11]
        - m[12] * m[7] * m[9];
    c[12] = -m[4] * m[9] * m[14] + m[4] * m[10] * m[13] + m[8] * m[5] * m[14]
        - m[8] * m[6] * m[13]
        - m[12] * m[5] * m[10]
        + m[12] * m[6] * m[9];
    c[1] = -m[1] * m[10] * m[15] + m[1] * m[11] * m[14] + m[9] * m[2] * m[15]
        - m[9] * m[3] * m[14]
        - m[13] * m[2] * m[11]
        + m[13] * m[3] * m[10];
    c[5] = m[0] * m[10] * m[15] - m[0] * m[11] * m[14] - m[8] * m[2] * m[15]
        + m[8] * m[3] * m[14]
        + m[12] * m[2] * m[11]
        - m[12] * m[3] * m[10];
    c[9] = -m[0] * m[9] * m[15] + m[0] * m[11] * m[13] + m[8] * m[1] * m[15]
        - m[8] * m[3] * m[13]
        - m[12] * m[1] * m[11]
        + m[12] * m[3] * m[9];
    c[13] = m[0] * m[9] * m[14] - m[0] * m[10] * m[13] - m[8] * m[1] * m[14]
        + m[8] * m[2] * m[13]
        + m[12] * m[1] * m[10]
        - m[12] * m[2] * m[9];
    c[2] = m[1] * m[6] * m[15] - m[1] * m[7] * m[14] - m[5] * m[2] * m[15]
        + m[5] * m[3] * m[14]
        + m[13] * m[2] * m[7]
        - m[13] * m[3] * m[6];
    c[6] = -m[0] * m[6] * m[15] + m[0] * m[7] * m[14] + m[4] * m[2] * m[15]
        - m[4] * m[3] * m[14]
        - m[12] * m[2] * m[7]
        + m[12] * m[3] * m[6];
    c[10] = m[0] * m[5] * m[15] - m[0] * m[7] * m[13] - m[4] * m[1] * m[15]
        + m[4] * m[3] * m[13]
        + m[12] * m[1] * m[7]
        - m[12] * m[3] * m[5];
    c[14] = -m[0] * m[5] * m[14] + m[0] * m[6] * m[13] + m[4] * m[1] * m[14]
        - m[4] * m[2] * m[13]
        - m[12] * m[1] * m[6]
        + m[12] * m[2] * m[5];
    c[3] = -m[1] * m[6] * m[11] + m[1] * m[7] * m[10] + m[5] * m[2] * m[11]
        - m[5] * m[3] * m[10]
        - m[9] * m[2] * m[7]
        + m[9] * m[3] * m[6];
    c[7] = m[0] * m[6] * m[11] - m[0] * m[7] * m[10] - m[4] * m[2] * m[11]
        + m[4] * m[3] * m[10]
        + m[8] * m[2] * m[7]
        - m[8] * m[3] * m[6];
    c[11] = -m[0] * m[5] * m[11] + m[0] * m[7] * m[9] + m[4] * m[1] * m[11]
        - m[4] * m[3] * m[9]
        - m[8] * m[1] * m[7]
        + m[8] * m[3] * m[5];
    c[15] = m[0] * m[5] * m[10] - m[0] * m[6] * m[9] - m[4] * m[1] * m[10]
        + m[4] * m[2] * m[9]
        + m[8] * m[1] * m[6]
        - m[8] * m[2] * m[5];
    c
}

fn invert4(m: &mut [f64]) -> bool {
    let adjugate = cofactors4(m);
    let determinant =
        m[0] * adjugate[0] + m[1] * adjugate[4] + m[2] * adjugate[8] + m[3] * adjugate[12];
    if determinant == 0.0 {
        return false;
    }
    let reciprocal = 1.0 / determinant;
    for (target, value) in m.iter_mut().zip(adjugate) {
        *target = value * reciprocal;
    }
    true
}

/// Determinant of the square matrix `m`.
///
/// Orders up to 3 use the cofactor expansion along the first row; larger
/// orders the LU factorization ([`Lu::determinant`]). Returns an error for a
/// non-square matrix.
pub fn determinant(m: &DenseMatrix) -> Result<f64, Error> {
    let n = m.nrows();
    if m.ncols() != n {
        return Err(Error::MatrixShapeMismatch {
            context: "determinant of a non-square matrix".to_string(),
            expected: (n, n),
            actual: m.shape(),
        });
    }
    let e = |i: usize, j: usize| m[(i, j)];
    Ok(match n {
        0 => 1.0,
        1 => e(0, 0),
        2 => e(0, 0) * e(1, 1) - e(1, 0) * e(0, 1),
        3 => {
            let minor_00 = e(1, 1) * e(2, 2) - e(2, 1) * e(1, 2);
            let minor_01 = e(1, 0) * e(2, 2) - e(2, 0) * e(1, 2);
            let minor_02 = e(1, 0) * e(2, 1) - e(2, 0) * e(1, 1);
            e(0, 0) * minor_00 - e(0, 1) * minor_01 + e(0, 2) * minor_02
        }
        _ => Lu::new(m.clone())?.determinant(),
    })
}

#[cfg(test)]
mod tests {
    use super::super::matmul;
    use super::super::test_support::{Lcg, max_abs_diff};
    use super::*;

    #[test]
    fn inverses_of_random_matrices_have_small_residuals() {
        let mut rng = Lcg(5);
        for n in 0..=8 {
            for _ in 0..4 {
                let a = rng.well_conditioned(n);
                let mut inverse = a.clone();
                assert!(try_inverse_in_place(&mut inverse).unwrap());
                let product = matmul(&a, &inverse).unwrap();
                let identity = DenseMatrix::identity(n).unwrap();
                assert!(max_abs_diff(product.as_slice(), identity.as_slice()) < 1e-13);
                if n > 0 {
                    let lu = Lu::new(a.clone()).unwrap().determinant();
                    let closed = determinant(&a).unwrap();
                    assert!((lu - closed).abs() <= 1e-12 * lu.abs());
                }
            }
        }
    }

    #[test]
    fn known_inverses_and_determinants() {
        let mut a = DenseMatrix::from_row_major(2, 2, &[4.0, 7.0, 2.0, 6.0]).unwrap();
        assert_eq!(determinant(&a).unwrap(), 10.0);
        assert!(try_inverse_in_place(&mut a).unwrap());
        assert_eq!(a.as_slice(), &[0.6, -0.2, -0.7, 0.4]);
        let mut d = DenseMatrix::zeros(4, 4).unwrap();
        for (i, v) in [2.0, 4.0, 0.5, -8.0].into_iter().enumerate() {
            d[(i, i)] = v;
        }
        assert_eq!(determinant(&d).unwrap(), -32.0);
        assert!(try_inverse_in_place(&mut d).unwrap());
        assert_eq!(
            (d[(0, 0)], d[(1, 1)], d[(2, 2)], d[(3, 3)], d[(0, 1)]),
            (0.5, 0.25, 2.0, -0.125, 0.0)
        );
    }

    #[test]
    fn singular_matrices_are_detected() {
        for n in 1..=6 {
            let mut a = DenseMatrix::from_element(n, n, 1.0).unwrap();
            if n == 1 {
                a[(0, 0)] = 0.0;
            }
            let original = a.clone();
            assert!(!try_inverse_in_place(&mut a).unwrap(), "order {n}");
            if n <= 4 {
                assert_eq!(a, original);
            }
            assert_eq!(determinant(&original).unwrap(), 0.0);
        }
        assert!(try_inverse_in_place(&mut DenseMatrix::zeros(2, 3).unwrap()).is_err());
        assert!(determinant(&DenseMatrix::zeros(3, 2).unwrap()).is_err());
    }
}
