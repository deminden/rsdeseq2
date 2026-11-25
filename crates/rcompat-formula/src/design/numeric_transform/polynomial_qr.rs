//! Householder QR basis for `poly()` and `contr.poly()`.
//!
//! Both bases are the columns of R's `qr.qy(QR, diag(diag(QR$qr), n, p))`
//! for the default `qr()` with tolerance 1e-7, computed by the shared
//! factorization in `rcompat-numeric` (`docs/provenance/specs/householder_qr.md`).
use crate::errors::Error;
use crate::math::householder::HouseholderQr;
use crate::math::rank::R_QR_TOLERANCE;
use crate::matrix::try_alloc_matrix;

/// Allocate `width` columns of length `rows`, entry `(row, column)` from `value`.
pub(in crate::design) fn try_columns(
    rows: usize,
    width: usize,
    value: impl Fn(usize, usize) -> f64,
) -> Result<Vec<Vec<f64>>, Error> {
    let mut columns = try_alloc_matrix(width, 1, "polynomial QR columns")?;
    for column in 0..width {
        let mut values = try_alloc_matrix(rows, 1, "polynomial QR column")?;
        values.extend((0..rows).map(|row| value(row, column)));
        columns.push(values);
    }
    Ok(columns)
}

/// Factorize `columns` and visit the basis vector of every pivoted column.
///
/// The caller supplies a finite, nonempty column-major problem with
/// `width <= rows`. `visit(j, basis)` is called for `j = 0..width` in pivoted
/// order with `basis = Q (d_j e_j)`, where `d_j` is the `j`-th diagonal entry
/// of the compact factorization. Returns the numerical rank.
pub(in crate::design) fn visit_basis(
    columns: Vec<Vec<f64>>,
    visit: impl FnMut(usize, &[f64]),
) -> Result<usize, Error> {
    let rows = columns.first().map_or(0, Vec::len);
    let mut qr = HouseholderQr::factorize(columns, rows, R_QR_TOLERANCE)?;
    qr.visit_basis(visit)?;
    Ok(qr.rank())
}
