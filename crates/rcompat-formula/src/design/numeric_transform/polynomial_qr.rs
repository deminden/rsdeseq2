//! Householder QR basis for `poly()` and `contr.poly()`.
//!
//! Not available in this revision: the previous implementation was removed
//! and is being replaced by one written from a specification.
use crate::errors::Error;

pub(in crate::design) fn try_columns(
    rows: usize,
    width: usize,
    value: impl Fn(usize, usize) -> f64,
) -> Result<Vec<Vec<f64>>, Error> {
    let _ = (rows, width, value);
    Err(unavailable())
}

// Caller supplies a finite, nonempty column-major problem with width <= rows.
pub(in crate::design) fn visit_basis(
    columns: Vec<Vec<f64>>,
    visit: impl FnMut(usize, &[f64]),
) -> Result<usize, Error> {
    let _ = (columns, visit);
    Err(unavailable())
}

fn unavailable() -> Error {
    Error::UnsupportedFeature {
        feature: "polynomial QR basis (being rewritten)".to_string(),
    }
}
