//! Polynomial contrast construction with relative-norm column deferral.
//!
//! This is separate from fitting rank checks: changing the order/number of
//! applied reflections changes the basis even when its span is unchanged.
use super::super::numeric_transform::{column_sum, polynomial_qr};
use crate::errors::Error;
use crate::numeric::math::ProfileMath;
pub(super) fn contrasts(n: usize, output: &mut [f64]) -> Result<(), Error> {
    let center = (n as f64 + 1.0) / 2.0;
    let columns = polynomial_qr::try_columns(n, n, |row, degree| {
        let x = row as f64 + 1.0 - center;
        match degree {
            0 => 1.0,
            1 => x,
            2 => x * x,
            _ => x.powf_p(degree as f64),
        }
    })?;
    polynomial_qr::visit_basis(columns, |column, basis| {
        if column > 0 {
            let length = column_sum(basis.iter().map(|v| v * v)).sqrt();
            for (row, &value) in basis.iter().enumerate() {
                output[row * (n - 1) + column - 1] = value / length;
            }
        }
    })?;
    Ok(())
}
