//! Polynomial contrast construction with relative-norm column deferral.
//!
//! This is separate from fitting rank checks: changing the order/number of
//! applied reflections changes the basis even when its span is unchanged.
use super::super::numeric_transform::{column_sum, polynomial_qr};
pub(super) fn contrasts(n: usize) -> Vec<f64> {
    let center = (n as f64 + 1.0) / 2.0;
    let columns = (0..n)
        .map(|degree| {
            (0..n)
                .map(|r| {
                    let x = r as f64 + 1.0 - center;
                    match degree {
                        0 => 1.0,
                        1 => x,
                        2 => x * x,
                        _ => x.powf(degree as f64),
                    }
                })
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let (columns, _) = polynomial_qr::basis(columns);
    let mut output = vec![0.0; n * (n - 1)];
    for column in 1..n {
        let length = column_sum(columns[column].iter().map(|v| v * v)).sqrt();
        for row in 0..n {
            output[row * (n - 1) + column - 1] = columns[column][row] / length;
        }
    }
    output
}
