use super::{invalid, polynomial_qr, scalar_mean};
use crate::{errors::DeseqError, matrix::RowMajorMatrix};

/// Recurrence metadata for evaluating the same polynomial basis on new values.
#[derive(Clone, Debug, PartialEq)]
pub struct PolynomialCoefficients {
    /// Centering recurrence coefficients, one per nonconstant degree.
    pub alpha: Vec<f64>,
    /// Leading one followed by squared norms of constant through final degree.
    pub norm2: Vec<f64>,
}

/// Polynomial values and reusable recurrence metadata.
#[derive(Clone, Debug, PartialEq)]
pub struct PolynomialBasis {
    /// Sample-by-degree values; the constant column is omitted.
    pub values: RowMajorMatrix<f64>,
    /// Absent for raw powers.
    pub coefficients: Option<PolynomialCoefficients>,
}

/// Construct raw or orthogonal polynomial model columns. Orthogonal fitting
/// centers inputs before QR and uses relative column deferral at 1e-7, matching
/// the polynomial operation rather than changing the fitting rank policy.
/// Supplied coefficients reuse the recurrence without refitting a QR basis.
pub fn polynomial_basis(
    input: &[f64],
    degree: usize,
    raw: bool,
    coefficients: Option<&PolynomialCoefficients>,
) -> Result<PolynomialBasis, DeseqError> {
    if input.is_empty() || degree == 0 || input.iter().any(|v| !v.is_finite()) {
        return Err(invalid(
            "polynomials require finite values and positive degree",
        ));
    }
    let result = polynomial_expression(input, degree, raw, coefficients)?;
    if result.values.as_slice().iter().any(|v| !v.is_finite()) {
        return Err(invalid("polynomial transform produced non-finite values"));
    }
    Ok(result)
}

pub(in crate::design) fn polynomial_expression(
    input: &[f64],
    degree: usize,
    raw: bool,
    coefficients: Option<&PolynomialCoefficients>,
) -> Result<PolynomialBasis, DeseqError> {
    if degree == 0 || (!raw && coefficients.is_none() && input.iter().any(|v| !v.is_finite())) {
        return Err(invalid(
            "orthogonal polynomials require nonmissing finite input and positive degree",
        ));
    }
    let width = degree
        .checked_add(1)
        .ok_or_else(|| invalid("polynomial degree overflow"))?;
    let size = input
        .len()
        .checked_mul(degree)
        .ok_or_else(|| invalid("polynomial shape overflow"))?;
    if size > isize::MAX as usize / std::mem::size_of::<f64>() {
        return Err(invalid("polynomial allocation shape overflow"));
    }
    if !raw && coefficients.is_none() && degree >= input.len() {
        return Err(invalid(
            "polynomial degree must be less than the number of unique values",
        ));
    }
    let mut output = vec![0.0; size];
    if raw {
        for (row, &x) in input.iter().enumerate() {
            for column in 0..degree {
                output[row * degree + column] = power(x, column + 1);
            }
        }
        return finish(input.len(), degree, output, None);
    }
    if let Some(coefs) = coefficients {
        if coefs.alpha.len() < degree
            || coefs.norm2.len() != width + 1
            || coefs.alpha[..degree].iter().any(|x| !x.is_finite())
            || coefs.norm2[1..=width]
                .iter()
                .any(|x| !x.is_finite() || *x <= 0.0)
        {
            return Err(invalid(
                "polynomial recurrence coefficients have invalid dimensions or values",
            ));
        }
        for (row, &x) in input.iter().enumerate() {
            let mut previous = 1.0;
            let mut current = x - coefs.alpha[0];
            output[row * degree] = current / coefs.norm2[2].sqrt();
            for column in 1..degree {
                let next = (x - coefs.alpha[column]) * current
                    - (coefs.norm2[column + 1] / coefs.norm2[column]) * previous;
                output[row * degree + column] = next / coefs.norm2[column + 2].sqrt();
                previous = current;
                current = next;
            }
        }
        return finish(input.len(), degree, output, Some(coefs.clone()));
    }
    let mut unique = input.to_vec();
    unique.sort_unstable_by(f64::total_cmp);
    unique.dedup_by(|a, b| *a == *b);
    if degree >= unique.len() {
        return Err(invalid(
            "polynomial degree must be less than the number of unique values",
        ));
    }
    let mean = scalar_mean(input);
    let centered = input.iter().map(|x| x - mean).collect::<Vec<_>>();
    let columns = (0..width)
        .map(|d| centered.iter().map(|&x| power(x, d)).collect::<Vec<_>>())
        .collect::<Vec<_>>();
    if columns.iter().flatten().any(|v| !v.is_finite()) {
        return Err(invalid("polynomial powers produced non-finite values"));
    }
    let (columns, rank) = polynomial_qr::basis(columns);
    if rank < degree {
        return Err(invalid(
            "polynomial degree is too high for the numerical rank",
        ));
    }
    let norms = columns
        .iter()
        .map(|c| c.iter().map(|x| x * x).sum::<f64>())
        .collect::<Vec<_>>();
    let alpha = (0..degree)
        .map(|c| {
            centered
                .iter()
                .zip(&columns[c])
                .map(|(x, z)| x * (z * z))
                .sum::<f64>()
                / norms[c]
                + mean
        })
        .collect();
    let mut norm2 = Vec::with_capacity(width + 1);
    norm2.push(1.0);
    norm2.extend_from_slice(&norms);
    for column in 1..width {
        let divisor = norms[column].sqrt();
        for row in 0..input.len() {
            output[row * degree + column - 1] = columns[column][row] / divisor;
        }
    }
    finish(
        input.len(),
        degree,
        output,
        Some(PolynomialCoefficients { alpha, norm2 }),
    )
}

fn power(x: f64, exponent: usize) -> f64 {
    if x.is_nan() {
        return if exponent == 0 {
            1.0
        } else {
            x + exponent as f64
        };
    }
    if x == 0.0 && exponent != 0 {
        return 0.0;
    }
    if (-11.0..=11.0).contains(&x) {
        match exponent {
            3 => return x * x * x,
            4 => return x * x * x * x,
            _ => {}
        }
    }
    match exponent {
        0 => 1.0,
        1 => x,
        2 => x * x,
        _ => x.powf(exponent as f64),
    }
}
fn finish(
    rows: usize,
    columns: usize,
    values: Vec<f64>,
    coefficients: Option<PolynomialCoefficients>,
) -> Result<PolynomialBasis, DeseqError> {
    Ok(PolynomialBasis {
        values: RowMajorMatrix::from_row_major(rows, columns, values)?,
        coefficients,
    })
}
