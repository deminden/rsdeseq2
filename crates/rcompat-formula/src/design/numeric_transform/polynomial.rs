use super::{column_sum, invalid_options, polynomial_qr, scalar_mean};
use crate::errors::{ErrorKind, missing_prerequisite};
use crate::numeric::math::ProfileMath;
use crate::{
    errors::Error,
    matrix::{RowMajorMatrix, try_alloc_matrix},
};

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
) -> Result<PolynomialBasis, Error> {
    if input.is_empty() || degree == 0 || input.iter().any(|v| !v.is_finite()) {
        return Err(missing_prerequisite(
            ErrorKind::Options,
            "poly()",
            "finite values and positive degree",
        ));
    }
    let result = polynomial_expression(input, degree, raw, coefficients)?;
    if result.values.as_slice().iter().any(|v| !v.is_finite()) {
        return Err(invalid_options(
            "polynomial transform produced non-finite values",
        ));
    }
    Ok(result)
}

pub(in crate::design) fn polynomial_expression(
    input: &[f64],
    degree: usize,
    raw: bool,
    coefficients: Option<&PolynomialCoefficients>,
) -> Result<PolynomialBasis, Error> {
    if degree == 0 || (!raw && coefficients.is_none() && input.iter().any(|v| !v.is_finite())) {
        return Err(missing_prerequisite(
            ErrorKind::Options,
            "an orthogonal polynomial",
            "nonmissing finite input and positive degree",
        ));
    }
    let width = degree
        .checked_add(1)
        .ok_or_else(|| invalid_options("polynomial degree overflow"))?;
    let size = input
        .len()
        .checked_mul(degree)
        .ok_or_else(|| invalid_options("polynomial shape overflow"))?;
    if size > isize::MAX as usize / std::mem::size_of::<f64>() {
        return Err(invalid_options("polynomial allocation shape overflow"));
    }
    if !raw && coefficients.is_none() && degree >= input.len() {
        return Err(invalid_options(
            "polynomial degree must be less than the number of unique values",
        ));
    }
    let mut output = try_alloc_matrix(input.len(), degree, "polynomial basis values")?;
    output.resize(size, 0.0);
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
            return Err(invalid_options(
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
    let mut centered = try_alloc_matrix(input.len(), 1, "polynomial input workspace")?;
    centered.extend_from_slice(input);
    centered.sort_unstable_by(f64::total_cmp);
    centered.dedup_by(|a, b| *a == *b);
    if degree >= centered.len() {
        return Err(invalid_options(
            "polynomial degree must be less than the number of unique values",
        ));
    }
    let mean = scalar_mean(input);
    centered.clear();
    centered.extend(input.iter().map(|x| x - mean));
    let columns = polynomial_qr::try_columns(input.len(), width, |row, column| {
        power(centered[row], column)
    })?;
    if columns.iter().flatten().any(|v| !v.is_finite()) {
        return Err(invalid_options(
            "polynomial powers produced non-finite values",
        ));
    }
    let mut alpha = try_alloc_matrix(degree, 1, "polynomial recurrence centers")?;
    let mut norm2 = try_alloc_matrix(width + 1, 1, "polynomial recurrence squared norms")?;
    norm2.push(1.0);
    let rank = polynomial_qr::visit_basis(columns, |column, basis| {
        let norm = column_sum(basis.iter().map(|x| x * x));
        norm2.push(norm);
        if column < degree {
            let weighted = centered.iter().zip(basis).map(|(x, z)| x * (z * z));
            alpha.push(column_sum(weighted) / norm + mean);
        }
        if column > 0 {
            let divisor = norm.sqrt();
            for (row, &value) in basis.iter().enumerate() {
                output[row * degree + column - 1] = value / divisor;
            }
        }
    })?;
    if rank < degree {
        return Err(invalid_options(
            "polynomial degree is too high for the numerical rank",
        ));
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
        _ => x.powf_p(exponent as f64),
    }
}
fn finish(
    rows: usize,
    columns: usize,
    values: Vec<f64>,
    coefficients: Option<PolynomialCoefficients>,
) -> Result<PolynomialBasis, Error> {
    Ok(PolynomialBasis {
        values: RowMajorMatrix::from_row_major(rows, columns, values)?,
        coefficients,
    })
}
