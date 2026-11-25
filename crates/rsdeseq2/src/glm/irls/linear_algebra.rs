use super::*;
use rcompat_numeric::math::dense::householder::{reduce_column, reflect_with_sign};
use rcompat_numeric::math::dense::triangular::{solve_lower_unit_in_place, solve_upper_in_place};
use rcompat_numeric::math::dense::{lu_solve_in_place, matmul_to, try_inverse_in_place};
// Weighted least squares, QR, covariance and design crossproducts.

#[derive(Default)]
pub(super) struct LeastSquaresScratch {
    qr: QrScratch,
    normal: NormalEquationsScratch,
}

#[derive(Default)]
pub(super) struct QrScratch {
    augmented_design: DenseMatrix,
    augmented_response: Vec<f64>,
    sqrt_weights: Vec<f64>,
    diagonal: Vec<f64>,
    basis_column: Vec<f64>,
    triangular: DenseMatrix,
    coefficients: Vec<f64>,
}

pub(super) fn solve_weighted_least_squares(
    scratch: &mut LeastSquaresScratch,
    x: &DenseMatrix,
    w: &[f64],
    z: &[f64],
    ridge_lambda: &[f64],
    solver: IrlsSolver,
    beta: &mut Vec<f64>,
) -> Result<bool, DeseqError> {
    match solver {
        IrlsSolver::NormalEquations => scratch.normal.solve(x, w, z, ridge_lambda, beta),
        IrlsSolver::Qr => {
            solve_weighted_least_squares_qr(&mut scratch.qr, x, w, z, ridge_lambda, beta)
        }
    }
}

pub(super) fn solve_weighted_least_squares_qr(
    scratch: &mut QrScratch,
    x: &DenseMatrix,
    w: &[f64],
    z: &[f64],
    ridge_lambda: &[f64],
    beta: &mut Vec<f64>,
) -> Result<bool, DeseqError> {
    let n = x.nrows();
    let p = x.ncols();
    let m = n
        .checked_add(p)
        .ok_or_else(|| invalid_dimensions("IRLS QR augmented rows", usize::MAX - p, n))?;
    // Preserve the R-host reduction order for narrow GLM designs.
    let column_basis = p <= 4 && crate::numeric_profile() != crate::NumericProfile::Portable;
    if scratch.augmented_design.shape() != (m, p) {
        *scratch = QrScratch {
            augmented_design: try_zero_matrix(m, p, "IRLS QR augmented design")?,
            augmented_response: try_zero_vector(m, "IRLS QR augmented response")?,
            sqrt_weights: try_filled_vec(n, 0.0, "IRLS QR square root weights")?,
            diagonal: try_filled_vec(p, 0.0, "IRLS QR diagonal")?,
            basis_column: if column_basis {
                try_zero_vector(m, "IRLS QR basis column")?
            } else {
                Vec::new()
            },
            triangular: try_zero_matrix(p, p, "IRLS QR triangular matrix")?,
            coefficients: try_zero_vector(p, "IRLS QR coefficients")?,
        };
    }
    if fill_augmented_qr_system(scratch, x, w, z, ridge_lambda).is_none() {
        return Ok(false);
    }
    let QrScratch {
        augmented_design,
        augmented_response,
        diagonal,
        basis_column,
        triangular,
        coefficients,
        ..
    } = scratch;

    for (col, signed_norm) in diagonal.iter_mut().enumerate() {
        *signed_norm = reduce_column(augmented_design, col);
    }
    if column_basis {
        // Each coefficient uses the original Q-transpose reduction order.
        for (col, coefficient) in coefficients.iter_mut().enumerate() {
            basis_column.fill(0.0);
            basis_column[col] = 1.0;
            for (step, signed_norm) in diagonal[..=col].iter().copied().enumerate().rev() {
                reflect_with_sign(
                    &augmented_design.column(step)[step..],
                    &mut basis_column[step..],
                    signed_norm.signum(),
                );
            }
            *coefficient = basis_column
                .iter()
                .zip(augmented_response.iter())
                .map(|(basis, response)| basis * response)
                .reduce(|sum, product| product + sum)
                .unwrap_or(0.0);
        }
    } else {
        // The augmented ridge rows contribute zero to the response.
        augmented_response[n..].fill(0.0);
        // Apply Q^T to the response in reflector order without forming Q.
        for (col, signed_norm) in diagonal.iter().copied().enumerate() {
            reflect_with_sign(
                &augmented_design.column(col)[col..],
                &mut augmented_response[col..],
                signed_norm.signum(),
            );
        }
        coefficients.copy_from_slice(&augmented_response[..p]);
    }
    if !solve_qr_triangle(augmented_design, diagonal, triangular, coefficients) {
        return Ok(false);
    }
    std::mem::swap(coefficients, beta);
    Ok(true)
}

fn fill_augmented_qr_system(
    scratch: &mut QrScratch,
    x: &DenseMatrix,
    w: &[f64],
    z: &[f64],
    ridge_lambda: &[f64],
) -> Option<()> {
    let n = x.nrows();
    let QrScratch {
        augmented_design,
        augmented_response,
        sqrt_weights,
        ..
    } = scratch;
    for (row, sqrt_weight) in sqrt_weights.iter_mut().enumerate() {
        *sqrt_weight = w[row].sqrt();
        if !sqrt_weight.is_finite() {
            return None;
        }
        augmented_response[row] = checked::mul(z[row], *sqrt_weight)?;
    }

    for col in 0..x.ncols() {
        let original = x.column(col);
        let column = augmented_design.column_mut(col);
        for ((target, design_value), sqrt_weight) in
            column[..n].iter_mut().zip(original).zip(&*sqrt_weights)
        {
            *target = checked::mul(*design_value, *sqrt_weight)?;
        }
        let sqrt_ridge = ridge_lambda[col].sqrt();
        if !sqrt_ridge.is_finite() {
            return None;
        }
        column[n..].fill(0.0);
        column[n + col] = sqrt_ridge;
    }

    Some(())
}

pub(super) fn solve_qr_triangle(
    augmented_design: &DenseMatrix,
    diagonal: &[f64],
    triangular: &mut DenseMatrix,
    coefficients: &mut [f64],
) -> bool {
    let p = diagonal.len();
    let mut regular = true;
    for col in 0..p {
        let mut negative_above = false;
        for row in 0..col {
            let value = augmented_design[(row, col)];
            // Elimination with a zero lower factor can only clear the sign
            // of a zero entry after a negative entry in the same column.
            triangular[(row, col)] = if value == 0.0 && negative_above {
                0.0
            } else {
                value
            };
            negative_above |= value.is_sign_negative();
            regular &= value.is_finite();
        }
        triangular[(col, col)] = diagonal[col].abs();
        let pivot = triangular[(col, col)];
        // LU scales each lower entry by 1/pivot. A subnormal pivot can
        // overflow that reciprocal even when the pivot itself is finite.
        regular &= pivot == 0.0 || pivot.is_normal();
    }

    if regular {
        let needs_lower_solve = coefficients
            .iter()
            .any(|value| !value.is_finite() || value.to_bits() == (-0.0_f64).to_bits());
        if needs_lower_solve {
            solve_lower_unit_in_place(triangular, coefficients);
        }
        solve_upper_in_place(triangular, coefficients)
    } else {
        for col in 0..p {
            for row in 0..col {
                triangular[(row, col)] = augmented_design[(row, col)];
            }
        }
        let solved = lu_solve_in_place(triangular, coefficients);
        // The next QR solve overwrites the upper triangle only.
        for col in 0..p {
            triangular.column_mut(col)[col + 1..].fill(0.0);
        }
        solved
    }
}

#[derive(Default)]
pub(super) struct CovarianceScratch {
    crossproduct: DenseMatrix,
    inverse: DenseMatrix,
    product: DenseMatrix,
    covariance: DenseMatrix,
    active: Vec<(usize, f64)>,
}

pub(super) enum FitCovariance<'a> {
    Full(&'a DenseMatrix),
    PositiveDiagonal,
}

pub(super) fn covariance_and_hat_diagonal<'a>(
    scratch: &'a mut CovarianceScratch,
    x: &DenseMatrix,
    w: &[f64],
    ridge_lambda: &[f64],
    hat: Option<&mut [f64]>,
    max_design: Option<f64>,
    retain_covariance: bool,
) -> Result<Option<FitCovariance<'a>>, DeseqError> {
    let p = x.ncols();
    if scratch.crossproduct.shape() != (p, p) {
        *scratch = CovarianceScratch {
            crossproduct: try_zero_matrix(p, p, "IRLS covariance crossproduct")?,
            inverse: try_zero_matrix(p, p, "IRLS covariance inverse")?,
            product: try_zero_matrix(p, p, "IRLS covariance intermediate")?,
            active: try_alloc_matrix(p, 1, "IRLS covariance active columns")?,
            ..CovarianceScratch::default()
        };
    } else {
        scratch.crossproduct.fill(0.0);
    }
    let CovarianceScratch {
        crossproduct,
        inverse,
        product,
        covariance,
        active,
    } = scratch;
    if fill_xtwx(x, w, crossproduct, active).is_none()
        || invert_regularized_crossproduct(crossproduct, ridge_lambda, inverse, product).is_none()
    {
        return Ok(None);
    }
    let max_inverse = if !retain_covariance || hat.is_none() {
        finite_abs_max(inverse.as_slice())
    } else {
        None
    };
    if matmul_to(inverse, crossproduct, product).is_err() {
        return Ok(None);
    }
    let fitted = if !retain_covariance
        && max_inverse.is_some_and(|max| covariance_diagonal_is_positive(product, inverse, max))
    {
        FitCovariance::PositiveDiagonal
    } else {
        if covariance.shape() != (p, p) {
            *covariance = try_zero_matrix(p, p, "IRLS covariance")?;
        }
        matmul_to(product, inverse, covariance)?;
        if !covariance.iter().all(|value| value.is_finite()) {
            return Ok(None);
        }
        FitCovariance::Full(covariance)
    };
    let valid_hat = if let Some(hat) = hat {
        visit_hat_diagonal(x, w, inverse, active, |sample, value| hat[sample] = value)
    } else {
        validate_unretained_hat(x, w, inverse, max_design, max_inverse, active)
    };
    Ok(valid_hat.map(|()| fitted))
}

fn invert_regularized_crossproduct(
    crossproduct: &DenseMatrix,
    ridge_lambda: &[f64],
    inverse: &mut DenseMatrix,
    product: &mut DenseMatrix,
) -> Option<()> {
    let p = crossproduct.ncols();
    // Once the inverse is ready, the covariance multiply overwrites product.
    // Reuse it for the wide pivoted LU that otherwise clones the input.
    let regularized = if p > 4 { &mut *product } else { &mut *inverse };
    regularized.copy_from(crossproduct);
    for (diagonal, &ridge) in ridge_lambda.iter().enumerate() {
        let entry = &mut regularized[(diagonal, diagonal)];
        *entry = checked::add(*entry, ridge)?;
    }
    let invertible = match p {
        // Closed-form adjugate inverse up to four coefficients.
        0..=4 => try_inverse_in_place(inverse).ok()?,
        _ => {
            inverse.fill_with_identity();
            lu_solve_in_place(product, inverse.as_mut_slice())
        }
    };
    invertible.then_some(())
}

pub(super) fn covariance_diagonal_is_positive(
    product: &DenseMatrix,
    inverse: &DenseMatrix,
    max_inverse: f64,
) -> bool {
    let p = product.ncols() as f64;
    let Some(max_product) = finite_abs_max(product.as_slice()) else {
        return false;
    };
    let term = max_product * max_inverse;
    let magnitude = term * p;
    if !term.is_normal() || magnitude > f64::MAX / 16.0 {
        return false;
    }
    // The magnitude bound leaves room for rounding and log2 scaling in all
    // entries. A dot product uses at most p products and p additions. Allow
    // for both reduction orders, block accumulation, and subnormal terms.
    let error_bound = magnitude * (32.0 * p * f64::EPSILON) + p * f64::MIN_POSITIVE;
    (0..product.nrows()).all(|row| {
        let diagonal: f64 = (0..product.ncols())
            .map(|col| product[(row, col)] * inverse[(col, row)])
            .sum();
        diagonal > error_bound
    })
}

pub(super) fn finite_abs_max(values: &[f64]) -> Option<f64> {
    values.iter().try_fold(0.0_f64, |max, value| {
        value.is_finite().then_some(max.max(value.abs()))
    })
}

pub(super) fn validate_unretained_hat(
    x: &DenseMatrix,
    w: &[f64],
    inverse: &DenseMatrix,
    max_design: Option<f64>,
    max_inverse: Option<f64>,
    active: &mut Vec<(usize, f64)>,
) -> Option<()> {
    let bounds = || {
        let max_design = max_design?;
        let left = checked::mul(max_design, finite_abs_max(w)?)?;
        let right = checked::mul(max_design, max_inverse?)?;
        let term = checked::mul(left, right)?;
        checked::mul(term, (x.ncols() as f64).powi(2))
    };
    // Each rounded addition can increase the magnitude by at most twice
    // the term's magnitude. Leave ample room for rounding and fall back
    // to the element-wise checks for extreme inputs.
    if bounds().is_some_and(|bound| bound <= f64::MAX / 4.0) {
        Some(())
    } else {
        visit_hat_diagonal(x, w, inverse, active, |_, _| {})
    }
}

fn visit_hat_diagonal(
    x: &DenseMatrix,
    w: &[f64],
    inverse: &DenseMatrix,
    active: &mut Vec<(usize, f64)>,
    mut use_value: impl FnMut(usize, f64),
) -> Option<()> {
    for (sample, &weight) in w.iter().enumerate().take(x.nrows()) {
        let mut value = 0.0;
        if !weight.is_finite() {
            return None;
        }
        fill_active_design_columns(x, sample, active);
        for (left, left_value) in active.iter().copied() {
            for (right, right_value) in active.iter().copied() {
                value = checked::add(
                    value,
                    checked_product4(left_value, weight, right_value, inverse[(right, left)])?,
                )?;
            }
        }
        use_value(sample, value);
    }
    Some(())
}

#[derive(Default)]
struct NormalEquationsScratch {
    crossproduct: DenseMatrix,
    response: Vec<f64>,
    active: Vec<(usize, f64)>,
}

impl NormalEquationsScratch {
    fn solve(
        &mut self,
        x: &DenseMatrix,
        w: &[f64],
        z: &[f64],
        ridge_lambda: &[f64],
        beta: &mut Vec<f64>,
    ) -> Result<bool, DeseqError> {
        let p = x.ncols();
        if self.crossproduct.shape() != (p, p) {
            *self = Self {
                crossproduct: try_zero_matrix(p, p, "IRLS normal equations crossproduct")?,
                response: try_zero_vector(p, "IRLS normal equations response")?,
                active: try_alloc_matrix(p, 1, "IRLS normal equations active columns")?,
            };
        } else {
            self.crossproduct.fill(0.0);
            self.response.fill(0.0);
        }
        if fill_normal_equations(x, w, z, ridge_lambda, self).is_none()
            || !lu_solve_in_place(&mut self.crossproduct, &mut self.response)
        {
            return Ok(false);
        }
        std::mem::swap(&mut self.response, beta);
        Ok(true)
    }
}

fn fill_normal_equations(
    x: &DenseMatrix,
    w: &[f64],
    z: &[f64],
    ridge_lambda: &[f64],
    scratch: &mut NormalEquationsScratch,
) -> Option<()> {
    let NormalEquationsScratch {
        crossproduct: xtwx,
        response: xtwz,
        active,
    } = scratch;
    for (sample, &weight) in w.iter().enumerate().take(x.nrows()) {
        fill_active_design_columns(x, sample, active);
        for (col, col_value) in active.iter().copied() {
            xtwz[col] = checked::add(xtwz[col], checked_product3(col_value, weight, z[sample])?)?;
        }
        for (left_pos, (row, row_value)) in active.iter().copied().enumerate() {
            for (col, col_value) in active.iter().copied().skip(left_pos) {
                xtwx[(row, col)] = checked::add(
                    xtwx[(row, col)],
                    checked_product3(row_value, weight, col_value)?,
                )?;
                if row != col {
                    xtwx[(col, row)] = xtwx[(row, col)];
                }
            }
        }
    }
    for diagonal in 0..x.ncols() {
        xtwx[(diagonal, diagonal)] =
            checked::add(xtwx[(diagonal, diagonal)], ridge_lambda[diagonal])?;
    }
    Some(())
}

#[cfg(test)]
pub(super) fn xtwx_without_ridge(x: &DenseMatrix, w: &[f64]) -> Option<DenseMatrix> {
    let mut xtwx = DenseMatrix::zeros(x.ncols(), x.ncols()).ok()?;
    let mut active = Vec::with_capacity(x.ncols());
    fill_xtwx(x, w, &mut xtwx, &mut active)?;
    Some(xtwx)
}

fn fill_xtwx(
    x: &DenseMatrix,
    w: &[f64],
    xtwx: &mut DenseMatrix,
    active: &mut Vec<(usize, f64)>,
) -> Option<()> {
    for (sample, &weight) in w.iter().enumerate().take(x.nrows()) {
        fill_active_design_columns(x, sample, active);
        for (left_pos, (row, row_value)) in active.iter().copied().enumerate() {
            for (col, col_value) in active.iter().copied().skip(left_pos) {
                xtwx[(row, col)] = checked::add(
                    xtwx[(row, col)],
                    checked_product3(row_value, weight, col_value)?,
                )?;
                if row != col {
                    xtwx[(col, row)] = xtwx[(row, col)];
                }
            }
        }
    }
    Some(())
}

pub(super) fn fill_active_design_columns(
    x: &DenseMatrix,
    row: usize,
    active: &mut Vec<(usize, f64)>,
) {
    active.clear();
    active.extend((0..x.ncols()).filter_map(|col| {
        let value = x[(row, col)];
        (value != 0.0).then_some((col, value))
    }));
}
