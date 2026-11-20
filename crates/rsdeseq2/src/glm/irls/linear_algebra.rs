use super::*;
// Weighted least squares, QR, covariance and design crossproducts.

#[derive(Default)]
pub(super) struct LeastSquaresScratch {
    qr: QrScratch,
    normal: NormalEquationsScratch,
}

#[derive(Default)]
pub(super) struct QrScratch {
    augmented_design: DMatrix<f64>,
    augmented_response: DVector<f64>,
    sqrt_weights: Vec<f64>,
    diagonal: Vec<f64>,
    basis_column: DVector<f64>,
    triangular: DMatrix<f64>,
    coefficients: DVector<f64>,
}

pub(super) fn solve_weighted_least_squares(
    scratch: &mut LeastSquaresScratch,
    x: &DMatrix<f64>,
    w: &[f64],
    z: &[f64],
    ridge_lambda: &[f64],
    solver: IrlsSolver,
    beta: &mut DVector<f64>,
) -> Option<()> {
    match solver {
        IrlsSolver::NormalEquations => scratch.normal.solve(x, w, z, ridge_lambda, beta),
        IrlsSolver::Qr => {
            solve_weighted_least_squares_qr(&mut scratch.qr, x, w, z, ridge_lambda, beta)
        }
    }
}

pub(super) fn solve_weighted_least_squares_qr(
    scratch: &mut QrScratch,
    x: &DMatrix<f64>,
    w: &[f64],
    z: &[f64],
    ridge_lambda: &[f64],
    beta: &mut DVector<f64>,
) -> Option<()> {
    use nalgebra::{
        ComplexField, Unit, geometry::Reflection, linalg::householder::clear_column_unchecked,
    };

    let n = x.nrows();
    let p = x.ncols();
    let m = n + p;
    // Preserve the R-host reduction order for narrow GLM designs.
    let column_basis = p <= 4 && crate::numeric_profile() != crate::NumericProfile::Portable;
    if scratch.augmented_design.shape() != (m, p) {
        scratch.augmented_design = DMatrix::zeros(m, p);
        scratch.augmented_response = DVector::zeros(m);
        scratch.sqrt_weights.resize(n, 0.0);
        scratch.diagonal.resize(p, 0.0);
        if column_basis {
            scratch.basis_column = DVector::zeros(m);
        }
        scratch.triangular = DMatrix::zeros(p, p);
        scratch.coefficients = DVector::zeros(p);
    }

    let QrScratch {
        augmented_design,
        augmented_response,
        sqrt_weights,
        diagonal,
        basis_column,
        triangular,
        coefficients,
    } = scratch;

    for (row, sqrt_weight) in sqrt_weights.iter_mut().enumerate() {
        *sqrt_weight = w[row].sqrt();
        if !sqrt_weight.is_finite() {
            return None;
        }
        augmented_response[row] = checked::mul(z[row], *sqrt_weight)?;
    }

    for (col, (mut column, original)) in augmented_design
        .column_iter_mut()
        .zip(x.column_iter())
        .enumerate()
    {
        let column = column.as_mut_slice();
        for ((target, design_value), sqrt_weight) in column[..n]
            .iter_mut()
            .zip(original.as_slice())
            .zip(&*sqrt_weights)
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

    for (col, signed_norm) in diagonal.iter_mut().enumerate() {
        *signed_norm = clear_column_unchecked(augmented_design, col, 0, None);
    }
    if column_basis {
        // Each coefficient uses the original Q-transpose reduction order.
        for (col, coefficient) in coefficients.iter_mut().enumerate() {
            basis_column.fill(0.0);
            basis_column[col] = 1.0;
            for (step, signed_norm) in diagonal[..=col].iter().copied().enumerate().rev() {
                let axis = augmented_design.view_range(step.., step);
                let reflection = Reflection::new(Unit::new_unchecked(axis), 0.0);
                reflection.reflect_with_sign(
                    &mut basis_column.rows_mut(step, m - step),
                    ComplexField::signum(signed_norm),
                );
            }
            *coefficient = basis_column
                .as_slice()
                .iter()
                .zip(augmented_response.as_slice())
                .map(|(basis, response)| basis * response)
                .reduce(|sum, product| product + sum)
                .unwrap_or(0.0);
        }
    } else {
        // The augmented ridge rows contribute zero to the response.
        augmented_response.as_mut_slice()[n..].fill(0.0);
        // Apply Q^T to the response in reflector order without forming Q.
        for (col, signed_norm) in diagonal.iter().copied().enumerate() {
            let axis = augmented_design.view_range(col.., col);
            let reflection = Reflection::new(Unit::new_unchecked(axis), 0.0);
            reflection.reflect_with_sign(
                &mut augmented_response.rows_mut(col, m - col),
                ComplexField::signum(signed_norm),
            );
        }
        coefficients.copy_from(&augmented_response.rows(0, p));
    }
    if !solve_qr_triangle(augmented_design, diagonal, triangular, coefficients) {
        return None;
    }
    std::mem::swap(coefficients, beta);
    Some(())
}

pub(super) fn solve_qr_triangle(
    augmented_design: &DMatrix<f64>,
    diagonal: &[f64],
    triangular: &mut DMatrix<f64>,
    coefficients: &mut DVector<f64>,
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
        (!needs_lower_solve || triangular.solve_lower_triangular_with_diag_mut(coefficients, 1.0))
            && triangular.solve_upper_triangular_mut(coefficients)
    } else {
        for col in 0..p {
            for row in 0..col {
                triangular[(row, col)] = augmented_design[(row, col)];
            }
        }
        triangular.clone().lu().solve_mut(coefficients)
    }
}

#[derive(Default)]
pub(super) struct CovarianceScratch {
    crossproduct: DMatrix<f64>,
    inverse: DMatrix<f64>,
    product: DMatrix<f64>,
    covariance: DMatrix<f64>,
    active: Vec<(usize, f64)>,
}

pub(super) fn covariance_and_hat_diagonal<'a>(
    scratch: &'a mut CovarianceScratch,
    x: &DMatrix<f64>,
    w: &[f64],
    ridge_lambda: &[f64],
    hat: Option<&mut [f64]>,
    max_design: Option<f64>,
) -> Option<&'a DMatrix<f64>> {
    let p = x.ncols();
    if scratch.crossproduct.shape() != (p, p) {
        scratch.crossproduct = DMatrix::zeros(p, p);
        scratch.inverse = DMatrix::zeros(p, p);
        scratch.product = DMatrix::zeros(p, p);
        scratch.covariance = DMatrix::zeros(p, p);
        scratch.active = Vec::with_capacity(p);
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
    fill_xtwx(x, w, crossproduct, active)?;
    inverse.copy_from(crossproduct);
    for diagonal in 0..p {
        inverse[(diagonal, diagonal)] =
            checked::add(inverse[(diagonal, diagonal)], ridge_lambda[diagonal])?;
    }
    if !inverse.try_inverse_mut() {
        return None;
    }
    product.gemm(1.0, inverse, crossproduct, 0.0);
    covariance.gemm(1.0, product, inverse, 0.0);
    if !covariance.iter().all(|value| value.is_finite()) {
        return None;
    }
    if let Some(hat) = hat {
        visit_hat_diagonal(x, w, inverse, active, |sample, value| hat[sample] = value)?;
    } else {
        validate_unretained_hat(x, w, inverse, max_design, active)?;
    }
    Some(covariance)
}

pub(super) fn finite_abs_max(values: &[f64]) -> Option<f64> {
    values.iter().try_fold(0.0_f64, |max, value| {
        value.is_finite().then_some(max.max(value.abs()))
    })
}

pub(super) fn validate_unretained_hat(
    x: &DMatrix<f64>,
    w: &[f64],
    inverse: &DMatrix<f64>,
    max_design: Option<f64>,
    active: &mut Vec<(usize, f64)>,
) -> Option<()> {
    let bounds = || {
        let max_design = max_design?;
        let left = checked::mul(max_design, finite_abs_max(w)?)?;
        let right = checked::mul(max_design, finite_abs_max(inverse.as_slice())?)?;
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
    x: &DMatrix<f64>,
    w: &[f64],
    inverse: &DMatrix<f64>,
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
    crossproduct: DMatrix<f64>,
    response: DVector<f64>,
    active: Vec<(usize, f64)>,
}

impl NormalEquationsScratch {
    fn solve(
        &mut self,
        x: &DMatrix<f64>,
        w: &[f64],
        z: &[f64],
        ridge_lambda: &[f64],
        beta: &mut DVector<f64>,
    ) -> Option<()> {
        let p = x.ncols();
        if self.crossproduct.shape() != (p, p) {
            self.crossproduct = DMatrix::zeros(p, p);
            self.response = DVector::zeros(p);
            self.active = Vec::with_capacity(p);
        } else {
            self.crossproduct.fill(0.0);
            self.response.fill(0.0);
        }
        fill_normal_equations(x, w, z, ridge_lambda, self)?;
        if !solve_lu_in_place(&mut self.crossproduct, &mut self.response) {
            return None;
        }
        std::mem::swap(&mut self.response, beta);
        Some(())
    }
}

fn fill_normal_equations(
    x: &DMatrix<f64>,
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

fn solve_lu_in_place(matrix: &mut DMatrix<f64>, response: &mut DVector<f64>) -> bool {
    use nalgebra::linalg::{gauss_step, gauss_step_swap};

    let p = matrix.ncols();
    for col in 0..p {
        let pivot = matrix.view_range(col.., col).icamax() + col;
        let value = matrix[(pivot, col)];
        if value == 0.0 {
            continue;
        }
        if pivot == col {
            gauss_step(matrix, value, col);
        } else {
            // Apply each LU permutation to the right-hand side immediately.
            response.swap_rows(col, pivot);
            matrix.columns_range_mut(..col).swap_rows(col, pivot);
            gauss_step_swap(matrix, value, col, pivot);
        }
    }
    p == 0
        || (matrix.solve_lower_triangular_with_diag_mut(response, 1.0)
            && matrix.solve_upper_triangular_mut(response))
}

#[cfg(test)]
pub(super) fn xtwx_without_ridge(x: &DMatrix<f64>, w: &[f64]) -> Option<DMatrix<f64>> {
    let mut xtwx = DMatrix::zeros(x.ncols(), x.ncols());
    let mut active = Vec::with_capacity(x.ncols());
    fill_xtwx(x, w, &mut xtwx, &mut active)?;
    Some(xtwx)
}

fn fill_xtwx(
    x: &DMatrix<f64>,
    w: &[f64],
    xtwx: &mut DMatrix<f64>,
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
    x: &DMatrix<f64>,
    row: usize,
    active: &mut Vec<(usize, f64)>,
) {
    active.clear();
    active.extend((0..x.ncols()).filter_map(|col| {
        let value = x[(row, col)];
        (value != 0.0).then_some((col, value))
    }));
}
