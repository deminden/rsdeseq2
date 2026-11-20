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
    basis: DMatrix<f64>,
    transposed_basis: DMatrix<f64>,
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
    if scratch.augmented_design.shape() != (m, p) {
        scratch.augmented_design = DMatrix::zeros(m, p);
        scratch.augmented_response = DVector::zeros(m);
        scratch.sqrt_weights.resize(n, 0.0);
        scratch.diagonal.resize(p, 0.0);
        scratch.basis = DMatrix::zeros(m, p);
        scratch.transposed_basis = DMatrix::zeros(p, m);
        scratch.triangular = DMatrix::zeros(p, p);
        scratch.coefficients = DVector::zeros(p);
    }

    let QrScratch {
        augmented_design,
        augmented_response,
        sqrt_weights,
        diagonal,
        basis,
        transposed_basis,
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
    // Replay the economic QR basis in nalgebra's original reflector order.
    basis.fill(0.0);
    basis.fill_diagonal(1.0);
    for (col, signed_norm) in diagonal.iter().copied().enumerate().rev() {
        let axis = augmented_design.view_range(col.., col);
        let reflection = Reflection::new(Unit::new_unchecked(axis), 0.0);
        reflection.reflect_with_sign(
            &mut basis.view_range_mut(col.., col..),
            ComplexField::signum(signed_norm),
        );
    }
    basis.transpose_to(transposed_basis);
    coefficients.gemv(1.0, transposed_basis, augmented_response, 0.0);
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

pub(super) fn covariance_and_hat_diagonal(
    x: &DMatrix<f64>,
    w: &[f64],
    ridge_lambda: &[f64],
    hat: &mut [f64],
) -> Option<Vec<f64>> {
    let (inverse, beta_covariance) = covariance_with_inverse(x, w, ridge_lambda)?;
    visit_hat_diagonal(x, w, &inverse, |sample, value| hat[sample] = value)?;
    Some(beta_covariance)
}

fn visit_hat_diagonal(
    x: &DMatrix<f64>,
    w: &[f64],
    inverse: &DMatrix<f64>,
    mut use_value: impl FnMut(usize, f64),
) -> Option<()> {
    let mut active = Vec::with_capacity(x.ncols());
    for (sample, &weight) in w.iter().enumerate().take(x.nrows()) {
        let mut value = 0.0;
        if !weight.is_finite() {
            return None;
        }
        fill_active_design_columns(x, sample, &mut active);
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

// Optim fallback preserves the original IRLS hat diagonal.
pub(super) fn fallback_covariance(
    x: &DMatrix<f64>,
    w: &[f64],
    ridge_lambda: &[f64],
) -> Option<Vec<f64>> {
    let (inverse, covariance) = covariance_with_inverse(x, w, ridge_lambda)?;
    visit_hat_diagonal(x, w, &inverse, |_, _| {})?;
    Some(covariance)
}

fn covariance_with_inverse(
    x: &DMatrix<f64>,
    w: &[f64],
    ridge_lambda: &[f64],
) -> Option<(DMatrix<f64>, Vec<f64>)> {
    let xtwx = xtwx_without_ridge(x, w)?;
    let mut xtwx_ridge = xtwx.clone();
    for diagonal in 0..x.ncols() {
        xtwx_ridge[(diagonal, diagonal)] =
            checked::add(xtwx_ridge[(diagonal, diagonal)], ridge_lambda[diagonal])?;
    }
    let inverse = xtwx_ridge.try_inverse()?;
    let sigma = &inverse * xtwx * &inverse;
    let mut beta_covariance = Vec::with_capacity(x.ncols() * x.ncols());
    for row in 0..x.ncols() {
        for col in 0..x.ncols() {
            let value = sigma[(row, col)];
            if !value.is_finite() {
                return None;
            }
            beta_covariance.push(value);
        }
    }
    Some((inverse, beta_covariance))
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

pub(super) fn xtwx_without_ridge(x: &DMatrix<f64>, w: &[f64]) -> Option<DMatrix<f64>> {
    let mut xtwx = DMatrix::zeros(x.ncols(), x.ncols());
    let mut active = Vec::with_capacity(x.ncols());
    for (sample, &weight) in w.iter().enumerate().take(x.nrows()) {
        fill_active_design_columns(x, sample, &mut active);
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
    Some(xtwx)
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
