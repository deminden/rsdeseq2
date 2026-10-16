// Weighted least squares, QR, covariance and design crossproducts.

fn solve_weighted_least_squares(
    x: &DMatrix<f64>,
    w: &[f64],
    z: &[f64],
    ridge_lambda: &[f64],
    solver: IrlsSolver,
) -> Option<DVector<f64>> {
    match solver {
        IrlsSolver::NormalEquations => {
            let (xtwx, xtwz) = xtwx_and_xtwz(x, w, z, ridge_lambda)?;
            xtwx.lu().solve(&xtwz)
        }
        IrlsSolver::Qr => solve_weighted_least_squares_qr(x, w, z, ridge_lambda),
    }
}

fn solve_weighted_least_squares_qr(
    x: &DMatrix<f64>,
    w: &[f64],
    z: &[f64],
    ridge_lambda: &[f64],
) -> Option<DVector<f64>> {
    let n = x.nrows();
    let p = x.ncols();
    let mut augmented_x = DMatrix::zeros(n + p, p);
    let mut augmented_z = DVector::zeros(n + p);
    for row in 0..n {
        let sqrt_weight = w[row].sqrt();
        if !sqrt_weight.is_finite() {
            return None;
        }
        augmented_z[row] = checked_product2(z[row], sqrt_weight)?;
        for col in 0..p {
            augmented_x[(row, col)] = checked_product2(x[(row, col)], sqrt_weight)?;
        }
    }
    for col in 0..p {
        let sqrt_ridge = ridge_lambda[col].sqrt();
        if !sqrt_ridge.is_finite() {
            return None;
        }
        augmented_x[(n + col, col)] = sqrt_ridge;
    }
    let (q, r) = augmented_x.qr().unpack();
    let rhs = q.transpose() * augmented_z;
    // The augmented matrix has n+p >= p rows; unpack already returns an
    // economic p-by-p R and a p-element transformed response.
    r.lu().solve(&rhs)
}

fn covariance_and_hat_diagonal(
    x: &DMatrix<f64>,
    w: &[f64],
    ridge_lambda: &[f64],
) -> Option<(Vec<f64>, Vec<f64>)> {
    let xtwx = xtwx_without_ridge(x, w)?;
    let mut xtwx_ridge = xtwx.clone();
    for diagonal in 0..x.ncols() {
        xtwx_ridge[(diagonal, diagonal)] =
            checked_sum2(xtwx_ridge[(diagonal, diagonal)], ridge_lambda[diagonal])?;
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
    let mut hat = Vec::with_capacity(x.nrows());
    let mut active = Vec::with_capacity(x.ncols());
    for (sample, &weight) in w.iter().enumerate().take(x.nrows()) {
        let mut value = 0.0;
        if !weight.is_finite() {
            return None;
        }
        fill_active_design_columns(x, sample, &mut active);
        for (left, left_value) in active.iter().copied() {
            for (right, right_value) in active.iter().copied() {
                value = checked_sum2(
                    value,
                    checked_product4(left_value, weight, right_value, inverse[(right, left)])?,
                )?;
            }
        }
        hat.push(value);
    }
    Some((beta_covariance, hat))
}

fn xtwx_and_xtwz(
    x: &DMatrix<f64>,
    w: &[f64],
    z: &[f64],
    ridge_lambda: &[f64],
) -> Option<(DMatrix<f64>, DVector<f64>)> {
    let mut xtwx = DMatrix::zeros(x.ncols(), x.ncols());
    let mut xtwz = DVector::zeros(x.ncols());
    let mut active = Vec::with_capacity(x.ncols());
    for (sample, &weight) in w.iter().enumerate().take(x.nrows()) {
        fill_active_design_columns(x, sample, &mut active);
        for (col, col_value) in active.iter().copied() {
            xtwz[col] = checked_sum2(xtwz[col], checked_product3(col_value, weight, z[sample])?)?;
        }
        for (left_pos, (row, row_value)) in active.iter().copied().enumerate() {
            for (col, col_value) in active.iter().copied().skip(left_pos) {
                xtwx[(row, col)] = checked_sum2(
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
            checked_sum2(xtwx[(diagonal, diagonal)], ridge_lambda[diagonal])?;
    }
    Some((xtwx, xtwz))
}

fn xtwx_without_ridge(x: &DMatrix<f64>, w: &[f64]) -> Option<DMatrix<f64>> {
    let mut xtwx = DMatrix::zeros(x.ncols(), x.ncols());
    let mut active = Vec::with_capacity(x.ncols());
    for (sample, &weight) in w.iter().enumerate().take(x.nrows()) {
        fill_active_design_columns(x, sample, &mut active);
        for (left_pos, (row, row_value)) in active.iter().copied().enumerate() {
            for (col, col_value) in active.iter().copied().skip(left_pos) {
                xtwx[(row, col)] = checked_sum2(
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

fn fill_active_design_columns(x: &DMatrix<f64>, row: usize, active: &mut Vec<(usize, f64)>) {
    active.clear();
    active.extend((0..x.ncols()).filter_map(|col| {
        let value = x[(row, col)];
        (value != 0.0).then_some((col, value))
    }));
}
