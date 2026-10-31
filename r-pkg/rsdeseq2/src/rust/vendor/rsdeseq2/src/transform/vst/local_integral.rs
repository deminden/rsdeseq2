use super::*;

#[derive(Clone, Debug)]
pub(super) struct LocalVstIntegral {
    x: Vec<f64>,
    y: Vec<f64>,
    second_derivatives: Vec<f64>,
}

impl LocalVstIntegral {
    pub(super) fn fit(
        max_count: f64,
        trend: &LocalDispersionTrend,
        inverse_size_factor_mean: f64,
    ) -> Result<Self, DeseqError> {
        let grid_len = 1000usize;
        let max_asinh = max_count.asinh();
        let mut grid = Vec::with_capacity(grid_len - 1);
        for step in 1..grid_len {
            let value = (max_asinh * step as f64 / (grid_len as f64 - 1.0)).sinh();
            grid.push(value);
        }
        let dispersions = trend.evaluate_many_allow_missing(&grid)?;
        Self::fit_from_grid(&grid, &dispersions, inverse_size_factor_mean)
    }

    pub(super) fn fit_from_grid(
        grid: &[f64],
        dispersions: &[f64],
        inverse_size_factor_mean: f64,
    ) -> Result<Self, DeseqError> {
        if grid.len() != dispersions.len() || grid.len() < 2 {
            return Err(DeseqError::InvalidDimensions {
                context: "local VST frozen dispersion grid".to_string(),
                expected: grid.len().max(2),
                actual: dispersions.len(),
            });
        }
        for (idx, (q, dispersion)) in grid.iter().zip(dispersions).enumerate() {
            if !q.is_finite() || *q <= 0.0 || (idx > 0 && *q <= grid[idx - 1]) {
                return Err(DeseqError::InvalidDispersion {
                    reason: format!(
                        "local VST frozen grid value at index {idx} must be finite, positive, and increasing"
                    ),
                });
            }
            if !dispersion.is_finite() || *dispersion <= 0.0 {
                return Err(DeseqError::InvalidDispersion {
                    reason: format!(
                        "local VST frozen dispersion at index {idx} must be finite and positive"
                    ),
                });
            }
        }
        let integrand = grid
            .iter()
            .copied()
            .zip(dispersions.iter().copied())
            .map(|(q, dispersion)| {
                let dispersion_q = checked_mul(dispersion, q, 0, "local VST dispersion q")?;
                let dispersion_q2 = checked_mul(dispersion_q, q, 0, "local VST dispersion q2")?;
                let poisson_q = checked_mul(inverse_size_factor_mean, q, 0, "local VST Poisson q")?;
                let variance =
                    checked_add(dispersion_q2, poisson_q, 0, "local VST variance curve")?;
                if !variance.is_finite() || variance <= 0.0 {
                    return Err(DeseqError::InvalidDispersion {
                        reason: "local VST variance curve must be finite and positive".to_string(),
                    });
                }
                Ok(variance.sqrt().recip())
            })
            .collect::<Result<Vec<_>, _>>()?;
        let mut x = Vec::with_capacity(grid.len() - 1);
        let mut y = Vec::with_capacity(grid.len() - 1);
        let mut cumulative = 0.0;
        for idx in 1..grid.len() {
            let width = checked_sub(grid[idx], grid[idx - 1], idx, "local VST integration width")?;
            let height_sum = checked_add(
                integrand[idx],
                integrand[idx - 1],
                idx,
                "local VST integration height sum",
            )?;
            let area = checked_div(
                checked_mul(width, height_sum, idx, "local VST integration area")?,
                2.0,
                idx,
                "local VST integration area",
            )?;
            cumulative = checked_add(cumulative, area, idx, "local VST integration cumulative")?;
            let midpoint_sum = checked_add(
                grid[idx],
                grid[idx - 1],
                idx,
                "local VST integration midpoint sum",
            )?;
            let midpoint = checked_div(midpoint_sum, 2.0, idx, "local VST integration midpoint")?;
            x.push(midpoint.asinh());
            y.push(cumulative);
        }
        if x.is_empty() {
            return Err(DeseqError::InvalidDimensions {
                context: "local VST integration grid".to_string(),
                expected: 2,
                actual: x.len(),
            });
        }
        let second_derivatives = fmm_spline_second_derivatives(&x, &y)?;
        Ok(Self {
            x,
            y,
            second_derivatives,
        })
    }

    pub(super) fn evaluate(&self, target: f64) -> Result<f64, DeseqError> {
        if !target.is_finite() || target < 0.0 {
            return Err(DeseqError::NonFiniteValue {
                context: "local VST interpolation target".to_string(),
                index: None,
                value: target,
            });
        }
        let last = self.x.len() - 1;
        let lower = if target <= self.x[0] {
            0
        } else if target >= self.x[last] {
            last - 1
        } else {
            self.x.partition_point(|value| *value < target) - 1
        };
        let upper = lower + 1;
        let width = checked_sub(
            self.x[upper],
            self.x[lower],
            upper,
            "local VST spline width",
        )?;
        let left = checked_sub(
            self.x[upper],
            target,
            lower,
            "local VST spline left distance",
        )?;
        let right = checked_sub(
            target,
            self.x[lower],
            lower,
            "local VST spline right distance",
        )?;
        let left_cube = checked_mul(
            checked_mul(left, left, lower, "local VST spline left cube")?,
            left,
            lower,
            "local VST spline left cube",
        )?;
        let right_cube = checked_mul(
            checked_mul(right, right, lower, "local VST spline right cube")?,
            right,
            lower,
            "local VST spline right cube",
        )?;
        let width_six = checked_mul(6.0, width, lower, "local VST spline denominator")?;
        let curvature_left = checked_div(
            checked_mul(
                self.second_derivatives[lower],
                left_cube,
                lower,
                "local VST spline curvature",
            )?,
            width_six,
            lower,
            "local VST spline curvature",
        )?;
        let curvature_right = checked_div(
            checked_mul(
                self.second_derivatives[upper],
                right_cube,
                upper,
                "local VST spline curvature",
            )?,
            width_six,
            upper,
            "local VST spline curvature",
        )?;
        let width_squared = checked_mul(width, width, lower, "local VST spline width squared")?;
        let correction_left = checked_div(
            checked_mul(
                self.second_derivatives[lower],
                width_squared,
                lower,
                "local VST spline correction",
            )?,
            6.0,
            lower,
            "local VST spline correction",
        )?;
        let correction_right = checked_div(
            checked_mul(
                self.second_derivatives[upper],
                width_squared,
                upper,
                "local VST spline correction",
            )?,
            6.0,
            upper,
            "local VST spline correction",
        )?;
        let linear_left = checked_div(
            checked_mul(
                checked_sub(
                    self.y[lower],
                    correction_left,
                    lower,
                    "local VST spline linear term",
                )?,
                left,
                lower,
                "local VST spline linear term",
            )?,
            width,
            lower,
            "local VST spline linear term",
        )?;
        let linear_right = checked_div(
            checked_mul(
                checked_sub(
                    self.y[upper],
                    correction_right,
                    upper,
                    "local VST spline linear term",
                )?,
                right,
                upper,
                "local VST spline linear term",
            )?,
            width,
            upper,
            "local VST spline linear term",
        )?;
        finite_value(
            checked_add(
                checked_add(
                    curvature_left,
                    curvature_right,
                    lower,
                    "local VST spline value",
                )?,
                checked_add(linear_left, linear_right, lower, "local VST spline value")?,
                lower,
                "local VST spline value",
            )?,
            None,
            "local VST spline value",
        )
    }
}

/// Compute the second derivatives for R's `splinefun(method = "fmm")` cubic
/// spline. The end equations match the spline's third derivatives to cubic
/// fits through the first and last four knots, as in R's `SplineCoef`.
fn fmm_spline_second_derivatives(x: &[f64], y: &[f64]) -> Result<Vec<f64>, DeseqError> {
    if x.len() != y.len() || x.len() < 4 {
        return Err(DeseqError::InvalidDimensions {
            context: "local VST FMM spline knots".to_string(),
            expected: 4,
            actual: x.len().min(y.len()),
        });
    }
    let n = x.len();
    let mut h = Vec::with_capacity(n - 1);
    let mut slopes = Vec::with_capacity(n - 1);
    for idx in 0..n - 1 {
        let width = checked_sub(x[idx + 1], x[idx], idx, "local VST spline knot width")?;
        if width <= 0.0 {
            return Err(DeseqError::InvalidOptions {
                reason: "local VST spline knots must be strictly increasing".to_string(),
            });
        }
        h.push(width);
        slopes.push(checked_div(
            checked_sub(y[idx + 1], y[idx], idx, "local VST spline secant")?,
            width,
            idx,
            "local VST spline secant",
        )?);
    }
    // R solves for sigma = second_derivative / 6. The endpoint right-hand
    // sides are divided-difference estimates of the cubic third derivatives.
    let mut diagonal = vec![0.0; n];
    let off_diagonal = h.clone();
    let mut rhs = vec![0.0; n];
    diagonal[0] = -h[0];
    diagonal[n - 1] = -h[n - 2];
    for knot in 1..n - 1 {
        diagonal[knot] = 2.0 * (h[knot - 1] + h[knot]);
        rhs[knot] = slopes[knot] - slopes[knot - 1];
    }
    let h0 = h[0];
    let h1 = h[1];
    let h2 = h[2];
    rhs[0] = ((slopes[2] - slopes[1]) / (h1 + h2) - (slopes[1] - slopes[0]) / (h0 + h1)) * h0 * h0
        / (h0 + h1 + h2);
    let last = n - 2;
    let before_last = n - 3;
    let two_before_last = n - 4;
    rhs[n - 1] = -((slopes[last] - slopes[before_last]) / (h[before_last] + h[last])
        - (slopes[before_last] - slopes[two_before_last]) / (h[two_before_last] + h[before_last]))
        * h[last]
        * h[last]
        / (h[two_before_last] + h[before_last] + h[last]);

    for row in 1..n {
        let factor = off_diagonal[row - 1] / diagonal[row - 1];
        diagonal[row] -= factor * off_diagonal[row - 1];
        rhs[row] -= factor * rhs[row - 1];
    }
    if diagonal
        .iter()
        .any(|value| !value.is_finite() || *value == 0.0)
    {
        return Err(DeseqError::InvalidDispersion {
            reason: "local VST FMM spline system is singular".to_string(),
        });
    }
    let mut sigma = vec![0.0; n];
    sigma[n - 1] = rhs[n - 1] / diagonal[n - 1];
    for row in (0..n - 1).rev() {
        sigma[row] = (rhs[row] - off_diagonal[row] * sigma[row + 1]) / diagonal[row];
    }
    let second: Vec<f64> = sigma.into_iter().map(|value| 6.0 * value).collect();
    if second.iter().any(|value| !value.is_finite()) {
        return Err(DeseqError::InvalidDispersion {
            reason: "local VST FMM spline coefficients must be finite".to_string(),
        });
    }
    Ok(second)
}
