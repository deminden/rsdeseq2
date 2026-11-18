use super::*;
use crate::numeric::checked::{self, CheckedExt};
use crate::numeric::math::ProfileMath;

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
        let max_asinh = max_count.asinh_p();
        let mut grid = Vec::with_capacity(grid_len - 1);
        for step in 1..grid_len {
            let value = (max_asinh * step as f64 / (grid_len as f64 - 1.0)).sinh_p();
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
                let dispersion_q = checked::try_mul(dispersion, q)
                    .nonfinite_at("local VST dispersion q", Some(0))?;
                let dispersion_q2 = checked::try_mul(dispersion_q, q)
                    .nonfinite_at("local VST dispersion q2", Some(0))?;
                let poisson_q = checked::try_mul(inverse_size_factor_mean, q)
                    .nonfinite_at("local VST Poisson q", Some(0))?;
                let variance = checked::try_add(dispersion_q2, poisson_q)
                    .nonfinite_at("local VST variance curve", Some(0))?;
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
            let width = checked::try_sub(grid[idx], grid[idx - 1])
                .nonfinite_at("local VST integration width", Some(idx))?;
            let height_sum = checked::try_add(integrand[idx], integrand[idx - 1])
                .nonfinite_at("local VST integration height sum", Some(idx))?;
            let area = checked::try_div_finite_result(
                checked::try_mul(width, height_sum)
                    .nonfinite_at("local VST integration area", Some(idx))?,
                2.0,
            )
            .nonfinite_at("local VST integration area", Some(idx))?;
            cumulative = checked::try_add(cumulative, area)
                .nonfinite_at("local VST integration cumulative", Some(idx))?;
            let midpoint_sum = checked::try_add(grid[idx], grid[idx - 1])
                .nonfinite_at("local VST integration midpoint sum", Some(idx))?;
            let midpoint = checked::try_div_finite_result(midpoint_sum, 2.0)
                .nonfinite_at("local VST integration midpoint", Some(idx))?;
            x.push(midpoint.asinh_p());
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
        let width = checked::try_sub(self.x[upper], self.x[lower])
            .nonfinite_at("local VST spline width", Some(upper))?;
        let left = checked::try_sub(self.x[upper], target)
            .nonfinite_at("local VST spline left distance", Some(lower))?;
        let right = checked::try_sub(target, self.x[lower])
            .nonfinite_at("local VST spline right distance", Some(lower))?;
        let left_cube = checked::try_mul(
            checked::try_mul(left, left).nonfinite_at("local VST spline left cube", Some(lower))?,
            left,
        )
        .nonfinite_at("local VST spline left cube", Some(lower))?;
        let right_cube = checked::try_mul(
            checked::try_mul(right, right)
                .nonfinite_at("local VST spline right cube", Some(lower))?,
            right,
        )
        .nonfinite_at("local VST spline right cube", Some(lower))?;
        let width_six = checked::try_mul(6.0, width)
            .nonfinite_at("local VST spline denominator", Some(lower))?;
        let curvature_left = checked::try_div_finite_result(
            checked::try_mul(self.second_derivatives[lower], left_cube)
                .nonfinite_at("local VST spline curvature", Some(lower))?,
            width_six,
        )
        .nonfinite_at("local VST spline curvature", Some(lower))?;
        let curvature_right = checked::try_div_finite_result(
            checked::try_mul(self.second_derivatives[upper], right_cube)
                .nonfinite_at("local VST spline curvature", Some(upper))?,
            width_six,
        )
        .nonfinite_at("local VST spline curvature", Some(upper))?;
        let width_squared = checked::try_mul(width, width)
            .nonfinite_at("local VST spline width squared", Some(lower))?;
        let correction_left = checked::try_div_finite_result(
            checked::try_mul(self.second_derivatives[lower], width_squared)
                .nonfinite_at("local VST spline correction", Some(lower))?,
            6.0,
        )
        .nonfinite_at("local VST spline correction", Some(lower))?;
        let correction_right = checked::try_div_finite_result(
            checked::try_mul(self.second_derivatives[upper], width_squared)
                .nonfinite_at("local VST spline correction", Some(upper))?,
            6.0,
        )
        .nonfinite_at("local VST spline correction", Some(upper))?;
        let linear_left = checked::try_div_finite_result(
            checked::try_mul(
                checked::try_sub(self.y[lower], correction_left)
                    .nonfinite_at("local VST spline linear term", Some(lower))?,
                left,
            )
            .nonfinite_at("local VST spline linear term", Some(lower))?,
            width,
        )
        .nonfinite_at("local VST spline linear term", Some(lower))?;
        let linear_right = checked::try_div_finite_result(
            checked::try_mul(
                checked::try_sub(self.y[upper], correction_right)
                    .nonfinite_at("local VST spline linear term", Some(upper))?,
                right,
            )
            .nonfinite_at("local VST spline linear term", Some(upper))?,
            width,
        )
        .nonfinite_at("local VST spline linear term", Some(upper))?;
        finite_value(
            checked::try_add(
                checked::try_add(curvature_left, curvature_right)
                    .nonfinite_at("local VST spline value", Some(lower))?,
                checked::try_add(linear_left, linear_right)
                    .nonfinite_at("local VST spline value", Some(lower))?,
            )
            .nonfinite_at("local VST spline value", Some(lower))?,
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
        let width = checked::try_sub(x[idx + 1], x[idx])
            .nonfinite_at("local VST spline knot width", Some(idx))?;
        if width <= 0.0 {
            return Err(DeseqError::InvalidOptions {
                reason: "local VST spline knots must be strictly increasing".to_string(),
            });
        }
        h.push(width);
        slopes.push(
            checked::try_div_finite_result(
                checked::try_sub(y[idx + 1], y[idx])
                    .nonfinite_at("local VST spline secant", Some(idx))?,
                width,
            )
            .nonfinite_at("local VST spline secant", Some(idx))?,
        );
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
