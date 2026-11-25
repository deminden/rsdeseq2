//! Standalone inverse-BFGS with safeguarded Armijo steps and metric restarts.
//! The update is the standard rank-two inverse-Hessian formula; no R code runs.
//! Reference stopping policy includes periodic metric resets, a fresh-metric
//! convergence probe retaining the previous score, and an iteration budget
//! counted by accepted gradient updates, including the initial gradient.
//! These policies are independently exercised by the saturated BB profiles.
use super::try_dense_matrix;
use crate::{Error, matrix::try_alloc_matrix};
use rcompat_lbfgsb::{OptimCounts, OptimResult};
use rcompat_numeric::math::dense::{DenseMatrix, dot, gemv_to};
pub(super) fn minimize(
    mut evaluate: impl FnMut(&[f64]) -> (f64, Vec<f64>),
    mut x: Vec<f64>,
) -> Result<OptimResult, Error> {
    let p = x.len();
    let (mut value, mut gradient) = evaluate(&x);
    let mut count = 1;
    let mut inverse = try_dense_matrix(p, p, "posterior BFGS metric", |i, j| {
        if i == j { 1.0 } else { 0.0 }
    })?;
    let workspace = |context| -> Result<Vec<f64>, Error> {
        let mut values = try_alloc_matrix(p, 1, context)?;
        values.resize(p, 0.0);
        Ok(values)
    };
    let mut direction = workspace("posterior BFGS direction")?;
    let mut s = workspace("posterior BFGS step")?;
    let mut y = workspace("posterior BFGS gradient change")?;
    let mut hy = workspace("posterior BFGS metric product")?;
    let mut fresh = true;
    let mut metric_age = 0;
    let mut convergence = 0;
    let mut message = "BFGS relative objective convergence";
    let mut accepted_gradients = 1;
    while accepted_gradients < 100 {
        if metric_age > 2 * p {
            inverse.fill_with_identity();
            fresh = true;
            metric_age = 0;
        }
        // direction = -H g; negating the column gaxpy rounds identically.
        gemv_to(&inverse, &gradient, &mut direction)?;
        for d in direction.iter_mut() {
            *d = -*d;
        }
        let mut slope = dot(&direction, &gradient);
        if !slope.is_finite() || slope >= 0.0 {
            inverse.fill_with_identity();
            for (d, &value) in direction.iter_mut().zip(gradient.iter()) {
                *d = -value;
            }
            slope = -dot(&gradient, &gradient);
            fresh = true;
            metric_age = 0;
        }
        if slope == 0.0 {
            break;
        }
        let mut step = 1.0;
        let mut accepted = None;
        for _ in 0..1024 {
            let next = x
                .iter()
                .zip(direction.iter())
                .map(|(x, d)| x + step * d)
                .collect::<Vec<_>>();
            if next == x {
                break;
            }
            let (v, gr) = evaluate(&next);
            count += 1;
            if v.is_finite() && gr.iter().all(|v| v.is_finite()) && v <= value + 1e-4 * step * slope
            {
                accepted = Some((next, v, gr));
                break;
            }
            step *= 0.2;
        }
        let Some((next, next_value, next_gradient)) = accepted else {
            if fresh {
                message = "BFGS no representable descent step";
                break;
            }
            inverse.fill_with_identity();
            fresh = true;
            metric_age = 0;
            continue;
        };
        let small =
            (value - next_value).abs() <= f64::EPSILON.sqrt() * (value.abs() + f64::EPSILON.sqrt());
        for (step, (&a, &b)) in s.iter_mut().zip(next.iter().zip(&x)) {
            *step = a - b;
        }
        for (change, (&a, &b)) in y.iter_mut().zip(next_gradient.iter().zip(&gradient)) {
            *change = a - b;
        }
        x = next;
        value = next_value;
        if small {
            if fresh {
                break;
            }
            inverse.fill_with_identity();
            fresh = true;
            metric_age = 0;
            continue;
        }
        gradient = next_gradient;
        accepted_gradients += 1;
        metric_age += 1;
        let curvature = dot(&s, &y);
        if curvature > 0.0 && curvature.is_finite() {
            gemv_to(&inverse, &y, &mut hy)?;
            let ratio = (1.0 + dot(&y, &hy) / curvature) / curvature;
            update_inverse_metric(&mut inverse, &s, &hy, ratio, curvature);
            fresh = false;
        } else {
            inverse.fill_with_identity();
            fresh = true;
            metric_age = 0;
        }
        if accepted_gradients == 100 {
            convergence = 1;
            message = "BFGS iteration limit";
        }
    }
    Ok(OptimResult {
        par: x,
        value,
        counts: OptimCounts {
            function: count,
            gradient: count,
        },
        convergence,
        message: message.into(),
    })
}

fn update_inverse_metric(
    inverse: &mut DenseMatrix,
    s: &[f64],
    hy: &[f64],
    ratio: f64,
    curvature: f64,
) {
    for j in 0..inverse.ncols() {
        for (i, entry) in inverse.column_mut(j).iter_mut().enumerate() {
            let outer = s[i] * s[j];
            let mixed = s[i] * hy[j] + hy[i] * s[j];
            *entry += ratio * outer - mixed / curvature;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bfgs_metric_update_matches_the_rank_two_expression() {
        let s = [2.0 / 7.0, -0.0, -3.0 / 7.0, 4.0 / 7.0, 0.0, 1.0];
        let hy = [-0.0, -0.4, 0.0, 0.3, 5.0, -6.0];
        let n = s.len();
        let mut metric = DenseMatrix::identity(n).unwrap();
        let ratio = 0.53;
        let curvature = 0.27;
        // I + (ratio * s sᵀ - (s hyᵀ + hy sᵀ) / curvature), entry by entry.
        let expected = DenseMatrix::try_from_fn(n, n, "test", |i, j| {
            let identity = if i == j { 1.0 } else { 0.0 };
            identity + (ratio * (s[i] * s[j]) - (s[i] * hy[j] + hy[i] * s[j]) / curvature)
        })
        .unwrap();

        update_inverse_metric(&mut metric, &s, &hy, ratio, curvature);
        for (actual, reference) in metric.iter().zip(expected.iter()) {
            assert_eq!(actual.to_bits(), reference.to_bits());
        }
    }
}
