//! Standalone inverse-BFGS with safeguarded Armijo steps and metric restarts.
//! The update is the standard rank-two inverse-Hessian formula; no R code runs.
use nalgebra::{DMatrix, DVector};
use rcompat_lbfgsb::{OptimCounts, OptimResult};
pub(super) fn minimize(
    mut evaluate: impl FnMut(&[f64]) -> (f64, Vec<f64>),
    mut x: Vec<f64>,
) -> OptimResult {
    let p = x.len();
    let (mut value, mut gradient) = evaluate(&x);
    let mut count = 1;
    let mut inverse = DMatrix::identity(p, p);
    let mut fresh = true;
    let mut convergence = 0;
    let mut message = "BFGS relative objective convergence";
    for iteration in 0..100 {
        let g = DVector::from_column_slice(&gradient);
        let mut direction = -&inverse * &g;
        let mut slope = direction.dot(&g);
        if !slope.is_finite() || slope >= 0.0 {
            inverse.fill(0.0);
            inverse.set_diagonal(&DVector::from_element(p, 1.0));
            direction = -&g;
            slope = -g.dot(&g);
            fresh = true;
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
            inverse = DMatrix::identity(p, p);
            fresh = true;
            continue;
        };
        let small =
            (value - next_value).abs() <= f64::EPSILON.sqrt() * (value.abs() + f64::EPSILON.sqrt());
        let s = DVector::from_iterator(p, next.iter().zip(&x).map(|(a, b)| a - b));
        let y = DVector::from_iterator(p, next_gradient.iter().zip(&gradient).map(|(a, b)| a - b));
        x = next;
        value = next_value;
        gradient = next_gradient;
        if small {
            if fresh {
                break;
            }
            inverse = DMatrix::identity(p, p);
            fresh = true;
            continue;
        }
        let curvature = s.dot(&y);
        if curvature > 0.0 && curvature.is_finite() {
            let hy = &inverse * &y;
            let ratio = (1.0 + y.dot(&hy) / curvature) / curvature;
            inverse += ratio * (&s * s.transpose())
                - (&s * hy.transpose() + &hy * s.transpose()) / curvature;
            fresh = false;
        } else {
            inverse = DMatrix::identity(p, p);
            fresh = true;
        }
        if iteration == 99 {
            convergence = 1;
            message = "BFGS iteration limit";
        }
    }
    OptimResult {
        par: x,
        value,
        counts: OptimCounts {
            function: count,
            gradient: count,
        },
        convergence,
        message: message.into(),
    }
}
