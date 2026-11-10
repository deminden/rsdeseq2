//! Unbounded L-BFGS with the reference strong-Wolfe search and stopping rules.
//! Rust adaptation of the MIT-licensed LBFGSpp iteration/search algorithms:
//! Copyright (C) 2016-2026 Yixuan Qiu and Dirk Toewe. See adjacent LBFGSPP-LICENSE
//! and the repository's THIRD_PARTY_NOTICES.md.
//! No Rcpp, Eigen, R runtime, or foreign callback is used.
use rcompat_lbfgsb::{OptimCounts, OptimResult};
use std::collections::VecDeque;
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}
fn norm(a: &[f64]) -> f64 {
    dot(a, a).sqrt()
}
#[derive(Clone)]
struct Point {
    x: Vec<f64>,
    value: f64,
    gradient: Vec<f64>,
}
impl Point {
    fn finite(&self) -> bool {
        self.value.is_finite() && self.gradient.iter().all(|v| v.is_finite())
    }
}
struct Correction {
    s: Vec<f64>,
    y: Vec<f64>,
    ys: f64,
}
fn quadratic(lo: f64, hi: f64, flo: f64, fhi: f64, glo: f64) -> f64 {
    let fd = fhi - flo;
    let sd = hi - lo;
    let middle = (hi + lo) / 2.0;
    let candidate = (fd * lo - middle * sd * glo) / (fd - sd * glo);
    let edge_distance = (candidate - lo).abs().min((candidate - hi).abs());
    if !candidate.is_finite()
        || candidate <= lo.min(hi)
        || candidate >= lo.max(hi)
        || edge_distance < 0.01 * sd.abs()
    {
        middle
    } else {
        candidate
    }
}
fn search(
    f: &mut impl FnMut(&[f64]) -> (f64, Vec<f64>),
    initial: &Point,
    direction: &[f64],
    mut step: f64,
    evaluations: &mut usize,
) -> (Point, f64, bool) {
    let initial_slope = dot(&initial.gradient, direction);
    if step <= 0.0 || initial_slope > 0.0 || !initial_slope.is_finite() {
        return (initial.clone(), 0.0, false);
    }
    let decrease = 1e-4 * initial_slope;
    let curvature = -0.9 * initial_slope;
    let mut low = 0.0;
    let mut low_slope = initial_slope;
    let mut low_point = initial.clone();
    let mut high;
    let mut high_value;
    let mut count = 0;
    let mut at = |step: f64| {
        let x = initial
            .x
            .iter()
            .zip(direction)
            .map(|(x, d)| x + step * d)
            .collect::<Vec<_>>();
        let (value, gradient) = f(&x);
        *evaluations += 1;
        Point { x, value, gradient }
    };
    loop {
        let point = at(step);
        let slope = dot(&point.gradient, direction);
        if !point.finite()
            || point.value - initial.value > step * decrease
            || (low > 0.0 && point.value >= low_point.value)
        {
            high = step;
            high_value = point.value;
            break;
        }
        if slope.abs() <= curvature {
            return (point, step, true);
        }
        high = low;
        high_value = low_point.value;
        low = step;
        low_slope = slope;
        low_point = point;
        if slope >= 0.0 {
            break;
        }
        count += 1;
        if count >= 100 {
            return (low_point, low, true);
        }
        step *= 2.0;
    }
    loop {
        step = quadratic(low, high, low_point.value, high_value, low_slope);
        let point = at(step);
        let slope = dot(&point.gradient, direction);
        if !point.finite()
            || point.value - initial.value > step * decrease
            || point.value >= low_point.value
        {
            if step == high {
                return (low_point, low, false);
            }
            high = step;
            high_value = point.value;
        } else {
            if slope.abs() <= curvature {
                return (point, step, true);
            }
            if slope * (high - low) >= 0.0 {
                high = low;
                high_value = low_point.value;
            }
            if step == low {
                return (low_point, low, false);
            }
            low = step;
            low_slope = slope;
            low_point = point;
        }
        count += 1;
        if count >= 100 {
            let valid = low > 0.0;
            return (low_point, low, valid);
        }
    }
}
/// Reference nbinom prefit settings: six corrections, 300 iterations,
/// gradient tolerances 1e-8 and one-step relative function tolerance 1e-8.
pub(super) fn minimize(
    mut f: impl FnMut(&[f64]) -> (f64, Vec<f64>),
    start: Vec<f64>,
) -> OptimResult {
    let (value, gradient) = f(&start);
    let mut point = Point {
        x: start,
        value,
        gradient,
    };
    let mut evaluations = 1;
    let mut history: VecDeque<Correction> = VecDeque::new();
    let mut direction = point.gradient.iter().map(|g| -g).collect::<Vec<_>>();
    let mut step = 1.0 / norm(&direction);
    let mut convergence = 0;
    let mut message = "initial gradient convergence";
    for iteration in 0..300 {
        let gradient_norm = norm(&point.gradient);
        if gradient_norm <= 1e-8 || gradient_norm <= 1e-8 * norm(&point.x) {
            message = "gradient convergence";
            break;
        }
        let previous = point.clone();
        let (next, _, success) = search(&mut f, &previous, &direction, step, &mut evaluations);
        point = next;
        if !success {
            convergence = -1;
            message = "strong-Wolfe line search failed";
            break;
        }
        let gradient_norm = norm(&point.gradient);
        if gradient_norm <= 1e-8 || gradient_norm <= 1e-8 * norm(&point.x) {
            message = "gradient convergence";
            break;
        }
        if (previous.value - point.value).abs()
            <= 1e-8 * point.value.abs().max(previous.value.abs()).max(1.0)
        {
            message = "relative objective convergence";
            break;
        }
        if iteration == 299 {
            convergence = 1;
            message = "iteration limit";
            break;
        }
        let s = point
            .x
            .iter()
            .zip(&previous.x)
            .map(|(x, p)| x - p)
            .collect::<Vec<_>>();
        let y = point
            .gradient
            .iter()
            .zip(&previous.gradient)
            .map(|(g, p)| g - p)
            .collect::<Vec<_>>();
        let ys = dot(&s, &y);
        let yy = dot(&y, &y);
        if ys > f64::EPSILON * yy {
            if history.len() == 6 {
                history.pop_front();
            }
            history.push_back(Correction { s, y, ys });
        }
        direction = point.gradient.iter().map(|g| -g).collect();
        let mut alpha = Vec::with_capacity(history.len());
        for c in history.iter().rev() {
            let a = dot(&c.s, &direction) / c.ys;
            alpha.push(a);
            for (d, y) in direction.iter_mut().zip(&c.y) {
                *d -= a * y;
            }
        }
        if let Some(last) = history.back() {
            let theta = dot(&last.y, &last.y) / last.ys;
            for d in &mut direction {
                *d /= theta;
            }
        }
        for (c, a) in history.iter().zip(alpha.into_iter().rev()) {
            let b = dot(&c.y, &direction) / c.ys;
            for (d, s) in direction.iter_mut().zip(&c.s) {
                *d += (a - b) * s;
            }
        }
        step = 1.0;
    }
    OptimResult {
        par: point.x,
        value: point.value,
        counts: OptimCounts {
            function: evaluations,
            gradient: evaluations,
        },
        convergence,
        message: message.into(),
    }
}
