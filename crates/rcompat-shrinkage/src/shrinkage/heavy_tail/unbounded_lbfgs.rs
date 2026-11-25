//! Unbounded limited-memory BFGS for posterior modes.
//!
//! The iteration is the L-BFGS method of Liu and Nocedal (1989), "On the
//! limited memory BFGS method for large scale optimization", Math.
//! Programming 45:503-528, in the form of Nocedal and Wright (2006),
//! *Numerical Optimization*, 2nd ed., Algorithm 7.5: the search direction is
//! the two-loop recursion of Algorithm 7.4 over the newest six correction
//! pairs, with the initial inverse Hessian scaled by `s'y / y'y` of the
//! newest pair (equation 7.20). The first trial step has unit length; later
//! trial steps are 1.
//!
//! Steps satisfy the strong Wolfe conditions (`c1 = 1e-4`, `c2 = 0.9`) and
//! are found in two phases, after Algorithms 3.5 and 3.6 of the same book.
//! The bracketing phase doubles the step until it fails sufficient decrease,
//! stops decreasing, or meets a non-negative slope; after 100 doublings
//! without a bracket the last step is taken. The zoom phase picks each trial
//! as the minimiser of the quadratic through the value and slope at the
//! better bracket end and the value at the other end (section 3.5), and
//! bisects instead when that minimiser lies within 1% of the bracket width
//! of either end or is not finite. A zoom that finds no acceptable step in
//! 100 trials is a line-search failure. Non-finite objective or gradient
//! values count as failed sufficient decrease.
//!
//! The run stops with `convergence` 0 when `||g|| <= gradient_tolerance *
//! max(1, ||x||)`, when one iteration lowers the objective by less than
//! `function_tolerance * max(1, |f_previous|)`, or when the 300-iteration
//! budget is spent (also reported as 0, with its own message). A
//! line-search failure returns the last finite trial point with
//! `convergence` -1; a non-finite start returns -3.
//!
//! The trial-step formula, the safeguard width, the trial limits, the form
//! of the objective test and the iteration-limit status were fixed by
//! black-box runs (`scripts/provenance/probe_unbounded_lbfgs.R`, output in
//! `docs/provenance/observations/unbounded_lbfgs_probe_r462_macos.txt`).
use rcompat_lbfgsb::{OptimCounts, OptimResult};
use std::collections::VecDeque;

/// Number of stored correction pairs.
const MEMORY: usize = 6;
/// Iteration budget (accepted steps).
const MAX_ITERATIONS: usize = 300;
/// Sufficient-decrease constant `c1` of the Wolfe conditions.
const SUFFICIENT_DECREASE: f64 = 1e-4;
/// Curvature constant `c2` of the strong Wolfe conditions.
const CURVATURE: f64 = 0.9;
/// Trials allowed in each line-search phase.
const MAX_PHASE_TRIALS: usize = 100;
/// Interpolated trials closer than this fraction of the bracket width to
/// either end are replaced by the midpoint.
const SAFEGUARD: f64 = 0.01;

const LINE_SEARCH_FAILED: &str =
    "the line search routine failed, unable to sufficiently decrease the function value";

/// Minimise with six corrections, 300 iterations and tolerances 1e-8.
pub(super) fn minimize(f: impl FnMut(&[f64]) -> (f64, Vec<f64>), start: Vec<f64>) -> OptimResult {
    minimize_with_tolerances(f, start, 1e-8, 1e-8)
}

/// Minimise with the given function and gradient tolerances (see the module
/// documentation for the stopping rules).
pub(super) fn minimize_with_tolerances(
    f: impl FnMut(&[f64]) -> (f64, Vec<f64>),
    start: Vec<f64>,
    function_tolerance: f64,
    gradient_tolerance: f64,
) -> OptimResult {
    let mut objective = Objective { f, evaluations: 0 };
    let mut x = start;
    let Some((mut value, mut gradient)) = objective.finite_at(&x) else {
        return objective.result(
            x,
            f64::NAN,
            -3,
            "the objective or its gradient is not finite at the start",
        );
    };
    let gradient_converged =
        |g: &[f64], x: &[f64]| norm(g) <= gradient_tolerance * norm(x).max(1.0);
    if gradient_converged(&gradient, &x) {
        return objective.result(x, value, 0, "gradient norm below tolerance");
    }
    let mut memory = Memory::default();
    let mut direction: Vec<f64> = gradient.iter().map(|g| -g).collect();
    let mut step = 1.0 / norm(&gradient);
    for _ in 0..MAX_ITERATIONS {
        let mut slope = dot(&gradient, &direction);
        if slope.is_nan() || slope >= 0.0 {
            // Not a descent direction: drop the curvature pairs and restart
            // from steepest descent with a unit-length first move.
            memory.clear();
            direction = gradient.iter().map(|g| -g).collect();
            slope = -dot(&gradient, &gradient);
            step = 1.0 / norm(&gradient);
        }
        let trial = match line_search(&mut objective, &x, value, slope, &direction, step) {
            Ok(trial) => trial,
            Err(last) => {
                let (x, value) = last.unwrap_or((x, value));
                return objective.result(x, value, -1, LINE_SEARCH_FAILED);
            }
        };
        let s: Vec<f64> = trial.x.iter().zip(&x).map(|(a, b)| a - b).collect();
        let y: Vec<f64> = trial
            .gradient
            .iter()
            .zip(&gradient)
            .map(|(a, b)| a - b)
            .collect();
        memory.push(s, y);
        let previous = value;
        x = trial.x;
        value = trial.value;
        gradient = trial.gradient;
        if gradient_converged(&gradient, &x) {
            return objective.result(x, value, 0, "gradient norm below tolerance");
        }
        if (previous - value).abs() < function_tolerance * previous.abs().max(1.0) {
            return objective.result(x, value, 0, "relative objective decrease below tolerance");
        }
        direction = memory.direction(&gradient);
        step = 1.0;
    }
    objective.result(x, value, 0, "iteration limit reached")
}

struct Objective<F> {
    f: F,
    evaluations: usize,
}

impl<F: FnMut(&[f64]) -> (f64, Vec<f64>)> Objective<F> {
    /// Value and gradient, or `None` when either is not finite.
    fn finite_at(&mut self, x: &[f64]) -> Option<(f64, Vec<f64>)> {
        self.evaluations = self.evaluations.saturating_add(1);
        let (value, gradient) = (self.f)(x);
        (value.is_finite() && gradient.len() == x.len() && gradient.iter().all(|g| g.is_finite()))
            .then_some((value, gradient))
    }

    fn result(&self, par: Vec<f64>, value: f64, convergence: i32, message: &str) -> OptimResult {
        OptimResult {
            par,
            value,
            counts: OptimCounts {
                function: self.evaluations,
                gradient: self.evaluations,
            },
            convergence,
            message: message.to_string(),
        }
    }
}

/// The newest correction pairs `(s, y, 1 / y's)`, oldest first.
#[derive(Default)]
struct Memory {
    pairs: VecDeque<(Vec<f64>, Vec<f64>, f64)>,
}

impl Memory {
    fn clear(&mut self) {
        self.pairs.clear();
    }

    /// Store a pair when its curvature `y's` is positive (Nocedal and
    /// Wright, section 7.2), dropping the oldest beyond the memory.
    fn push(&mut self, s: Vec<f64>, y: Vec<f64>) {
        let sy = dot(&s, &y);
        if !(sy > 0.0 && sy.is_finite()) {
            return;
        }
        if self.pairs.len() == MEMORY {
            self.pairs.pop_front();
        }
        self.pairs.push_back((s, y, 1.0 / sy));
    }

    /// `-H g` by the two-loop recursion (Nocedal and Wright, Algorithm 7.4).
    fn direction(&self, gradient: &[f64]) -> Vec<f64> {
        let mut q = gradient.to_vec();
        let mut alphas = [0.0; MEMORY];
        for (k, (s, y, rho)) in self.pairs.iter().enumerate().rev() {
            let alpha = rho * dot(s, &q);
            alphas[k] = alpha;
            axpy(-alpha, y, &mut q);
        }
        if let Some((s, y, _)) = self.pairs.back() {
            let gamma = dot(s, y) / dot(y, y);
            q.iter_mut().for_each(|v| *v *= gamma);
        }
        for (k, (s, y, rho)) in self.pairs.iter().enumerate() {
            let beta = rho * dot(y, &q);
            axpy(alphas[k] - beta, s, &mut q);
        }
        q.iter_mut().for_each(|v| *v = -*v);
        q
    }
}

/// One evaluated point on the search line.
struct Point {
    step: f64,
    /// `+inf` when the objective or gradient was not finite.
    value: f64,
    /// Directional derivative; NaN when not available.
    slope: f64,
    x: Vec<f64>,
    gradient: Vec<f64>,
}

/// Strong Wolfe line search. On failure returns the last finite zoom trial
/// (point and value), if any.
fn line_search<F: FnMut(&[f64]) -> (f64, Vec<f64>)>(
    objective: &mut Objective<F>,
    x: &[f64],
    value: f64,
    slope: f64,
    direction: &[f64],
    initial_step: f64,
) -> Result<Point, Option<(Vec<f64>, f64)>> {
    let mut evaluate = |step: f64| {
        let point: Vec<f64> = x.iter().zip(direction).map(|(x, d)| x + step * d).collect();
        match objective.finite_at(&point) {
            Some((v, g)) => Point {
                step,
                value: v,
                slope: dot(&g, direction),
                x: point,
                gradient: g,
            },
            None => Point {
                step,
                value: f64::INFINITY,
                slope: f64::NAN,
                x: point,
                gradient: vec![],
            },
        }
    };
    // Written as a difference so that rounding of `value + ...` cannot
    // accept a trial that does not lower the objective. Non-finite trials
    // carry an infinite value and always fail.
    let insufficient = |p: &Point| p.value - value > SUFFICIENT_DECREASE * p.step * slope;
    let flat_enough = |p: &Point| p.slope.abs() <= -CURVATURE * slope;
    let mut low = Point {
        step: 0.0,
        value,
        slope,
        x: x.to_vec(),
        gradient: vec![],
    };
    let mut step = initial_step;
    // Bracketing phase (Algorithm 3.5).
    let mut trials = 0;
    let (mut low, mut high) = loop {
        if trials == MAX_PHASE_TRIALS {
            // No bracket after the allowed doublings: take the last step.
            return Ok(low);
        }
        trials += 1;
        let trial = evaluate(step);
        if insufficient(&trial) || (low.step > 0.0 && trial.value >= low.value) {
            break (low, trial);
        }
        if flat_enough(&trial) {
            return Ok(trial);
        }
        if trial.slope >= 0.0 {
            break (trial, low);
        }
        step = 2.0 * trial.step;
        low = trial;
    };
    // Zoom phase (Algorithm 3.6): `low` has sufficient decrease and the
    // lowest value seen; acceptable steps lie between `low` and `high`.
    let mut last_finite = None;
    for _ in 0..MAX_PHASE_TRIALS {
        let width = high.step - low.step;
        let mut step = low.step
            - low.slope * width * width / (2.0 * (high.value - low.value - low.slope * width));
        let (a, b) = (low.step.min(high.step), low.step.max(high.step));
        let margin = SAFEGUARD * (b - a);
        if !(step >= a + margin && step <= b - margin) {
            step = 0.5 * (low.step + high.step);
        }
        let trial = evaluate(step);
        if trial.value.is_finite() {
            last_finite = Some((trial.x.clone(), trial.value));
        }
        if insufficient(&trial) || trial.value >= low.value {
            high = trial;
        } else {
            if flat_enough(&trial) {
                return Ok(trial);
            }
            if trial.slope * (high.step - low.step) >= 0.0 {
                high = low;
            }
            low = trial;
        }
    }
    Err(last_finite)
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}

fn norm(a: &[f64]) -> f64 {
    dot(a, a).sqrt()
}

fn axpy(alpha: f64, x: &[f64], y: &mut [f64]) {
    y.iter_mut().zip(x).for_each(|(y, x)| *y += alpha * x);
}

#[cfg(test)]
// Observed values are copied with the 17 digits the probe printed.
#[allow(clippy::excessive_precision, clippy::needless_range_loop)]
mod tests {
    use super::*;

    fn rosenbrock(x: &[f64]) -> (f64, Vec<f64>) {
        let t1 = x[1] - x[0] * x[0];
        let t2 = 1.0 - x[0];
        (
            100.0 * t1 * t1 + t2 * t2,
            vec![-400.0 * x[0] * t1 - 2.0 * t2, 200.0 * t1],
        )
    }

    /// Runs the optimiser and records every evaluation point.
    fn traced(
        f: impl Fn(&[f64]) -> (f64, Vec<f64>),
        start: Vec<f64>,
        function_tolerance: f64,
        gradient_tolerance: f64,
    ) -> (OptimResult, Vec<Vec<f64>>) {
        let mut points = vec![];
        let result = minimize_with_tolerances(
            |x| {
                points.push(x.to_vec());
                f(x)
            },
            start,
            function_tolerance,
            gradient_tolerance,
        );
        (result, points)
    }

    fn close(a: f64, b: f64, tolerance: f64) -> bool {
        (a - b).abs() <= tolerance * b.abs().max(1e-300)
    }

    #[test]
    fn rosenbrock_matches_the_documented_and_observed_run() {
        // Documented example and probe run with tolerances 1e-6 and 1e-5:
        // 54 evaluations ending at (0.99996829440208279, 0.9999353850282976).
        let (r, points) = traced(rosenbrock, vec![-1.2, 1.0], 1e-6, 1e-5);
        assert_eq!((r.convergence, points.len()), (0, 54), "{}", r.message);
        assert_eq!(r.counts.function, 54);
        assert!(close(r.par[0], 0.99996829440208279, 1e-8), "{:?}", r.par);
        assert!(close(r.par[1], 0.9999353850282976, 1e-8), "{:?}", r.par);
        assert!(close(r.value, 1.1503946922813757e-09, 1e-6));
        let r = minimize(rosenbrock, vec![-1.2, 1.0]);
        assert_eq!(r.convergence, 0);
        assert!(r.par.iter().all(|v| (v - 1.0).abs() < 1e-6), "{:?}", r.par);
    }

    #[test]
    fn convex_quadratic_in_ten_dimensions_reaches_its_minimiser() {
        // f(x) = x'Ax/2 - b'x with tridiagonal A and minimiser x* = (1, -2, 3, ...).
        let n = 10;
        let diagonal: Vec<f64> = (0..n).map(|i| 3.0 + i as f64).collect();
        let product = |x: &[f64]| -> Vec<f64> {
            (0..n)
                .map(|i| {
                    let mut v = diagonal[i] * x[i];
                    if i > 0 {
                        v -= x[i - 1];
                    }
                    if i + 1 < n {
                        v -= x[i + 1];
                    }
                    v
                })
                .collect()
        };
        let target: Vec<f64> = (0..n)
            .map(|i| {
                if i % 2 == 0 {
                    1.0 + i as f64
                } else {
                    -1.0 - i as f64
                }
            })
            .collect();
        let b = product(&target);
        let f = |x: &[f64]| {
            let ax = product(x);
            let value = 0.5 * dot(x, &ax) - dot(&b, x);
            (value, ax.iter().zip(&b).map(|(a, b)| a - b).collect())
        };
        // Gradient test only: the objective test would stop at a relative
        // decrease of 1e-12, before the minimiser is resolved to 1e-6.
        let r = minimize_with_tolerances(f, vec![0.0; n], 0.0, 1e-10);
        assert!(r.convergence == 0 || r.convergence == -1, "{}", r.message);
        for (x, t) in r.par.iter().zip(&target) {
            assert!((x - t).abs() < 1e-6, "{:?}", r.par);
        }
    }

    #[test]
    fn logistic_regression_matches_newton_iterations() {
        // Overlapping two-covariate data with an intercept.
        let n = 40;
        let rows: Vec<[f64; 3]> = (0..n)
            .map(|i| {
                let t = i as f64;
                [1.0, (0.37 * t).sin() * 2.0, ((t * 0.11).cos() - 0.3) * 1.5]
            })
            .collect();
        let y: Vec<f64> = (0..n)
            .map(|i| {
                let r = &rows[i];
                let score = 0.4 + 1.3 * r[1] - 0.8 * r[2] + (1.7 * i as f64).sin();
                if score > 0.0 { 1.0 } else { 0.0 }
            })
            .collect();
        let negative_log_likelihood = |beta: &[f64]| {
            let mut value = 0.0;
            let mut gradient = vec![0.0; 3];
            for (r, &yi) in rows.iter().zip(&y) {
                let eta: f64 = r.iter().zip(beta).map(|(a, b)| a * b).sum();
                let p = 1.0 / (1.0 + (-eta).exp());
                value += eta.max(0.0) + (-eta.abs()).exp().ln_1p() - yi * eta;
                for j in 0..3 {
                    gradient[j] += (p - yi) * r[j];
                }
            }
            (value, gradient)
        };
        // Newton-Raphson reference.
        let mut newton = [0.0; 3];
        for _ in 0..50 {
            let mut hessian = [[0.0; 3]; 3];
            let (_, gradient) = negative_log_likelihood(&newton);
            for r in &rows {
                let eta: f64 = r.iter().zip(&newton).map(|(a, b)| a * b).sum();
                let p = 1.0 / (1.0 + (-eta).exp());
                for j in 0..3 {
                    for k in 0..3 {
                        hessian[j][k] += p * (1.0 - p) * r[j] * r[k];
                    }
                }
            }
            // Solve the 3x3 system by Gaussian elimination.
            let mut a = hessian;
            let mut g = gradient.clone();
            for c in 0..3 {
                for r in c + 1..3 {
                    let m = a[r][c] / a[c][c];
                    for k in c..3 {
                        a[r][k] -= m * a[c][k];
                    }
                    g[r] -= m * g[c];
                }
            }
            for c in (0..3).rev() {
                let s: f64 = (c + 1..3).map(|k| a[c][k] * g[k]).sum();
                g[c] = (g[c] - s) / a[c][c];
            }
            for j in 0..3 {
                newton[j] -= g[j];
            }
        }
        let r = minimize_with_tolerances(negative_log_likelihood, vec![0.0; 3], 1e-12, 1e-10);
        assert!(r.convergence == 0 || r.convergence == -1, "{}", r.message);
        for (b, t) in r.par.iter().zip(&newton) {
            assert!((b - t).abs() < 1e-6, "{:?} vs {newton:?}", r.par);
        }
    }

    #[test]
    fn non_finite_start_is_rejected_without_iterating() {
        for bad in [f64::NAN, f64::INFINITY] {
            let mut calls = 0;
            let r = minimize(
                |_| {
                    calls += 1;
                    (bad, vec![0.0])
                },
                vec![0.5],
            );
            assert_eq!((r.convergence, calls, r.par.clone()), (-3, 1, vec![0.5]));
            let r = minimize(|_| (1.0, vec![bad]), vec![0.5]);
            assert_eq!(r.convergence, -3);
        }
    }

    #[test]
    fn infinite_trials_are_shrunk_and_end_in_a_reported_failure() {
        // Observation `inf_after_unit_step`: bisection away from the infinite
        // region, one accepted step, then a failed search of 101 trials.
        let f = |x: &[f64]| {
            if x[0] > 0.5 {
                (f64::INFINITY, vec![f64::INFINITY])
            } else {
                ((x[0] - 2.0).powi(2) + 10.0, vec![2.0 * (x[0] - 2.0)])
            }
        };
        let (r, points) = traced(f, vec![0.0], 1e-8, 1e-8);
        assert_eq!(points.len(), 104);
        assert_eq!(&points[..4], &[vec![0.0], vec![1.0], vec![0.5], vec![2.0]]);
        assert_eq!(
            (r.convergence, r.par.clone(), r.value),
            (-1, vec![0.5], 12.25)
        );
        assert!(r.message.starts_with("the line search routine failed"));
    }

    #[test]
    fn trial_steps_follow_the_observed_interpolation_and_safeguards() {
        let expo = |x: &[f64]| (x[0].exp() - 2.0 * x[0] + 10.0, vec![x[0].exp() - 2.0]);
        // Observation `exp_1d_far`: doubling, bisection while the quadratic
        // trial lies within 1% of a bracket end, then interpolation.
        let observed = [
            -5.0,
            -4.0,
            -3.0,
            -1.0,
            17.077352895368133,
            8.0386764476840664,
            3.5193382238420332,
            -0.47473364008818231,
            2.3729145798329236,
            0.19555325837005355,
            1.0803603038228453,
            0.65807055853536989,
            0.70265948328706462,
            0.69297964469545126,
            0.69314638497514414,
        ];
        let (r, points) = traced(expo, vec![-5.0], 1e-8, 1e-8);
        assert_eq!((r.convergence, points.len()), (0, observed.len()));
        for (p, o) in points.iter().zip(observed) {
            assert!(close(p[0], o, 1e-12), "{} vs {o}", p[0]);
        }
        // Observations `quadratic_minimum_at_0.01` and `..._0.0105`.
        for (m, second) in [(0.01, 0.5), (0.0105, 0.0105)] {
            let f = move |x: &[f64]| (3.0 * (x[0] - m).powi(2) + 10.0, vec![6.0 * (x[0] - m)]);
            let (_, points) = traced(f, vec![0.0], 0.0, 1e-6);
            assert!(close(points[2][0], second, 1e-12), "{m}: {:?}", points);
        }
        // Observation `steep_return_f0_60`: a too steep acceptable step is
        // zoomed from its own end; a trial within 1% of it is bisected.
        let f = |x: &[f64]| {
            if x[0] == 0.0 {
                (60.0, vec![-1.0])
            } else if x[0] == 1.0 {
                (10.0, vec![0.95])
            } else {
                (9.5 + (x[0] - 0.99).powi(2), vec![2.0 * (x[0] - 0.99)])
            }
        };
        let (_, points) = traced(f, vec![0.0], 0.0, 0.0);
        assert_eq!(&points[..3], &[vec![0.0], vec![1.0], vec![0.5]]);
        assert!(close(points[3][0], 0.66335511401520197, 1e-14));
        // Observation `linear_unbounded`: 100 doublings, then the last step
        // is taken and the relative gradient test stops the run.
        let (r, points) = traced(|x| (-x[0], vec![-1.0]), vec![0.0], 1e-8, 1e-8);
        assert_eq!((r.convergence, points.len()), (0, 101));
        assert_eq!(r.par, vec![2f64.powi(99)]);
    }

    #[test]
    fn objective_test_scales_by_the_previous_value_floored_at_one() {
        // Observations `logcosh_1d_*`: the run stops once a step lowers the
        // objective by less than eps * max(1, |f_previous|).
        let f = |x: &[f64]| (x[0].cosh().ln(), vec![x[0].tanh()]);
        for (start, tolerance, evaluations) in
            [(1.5, 0.05, 5), (1.5, 0.5, 4), (3.0, 0.7, 5), (3.0, 0.9, 3)]
        {
            let (r, points) = traced(f, vec![start], tolerance, 0.0);
            assert_eq!(
                (r.convergence, points.len()),
                (0, evaluations),
                "{start} {tolerance}"
            );
        }
    }
}
