//! Native log-barrier path following for the convex mixture problem.
//! Scaled Newton equations avoid inverse-square weight conditioning. The
//! accepted SQP path is neither called nor changed by this solver.
use super::{DMatrix, DVector, Problem, invalid};
use crate::DeseqError;
#[derive(Clone, Copy, Debug)]
pub struct InteriorPointOptions {
    /// Upper bound on the normalized convex dual gap. A mass-scaled
    /// accuracy target may tighten it to preserve log-likelihood precision.
    pub relative_tolerance: f64,
    /// Total Newton and barrier-centering iterations, including early exits.
    pub max_iterations: usize,
}
impl Default for InteriorPointOptions {
    fn default() -> Self {
        Self {
            relative_tolerance: 1e-6,
            max_iterations: 300,
        }
    }
}
pub(super) struct InteriorPointFit {
    pub weights: Vec<f64>,
    pub converged: bool,
    pub residual: f64,
    pub iterations: usize,
    pub evaluations: usize,
    pub message: String,
}
fn barrier(problem: &Problem, w: &[f64], mu: f64, evaluations: &mut usize) -> f64 {
    if w.iter().any(|v| *v <= 0.0 || !v.is_finite()) {
        return f64::INFINITY;
    }
    *evaluations += 1;
    problem.objective_gradient(w).0 - mu * w.iter().map(|v| v.ln()).sum::<f64>()
}
impl Problem {
    pub(super) fn solve_interior_point(
        &self,
        start: Vec<f64>,
        o: InteriorPointOptions,
    ) -> Result<InteriorPointFit, DeseqError> {
        if !o.relative_tolerance.is_finite() || o.relative_tolerance <= 0.0 {
            return Err(invalid(
                "interior-point relative tolerance must be positive finite",
            ));
        }
        let k = self.components;
        let mut x = start
            .iter()
            .map(|v| (v + 1.0 / k as f64) / 2.0)
            .collect::<Vec<_>>();
        let mut mu = 0.1 / k as f64;
        let accuracy = o
            .relative_tolerance
            .min((1e-8 / self.mass).max(64.0 * f64::EPSILON));
        let mut iterations = 0;
        let mut evaluations = 0;
        while iterations < o.max_iterations {
            let sum = x.iter().sum::<f64>();
            let normalized = x.iter().map(|v| v / sum).collect::<Vec<_>>();
            evaluations += 1;
            let residual = self
                .objective_gradient(&normalized)
                .1
                .iter()
                .map(|v| (-v).max(0.0))
                .fold(0.0, f64::max);
            if residual <= accuracy {
                return Ok(InteriorPointFit {
                    weights: x,
                    converged: true,
                    residual,
                    iterations,
                    evaluations,
                    message: "interior-point dual gap converged".into(),
                });
            }

            evaluations += 1;
            let (f, g) = self.objective_gradient(&x);
            let h = self.hessian(&x);
            let mut scaled = DMatrix::from_fn(k, k, |i, j| x[i] * h[(i, j)] * x[j]);
            for j in 0..k {
                scaled[(j, j)] += mu;
            }
            let rhs = DVector::from_iterator(k, g.iter().zip(&x).map(|(g, x)| g * x - mu));
            let direction = scaled
                .cholesky()
                .and_then(|h| {
                    let v = h.solve(&rhs);
                    v.iter().all(|v| v.is_finite()).then_some(v)
                })
                .ok_or_else(|| invalid("interior-point scaled Hessian is not positive definite"))?;
            let decrement = rhs.dot(&direction);
            if decrement <= (mu * 0.1).max(1e-25) {
                mu *= 0.1;
                iterations += 1;
                continue;
            }
            let delta = x
                .iter()
                .zip(direction.iter())
                .map(|(x, d)| -x * d)
                .collect::<Vec<_>>();
            let slope = g
                .iter()
                .zip(&x)
                .zip(&delta)
                .map(|((g, x), d)| (g - mu / x) * d)
                .sum::<f64>();
            let before = f - mu * x.iter().map(|v| v.ln()).sum::<f64>();
            let mut step = 1.0_f64;
            for (&w, &d) in x.iter().zip(&delta) {
                if d < 0.0 {
                    step = step.min(-0.99 * w / d);
                }
            }
            let mut candidate = Vec::new();
            let mut accepted = false;
            for _ in 0..60 {
                candidate = x.iter().zip(&delta).map(|(x, d)| x + step * d).collect();
                if barrier(self, &candidate, mu, &mut evaluations) <= before + 0.01 * step * slope {
                    accepted = true;
                    break;
                }
                step *= 0.5;
            }
            iterations += 1;
            if accepted {
                x = candidate;
            } else if decrement <= 1e-12 {
                mu *= 0.1;
            } else {
                return Err(invalid(
                    "interior-point line search did not decrease the barrier",
                ));
            }
        }
        let total = x.iter().sum::<f64>();
        let normalized = x.iter().map(|v| v / total).collect::<Vec<_>>();
        evaluations += 1;
        let residual = self
            .objective_gradient(&normalized)
            .1
            .iter()
            .map(|g| (-g).max(0.0))
            .fold(0.0, f64::max);
        Ok(InteriorPointFit {
            weights: x,
            converged: false,
            residual,
            iterations,
            evaluations,
            message: "interior-point iteration budget exhausted".into(),
        })
    }
}
