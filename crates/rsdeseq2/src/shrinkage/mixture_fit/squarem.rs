//! Independent squared extrapolation (Varadhan/Roland 2008), with explicit
//! frontend policy controls. No SQUAREM/ashr source implementation is imported.
use super::invalid;
use crate::DeseqError;
#[derive(Clone, Copy, Debug)]
pub struct SquaremOptions {
    /// None selects min(1e-7, 0.1/n), the source observation-count rule.
    pub tolerance: Option<f64>,
    pub max_evaluations: usize,
    pub method: u8,
    pub minimum_step: f64,
    pub maximum_step: f64,
    pub step_multiplier: f64,
    pub objective_increase: f64,
    /// Retained for source control transport; both K=1 policies ignore it.
    pub square: bool,
}
impl Default for SquaremOptions {
    fn default() -> Self {
        Self {
            tolerance: None,
            max_evaluations: 5000,
            method: 3,
            minimum_step: 1.0,
            maximum_step: 1.0,
            step_multiplier: 4.0,
            objective_increase: 1.0,
            square: true,
        }
    }
}
impl SquaremOptions {
    pub(super) fn validate(self) -> Result<(), DeseqError> {
        if self.tolerance.is_some_and(|v| !v.is_finite() || v < 0.0)
            || !(1..=3).contains(&self.method)
            || !self.minimum_step.is_finite()
            || !self.maximum_step.is_finite()
            || self.minimum_step > self.maximum_step
            || !self.step_multiplier.is_finite()
            || self.step_multiplier <= 1.0
            || self.objective_increase.is_nan()
            || self.objective_increase < 0.0
        {
            return Err(invalid("invalid squared-EM controls"));
        }
        Ok(())
    }
}
pub(super) struct Iteration {
    pub state: Vec<f64>,
    pub converged: bool,
    pub evaluations: usize,
    pub iterations: usize,
    pub objective: f64,
    pub residual: f64,
}
fn norm2(v: &[f64]) -> f64 {
    v.iter().map(|v| v * v).sum::<f64>().sqrt()
}
pub(super) fn run(
    mut state: Vec<f64>,
    options: SquaremOptions,
    tolerance: f64,
    compiled: bool,
    update: impl Fn(&[f64]) -> Result<Vec<f64>, DeseqError>,
    objective: impl Fn(&[f64]) -> f64,
) -> Result<Iteration, DeseqError> {
    options.validate()?;
    let mut value = objective(&state);
    let mut calls = 0_usize;
    let mut iterations = 1_usize;
    let mut residual = f64::INFINITY;
    let mut converged = false;
    let baseline = options.maximum_step;
    let mut low = options.minimum_step;
    let mut high = baseline;
    while if compiled {
        calls.saturating_add(1) < options.max_evaluations
    } else {
        calls < options.max_evaluations
    } {
        let first = update(&state)?;
        calls += 1;
        let delta = first
            .iter()
            .zip(&state)
            .map(|(a, b)| a - b)
            .collect::<Vec<_>>();
        residual = norm2(&delta);
        if residual < tolerance {
            if !compiled {
                state = first;
            }
            converged = calls < options.max_evaluations;
            break;
        }
        let second = update(&first)?;
        calls += 1;
        let next_delta = second
            .iter()
            .zip(&first)
            .map(|(a, b)| a - b)
            .collect::<Vec<_>>();
        residual = norm2(&next_delta);
        if residual < tolerance {
            // The source compiled routine reports convergence here but
            // returns the old state, even if the first update moved far.
            // Commit the stationary second update to avoid false convergence.
            state = second;
            converged = calls < options.max_evaluations;
            break;
        }
        let curvature = next_delta
            .iter()
            .zip(&delta)
            .map(|(a, b)| a - b)
            .collect::<Vec<_>>();
        let rr = delta.iter().map(|r| r * r).sum::<f64>();
        let vv = curvature.iter().map(|v| v * v).sum::<f64>();
        let rv = delta
            .iter()
            .zip(&curvature)
            .map(|(r, v)| r * v)
            .sum::<f64>();
        let raw = match if compiled { 3 } else { options.method } {
            1 => -rv / vv,
            2 => -rr / rv,
            _ => (rr / vv).sqrt(),
        };
        let mut step = raw.min(high).max(low);
        let candidate = state
            .iter()
            .zip(&delta)
            .zip(&curvature)
            .map(|((p, r), v)| p + 2.0 * step * r + step * step * v)
            .collect::<Vec<_>>();
        let candidate = if (step - 1.0).abs() > 0.01 {
            calls += 1;
            update(&candidate)
        } else {
            Ok(candidate)
        };
        let (candidate, score, rejected) = match candidate {
            Ok(candidate) if candidate.iter().all(|v| v.is_finite()) => {
                let score = if options.objective_increase.is_infinite() {
                    value
                } else {
                    objective(&candidate)
                };
                if score.is_nan() || score > value + options.objective_increase {
                    (second.clone(), objective(&second), true)
                } else {
                    (candidate, score, false)
                }
            }
            _ => (second.clone(), objective(&second), true),
        };
        if rejected {
            if step == high {
                high = baseline.max(high / options.step_multiplier);
            }
            step = 1.0;
        }
        if step == high {
            high *= options.step_multiplier;
        }
        if low < 0.0 && step == low {
            low *= options.step_multiplier;
        }
        state = candidate;
        if !score.is_nan() {
            value = score;
        }
        iterations += 1;
    }
    if compiled || options.objective_increase.is_infinite() {
        value = objective(&state);
    }
    Ok(Iteration {
        state,
        converged,
        evaluations: calls,
        iterations,
        objective: value,
        residual,
    })
}
