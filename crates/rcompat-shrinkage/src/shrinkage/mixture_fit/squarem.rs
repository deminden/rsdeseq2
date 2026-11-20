//! Independent squared extrapolation (Varadhan/Roland 2008), with explicit
//! frontend policy controls. No SQUAREM/ashr source implementation is imported.
use super::invalid_mixture_value;
use crate::Error;
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
/// Polynomial extrapolation scheme of higher-order SQUAREM.
pub enum PolynomialExtrapolation {
    #[default]
    /// Reduced-rank extrapolation.
    ReducedRank,
    /// Minimal-polynomial extrapolation.
    MinimalPolynomial,
}
#[derive(Clone, Copy, Debug)]
/// SQUAREM acceleration controls (as in the R `SQUAREM` package).
pub struct SquaremOptions {
    /// None selects min(1e-7, 0.1/n), the source observation-count rule.
    pub tolerance: Option<f64>,
    /// Maximum fixed-point evaluations.
    pub max_evaluations: usize,
    /// Step-length scheme (1, 2 or 3).
    pub method: u8,
    /// Smallest extrapolation step.
    pub minimum_step: f64,
    /// Largest initial extrapolation step.
    pub maximum_step: f64,
    /// Factor by which the maximum step grows.
    pub step_multiplier: f64,
    /// Allowed objective increase before falling back to a plain EM step.
    pub objective_increase: f64,
    /// Retained for source control transport; both K=1 policies ignore it.
    pub square: bool,
    /// Extrapolation order.
    pub polynomial_order: usize,
    /// Higher-order extrapolation scheme.
    pub polynomial_method: PolynomialExtrapolation,
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
            polynomial_order: 1,
            polynomial_method: Default::default(),
        }
    }
}
impl SquaremOptions {
    pub(super) fn validate(self) -> Result<(), Error> {
        let (operation, value) = if self.polynomial_order == 0 {
            ("squared-EM polynomial order", 0.0)
        } else if let Some(value) = self.tolerance
            && (!value.is_finite() || value < 0.0)
        {
            ("squared-EM tolerance", value)
        } else if !(1..=3).contains(&self.method) {
            ("squared-EM method", f64::from(self.method))
        } else if !self.minimum_step.is_finite() {
            ("squared-EM minimum step", self.minimum_step)
        } else if !self.maximum_step.is_finite() || self.minimum_step > self.maximum_step {
            ("squared-EM maximum step", self.maximum_step)
        } else if !self.step_multiplier.is_finite() || self.step_multiplier <= 1.0 {
            ("squared-EM step multiplier", self.step_multiplier)
        } else if self.objective_increase.is_nan() || self.objective_increase < 0.0 {
            ("squared-EM objective increase", self.objective_increase)
        } else {
            return Ok(());
        };
        Err(invalid_mixture_value(
            operation,
            None,
            value,
            "invalid squared-EM controls",
        ))
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
    update: impl Fn(&[f64]) -> Result<Vec<f64>, Error>,
    mut objective: impl FnMut(&[f64]) -> f64,
) -> Result<Iteration, Error> {
    options.validate()?;
    if !compiled && options.polynomial_order > 1 {
        return super::polynomial::run(state, options, tolerance, update, objective);
    }
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
