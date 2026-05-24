use crate::errors::DeseqError;

/// Options for a small deterministic bounded optimizer.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BoundedOptimizerOptions {
    /// Maximum number of accepted or attempted line-search iterations.
    pub maxit: usize,
    /// Projected-gradient convergence tolerance.
    pub gradient_tol: f64,
    /// Stop when the accepted parameter move is smaller than this value.
    pub step_tol: f64,
    /// Initial backtracking step size.
    pub initial_step: f64,
    /// Smallest step size considered during backtracking.
    pub min_step: f64,
    /// Armijo sufficient-decrease factor.
    pub armijo: f64,
}

impl Default for BoundedOptimizerOptions {
    fn default() -> Self {
        Self {
            maxit: 200,
            gradient_tol: 1e-8,
            step_tol: 1e-10,
            initial_step: 1.0,
            min_step: 1e-12,
            armijo: 1e-4,
        }
    }
}

/// Output from bounded minimization.
#[derive(Clone, Debug, PartialEq)]
pub struct BoundedOptimizationOutput {
    /// Optimized parameters.
    pub parameters: Vec<f64>,
    /// Final objective value.
    pub value: f64,
    /// Whether the projected-gradient or small-step convergence rule was met.
    pub converged: bool,
    /// Number of outer iterations attempted.
    pub iterations: usize,
}

/// Minimize a differentiable objective over shared lower/upper bounds.
///
/// This is intentionally compact and dependency-free. It uses projected
/// gradients with Armijo backtracking, which is enough for the low-dimensional
/// GLM fallback rows that need a bounded pure-Rust optimizer.
pub fn minimize_bounded_projected_gradient<F>(
    start: &[f64],
    lower: f64,
    upper: f64,
    options: BoundedOptimizerOptions,
    mut objective_and_gradient: F,
) -> Result<BoundedOptimizationOutput, DeseqError>
where
    F: FnMut(&[f64]) -> Result<(f64, Vec<f64>), DeseqError>,
{
    validate_options(lower, upper, options)?;
    let mut parameters = start
        .iter()
        .copied()
        .map(|value| value.clamp(lower, upper))
        .collect::<Vec<_>>();
    let (mut value, mut gradient) = objective_and_gradient(&parameters)?;
    validate_state(value, &gradient, parameters.len())?;

    for iter in 0..options.maxit {
        let direction = projected_descent_direction(&parameters, &gradient, lower, upper);
        let direction_norm = l2_norm(&direction);
        if direction_norm <= options.gradient_tol {
            return Ok(BoundedOptimizationOutput {
                parameters,
                value,
                converged: true,
                iterations: iter,
            });
        }

        let directional_derivative = dot(&gradient, &direction);
        if directional_derivative >= 0.0 || !directional_derivative.is_finite() {
            return Ok(BoundedOptimizationOutput {
                parameters,
                value,
                converged: false,
                iterations: iter,
            });
        }

        let mut step = options.initial_step;
        let mut accepted = None;
        while step >= options.min_step {
            let candidate = parameters
                .iter()
                .copied()
                .zip(direction.iter().copied())
                .map(|(value, delta)| (value + step * delta).clamp(lower, upper))
                .collect::<Vec<_>>();
            let movement = max_abs_difference(&parameters, &candidate);
            if movement <= options.step_tol {
                return Ok(BoundedOptimizationOutput {
                    parameters,
                    value,
                    converged: true,
                    iterations: iter + 1,
                });
            }
            let actual_directional_derivative =
                actual_directional_derivative(&parameters, &candidate, &gradient);
            if !actual_directional_derivative.is_finite() || actual_directional_derivative >= 0.0 {
                step *= 0.5;
                continue;
            }
            let (candidate_value, candidate_gradient) = objective_and_gradient(&candidate)?;
            if candidate_value.is_finite()
                && candidate_value <= value + options.armijo * actual_directional_derivative
            {
                validate_state(candidate_value, &candidate_gradient, parameters.len())?;
                accepted = Some((candidate, candidate_value, candidate_gradient));
                break;
            }
            step *= 0.5;
        }

        let Some((candidate, candidate_value, candidate_gradient)) = accepted else {
            return Ok(BoundedOptimizationOutput {
                parameters,
                value,
                converged: false,
                iterations: iter + 1,
            });
        };
        parameters = candidate;
        value = candidate_value;
        gradient = candidate_gradient;
    }

    Ok(BoundedOptimizationOutput {
        parameters,
        value,
        converged: false,
        iterations: options.maxit,
    })
}

fn validate_options(
    lower: f64,
    upper: f64,
    options: BoundedOptimizerOptions,
) -> Result<(), DeseqError> {
    if !lower.is_finite()
        || !upper.is_finite()
        || lower >= upper
        || options.maxit == 0
        || !options.gradient_tol.is_finite()
        || options.gradient_tol <= 0.0
        || !options.step_tol.is_finite()
        || options.step_tol <= 0.0
        || !options.initial_step.is_finite()
        || options.initial_step <= 0.0
        || !options.min_step.is_finite()
        || options.min_step <= 0.0
        || options.min_step > options.initial_step
        || !options.armijo.is_finite()
        || options.armijo <= 0.0
        || options.armijo >= 1.0
    {
        return Err(DeseqError::InvalidOptions {
            reason: "invalid bounded optimizer options".to_string(),
        });
    }
    Ok(())
}

fn validate_state(value: f64, gradient: &[f64], expected_len: usize) -> Result<(), DeseqError> {
    if gradient.len() != expected_len {
        return Err(crate::errors::invalid_dimensions(
            "bounded optimizer gradient",
            expected_len,
            gradient.len(),
        ));
    }
    if !value.is_finite() {
        return Err(DeseqError::NonFiniteValue {
            context: "bounded optimizer objective".to_string(),
            index: None,
            value,
        });
    }
    for (idx, value) in gradient.iter().copied().enumerate() {
        if !value.is_finite() {
            return Err(DeseqError::NonFiniteValue {
                context: "bounded optimizer gradient".to_string(),
                index: Some(idx),
                value,
            });
        }
    }
    Ok(())
}

fn projected_descent_direction(
    parameters: &[f64],
    gradient: &[f64],
    lower: f64,
    upper: f64,
) -> Vec<f64> {
    parameters
        .iter()
        .copied()
        .zip(gradient.iter().copied())
        .map(|(parameter, gradient)| {
            if (parameter <= lower && gradient > 0.0) || (parameter >= upper && gradient < 0.0) {
                0.0
            } else {
                -gradient
            }
        })
        .collect()
}

fn l2_norm(values: &[f64]) -> f64 {
    values.iter().map(|value| value * value).sum::<f64>().sqrt()
}

fn dot(left: &[f64], right: &[f64]) -> f64 {
    left.iter()
        .copied()
        .zip(right.iter().copied())
        .map(|(left, right)| left * right)
        .sum()
}

fn max_abs_difference(left: &[f64], right: &[f64]) -> f64 {
    left.iter()
        .copied()
        .zip(right.iter().copied())
        .map(|(left, right)| (left - right).abs())
        .fold(0.0, f64::max)
}

fn actual_directional_derivative(parameters: &[f64], candidate: &[f64], gradient: &[f64]) -> f64 {
    gradient
        .iter()
        .copied()
        .zip(
            candidate
                .iter()
                .copied()
                .zip(parameters.iter().copied())
                .map(|(candidate, parameter)| candidate - parameter),
        )
        .map(|(gradient, direction)| gradient * direction)
        .sum()
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn bounded_optimizer_minimizes_quadratic_inside_bounds() {
        let output = minimize_bounded_projected_gradient(
            &[8.0],
            -10.0,
            10.0,
            BoundedOptimizerOptions::default(),
            |values| {
                let delta = values[0] - 2.0;
                Ok((delta * delta, vec![2.0 * delta]))
            },
        )
        .unwrap();

        assert!(output.converged);
        assert_relative_eq!(output.parameters[0], 2.0, epsilon = 1e-6);
    }

    #[test]
    fn bounded_optimizer_stops_at_active_bound() {
        let output = minimize_bounded_projected_gradient(
            &[0.0],
            -1.0,
            1.0,
            BoundedOptimizerOptions::default(),
            |values| {
                let delta = values[0] - 3.0;
                Ok((delta * delta, vec![2.0 * delta]))
            },
        )
        .unwrap();

        assert!(output.converged);
        assert_relative_eq!(output.parameters[0], 1.0, epsilon = 1e-10);
    }

    #[test]
    fn actual_directional_derivative_uses_clamped_candidate_movement() {
        let derivative = actual_directional_derivative(&[0.9], &[1.0], &[-2.0]);

        assert_relative_eq!(derivative, -0.2, epsilon = 1e-12);
    }
}
