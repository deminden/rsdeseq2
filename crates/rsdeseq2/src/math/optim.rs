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
        let Some(direction_norm) = l2_norm(&direction) else {
            return Ok(BoundedOptimizationOutput {
                parameters,
                value,
                converged: false,
                iterations: iter,
            });
        };
        if direction_norm <= options.gradient_tol {
            return Ok(BoundedOptimizationOutput {
                parameters,
                value,
                converged: true,
                iterations: iter,
            });
        }

        let Some(directional_derivative) = dot(&gradient, &direction) else {
            return Ok(BoundedOptimizationOutput {
                parameters,
                value,
                converged: false,
                iterations: iter,
            });
        };
        if directional_derivative >= 0.0 {
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
            let Some(movement) = max_abs_difference(&parameters, &candidate) else {
                step *= 0.5;
                continue;
            };
            if movement <= options.step_tol {
                return Ok(BoundedOptimizationOutput {
                    parameters,
                    value,
                    converged: true,
                    iterations: iter + 1,
                });
            }
            let Some(actual_directional_derivative) =
                actual_directional_derivative(&parameters, &candidate, &gradient)
            else {
                step *= 0.5;
                continue;
            };
            if actual_directional_derivative >= 0.0 {
                step *= 0.5;
                continue;
            }
            let (candidate_value, candidate_gradient) = objective_and_gradient(&candidate)?;
            let armijo_bound =
                checked_armijo_bound(value, options.armijo, actual_directional_derivative);
            if candidate_value.is_finite()
                && armijo_bound.is_some_and(|bound| candidate_value <= bound)
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

fn l2_norm(values: &[f64]) -> Option<f64> {
    let scale = values
        .iter()
        .copied()
        .map(f64::abs)
        .try_fold(0.0_f64, |scale, value| {
            value.is_finite().then_some(scale.max(value))
        })?;
    if scale == 0.0 {
        return Some(0.0);
    }
    let mut sum = 0.0;
    for value in values.iter().copied() {
        let scaled = value / scale;
        let term = scaled * scaled;
        let next = sum + term;
        if !term.is_finite() || !next.is_finite() {
            return None;
        }
        sum = next;
    }
    let norm = scale * sum.sqrt();
    norm.is_finite().then_some(norm)
}

fn dot(left: &[f64], right: &[f64]) -> Option<f64> {
    let mut terms = Vec::with_capacity(left.len().min(right.len()));
    for (left, right) in left.iter().copied().zip(right.iter().copied()) {
        let term = left * right;
        if !term.is_finite() {
            return None;
        }
        terms.push(term);
    }
    scaled_sum(&terms)
}

fn max_abs_difference(left: &[f64], right: &[f64]) -> Option<f64> {
    let mut max_difference = 0.0;
    for (left, right) in left.iter().copied().zip(right.iter().copied()) {
        let difference = (left - right).abs();
        if !difference.is_finite() {
            return None;
        }
        max_difference = f64::max(max_difference, difference);
    }
    Some(max_difference)
}

fn actual_directional_derivative(
    parameters: &[f64],
    candidate: &[f64],
    gradient: &[f64],
) -> Option<f64> {
    let mut terms = Vec::with_capacity(gradient.len().min(candidate.len()).min(parameters.len()));
    for (gradient, (candidate, parameter)) in gradient
        .iter()
        .copied()
        .zip(candidate.iter().copied().zip(parameters.iter().copied()))
    {
        let direction = candidate - parameter;
        let term = gradient * direction;
        if !direction.is_finite() || !term.is_finite() {
            return None;
        }
        terms.push(term);
    }
    scaled_sum(&terms)
}

fn scaled_sum(values: &[f64]) -> Option<f64> {
    let scale = values
        .iter()
        .copied()
        .map(f64::abs)
        .try_fold(0.0_f64, |scale, value| {
            value.is_finite().then_some(scale.max(value))
        })?;
    if scale == 0.0 {
        return Some(0.0);
    }
    let mut sum = 0.0;
    for value in values.iter().copied() {
        let next = sum + value / scale;
        if !next.is_finite() {
            return None;
        }
        sum = next;
    }
    let value = sum * scale;
    value.is_finite().then_some(value)
}

fn checked_armijo_bound(value: f64, armijo: f64, derivative: f64) -> Option<f64> {
    let term = armijo * derivative;
    let bound = value + term;
    (value.is_finite()
        && armijo.is_finite()
        && derivative.is_finite()
        && term.is_finite()
        && bound.is_finite())
    .then_some(bound)
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
        let derivative = actual_directional_derivative(&[0.9], &[1.0], &[-2.0]).unwrap();

        assert_relative_eq!(derivative, -0.2, epsilon = 1e-12);
    }

    #[test]
    fn scalar_helpers_reject_nonfinite_accumulation() {
        assert_relative_eq!(
            l2_norm(&[f64::MAX / 2.0, f64::MAX / 2.0]).unwrap(),
            f64::MAX / 2.0 * 2.0_f64.sqrt(),
            epsilon = 1e292
        );
        assert_eq!(dot(&[f64::MAX, f64::MAX], &[1.0, -1.0]).unwrap(), 0.0);
        assert_eq!(
            actual_directional_derivative(&[0.0, 0.0], &[1.0, 1.0], &[f64::MAX, -f64::MAX])
                .unwrap(),
            0.0
        );
        assert_eq!(l2_norm(&[f64::MAX, f64::MAX]), None);
        assert_eq!(dot(&[f64::MAX], &[2.0]), None);
        assert_eq!(max_abs_difference(&[-f64::MAX], &[f64::MAX]), None);
        assert_eq!(
            actual_directional_derivative(&[0.0], &[2.0], &[f64::MAX]),
            None
        );
        assert_eq!(checked_armijo_bound(f64::MAX, 1.0, f64::MAX), None);
    }

    #[test]
    fn bounded_optimizer_reports_nonconvergence_for_overflowed_gradient_norm() {
        let output = minimize_bounded_projected_gradient(
            &[0.0],
            -10.0,
            10.0,
            BoundedOptimizerOptions::default(),
            |_| Ok((1.0, vec![f64::MAX])),
        )
        .unwrap();

        assert!(!output.converged);
        assert_eq!(output.iterations, 0);
    }
}
