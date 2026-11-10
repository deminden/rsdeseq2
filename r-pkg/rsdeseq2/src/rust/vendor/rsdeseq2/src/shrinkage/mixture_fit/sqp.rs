use super::{DMatrix, DVector, DeseqError, MixtureFitOptions, Problem, invalid};

impl Problem {
    fn expectation_step(&self, weights: &[f64]) -> Result<Vec<f64>, DeseqError> {
        let mut next = self
            .prior_counts
            .iter()
            .map(|b| b / self.mass)
            .collect::<Vec<_>>();
        for (row, &a) in self
            .likelihoods
            .chunks_exact(self.components)
            .zip(&self.observation_weights)
        {
            if a == 0.0 {
                continue;
            }
            let marginal = row.iter().zip(weights).map(|(l, w)| l * w).sum::<f64>();
            if marginal <= 0.0 {
                return Err(invalid("EM update has no positive likelihood"));
            }
            for ((new, &likelihood), &w) in next.iter_mut().zip(row).zip(weights) {
                *new += (a / self.mass) * ((w * likelihood) / marginal);
            }
        }
        Ok(next)
    }
    fn hessian(&self, weights: &[f64]) -> DMatrix<f64> {
        let k = self.components;
        let mut h = DMatrix::zeros(k, k);
        for (row, &a) in self
            .likelihoods
            .chunks_exact(k)
            .zip(&self.observation_weights)
        {
            if a == 0.0 {
                continue;
            }
            let marginal = row.iter().zip(weights).map(|(l, w)| l * w).sum::<f64>() + self.offset;
            for i in 0..k {
                for j in 0..k {
                    h[(i, j)] += (a / self.mass) * (row[i] / marginal) * (row[j] / marginal);
                }
            }
        }
        for i in 0..k {
            if self.prior_counts[i] > 0.0 {
                h[(i, i)] +=
                    (self.prior_counts[i] / self.mass) / (weights[i] + self.offset).powi(2);
            }
        }
        h
    }
    pub(super) fn solve(
        &self,
        mut weights: Vec<f64>,
        options: MixtureFitOptions,
    ) -> Result<(Vec<f64>, bool, f64, usize, String), DeseqError> {
        for _ in 0..options.em_iterations {
            weights = self.expectation_step(&weights)?;
        }
        let mut evaluations = 0;
        for _ in 0..options.max_iterations {
            weights = self.expectation_step(&weights)?;
            for (w, &prior) in weights.iter_mut().zip(&self.prior_counts) {
                if *w <= options.zero_threshold && (prior == 0.0 || self.offset > 0.0) {
                    *w = 0.0;
                }
            }
            if !weights.iter().any(|w| *w > 0.0) {
                return Err(invalid("zero threshold removed every component"));
            }
            let (value, gradient) = self.objective_gradient(&weights);
            if !value.is_finite() || gradient.iter().any(|g| !g.is_finite()) {
                return Err(invalid(
                    "objective or gradient exceeds the representable range",
                ));
            }
            evaluations += 1;
            let residual = dual_residual(&gradient, &weights, options);
            if residual <= options.gradient_tolerance {
                return Ok((
                    weights,
                    true,
                    residual,
                    evaluations,
                    "SQP dual convergence reached".into(),
                ));
            }
            let hessian = self.hessian(&weights);
            if hessian.iter().any(|v| !v.is_finite()) {
                return Err(invalid("Hessian exceeds the representable range"));
            }
            let current = DVector::from_column_slice(&weights);
            let linear = DVector::from_column_slice(&gradient) - &hessian * &current;
            let candidate =
                nonnegative_quadratic(&hessian, &linear, &current, options.active_iterations)?;
            let direction = &candidate - &current;
            let slope = gradient
                .iter()
                .zip(direction.iter())
                .map(|(g, d)| g * d)
                .sum::<f64>();
            let mut step = 1.0;
            let mut accepted = None;
            while step >= 1e-8 {
                let proposed = weights
                    .iter()
                    .zip(direction.iter())
                    .map(|(w, d)| (w + step * d).max(0.0))
                    .collect::<Vec<_>>();
                let next = self.objective_gradient(&proposed).0;
                evaluations += 1;
                if next.is_finite() && next <= value + 0.01 * step * slope {
                    accepted = Some(proposed);
                    break;
                }
                step *= 0.75;
            }
            if let Some(next) = accepted {
                weights = next;
            } else {
                return Ok((
                    weights,
                    false,
                    residual,
                    evaluations,
                    "SQP line search failed".into(),
                ));
            }
        }
        let gradient = self.objective_gradient(&weights).1;
        let residual = dual_residual(&gradient, &weights, options);
        evaluations += 1;
        Ok((
            weights,
            false,
            residual,
            evaluations,
            "SQP iteration limit reached".into(),
        ))
    }
}

fn dual_residual(gradient: &[f64], weights: &[f64], options: MixtureFitOptions) -> f64 {
    (-gradient
        .iter()
        .zip(weights)
        .filter(|(_, w)| !options.positive_coordinate_convergence || **w > options.zero_threshold)
        .map(|(g, _)| *g)
        .fold(f64::INFINITY, f64::min))
    .max(0.0)
}

// Active-set nonnegative quadratic minimization. Each free subproblem solves
// for its target coordinates; a blocking bound shortens the move and removes
// the corresponding coordinate. Negative bound multipliers reopen coordinates.
fn nonnegative_quadratic(
    h: &DMatrix<f64>,
    c: &DVector<f64>,
    start: &DVector<f64>,
    max_iterations: usize,
) -> Result<DVector<f64>, DeseqError> {
    let k = c.len();
    let mut y = start.clone();
    let mut free = y.iter().map(|v| *v > 0.0).collect::<Vec<_>>();
    for _ in 0..max_iterations {
        let indices = (0..k).filter(|&i| free[i]).collect::<Vec<_>>();
        if indices.is_empty() {
            let gradient = h * &y + c;
            let i = (0..k)
                .min_by(|&i, &j| gradient[i].total_cmp(&gradient[j]))
                .unwrap();
            if gradient[i] >= -1e-10 {
                return Ok(y);
            }
            free[i] = true;
            continue;
        }
        let mut reduced = DMatrix::from_fn(indices.len(), indices.len(), |i, j| {
            h[(indices[i], indices[j])]
        });
        let rhs = DVector::from_iterator(indices.len(), indices.iter().map(|&i| -c[i]));
        let scale = reduced.diagonal().iter().copied().fold(1.0, f64::max);
        let mut ridge = 0.0;
        let mut solved = None;
        for _ in 0..20 {
            if let Some(cholesky) = reduced.clone().cholesky() {
                solved = Some(cholesky.solve(&rhs));
                break;
            }
            let next = if ridge == 0.0 {
                1e-15 * scale
            } else {
                10.0 * ridge
            };
            for i in 0..indices.len() {
                reduced[(i, i)] += next - ridge;
            }
            ridge = next;
        }
        let target = solved
            .ok_or_else(|| invalid("quadratic subproblem is not numerically positive definite"))?;
        if target.iter().any(|v| !v.is_finite()) {
            return Err(invalid(
                "quadratic solution exceeds the representable range",
            ));
        }
        let mut step = 1.0_f64;
        let mut blocked = None;
        for (n, &i) in indices.iter().enumerate() {
            let d = target[n] - y[i];
            if d < 0.0 && -y[i] / d < step {
                step = -y[i] / d;
                blocked = Some(i);
            }
        }
        for (n, &i) in indices.iter().enumerate() {
            y[i] = (y[i] + step * (target[n] - y[i])).max(0.0);
        }
        if let Some(i) = blocked {
            y[i] = 0.0;
            free[i] = false;
            continue;
        }
        let gradient = h * &y + c;
        let worst = (0..k)
            .filter(|&i| !free[i])
            .min_by(|&i, &j| gradient[i].total_cmp(&gradient[j]));
        if let Some(i) = worst
            && gradient[i] < -1e-10
        {
            free[i] = true;
            continue;
        }
        return Ok(y);
    }
    Ok(y)
}
