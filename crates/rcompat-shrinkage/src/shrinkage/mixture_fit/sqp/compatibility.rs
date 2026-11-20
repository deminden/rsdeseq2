//! Direction-based active-set and floored EM policy used by original mixSQP.
//! Algorithm reference: mixsqp 0.3-54, Kim/Carbonetto/Stephens/Anitescu (MIT).
//! See MIXSQP-LICENSE for the original copyright and permission notice.
use super::{DMatrix, DVector, ERROR_SCOPE, Error, Problem, try_mixture_matrix};
use crate::errors::scoped_invalid_options;
use crate::shrinkage::MixtureSqpOptions;

impl Problem {
    pub(super) fn ashr_hessian_into(&self, weights: &[f64], hessian: &mut DMatrix<f64>) {
        let k = self.components;
        hessian.fill(0.0);
        for i in 0..k {
            if self.prior_counts[i] > 0.0 {
                let mass = self.prior_counts[i] / self.mass;
                let marginal = weights[i] + self.offset;
                let z = mass.sqrt() / marginal;
                hessian[(i, i)] = z * z;
            }
        }
        for (row, &mass) in self
            .likelihoods
            .chunks_exact(k)
            .zip(&self.observation_weights)
        {
            if mass == 0.0 {
                continue;
            }
            let mass = mass / self.mass;
            let marginal = row.iter().zip(weights).map(|(l, w)| l * w).sum::<f64>() + self.offset;
            let z = mass.sqrt() / marginal;
            for i in 0..k {
                for j in 0..k {
                    hessian[(i, j)] += (z * row[i]) * (z * row[j]);
                }
            }
        }
    }
    pub(super) fn ashr_expectation_step(&self, weights: &[f64]) -> Result<Vec<f64>, Error> {
        const FLOOR: f64 = 1e-15;
        let mut next = vec![0.0; self.components];
        let mut probabilities = vec![0.0; self.components];
        let mut add = |row: &[f64], mass: f64| -> Result<(), Error> {
            for ((p, &likelihood), &weight) in probabilities.iter_mut().zip(row).zip(weights) {
                *p = likelihood * (weight + FLOOR);
            }
            let maximum = probabilities
                .iter()
                .copied()
                .fold(f64::NEG_INFINITY, f64::max);
            if maximum <= 0.0 || !maximum.is_finite() {
                return Err(scoped_invalid_options(
                    ERROR_SCOPE,
                    "EM update has no positive likelihood",
                ));
            }
            for p in &mut probabilities {
                *p = *p / maximum + FLOOR;
            }
            let total = probabilities.iter().sum::<f64>();
            for (value, &p) in next.iter_mut().zip(&probabilities) {
                *value += (mass / self.mass) * (p / total);
            }
            Ok(())
        };
        let mut prior_row = vec![0.0; self.components];
        for (i, &count) in self.prior_counts.iter().enumerate() {
            if count > 0.0 {
                prior_row[i] = 1.0;
                add(&prior_row, count)?;
                prior_row[i] = 0.0;
            }
        }
        for (row, &mass) in self
            .likelihoods
            .chunks_exact(self.components)
            .zip(&self.observation_weights)
        {
            if mass > 0.0 {
                add(row, mass)?;
            }
        }
        if self.low_rank {
            for x in &mut next {
                *x = x.max(0.0);
            }
        }
        Ok(next)
    }
}

pub(super) fn quadratic(
    h: &DMatrix<f64>,
    c: &DVector<f64>,
    start: &DVector<f64>,
    max_iterations: usize,
    options: MixtureSqpOptions,
) -> Result<DVector<f64>, Error> {
    let k = c.len();
    let mut y = start.clone();
    let mut free = y.iter().map(|&v| v > 0.0).collect::<Vec<_>>();
    for _ in 0..max_iterations {
        let indices = (0..k).filter(|&i| free[i]).collect::<Vec<_>>();
        if indices.is_empty() {
            let gradient = h * &y + c;
            let i = (0..k)
                .min_by(|&i, &j| gradient[i].total_cmp(&gradient[j]))
                .ok_or_else(|| crate::errors::invariant("value present"))?;
            if gradient[i] >= -options.active_tolerance {
                break;
            }
            free[i] = true;
            continue;
        }
        for i in 0..k {
            if !free[i] {
                y[i] = 0.0;
            }
        }
        let n = indices.len();
        let reduced = try_mixture_matrix(n, n, "mixture SQP active Hessian", |i, j| {
            h[(indices[i], indices[j])]
        })?;
        let current = DVector::from_iterator(indices.len(), indices.iter().map(|&i| y[i]));
        let rhs = -(DVector::from_iterator(indices.len(), indices.iter().map(|&i| c[i]))
            + &reduced * current);
        let direction = search_direction(&reduced, &rhs, options.ridge_multiplier)?;
        if direction.iter().any(|x| !x.is_finite()) {
            return Err(scoped_invalid_options(
                ERROR_SCOPE,
                "quadratic solution exceeds the representable range",
            ));
        }
        let threshold = options.search_direction_threshold.unwrap_or(1e-14);
        if direction.iter().all(|x| x.abs() <= threshold) {
            let gradient = h * &y + c;
            let worst = (0..k)
                .filter(|&i| !free[i])
                .min_by(|&i, &j| gradient[i].total_cmp(&gradient[j]));
            if let Some(i) = worst
                && gradient[i] < -options.active_tolerance
            {
                free[i] = true;
            } else {
                break;
            }
        } else {
            let mut step = 1.0_f64;
            let mut blocked = None;
            for (n, &i) in indices.iter().enumerate() {
                if direction[n] < 0.0 && -y[i] / direction[n] < step {
                    step = -y[i] / direction[n];
                    blocked = Some(i);
                }
            }
            for (n, &i) in indices.iter().enumerate() {
                y[i] = (y[i] + step * direction[n]).max(0.0);
            }
            if indices.len() > 1
                && let Some(i) = blocked
            {
                free[i] = false;
                y[i] = 0.0;
            }
        }
    }
    Ok(y)
}

// A successful Cholesky factorization alone is not the original solver's
// acceptance condition: near-singular directions use a truncated SVD. Compute
// the reciprocal 1-norm condition estimate without an R/LAPACK dependency.
fn search_direction(
    h: &DMatrix<f64>,
    rhs: &DVector<f64>,
    ridge_multiplier: f64,
) -> Result<DVector<f64>, Error> {
    let n = h.nrows();
    let minimum = (0..n).map(|i| h[(i, i)]).fold(f64::INFINITY, f64::min);
    let mut ridge = if minimum > 1e-15 {
        0.0
    } else {
        1e-15 - minimum
    };
    let mut solved = None;
    for _ in 0..128 {
        if ridge * ridge_multiplier > 1e15 {
            break;
        }
        let mut corrected =
            try_mixture_matrix(n, n, "mixture SQP ridge Hessian", |i, j| h[(i, j)])?;
        for i in 0..n {
            corrected[(i, i)] += ridge;
        }
        let factors = try_mixture_matrix(n, n, "mixture SQP Cholesky", |i, j| corrected[(i, j)])?;
        if let Some(cholesky) = factors.cholesky() {
            let norm = |m: &DMatrix<f64>| {
                m.column_iter()
                    .map(|c| c.iter().map(|x| x.abs()).sum::<f64>())
                    .fold(0.0_f64, f64::max)
            };
            let mut inverse = try_mixture_matrix(n, n, "mixture SQP Hessian inverse", |i, j| {
                if i == j { 1.0 } else { 0.0 }
            })?;
            cholesky.solve_mut(&mut inverse);
            let reciprocal = 1.0 / (norm(&corrected) * norm(&inverse));
            drop(inverse);
            solved = Some(if reciprocal < f64::EPSILON || !reciprocal.is_finite() {
                let svd = nalgebra::linalg::SVD::try_new(corrected, true, true, f64::EPSILON, 1000)
                    .ok_or_else(|| {
                        scoped_invalid_options(
                            ERROR_SCOPE,
                            "quadratic direction SVD did not converge",
                        )
                    })?;
                let cutoff = n as f64 * f64::EPSILON * svd.singular_values[0];
                svd.solve(rhs, cutoff).map_err(|e| {
                    scoped_invalid_options(ERROR_SCOPE, format!("quadratic direction SVD: {e}"))
                })?
            } else {
                cholesky.solve(rhs)
            });
            break;
        }
        let next = if ridge <= 0.0 {
            1e-15
        } else {
            ridge * ridge_multiplier
        };
        if next <= ridge {
            break;
        }
        ridge = next;
    }
    solved.ok_or_else(|| {
        scoped_invalid_options(
            ERROR_SCOPE,
            "quadratic subproblem is not numerically positive definite",
        )
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn directions_discard_only_numerically_singular_axes() {
        // H * [1, 0] = [1, 1] exactly. The nearly rank-one case admits a
        // Cholesky factorization but truncates its weak axis: the independent
        // minimum-norm rank-one solution is [1/2, 1/2]. Diagonal entries are
        // large enough that the initial ridge correction is zero.
        for (delta, expected) in [(1e-10, [1.0, 0.0]), (2.0 * f64::EPSILON, [0.5, 0.5])] {
            let h = DMatrix::from_row_slice(2, 2, &[1.0, 1.0, 1.0, 1.0 + delta]);
            let rhs = DVector::from_vec(vec![1.0, 1.0]);
            let actual = search_direction(&h, &rhs, 10.0).unwrap();
            for i in 0..2 {
                assert!((actual[i] - expected[i]).abs() < 2e-14);
            }
        }
    }
}
