//! ashr-compatible policy of the mixture SQP solver.
//!
//! Written from `docs/provenance/specs/mixture_sqp_compatibility.md`, which
//! follows Kim, Carbonetto, Stephens and Anitescu (2020, J. Comput. Graph.
//! Stat. 29:261-273) for the mixture problem and its SQP iteration, and
//! Nocedal and Wright (2006, Numerical Optimization, section 16.5) for the
//! primal active-set method of the convex quadratic subproblem. The EM floor
//! was fixed by black-box runs of R's mixsqp 0.3-54
//! (`scripts/provenance/probe_mixsqp_em_floor.R`, output in
//! `docs/provenance/observations/mixsqp_em_floor_r462_macos.txt`).
//!
//! The mixture problem is the one of the paper, section 2: rows of the
//! (row-scaled) likelihood matrix with weights `a_i / mass`, followed by one
//! pseudo-row `e_k` with weight `b_k / mass` for every positive Dirichlet
//! count `b_k`, so that the objective is
//! `f(x) = sum(x) - sum_i w_i log(L_i x + eps)` on the nonnegative orthant
//! (the relaxed problem of section 3).
use super::{DenseMatrix, ERROR_SCOPE, Error, Problem, affine, try_dense_matrix};
use crate::errors::scoped_invalid_options;
use crate::shrinkage::MixtureSqpOptions;
use rcompat_numeric::math::dense::{Cholesky, SymmetricEigen, dot, tr_gemv_to};

/// Absolute amount added to every proportion before an EM update, and, scaled
/// by the largest updated proportion, after it. Observed in mixsqp's EM
/// iterations: a proportion that starts at zero, or that the likelihood
/// cannot support, is kept strictly positive.
const EM_FLOOR: f64 = 1e-15;

impl Problem {
    /// Hessian of the relaxed objective, `sum_i w_i L_i L_i' / (L_i x + eps)^2`
    /// (paper, section 3), over the pseudo-rows and then the data rows. Only
    /// the upper triangle is accumulated, so the result is exactly symmetric.
    pub(super) fn ashr_hessian_into(&self, weights: &[f64], hessian: &mut DenseMatrix) {
        let k = self.components;
        hessian.fill(0.0);
        for (i, (&count, &w)) in self.prior_counts.iter().zip(weights).enumerate() {
            if count != 0.0 {
                let denominator = w + self.offset;
                hessian[(i, i)] += (count / self.mass) / (denominator * denominator);
            }
        }
        let mut scaled = vec![0.0; k];
        for (row, &a) in self
            .likelihoods
            .chunks_exact(k)
            .zip(&self.observation_weights)
        {
            if a == 0.0 {
                continue;
            }
            let marginal = row.iter().zip(weights).map(|(l, w)| l * w).sum::<f64>() + self.offset;
            let coefficient = (a / self.mass) / (marginal * marginal);
            for (s, &l) in scaled.iter_mut().zip(row) {
                *s = coefficient * l;
            }
            for j in 0..k {
                let lj = row[j];
                if lj == 0.0 {
                    continue;
                }
                for i in 0..=j {
                    hessian[(i, j)] += scaled[i] * lj;
                }
            }
        }
        for j in 0..k {
            for i in 0..j {
                hessian[(j, i)] = hessian[(i, j)];
            }
        }
    }

    /// One EM update of the proportions (paper, section 2): with `z = x + d`,
    /// `x_k <- z_k sum_i w_i L_ik / (L_i z)` over data and pseudo-rows, then
    /// `x <- x + d max(x)`, where `d` is [`EM_FLOOR`]. A pseudo-row
    /// contributes exactly `b_k / mass`. As in the black-box runs, the `eps`
    /// control does not enter the update.
    pub(super) fn ashr_expectation_step(&self, weights: &[f64]) -> Result<Vec<f64>, Error> {
        let k = self.components;
        let shifted = weights.iter().map(|w| w + EM_FLOOR).collect::<Vec<_>>();
        let mut factors = vec![0.0; k];
        for (row, &a) in self
            .likelihoods
            .chunks_exact(k)
            .zip(&self.observation_weights)
        {
            if a == 0.0 {
                continue;
            }
            let marginal = row.iter().zip(&shifted).map(|(l, z)| l * z).sum::<f64>();
            if marginal.is_nan() || marginal <= 0.0 || marginal.is_infinite() {
                return Err(scoped_invalid_options(
                    ERROR_SCOPE,
                    "EM update has no positive likelihood",
                ));
            }
            let coefficient = (a / self.mass) / marginal;
            for (f, &l) in factors.iter_mut().zip(row) {
                *f += coefficient * l;
            }
        }
        let mut next = shifted
            .iter()
            .zip(&factors)
            .zip(&self.prior_counts)
            .map(|((z, f), &count)| z * f + count / self.mass)
            .collect::<Vec<_>>();
        let largest = next.iter().copied().fold(0.0, f64::max);
        for value in &mut next {
            *value += EM_FLOOR * largest;
        }
        if next.iter().any(|v| !v.is_finite()) {
            return Err(scoped_invalid_options(
                ERROR_SCOPE,
                "EM update exceeds the representable range",
            ));
        }
        Ok(next)
    }
}

/// Primal active-set method for `min 1/2 y'Hy + c'y` subject to `y >= 0`
/// (paper, Algorithm 2; Nocedal and Wright, Algorithm 16.3), from `start`.
///
/// The working set holds the coordinates fixed at zero. Each iteration solves
/// the equality-constrained subproblem on the free coordinates for a search
/// direction `p`. A direction whose entries are all within
/// `search_direction_threshold` counts as zero; the multipliers of the
/// working set (the gradient `Hy + c` there) are then examined, and the
/// coordinate with the most negative multiplier below `-active_tolerance`
/// is released, or the current point is returned. A nonzero direction is
/// followed up to the first blocking bound, which joins the working set.
pub(super) fn quadratic(
    h: &DenseMatrix,
    c: &[f64],
    start: &[f64],
    max_iterations: usize,
    options: MixtureSqpOptions,
) -> Result<Vec<f64>, Error> {
    let k = c.len();
    if h.nrows() != k || h.ncols() != k || start.len() != k {
        return Err(scoped_invalid_options(
            ERROR_SCOPE,
            "quadratic subproblem dimensions disagree",
        ));
    }
    let mut y = start.iter().map(|v| v.max(0.0)).collect::<Vec<_>>();
    let mut free = y.iter().map(|v| *v > 0.0).collect::<Vec<_>>();
    // Without a direction threshold, a full step solves the subproblem on the
    // current free set exactly, so the next direction is taken to be zero.
    let threshold = options.search_direction_threshold;
    let mut solved_free_set = false;
    for _ in 0..max_iterations {
        let gradient = affine(h, &y, c)?;
        let indices = (0..k).filter(|&i| free[i]).collect::<Vec<_>>();
        let direction = if indices.is_empty() || (threshold.is_none() && solved_free_set) {
            None
        } else {
            let reduced = try_dense_matrix(
                indices.len(),
                indices.len(),
                "mixture SQP active Hessian",
                |i, j| h[(indices[i], indices[j])],
            )?;
            let rhs = indices.iter().map(|&i| -gradient[i]).collect::<Vec<_>>();
            let p = reduced_solve(reduced, rhs, options.ridge_multiplier)?;
            let limit = threshold.unwrap_or(0.0);
            if p.iter().all(|v| v.abs() <= limit) {
                None
            } else {
                Some(p)
            }
        };
        let Some(p) = direction else {
            let worst = (0..k)
                .filter(|&i| !free[i])
                .min_by(|&i, &j| gradient[i].total_cmp(&gradient[j]));
            match worst {
                Some(i) if gradient[i] < -options.active_tolerance => {
                    free[i] = true;
                    solved_free_set = false;
                    continue;
                }
                _ => return Ok(y),
            }
        };
        let mut step = 1.0_f64;
        let mut blocking = None;
        for (n, &i) in indices.iter().enumerate() {
            if p[n] < 0.0 {
                let limit = -y[i] / p[n];
                if limit < step {
                    step = limit;
                    blocking = Some(i);
                }
            }
        }
        for (n, &i) in indices.iter().enumerate() {
            y[i] = (y[i] + step * p[n]).max(0.0);
        }
        if let Some(i) = blocking {
            y[i] = 0.0;
            free[i] = false;
            solved_free_set = false;
        } else {
            solved_free_set = true;
        }
    }
    Ok(y)
}

/// Solves the reduced system `A p = r` for a symmetric positive
/// semidefinite `A`.
///
/// A Cholesky factor is used when its reciprocal condition estimate exceeds
/// machine precision. A numerically singular system instead gets the
/// minimum-norm solution restricted to the eigen-directions whose eigenvalues
/// exceed `n * eps * max|lambda|`; the others (including tiny negative ones
/// from rounding) are dropped. Black-box runs of mixsqp on an ill-conditioned
/// ashr problem (`scripts/provenance/probe_mixsqp_singular_direction.R`)
/// agree with this minimum-norm direction and not with a ridge-regularised
/// one, whose huge components send the active set to a different vertex.
/// Only if the eigen-decomposition fails is a ridge, growing by
/// `ridge_multiplier`, added until a well-conditioned factor exists. Any
/// remaining nonfinite entry is set to zero.
fn reduced_solve(a: DenseMatrix, rhs: Vec<f64>, ridge_multiplier: f64) -> Result<Vec<f64>, Error> {
    let n = a.nrows();
    let norm = one_norm(&a);
    if !norm.is_finite() || norm <= 0.0 {
        return Ok(vec![0.0; n]);
    }
    if let Some(p) = conditioned_cholesky_solve(&a, &rhs, norm, 0.0) {
        return Ok(p);
    }
    if let Some(eigen) = SymmetricEigen::new(a.clone())? {
        let largest = eigen
            .eigenvalues
            .iter()
            .copied()
            .fold(0.0_f64, |m, v| m.max(v.abs()));
        let cutoff = largest * n as f64 * f64::EPSILON;
        let mut projected = vec![0.0; n];
        tr_gemv_to(&eigen.eigenvectors, &rhs, &mut projected)?;
        let mut p = vec![0.0; n];
        for (j, &lambda) in eigen.eigenvalues.iter().enumerate() {
            if lambda > cutoff && lambda.is_finite() {
                let coefficient = projected[j] / lambda;
                for (value, &vector) in p.iter_mut().zip(eigen.eigenvectors.column(j)) {
                    *value += coefficient * vector;
                }
            }
        }
        if p.iter().all(|v| v.is_finite()) {
            return Ok(p);
        }
    }
    let scale = (0..n).map(|i| a[(i, i)].abs()).fold(0.0, f64::max);
    let mut ridge = f64::EPSILON * scale.max(f64::MIN_POSITIVE);
    for _ in 0..40 {
        if let Some(p) = conditioned_cholesky_solve(&a, &rhs, norm, ridge) {
            return Ok(p);
        }
        ridge *= ridge_multiplier;
        if !ridge.is_finite() {
            break;
        }
    }
    Err(scoped_invalid_options(
        ERROR_SCOPE,
        "quadratic subproblem is not numerically positive semidefinite",
    ))
}

/// Cholesky solve of `(A + ridge I) p = r`, accepted only when the
/// reciprocal condition estimate exceeds machine precision.
fn conditioned_cholesky_solve(
    a: &DenseMatrix,
    rhs: &[f64],
    norm: f64,
    ridge: f64,
) -> Option<Vec<f64>> {
    let mut shifted = a.clone();
    for i in 0..shifted.nrows() {
        shifted[(i, i)] += ridge;
    }
    let factor = Cholesky::new(shifted)?;
    let solve = |v: &[f64]| {
        let mut x = v.to_vec();
        factor.solve_in_place(&mut x);
        x
    };
    let inverse_norm = inverse_one_norm(solve, a.nrows());
    let rcond = 1.0 / ((norm + ridge) * inverse_norm);
    if !rcond.is_finite() || rcond <= f64::EPSILON {
        return None;
    }
    let p = solve(rhs);
    p.iter().all(|v| v.is_finite()).then_some(p)
}

fn one_norm(a: &DenseMatrix) -> f64 {
    (0..a.ncols())
        .map(|j| a.column(j).iter().map(|v| v.abs()).sum::<f64>())
        .fold(0.0, f64::max)
}

/// Estimate of `||A^{-1}||_1` for symmetric `A` from solves alone (Hager's
/// method as refined by Higham 1988, ACM TOMS 14:381-396, Algorithm 4.1,
/// including its alternative test vector).
fn inverse_one_norm(solve: impl Fn(&[f64]) -> Vec<f64>, n: usize) -> f64 {
    if n == 0 {
        return 0.0;
    }
    let mut x = vec![1.0 / n as f64; n];
    let mut estimate = 0.0_f64;
    let mut previous = None;
    for iteration in 0..5 {
        let y = solve(&x);
        let norm = y.iter().map(|v| v.abs()).sum::<f64>();
        if !norm.is_finite() {
            return f64::INFINITY;
        }
        if iteration > 0 && norm <= estimate {
            break;
        }
        estimate = norm;
        let sign = y
            .iter()
            .map(|&v| if v >= 0.0 { 1.0 } else { -1.0 })
            .collect::<Vec<_>>();
        let z = solve(&sign);
        let (j, largest) = z
            .iter()
            .enumerate()
            .fold((0, f64::NEG_INFINITY), |(bj, bv), (j, v)| {
                if v.abs() > bv { (j, v.abs()) } else { (bj, bv) }
            });
        if !largest.is_finite() {
            return f64::INFINITY;
        }
        if iteration > 0 && (previous == Some(j) || largest <= dot(&z, &x)) {
            break;
        }
        previous = Some(j);
        x = vec![0.0; n];
        x[j] = 1.0;
    }
    let alternative = (0..n)
        .map(|i| {
            let sign = if i % 2 == 0 { 1.0 } else { -1.0 };
            sign * (1.0 + i as f64 / (n.max(2) - 1) as f64)
        })
        .collect::<Vec<_>>();
    let tail = solve(&alternative).iter().map(|v| v.abs()).sum::<f64>();
    estimate.max(2.0 * tail / (3.0 * n as f64))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn problem(likelihoods: Vec<f64>, a: Vec<f64>, b: Vec<f64>, offset: f64) -> Problem {
        let components = b.len();
        let mass = a.iter().sum::<f64>() + b.iter().sum::<f64>();
        let rows = a.len();
        Problem {
            likelihoods,
            row_maxima: vec![0.0; rows],
            observation_weights: a,
            prior_counts: b,
            mass,
            components,
            offset,
            low_rank: false,
        }
    }

    #[test]
    fn hessian_sums_scaled_outer_products_and_prior_rows() {
        let p = problem(
            vec![1.0, 0.5, 0.25, 1.0],
            vec![1.0, 3.0],
            vec![2.0, 0.0],
            0.1,
        );
        let x = [0.4, 0.6];
        let mut h = DenseMatrix::zeros(2, 2).unwrap();
        p.ashr_hessian_into(&x, &mut h);
        let m1: f64 = 0.4 + 0.3 + 0.1;
        let m2: f64 = 0.1 + 0.6 + 0.1;
        let mass = 6.0;
        let expected = [
            (1.0 / mass) / m1.powi(2) + (3.0 / mass) * 0.0625 / m2.powi(2) + (2.0 / mass) / 0.25,
            (1.0 / mass) * 0.5 / m1.powi(2) + (3.0 / mass) * 0.25 / m2.powi(2),
            (1.0 / mass) * 0.25 / m1.powi(2) + (3.0 / mass) / m2.powi(2),
        ];
        assert!((h[(0, 0)] - expected[0]).abs() < 1e-14);
        assert!((h[(0, 1)] - expected[1]).abs() < 1e-14);
        assert_eq!(h[(0, 1)], h[(1, 0)]);
        assert!((h[(1, 1)] - expected[2]).abs() < 1e-14);
    }

    #[test]
    fn expectation_step_floors_and_keeps_prior_mass() {
        // One row L = (1, 1): the EM map is the identity apart from the floor.
        let p = problem(vec![1.0, 1.0], vec![1.0], vec![0.0, 0.0], 1e-6);
        let next = p.ashr_expectation_step(&[0.0, 1.0]).unwrap();
        assert!((next[0] - 2e-15).abs() < 1e-25);
        // (1 + d) / (1 + 2d) + d * max: unchanged to rounding.
        assert!((next[1] - 1.0).abs() < 1e-15);
        // Rows (2, 1) and (1, 3), weights 1 and 1, prior counts (2, 0).
        let p = problem(
            vec![2.0, 1.0, 1.0, 3.0],
            vec![1.0, 1.0],
            vec![2.0, 0.0],
            0.0,
        );
        let x = [0.5, 0.5];
        let next = p.ashr_expectation_step(&x).unwrap();
        let em = |k: usize| {
            let z = x[k] + EM_FLOOR;
            let m1 = 2.0 * (0.5 + EM_FLOOR) + (0.5 + EM_FLOOR);
            let m2 = (0.5 + EM_FLOOR) + 3.0 * (0.5 + EM_FLOOR);
            let l = [[2.0, 1.0], [1.0, 3.0]];
            z * (0.25 * l[0][k] / m1 + 0.25 * l[1][k] / m2) + [0.5, 0.0][k]
        };
        let largest = em(0).max(em(1));
        for (k, value) in next.iter().enumerate() {
            assert!((value - (em(k) + EM_FLOOR * largest)).abs() < 1e-15);
        }
        assert!((next.iter().sum::<f64>() - 1.0).abs() < 1e-14);
    }

    #[test]
    fn active_set_matches_hand_solutions() {
        let options = MixtureSqpOptions {
            search_direction_threshold: Some(1e-14),
            ..Default::default()
        };
        // Interior minimum: H = [[2, 0], [0, 4]], c = (-2, -4) gives y = (1, 1).
        let h = DenseMatrix::from_row_major(2, 2, &[2.0, 0.0, 0.0, 4.0]).unwrap();
        let c = vec![-2.0, -4.0];
        let y = quadratic(&h, &c, &[0.5, 0.5], 10, options).unwrap();
        assert!((y[0] - 1.0).abs() < 1e-15 && (y[1] - 1.0).abs() < 1e-15);
        // Bound active at the optimum: c = (-2, 4) gives y = (1, 0).
        let c = vec![-2.0, 4.0];
        let y = quadratic(&h, &c, &[0.5, 0.5], 10, options).unwrap();
        assert!((y[0] - 1.0).abs() < 1e-15 && y[1] == 0.0);
        // A zero coordinate with a negative multiplier is released.
        let h = DenseMatrix::from_row_major(2, 2, &[2.0, 1.0, 1.0, 2.0]).unwrap();
        let c = vec![-3.0, -3.0];
        let y = quadratic(&h, &c, &[1.0, 0.0], 10, options).unwrap();
        assert!((y[0] - 1.0).abs() < 1e-14 && (y[1] - 1.0).abs() < 1e-14);
        // The same without a direction threshold.
        let y = quadratic(&h, &c, &[1.0, 0.0], 10, MixtureSqpOptions::default()).unwrap();
        assert!((y[0] - 1.0).abs() < 1e-14 && (y[1] - 1.0).abs() < 1e-14);
    }

    #[test]
    fn singular_reduced_hessian_stays_finite() {
        let options = MixtureSqpOptions {
            search_direction_threshold: Some(1e-14),
            ..Default::default()
        };
        // Rank one: every minimiser has y0 + y1 = 1 when c = (-1, -1).
        let h = DenseMatrix::from_row_major(2, 2, &[1.0, 1.0, 1.0, 1.0]).unwrap();
        let c = vec![-1.0, -1.0];
        let y = quadratic(&h, &c, &[0.25, 0.25], 10, options).unwrap();
        assert!(y.iter().all(|v| v.is_finite() && *v >= 0.0));
        assert!((y[0] + y[1] - 1.0).abs() < 1e-6);
        let p = reduced_solve(DenseMatrix::zeros(2, 2).unwrap(), vec![1.0, 1.0], 10.0).unwrap();
        assert!(p.iter().all(|v| *v == 0.0));
    }

    #[test]
    fn inverse_norm_estimate_is_exact_for_diagonal_matrices() {
        let d = [4.0, 0.5, 2.0];
        let estimate = inverse_one_norm(|v| v.iter().zip(&d).map(|(x, d)| x / d).collect(), 3);
        assert!((estimate - 2.0).abs() < 1e-15);
    }
}
