//! Independent cycling RRE/MPE and squared polynomial extrapolation.
//! The polynomials minimize a residual or annihilate its last difference;
//! squaring convolves their coefficients (Roland/Varadhan/Frangakis 2007).
use super::{
    ERROR_SCOPE,
    squarem::{Iteration, PolynomialExtrapolation, SquaremOptions},
    try_dense_matrix,
};
use crate::errors::scoped_invalid_options;
use crate::{Error, matrix::try_alloc_matrix};
use nalgebra::DVector;
fn norm(a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .zip(b)
        .map(|(a, b)| (a - b) * (a - b))
        .sum::<f64>()
        .sqrt()
}
fn coefficients(states: &[Vec<f64>], o: SquaremOptions) -> Result<Option<Vec<f64>>, Error> {
    let k = o.polynomial_order;
    let p = states[0].len();
    let u = try_dense_matrix(p, k + 1, "mixture polynomial state differences", |i, j| {
        states[j + 1][i] - states[j][i]
    })?;
    let mut c = match o.polynomial_method {
        PolynomialExtrapolation::ReducedRank => {
            let transpose =
                try_dense_matrix(k + 1, p, "mixture polynomial transpose", |i, j| u[(j, i)])?;
            let mut gram =
                try_dense_matrix(k + 1, k + 1, "mixture polynomial Gram matrix", |_, _| 0.0)?;
            transpose.mul_to(&u, &mut gram);
            let Some(solution) = gram.col_piv_qr().solve(&DVector::from_element(k + 1, 1.0)) else {
                return Ok(None);
            };
            solution.as_slice().to_vec()
        }
        PolynomialExtrapolation::MinimalPolynomial => {
            let a = try_dense_matrix(p, k, "mixture polynomial QR design", |i, j| u[(i, j)])?;
            let b = -u.column(k).into_owned();
            let qr = a.col_piv_qr();
            let mut rhs = qr.q().transpose() * b;
            let r = qr.r();
            if !r.solve_upper_triangular_mut(&mut rhs) {
                return Ok(None);
            }
            qr.p().inv_permute_rows(&mut rhs);
            let mut c = rhs.as_slice().to_vec();
            c.push(1.0);
            c
        }
    };
    let total = c.iter().sum::<f64>();
    if !total.is_finite() || total.abs() < 1e-7 || c.iter().any(|v| !v.is_finite()) {
        return Ok(None);
    }
    for v in &mut c {
        *v /= total;
    }
    Ok(Some(c))
}
pub(super) fn run(
    initial: Vec<f64>,
    o: SquaremOptions,
    tolerance: f64,
    update: impl Fn(&[f64]) -> Result<Vec<f64>, Error>,
    mut objective: impl FnMut(&[f64]) -> f64,
) -> Result<Iteration, Error> {
    let k = o.polynomial_order;
    if k > initial.len() {
        return Err(Error::EmPolynomialOrderExceeded {
            order: k,
            components: initial.len(),
        });
    }
    let capacity = k
        .checked_mul(2)
        .and_then(|v| v.checked_add(2))
        .ok_or_else(|| scoped_invalid_options(ERROR_SCOPE, "polynomial dimensions overflow"))?;
    capacity.checked_mul(initial.len()).ok_or_else(|| {
        scoped_invalid_options(ERROR_SCOPE, "polynomial state dimensions overflow")
    })?;
    let mut score = objective(&initial);
    let mut state = update(&initial)?;
    let mut residual = norm(&state, &initial);
    let mut calls = 1;
    let mut iterations = 0;
    while calls < o.max_evaluations && residual > tolerance {
        iterations += 1;
        let mut states = try_alloc_matrix(capacity, 1, "mixture polynomial states")?;
        states.push(state);
        for _ in 0..=k {
            states.push(update(
                states
                    .last()
                    .ok_or_else(|| crate::errors::invariant("value present"))?,
            )?);
            calls += 1;
        }
        let mut fallback_residual = norm(&states[k + 1], &states[k]);
        let Some(c) = coefficients(&states, o)? else {
            state = states
                .pop()
                .ok_or_else(|| crate::errors::invariant("value present"))?;
            residual = fallback_residual;
            iterations += 1;
            continue;
        };
        let p = states[0].len();
        let mut proposal = try_alloc_matrix(p, 1, "mixture polynomial proposal")?;
        proposal.resize(p, 0.0);
        if o.square {
            for _ in 1..k {
                states.push(update(
                    states
                        .last()
                        .ok_or_else(|| crate::errors::invariant("value present"))?,
                )?);
                calls += 1;
            }
            for i in 0..=k {
                for j in 0..=k {
                    for (v, &p) in proposal.iter_mut().zip(&states[i + j]) {
                        *v += c[i] * c[j] * p;
                    }
                }
            }
            // Keep the established fallback probe one update behind the last
            // squared sequence point, independently of the extrapolated state.
            fallback_residual = norm(&states[2 * k - 1], &states[2 * k - 2]);
        } else {
            for i in 0..=k {
                for (v, &p) in proposal.iter_mut().zip(&states[i]) {
                    *v += c[i] * p;
                }
            }
        }
        calls += 1;
        let stabilized = update(&proposal);
        let accepted = match stabilized {
            Ok(candidate) if candidate.iter().all(|v| v.is_finite()) => {
                residual = norm(&candidate, &proposal);
                let value = objective(&candidate);
                if value.is_nan() || value > score + o.objective_increase {
                    None
                } else {
                    score = value;
                    Some(candidate)
                }
            }
            _ => None,
        };
        state = if let Some(candidate) = accepted {
            candidate
        } else {
            residual = fallback_residual;
            states
                .pop()
                .ok_or_else(|| crate::errors::invariant("value present"))?
        };
    }
    Ok(Iteration {
        state,
        converged: calls < o.max_evaluations,
        evaluations: calls,
        iterations,
        objective: score,
        residual,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn polynomial_extrapolation_solves_an_independent_affine_fixed_point() {
        let target = [0.2, -0.7, 1.4];
        let contraction = [0.2, 0.5, 0.8];
        for square in [false, true] {
            for method in [
                PolynomialExtrapolation::MinimalPolynomial,
                PolynomialExtrapolation::ReducedRank,
            ] {
                let update = |x: &[f64]| {
                    Ok(x.iter()
                        .zip(target)
                        .zip(contraction)
                        .map(|((x, t), a)| t + a * (x - t))
                        .collect())
                };
                let objective =
                    |x: &[f64]| x.iter().zip(target).map(|(x, t)| (x - t) * (x - t)).sum();
                let fit = run(
                    vec![-1.0, 2.0, 0.0],
                    SquaremOptions {
                        polynomial_order: if method == PolynomialExtrapolation::MinimalPolynomial {
                            3
                        } else {
                            2
                        },
                        polynomial_method: method,
                        square,
                        ..Default::default()
                    },
                    1e-10,
                    update,
                    objective,
                )
                .unwrap();
                assert!(fit.converged);
                assert!(
                    fit.state
                        .iter()
                        .zip(target)
                        .all(|(x, t)| (x - t).abs() < 1e-8),
                    "{method:?} {square}: {:?}",
                    fit.state
                );
            }
        }
    }
}
