//! Independent cycling RRE/MPE and squared polynomial extrapolation.
//! The polynomials minimize a residual or annihilate its last difference;
//! squaring convolves their coefficients (Roland/Varadhan/Frangakis 2007).
use super::{
    invalid,
    squarem::{Iteration, PolynomialExtrapolation, SquaremOptions},
};
use crate::DeseqError;
use nalgebra::{DMatrix, DVector};
fn norm(a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .zip(b)
        .map(|(a, b)| (a - b) * (a - b))
        .sum::<f64>()
        .sqrt()
}
fn coefficients(states: &[Vec<f64>], o: SquaremOptions) -> Option<Vec<f64>> {
    let k = o.polynomial_order;
    let p = states[0].len();
    let u = DMatrix::from_fn(p, k + 1, |i, j| states[j + 1][i] - states[j][i]);
    let mut c = match o.polynomial_method {
        PolynomialExtrapolation::ReducedRank => {
            let gram = u.transpose() * &u;
            gram.col_piv_qr()
                .solve(&DVector::from_element(k + 1, 1.0))?
                .as_slice()
                .to_vec()
        }
        PolynomialExtrapolation::MinimalPolynomial => {
            let a = u.columns(0, k).into_owned();
            let b = -u.column(k).into_owned();
            let qr = a.col_piv_qr();
            let mut rhs = qr.q().transpose() * b;
            let r = qr.r();
            r.solve_upper_triangular_mut(&mut rhs).then_some(())?;
            qr.p().inv_permute_rows(&mut rhs);
            let mut c = rhs.as_slice().to_vec();
            c.push(1.0);
            c
        }
    };
    let total = c.iter().sum::<f64>();
    if !total.is_finite() || total.abs() < 1e-7 || c.iter().any(|v| !v.is_finite()) {
        return None;
    }
    for v in &mut c {
        *v /= total;
    }
    Some(c)
}
pub(super) fn run(
    initial: Vec<f64>,
    o: SquaremOptions,
    tolerance: f64,
    update: impl Fn(&[f64]) -> Result<Vec<f64>, DeseqError>,
    objective: impl Fn(&[f64]) -> f64,
) -> Result<Iteration, DeseqError> {
    let k = o.polynomial_order;
    if k > initial.len() {
        return Err(invalid("polynomial order exceeds retained component count"));
    }
    let capacity = k
        .checked_mul(2)
        .and_then(|v| v.checked_add(2))
        .ok_or_else(|| invalid("polynomial dimensions overflow"))?;
    capacity
        .checked_mul(initial.len())
        .ok_or_else(|| invalid("polynomial state dimensions overflow"))?;
    let mut score = objective(&initial);
    let mut state = update(&initial)?;
    let mut residual = norm(&state, &initial);
    let mut calls = 1;
    let mut iterations = 0;
    while calls < o.max_evaluations && residual > tolerance {
        iterations += 1;
        let mut states = Vec::with_capacity(capacity);
        states.push(state);
        for _ in 0..=k {
            states.push(update(states.last().unwrap())?);
            calls += 1;
        }
        let mut fallback_residual = norm(&states[k + 1], &states[k]);
        let Some(c) = coefficients(&states, o) else {
            state = states.pop().unwrap();
            residual = fallback_residual;
            iterations += 1;
            continue;
        };
        let mut proposal = vec![0.0; states[0].len()];
        if o.square {
            for _ in 1..k {
                states.push(update(states.last().unwrap())?);
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
            states.pop().unwrap()
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
