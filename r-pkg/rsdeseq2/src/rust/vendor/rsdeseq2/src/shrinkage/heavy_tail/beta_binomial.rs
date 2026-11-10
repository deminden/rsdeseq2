//! Beta-binomial logit likelihood, analytic score and native fitting policies.
//! The parameter-dependent likelihood is log (a)_y + log (b)_(n-y), where
//! a=theta*p, b=theta*(1-p). Gamma/psi recurrences follow DLMF 5.5/5.11.
use super::{
    prior::{CauchyCoefficientPrior, invalid},
    regression::{RegressionOptimizer, RegressionPosteriorOptions},
};
use crate::probability::digamma_positive as psi;
use crate::{DeseqError, matrix::MatrixView, probability::stirling_correction};
use nalgebra::DMatrix;
use rcompat_lbfgsb::{Bounds, OptimControl, OptimResult, optim_lbfgsb_with_gradient};

fn psi_difference(x: f64, k: f64) -> f64 {
    if k == 0.0 {
        return 0.0;
    }
    if x < 32.0 {
        return psi(x + k) - psi(x);
    }
    let logarithm = (k / x).ln_1p();
    let inverse = 1.0 / x;
    let mut out = logarithm - 0.5 * inverse * (-logarithm).exp_m1();
    for (i, coefficient) in [
        1.0 / 12.0,
        -1.0 / 120.0,
        1.0 / 252.0,
        -1.0 / 240.0,
        1.0 / 132.0,
        -691.0 / 32760.0,
    ]
    .iter()
    .enumerate()
    {
        let power = 2 * (i + 1) as i32;
        out -= coefficient * inverse.powi(power) * (-f64::from(power) * logarithm).exp_m1();
    }
    out
}
pub(super) fn log_rising(x: f64, k: f64) -> f64 {
    if k == 0.0 {
        return 0.0;
    }
    if x < 32.0 {
        return libm::lgamma(x + k) - libm::lgamma(x);
    }
    // Algebraically subtract Stirling terms before evaluation, retaining small
    // increments even when x+k rounds to x or the two log-gammas nearly agree.
    let logarithm = (k / x).ln_1p();
    k * x.ln() + (x + k - 0.5) * logarithm - k + stirling_correction(x + k) - stirling_correction(x)
}
pub(super) struct Row<'a> {
    pub design: MatrixView<'a, f64>,
    pub successes: &'a [f64],
    pub trials: &'a [f64],
    pub weights: &'a [f64],
    pub concentration: f64,
    pub prior: &'a CauchyCoefficientPrior,
}
impl Row<'_> {
    pub fn objective_gradient(&self, beta: &[f64], capped: bool) -> (f64, Vec<f64>) {
        let mut gradient = vec![0.0; beta.len()];
        let mut value = self.prior.contribution(beta, &mut gradient);
        // Source nuisance repeats count in the objective, while its analytic
        // score assigns each coefficient once. Keep that observable policy.
        for (i, &j) in self.prior.nuisance.iter().enumerate() {
            if self.prior.nuisance[..i].contains(&j) {
                value += beta[j] * beta[j]
                    / (2.0 * self.prior.nuisance_scale * self.prior.nuisance_scale);
            }
        }
        for i in 0..self.successes.len() {
            let eta = beta
                .iter()
                .enumerate()
                .map(|(j, b)| self.design.at(i, j) * b)
                .sum::<f64>();
            let odds = (-eta).exp();
            let odds = if capped {
                odds.clamp(0.001, 999.0)
            } else {
                odds
            };
            let p = 1.0 / (1.0 + odds);
            let a = self.concentration * p;
            let b = self.concentration - a;
            if !(a > 0.0 && b > 0.0 && a.is_finite() && b.is_finite()) {
                return (f64::INFINITY, vec![f64::NAN; beta.len()]);
            }
            let y = self.successes[i];
            let failures = self.trials[i] - y;
            value -= self.weights[i] * (log_rising(a, y) + log_rising(b, failures));
            let score = self.weights[i]
                * self.concentration
                * p
                * (1.0 - p)
                * (psi_difference(a, y) - psi_difference(b, failures));
            for (j, g) in gradient.iter_mut().enumerate() {
                *g -= self.design.at(i, j) * score;
            }
        }
        (value, gradient)
    }
    pub fn hessian(&self, beta: &[f64], step: f64) -> DMatrix<f64> {
        let p = beta.len();
        let mut out = DMatrix::zeros(p, p);
        for j in 0..p {
            let mut plus = beta.to_vec();
            let mut minus = beta.to_vec();
            plus[j] += step;
            minus[j] -= step;
            let gp = self.objective_gradient(&plus, false).1;
            let gm = self.objective_gradient(&minus, false).1;
            for i in 0..p {
                out[(i, j)] = (gp[i] - gm[i]) / (2.0 * step);
            }
        }
        for j in 0..p {
            for i in 0..j {
                let mean = (out[(i, j)] + out[(j, i)]) / 2.0;
                out[(i, j)] = mean;
                out[(j, i)] = mean;
            }
        }
        out
    }
    pub fn optimize(
        &self,
        start: Vec<f64>,
        options: &RegressionPosteriorOptions,
    ) -> Result<OptimResult, DeseqError> {
        let constant = -self.objective_gradient(&start, false).0 - 1.0;
        let f = |b: &[f64]| {
            let (value, gradient) = self.objective_gradient(b, false);
            (value + constant, gradient)
        };
        match options.optimizer {
            RegressionOptimizer::Bfgs => Ok(super::bfgs::minimize(f, start)),
            RegressionOptimizer::Lbfgsb => {
                let p = start.len();
                let bounds = Bounds::new(vec![options.bounds.0; p], vec![options.bounds.1; p])
                    .map_err(|e| invalid(&e.to_string()))?;
                optim_lbfgsb_with_gradient(
                    start,
                    bounds,
                    |b| f(b).0,
                    |b| f(b).1,
                    OptimControl::default_for_dimension(p),
                )
                .map_err(|e| invalid(&e.to_string()))
            }
        }
    }
    pub fn prefit(&self, start: Vec<f64>) -> Result<OptimResult, DeseqError> {
        let mut anchor = vec![0.0; start.len()];
        anchor[0] = 3.0;
        let scale = self.objective_gradient(&anchor, false).0.abs().max(100.0);
        let (value, gradient) = self.objective_gradient(&start, true);
        if !scale.is_finite() || !value.is_finite() || gradient.iter().any(|v| !v.is_finite()) {
            return Err(invalid(
                "initial beta-binomial posterior exceeds representable range",
            ));
        }
        Ok(super::unbounded_lbfgs::minimize_with_tolerances(
            |beta| {
                let (v, g) = self.objective_gradient(beta, true);
                (v / scale - 10.0, g.into_iter().map(|v| v / scale).collect())
            },
            start,
            1e-12,
            1e-10,
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::matrix::MatrixOrder;
    fn number(v: &serde_json::Value) -> f64 {
        let s = v.as_str().unwrap();
        if let Some(h) = s.strip_prefix("bits:") {
            f64::from_bits(u64::from_str_radix(h, 16).unwrap())
        } else {
            s.parse().unwrap()
        }
    }
    fn numbers(v: &serde_json::Value) -> Vec<f64> {
        v.as_array().unwrap().iter().map(number).collect()
    }
    #[test]
    fn analytic_likelihood_and_score_match_independent_original_apeglm() {
        let reference: serde_json::Value = serde_json::from_str(include_str!(
            "../../../tests/data/beta_binomial_posterior_r462.json"
        ))
        .unwrap();
        let mut failures = vec![];
        for c in reference["analytic_leaves"].as_array().unwrap() {
            let x = numbers(&c["design"]["values"]);
            let y = numbers(&c["successes"]);
            let n = numbers(&c["trials"]);
            let w = numbers(&c["weights"]);
            let prior = CauchyCoefficientPrior {
                nuisance: vec![0],
                scale: number(&c["scale"]),
                nuisance_scale: number(&c["nuisance_scale"]),
                adapted_variance: None,
            };
            let row = Row {
                design: MatrixView::new(y.len(), 2, &x, MatrixOrder::ColumnMajor).unwrap(),
                successes: &y,
                trials: &n,
                weights: &w,
                concentration: number(&c["theta"]),
                prior: &prior,
            };
            let (value, gradient) = row.objective_gradient(&numbers(&c["beta"]), false);
            for (i, (a, b)) in std::iter::once((value, number(&c["value"])))
                .chain(gradient.into_iter().zip(numbers(&c["gradient"])))
                .enumerate()
            {
                if !a.is_finite() || (a - b).abs() > 1e-9 + 1e-12 * b.abs() {
                    failures.push(format!("{} component {i}: {a:.17e} != {b:.17e}", c["name"]));
                }
            }
        }
        assert!(failures.is_empty(), "{}", failures.join("\n"));
    }
}
