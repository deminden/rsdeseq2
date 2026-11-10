//! Negative-binomial log-link likelihood and finite-gradient Laplace curvature.
use super::prior::CauchyCoefficientPrior;
use crate::matrix::MatrixView;
use nalgebra::DMatrix;
pub(super) struct Row<'a> {
    pub design: MatrixView<'a, f64>,
    pub counts: Vec<f64>,
    pub weights: Vec<f64>,
    pub offsets: Vec<f64>,
    pub size: f64,
    pub prior: &'a CauchyCoefficientPrior,
}
impl Row<'_> {
    pub fn objective_gradient(&self, beta: &[f64]) -> (f64, Vec<f64>) {
        let mut gradient = vec![0.0; beta.len()];
        let mut value = self.prior.contribution(beta, &mut gradient);
        let log_size = self.size.ln();
        for (i, &y) in self.counts.iter().enumerate() {
            let linear = beta
                .iter()
                .enumerate()
                .map(|(j, b)| self.design.at(i, j) * b)
                .sum::<f64>();
            let eta = linear + self.offsets[i];
            let centered = eta - log_size;
            // log(size+exp(eta)) and its logistic derivative without overflow.
            let (log_sum, probability) = if centered > 0.0 {
                (
                    eta + (-centered).exp().ln_1p(),
                    1.0 / (1.0 + (-centered).exp()),
                )
            } else {
                let e = centered.exp();
                (log_size + e.ln_1p(), e / (1.0 + e))
            };
            let weight = self.weights[i];
            value += weight * ((self.size + y) * log_sum - y * linear);
            let score = weight * ((self.size + y) * probability - y);
            for (j, g) in gradient.iter_mut().enumerate() {
                *g += self.design.at(i, j) * score;
            }
        }
        (value, gradient)
    }
    pub fn hessian(&self, beta: &[f64], step: f64) -> DMatrix<f64> {
        let p = beta.len();
        let mut h = DMatrix::zeros(p, p);
        for j in 0..p {
            let mut plus = beta.to_vec();
            let mut minus = beta.to_vec();
            plus[j] += step;
            minus[j] -= step;
            let gplus = self.objective_gradient(&plus).1;
            let gminus = self.objective_gradient(&minus).1;
            for i in 0..p {
                h[(i, j)] = (gplus[i] - gminus[i]) / (2.0 * step);
            }
        }
        for j in 0..p {
            for i in 0..j {
                let a = (h[(i, j)] + h[(j, i)]) / 2.0;
                h[(i, j)] = a;
                h[(j, i)] = a;
            }
        }
        h
    }
}
