//! Negative-binomial log-link likelihood and finite-gradient Laplace curvature.
use super::prior::CauchyCoefficientPrior;
use crate::numeric::math::ProfileMath;
use crate::{Error, matrix::MatrixView};
use rcompat_numeric::math::dense::DenseMatrix;
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
        let log_size = self.size.ln_p();
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
                    eta + (-centered).exp_p().ln_1p_p(),
                    1.0 / (1.0 + (-centered).exp_p()),
                )
            } else {
                let e = centered.exp_p();
                (log_size + e.ln_1p_p(), e / (1.0 + e))
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
    pub fn hessian(&self, beta: &[f64], step: f64) -> Result<DenseMatrix, Error> {
        super::finite_difference_hessian(beta, step, |x| self.objective_gradient(x).1)
    }
}
