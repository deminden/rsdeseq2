//! Mixture EM and variational-Dirichlet updates with native acceleration.
//! Equations are likelihood responsibilities and the Dirichlet ELBO, not an
//! execution bridge or translation of foreign package source.
use super::{MixtureFit, MixtureOptimizer, invalid, squarem};
use crate::probability::digamma_positive;
use crate::{DeseqError, matrix::MatrixView};
use libm::lgamma as ln_gamma;
// Extrapolated coordinates may cross zero before the stabilizing update.
// The source fixed-point map continues digamma across noninteger negatives;
// rejecting those candidates changes sparse variational trajectories.
fn digamma(x: f64) -> f64 {
    if x > 0.0 {
        return digamma_positive(x);
    }
    if !x.is_finite() || x.fract() == 0.0 {
        return f64::NAN;
    }
    let angle = std::f64::consts::PI * x.rem_euclid(1.0);
    digamma_positive(1.0 - x) - std::f64::consts::PI * angle.cos() / angle.sin()
}
struct Problem {
    likelihood: Vec<f64>,
    row_maxima: Vec<f64>,
    prior: Vec<f64>,
    weights: Vec<f64>,
    k: usize,
    variational: bool,
    compiled: bool,
}
fn normalize(mut v: Vec<f64>) -> Vec<f64> {
    let total = v.iter().sum::<f64>();
    for x in &mut v {
        *x /= total;
    }
    v
}
impl Problem {
    fn effective(&self, state: &[f64]) -> Vec<f64> {
        if self.variational {
            let total = digamma(state.iter().sum());
            state.iter().map(|p| (digamma(*p) - total).exp()).collect()
        } else if self.compiled {
            state.to_vec()
        } else {
            normalize(state.iter().map(|v| v.max(0.0)).collect())
        }
    }
    fn update(&self, state: &[f64]) -> Result<Vec<f64>, DeseqError> {
        if self.variational && state.iter().any(|v| !v.is_finite()) {
            return Err(invalid(
                "variational concentrations must be positive finite",
            ));
        }
        let pi = self.effective(state);
        let mut counts = vec![0.0; self.k];
        for (row, &weight) in self.likelihood.chunks_exact(self.k).zip(&self.weights) {
            let marginal = row.iter().zip(&pi).map(|(l, p)| l * p).sum::<f64>();
            if !marginal.is_finite() || marginal <= 0.0 {
                return Err(invalid(
                    "EM responsibility has no finite positive likelihood",
                ));
            }
            for j in 0..self.k {
                counts[j] += weight * ((row[j] * pi[j]) / marginal);
            }
        }
        for (count, prior) in counts.iter_mut().zip(&self.prior) {
            *count += if self.variational {
                *prior
            } else {
                prior - 1.0
            };
        }
        let next = if self.variational {
            counts
        } else {
            normalize(counts)
        };
        if next.iter().any(|v| !v.is_finite()) {
            return Err(invalid("EM update returned nonfinite parameters"));
        }
        Ok(next)
    }
    fn log_likelihood(&self, pi: &[f64]) -> f64 {
        self.likelihood
            .chunks_exact(self.k)
            .zip(&self.weights)
            .zip(&self.row_maxima)
            .map(|((row, w), max)| {
                if *w == 0.0 {
                    0.0
                } else {
                    w * (max + row.iter().zip(pi).map(|(l, p)| l * p).sum::<f64>().ln())
                }
            })
            .sum()
    }
    fn penalty(&self, pi: &[f64]) -> f64 {
        self.prior
            .iter()
            .zip(pi)
            .filter(|(a, _)| **a != 1.0)
            .map(|(a, p)| (a - 1.0) * p.ln())
            .sum()
    }
    fn objective(&self, state: &[f64]) -> f64 {
        let pi = self.effective(state);
        if !self.variational {
            return -(self
                .likelihood
                .chunks_exact(self.k)
                .zip(&self.weights)
                .map(|(row, w)| {
                    if *w == 0.0 {
                        0.0
                    } else {
                        w * row.iter().zip(&pi).map(|(l, p)| l * p).sum::<f64>().ln()
                    }
                })
                .sum::<f64>()
                + self.penalty(&pi));
        }
        let mut energy = 0.0;
        let mut entropy = 0.0;
        let marginals = self
            .likelihood
            .chunks_exact(self.k)
            .map(|row| row.iter().zip(&pi).map(|(l, p)| l * p).sum::<f64>())
            .collect::<Vec<_>>();
        for (j, &p) in pi.iter().enumerate() {
            for (i, &marginal) in marginals.iter().enumerate() {
                let l = self.likelihood[i * self.k + j];
                let q = l * p / marginal;
                let term = q * (l * p).ln();
                if !term.is_nan() {
                    energy += term;
                }
                entropy += q * q.ln();
            }
        }
        let sum = state.iter().sum::<f64>();
        let prior_sum = self.prior.iter().sum::<f64>();
        let mut kl = ln_gamma(sum) - ln_gamma(prior_sum);
        kl += state
            .iter()
            .zip(&self.prior)
            .map(|(q, p)| (q - p) * (digamma(*q) - digamma(sum)))
            .sum::<f64>();
        kl += self
            .prior
            .iter()
            .zip(state)
            .map(|(p, q)| ln_gamma(*p) - ln_gamma(*q))
            .sum::<f64>();
        -(energy - kl - entropy)
    }
    fn null_optimal(&self) -> bool {
        if self.variational || self.prior.iter().skip(1).any(|a| *a > 1.0) {
            return false;
        }
        (0..self.k).all(|j| {
            let gradient = self.row_maxima.len() as f64
                - self
                    .likelihood
                    .chunks_exact(self.k)
                    .map(|row| row[j] / row[0])
                    .sum::<f64>()
                + self.prior[0]
                - 1.0;
            gradient.is_nan() || gradient >= 0.0
        })
    }
}
pub(super) fn fit(
    log_likelihoods: MatrixView<'_, f64>,
    weights: Option<&[f64]>,
    prior: &[f64],
    initial: Option<&[f64]>,
    optimizer: MixtureOptimizer,
) -> Result<MixtureFit, DeseqError> {
    let (controls, weighted, variational, compiled) = match optimizer {
        MixtureOptimizer::Em(o) => (o, false, false, false),
        MixtureOptimizer::WeightedEm(o) => (o, true, false, false),
        MixtureOptimizer::VariationalEm(o) => (o, false, true, false),
        MixtureOptimizer::CxxEm(o) => (o, false, false, true),
        _ => return Err(invalid("EM fitter requires an EM optimizer")),
    };
    controls.validate()?;
    let n = log_likelihoods.n_rows();
    let k = log_likelihoods.n_cols();
    if weights.is_some() && !weighted {
        return Err(invalid(
            "observation weights require the weighted EM optimizer",
        ));
    }
    if prior.len() != k
        || prior
            .iter()
            .any(|v| !v.is_finite() || *v <= 0.0 || (!variational && *v < 1.0))
    {
        return Err(invalid("invalid EM Dirichlet parameters"));
    }
    let weights = weights.map(<[f64]>::to_vec).unwrap_or_else(|| vec![1.0; n]);
    if weights.len() != n
        || weights.iter().any(|v| !v.is_finite() || *v < 0.0)
        || !weights.iter().sum::<f64>().is_finite()
        || weights.iter().sum::<f64>() + prior.iter().map(|p| p - 1.0).sum::<f64>() <= 0.0
    {
        return Err(invalid(
            "EM observation weights must be finite nonnegative and row aligned with positive total",
        ));
    }
    let mut likelihood = Vec::with_capacity(n * k);
    let mut row_maxima = Vec::with_capacity(n);
    for i in 0..n {
        let row = (0..k)
            .map(|j| *log_likelihoods.at(i, j))
            .collect::<Vec<_>>();
        if row.iter().any(|v| v.is_nan() || *v == f64::INFINITY) {
            return Err(invalid(
                "EM log likelihoods must be finite or negative infinity",
            ));
        }
        let max = row.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        if !max.is_finite() {
            return Err(invalid("EM row has no finite likelihood"));
        }
        row_maxima.push(max);
        likelihood.extend(row.iter().map(|v| (v - max).exp()));
    }
    let mut initial = initial
        .map(<[f64]>::to_vec)
        .unwrap_or_else(|| vec![1.0 / k as f64; k]);
    if initial.len() != k
        || initial.iter().any(|v| !v.is_finite() || *v < 0.0)
        || !initial.iter().any(|v| *v > 0.0)
    {
        return Err(invalid("invalid initial EM weights"));
    }
    if initial.contains(&0.0) {
        for v in &mut initial {
            *v = v.max(1e-6);
        }
        initial = normalize(initial);
    }
    let keep = (0..k)
        .filter(|&j| prior[j] > 1.0 || (0..n).any(|i| likelihood[i * k + j] > 0.0))
        .collect::<Vec<_>>();
    if keep.is_empty() {
        return Err(invalid("no nonzero likelihood components"));
    }
    let matrix = (0..n)
        .flat_map(|i| {
            keep.iter().map({
                let likelihood = &likelihood;
                move |&j| likelihood[i * k + j]
            })
        })
        .collect();
    let problem = Problem {
        likelihood: matrix,
        row_maxima,
        weights,
        prior: keep.iter().map(|&j| prior[j]).collect(),
        k: keep.len(),
        variational,
        compiled,
    };
    let start = if variational {
        vec![1.0; keep.len()]
    } else {
        keep.iter().map(|&j| initial[j]).collect()
    };
    let shortcut = keep.len() == 1 || (!weighted && problem.null_optimal());
    let iteration = if shortcut {
        squarem::Iteration {
            state: (0..keep.len())
                .map(|j| if j == 0 { 1.0 } else { 0.0 })
                .collect(),
            converged: true,
            evaluations: 0,
            iterations: 0,
            objective: 0.0,
            residual: 0.0,
        }
    } else {
        squarem::run(
            start,
            controls,
            controls.tolerance.unwrap_or((0.1 / n as f64).min(1e-7)),
            compiled,
            |p| problem.update(p),
            |p| problem.objective(p),
        )?
    };
    let sum = iteration.state.iter().sum::<f64>();
    let active = normalize(iteration.state.iter().map(|v| v.max(0.0)).collect());
    let mut weights = vec![0.0; k];
    for (&j, &v) in keep.iter().zip(&active) {
        weights[j] = v;
    }
    let mut concentrations = vec![0.0; k];
    for (&j, &v) in keep.iter().zip(&iteration.state) {
        concentrations[j] = v;
    }
    let ll = problem.log_likelihood(&active);
    let pll = ll + problem.penalty(&active);
    if !ll.is_finite() || !pll.is_finite() {
        return Err(invalid("EM fit returned nonfinite likelihood"));
    }
    Ok(MixtureFit {
        weights,
        log_likelihood: ll,
        penalized_log_likelihood: pll,
        converged: iteration.converged,
        dual_residual: iteration.residual,
        optimizer_weight_sum: sum,
        function_evaluations: iteration.evaluations,
        message: if shortcut {
            "gradient check selected first component".into()
        } else if iteration.converged {
            "squared EM converged".into()
        } else {
            "squared EM evaluation budget exhausted".into()
        },
        iterations: (!shortcut).then_some(iteration.iterations),
        posterior_concentrations: (variational && !shortcut).then_some(concentrations),
        optimizer_objective: (!shortcut).then_some(iteration.objective),
    })
}
