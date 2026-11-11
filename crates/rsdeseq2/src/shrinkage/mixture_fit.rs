use crate::{DeseqError, errors::invalid_dimensions, matrix::MatrixView};
use nalgebra::{DMatrix, DVector};

/// Controls for the native convex mixture-weight fitter, independent of any
/// R package's optimizer-specific control list.
#[derive(Clone, Copy, Debug)]
pub struct MixtureFitOptions {
    pub max_iterations: usize,
    pub gradient_tolerance: f64,
    pub em_iterations: usize,
    pub active_iterations: usize,
    pub zero_threshold: f64,
    /// Added to row-normalized component mixtures, including Dirichlet
    /// pseudo-observations. Zero fits the unregularized likelihood; ashr's
    /// mixSQP adapter uses 1e-6 by default.
    pub likelihood_offset: f64,
    /// Original mixSQP checks dual residuals only on positive coordinates.
    /// Keep false for the general convex fitter's stronger full KKT check.
    pub positive_coordinate_convergence: bool,
    pub optimizer: MixtureOptimizer,
    pub sqp: MixtureSqpOptions,
}
impl Default for MixtureFitOptions {
    fn default() -> Self {
        Self {
            max_iterations: 1000,
            gradient_tolerance: 1e-8,
            em_iterations: 20,
            active_iterations: 100,
            zero_threshold: 1e-8,
            likelihood_offset: 0.0,
            positive_coordinate_convergence: false,
            optimizer: MixtureOptimizer::SequentialQuadratic,
            sqp: MixtureSqpOptions::default(),
        }
    }
}

/// Normalized mixture proportions and explicit numerical convergence diagnostics.
#[derive(Clone, Debug)]
pub struct MixtureFit {
    pub weights: Vec<f64>,
    pub log_likelihood: f64,
    pub penalized_log_likelihood: f64,
    pub converged: bool,
    /// Convex dual/KKT residual, or last fixed-point residual for EM.
    pub dual_residual: f64,
    pub optimizer_weight_sum: f64,
    pub function_evaluations: usize,
    pub message: String,
    pub iterations: Option<usize>,
    pub posterior_concentrations: Option<Vec<f64>>,
    pub optimizer_objective: Option<f64>,
}

fn invalid(reason: impl Into<String>) -> DeseqError {
    DeseqError::InvalidOptions {
        reason: format!("mixture prior fit: {}", reason.into()),
    }
}

struct Problem {
    likelihoods: Vec<f64>,
    row_maxima: Vec<f64>,
    observation_weights: Vec<f64>,
    prior_counts: Vec<f64>,
    mass: f64,
    components: usize,
    offset: f64,
    low_rank: bool,
}
mod controls;
mod low_rank;
pub use controls::MixtureSqpOptions;
mod em;
mod interior_point;
mod polynomial;
mod sqp;
mod squarem;
pub use interior_point::InteriorPointOptions;
pub use squarem::{PolynomialExtrapolation, SquaremOptions};

/// Distinct native weight-fitting policies; the existing SQP path is the default.
#[derive(Clone, Copy, Debug, Default)]
pub enum MixtureOptimizer {
    #[default]
    SequentialQuadratic,
    Em(SquaremOptions),
    WeightedEm(SquaremOptions),
    VariationalEm(SquaremOptions),
    CxxEm(SquaremOptions),
    InteriorPoint(InteriorPointOptions),
}
impl MixtureOptimizer {
    pub fn is_variational(self) -> bool {
        matches!(self, Self::VariationalEm(_))
    }
    pub fn is_sqp(self) -> bool {
        matches!(self, Self::SequentialQuadratic)
    }
}

impl Problem {
    fn objective_gradient(&self, weights: &[f64]) -> (f64, Vec<f64>) {
        self.objective_gradient_with_policy(weights, false)
    }
    fn objective_gradient_with_policy(&self, weights: &[f64], ashr: bool) -> (f64, Vec<f64>) {
        // Without an equality constraint, minimize
        // sum(w) - (sum_i a_i log(L_i w) + sum_k b_k log(w_k))/N,
        // where N = sum(a) + sum(b). At zero offset the optimum has sum(w) == 1;
        // a nonzero offset is included in both logarithm and derivative below.
        let mut objective = weights.iter().sum::<f64>();
        let mut gradient = vec![if ashr { 0.0 } else { 1.0 }; self.components];
        if ashr {
            for ((&count, &w), gradient) in self.prior_counts.iter().zip(weights).zip(&mut gradient)
            {
                if count != 0.0 {
                    *gradient = -(count / self.mass) / (w + self.offset);
                }
            }
        }
        for (row, &a) in self
            .likelihoods
            .chunks_exact(self.components)
            .zip(&self.observation_weights)
        {
            if a == 0.0 {
                continue;
            }
            let marginal = row.iter().zip(weights).map(|(l, w)| l * w).sum::<f64>() + self.offset;
            if marginal <= 0.0 {
                // A line-search point outside the finite likelihood domain;
                // the final accepted point is validated separately.
                return (1e100, vec![-1e100; self.components]);
            }
            objective -= (a / self.mass) * marginal.ln();
            for (gradient, &likelihood) in gradient.iter_mut().zip(row) {
                *gradient -= if ashr {
                    ((a / self.mass) / marginal) * likelihood
                } else {
                    (a / self.mass) * (likelihood / marginal)
                };
            }
        }
        for ((&count, &w), gradient) in self.prior_counts.iter().zip(weights).zip(&mut gradient) {
            if count != 0.0 {
                objective -= (count / self.mass) * (w + self.offset).ln();
                if !ashr {
                    *gradient -= (count / self.mass) / (w + self.offset);
                }
            }
        }
        if ashr {
            for gradient in &mut gradient {
                *gradient += 1.0;
            }
        }
        (objective, gradient)
    }
    fn log_likelihood(&self, weights: &[f64]) -> f64 {
        self.likelihoods
            .chunks_exact(self.components)
            .zip(&self.observation_weights)
            .zip(&self.row_maxima)
            .filter(|((_, a), _)| **a > 0.0)
            .map(|((row, a), max)| {
                let marginal = row.iter().zip(weights).map(|(l, w)| l * w).sum::<f64>();
                a * (max + marginal.ln())
            })
            .sum()
    }
}

/// Fit mixture weights with a Dirichlet prior, using convex SQP by default.
/// Convex methods require parameters >= 1; variational EM accepts positive
/// parameters and returns their fitted posterior concentrations separately.
///
/// Each matrix row is one observation and each column one component. Log
/// likelihoods may be finite or negative infinity; positive-weight observations
/// require at least one finite component. Zero observation weights do not enter
/// fitting. This is a shared numerical kernel, not an ashr method adapter.
pub fn fit_mixture_proportions(
    log_likelihoods: MatrixView<'_, f64>,
    observation_weights: Option<&[f64]>,
    dirichlet: &[f64],
    initial_weights: Option<&[f64]>,
    options: MixtureFitOptions,
) -> Result<MixtureFit, DeseqError> {
    if !options.optimizer.is_sqp()
        && !matches!(options.optimizer, MixtureOptimizer::InteriorPoint(_))
    {
        return em::fit(
            log_likelihoods,
            observation_weights,
            dirichlet,
            initial_weights,
            options.optimizer,
        );
    }
    let n = log_likelihoods.n_rows();
    if options.optimizer.is_sqp() {
        options.sqp.validate()?;
    }
    let k = log_likelihoods.n_cols();
    if k == 0 {
        return Err(invalid("at least one component is required"));
    }
    if options.max_iterations == 0
        || !options.gradient_tolerance.is_finite()
        || options.gradient_tolerance < 0.0
        || options.active_iterations == 0
        || !options.zero_threshold.is_finite()
        || options.zero_threshold < 0.0
        || options.zero_threshold >= 1.0
        || !options.likelihood_offset.is_finite()
        || options.likelihood_offset < 0.0
    {
        return Err(invalid("invalid optimization controls"));
    }
    if dirichlet.len() != k {
        return Err(invalid_dimensions(
            "mixture Dirichlet parameters",
            k,
            dirichlet.len(),
        ));
    }
    if dirichlet.iter().any(|p| !p.is_finite() || *p < 1.0) {
        return Err(invalid("Dirichlet parameters must be finite and >= 1"));
    }
    let prior_counts = dirichlet.iter().map(|p| p - 1.0).collect::<Vec<_>>();
    let a = if let Some(a) = observation_weights {
        if a.len() != n {
            return Err(invalid_dimensions(
                "mixture observation weights",
                n,
                a.len(),
            ));
        }
        if a.iter().any(|a| !a.is_finite() || *a < 0.0) {
            return Err(invalid(
                "observation weights must be finite and nonnegative",
            ));
        }
        a.to_vec()
    } else {
        vec![1.0; n]
    };
    let mass = a.iter().sum::<f64>() + prior_counts.iter().sum::<f64>();
    if !mass.is_finite() || mass <= 0.0 {
        return Err(invalid("positive finite data or prior mass is required"));
    }
    let mut likelihoods = Vec::with_capacity(
        n.checked_mul(k)
            .ok_or_else(|| invalid("matrix size overflow"))?,
    );
    let mut row_maxima = Vec::with_capacity(n);
    for (i, &weight) in a.iter().enumerate() {
        let row = (0..k)
            .map(|j| log_likelihoods.get(i, j).copied())
            .collect::<Result<Vec<_>, _>>()?;
        if row.iter().any(|l| l.is_nan() || *l == f64::INFINITY) {
            return Err(invalid(
                "log likelihoods must be finite or negative infinity",
            ));
        }
        let max = row.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        if !max.is_finite() && weight > 0.0 {
            return Err(invalid(
                "positive-weight observation has no finite component",
            ));
        }
        row_maxima.push(if weight == 0.0 { 0.0 } else { max });
        likelihoods.extend(
            row.into_iter()
                .map(|l| if weight == 0.0 { 0.0 } else { (l - max).exp() }),
        );
    }
    let mut start = if let Some(initial) = initial_weights {
        if initial.len() != k {
            return Err(invalid_dimensions(
                "initial mixture weights",
                k,
                initial.len(),
            ));
        }
        if initial.iter().any(|v| !v.is_finite() || *v < 0.0) || !initial.iter().any(|v| *v > 0.0) {
            return Err(invalid(
                "initial weights must be finite, nonnegative and not all zero",
            ));
        }
        initial.to_vec()
    } else {
        vec![1.0; k]
    };
    let total = start.iter().sum::<f64>();
    if !total.is_finite() {
        return Err(invalid("initial weight sum overflow"));
    }
    for w in &mut start {
        *w /= total;
    }
    let problem = Problem {
        likelihoods,
        row_maxima,
        observation_weights: a,
        prior_counts,
        mass,
        components: k,
        offset: if matches!(options.optimizer, MixtureOptimizer::InteriorPoint(_)) {
            0.0
        } else {
            options.likelihood_offset
        },
        low_rank: false,
    };
    let approximate = if options.optimizer.is_sqp() {
        options
            .sqp
            .svd_tolerance
            .map(|tol| problem.truncated(tol))
            .transpose()?
            .flatten()
    } else {
        None
    };
    let optimization = approximate.as_ref().unwrap_or(&problem);
    let interior = match options.optimizer {
        MixtureOptimizer::InteriorPoint(o) => Some(o),
        _ => None,
    };
    let (mut raw_weights, converged, residual, evaluations, message, interior_iterations) =
        if let Some(o) = interior {
            let f = problem.solve_interior_point(start, o)?;
            (
                f.weights,
                f.converged,
                f.residual,
                f.evaluations,
                f.message,
                Some(f.iterations),
            )
        } else {
            let (w, c, r, e, m) = optimization.solve(start, options)?;
            (w, c, r, e, m, None)
        };
    let total = raw_weights.iter().sum::<f64>();
    if !total.is_finite() || total <= 0.0 {
        return Err(invalid("optimizer returned invalid total weight"));
    }
    for w in &mut raw_weights {
        *w /= total;
    }
    let weights = raw_weights;
    let dual_residual = residual;
    let objective = optimization.objective_gradient(&weights).0;
    let log_likelihood = problem.log_likelihood(&weights);
    let penalized_log_likelihood = log_likelihood
        + problem
            .prior_counts
            .iter()
            .zip(&weights)
            .filter(|(count, _)| **count > 0.0)
            .map(|(count, w)| count * w.ln())
            .sum::<f64>();
    if !objective.is_finite()
        || !log_likelihood.is_finite()
        || !penalized_log_likelihood.is_finite()
    {
        return Err(invalid(
            "optimizer returned nonfinite likelihood or penalty",
        ));
    }
    Ok(MixtureFit {
        weights,
        log_likelihood,
        penalized_log_likelihood,
        converged,
        dual_residual,
        optimizer_weight_sum: total,
        function_evaluations: evaluations,
        message,
        iterations: interior_iterations,
        posterior_concentrations: None,
        optimizer_objective: interior.map(|_| objective),
    })
}
