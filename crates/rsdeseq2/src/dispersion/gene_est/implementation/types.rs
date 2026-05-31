use nalgebra::DMatrix;
use statrs::function::gamma::{digamma, ln_gamma};

use crate::core::CountMatrix;
use crate::design::DesignMatrix;
use crate::errors::{invalid_dimensions, DeseqError};
use crate::glm::{
    fit_fixed_dispersion_irls_with_normalization_factors_and_weights,
    fit_fixed_dispersion_irls_with_weights, IrlsOptions,
};
use crate::math::trigamma;
use crate::matrix::RowMajorMatrix;

/// Options for the initial gene-wise dispersion estimator.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GeneWiseDispersionOptions {
    /// Minimum final dispersion estimate.
    pub min_disp: f64,
    /// Optional maximum final dispersion estimate. Defaults to `max(10, n_samples)`.
    pub max_disp: Option<f64>,
    /// Lower bound on fitted raw means during dispersion fitting.
    pub min_mu: f64,
    /// Number of points in each log-alpha grid pass.
    pub grid_points: usize,
    /// Apply DESeq2's Cox-Reid log determinant adjustment.
    pub use_cox_reid: bool,
    /// Threshold used by DESeq2 to choose samples for weighted Cox-Reid terms.
    pub weight_threshold: f64,
    /// Dispersion optimizer to use after rough/moments initialization.
    pub fit_method: GeneWiseDispersionFitMethod,
    /// DESeq2 Armijo line-search initial step size.
    pub kappa_0: f64,
    /// DESeq2 dispersion log-posterior convergence tolerance.
    pub disp_tol: f64,
    /// Maximum line-search iterations.
    pub maxit: usize,
    /// Number of mean/dispersion alternations for the GLM-mu branch.
    pub niter: usize,
}

impl Default for GeneWiseDispersionOptions {
    fn default() -> Self {
        Self {
            min_disp: 1e-8,
            max_disp: None,
            min_mu: 0.5,
            grid_points: 20,
            use_cox_reid: true,
            weight_threshold: 1e-2,
            fit_method: GeneWiseDispersionFitMethod::LineSearch,
            kappa_0: 1.0,
            disp_tol: 1e-6,
            maxit: 100,
            niter: 1,
        }
    }
}

/// Normal prior on `log(alpha)` used by DESeq2's MAP dispersion objective.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DispersionPrior {
    /// Prior mean on the log-dispersion scale.
    pub log_mean: f64,
    /// Prior variance on the log-dispersion scale.
    pub variance: f64,
}

impl DispersionPrior {
    /// Create and validate a log-dispersion prior.
    pub fn new(log_mean: f64, variance: f64) -> Result<Self, DeseqError> {
        let prior = Self { log_mean, variance };
        validate_dispersion_prior(Some(prior))?;
        Ok(prior)
    }
}

/// Optimizer used for fixed-mean gene-wise dispersion estimates.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum GeneWiseDispersionFitMethod {
    /// DESeq2-shaped Armijo line search, with grid fallback for non-converged rows.
    #[default]
    LineSearch,
    /// Two-pass log-alpha grid search.
    Grid,
}

/// Output from the current linear-mu gene-wise dispersion stage.
///
/// All-zero genes are expanded back with `NaN` numeric fields and
/// `converged=false`, mirroring DESeq2's missing-row expansion pattern.
#[derive(Clone, Debug, PartialEq)]
pub struct GeneWiseDispersionOutput {
    /// Gene-wise dispersion estimates, with `NaN` for all-zero rows.
    pub disp_gene_est: Vec<f64>,
    /// Number of objective evaluations used by the grid search.
    pub disp_iter: Vec<usize>,
    /// DESeq2-style rough dispersion starts.
    pub rough_disp: Vec<f64>,
    /// DESeq2-style moments dispersion starts.
    pub moments_disp: Vec<f64>,
    /// Bounded initial dispersion values.
    pub initial_disp: Vec<f64>,
    /// Fitted raw means used for dispersion estimation.
    pub mu: RowMajorMatrix<f64>,
    /// Convergence flags. The grid search is deterministic and marks fitted
    /// non-all-zero genes as converged; line search follows DESeq2's
    /// `dispIter < maxit & dispIter != 1` convergence shape.
    pub converged: Vec<bool>,
}

/// Diagnostics from one DESeq2-shaped dispersion line search.
#[derive(Clone, Debug, PartialEq)]
pub struct DispersionLineSearchOutput {
    /// Estimated dispersion on the alpha scale.
    pub dispersion: f64,
    /// Estimated dispersion on the log-alpha scale.
    pub log_alpha: f64,
    /// Number of line-search loop iterations.
    pub iter: usize,
    /// Number of accepted line-search proposals.
    pub iter_accept: usize,
    /// Initial objective value.
    pub initial_lp: f64,
    /// Initial derivative with respect to log alpha.
    pub initial_dlp: f64,
    /// Final objective value.
    pub last_lp: f64,
    /// Final derivative with respect to log alpha.
    pub last_dlp: f64,
    /// Final second derivative with respect to log alpha.
    pub last_d2lp: f64,
    /// Final accepted objective change, or `-1` if no step was accepted.
    pub last_change: f64,
    /// DESeq2-style convergence flag.
    pub converged: bool,
}

/// Borrowed inputs for weighted prior-aware dispersion optimizer helpers.
#[derive(Clone, Copy, Debug)]
pub struct WeightedDispersionFitInput<'a> {
    /// Raw counts for one gene.
    pub counts: &'a [u32],
    /// Fitted means for one gene.
    pub mu: &'a [f64],
    /// Design matrix generated by R or caller code.
    pub design: &'a DesignMatrix,
    /// Starting dispersion on the alpha scale.
    pub initial_dispersion: f64,
    /// Dispersion optimizer options.
    pub options: GeneWiseDispersionOptions,
    /// Sample count used for default max-dispersion bounds.
    pub n_samples: usize,
    /// Normal prior on log dispersion.
    pub prior: DispersionPrior,
    /// Row-normalized observation weights for one gene.
    pub weights: &'a [f64],
}

/// Borrowed inputs for gene-wise dispersion estimation.
#[derive(Clone, Copy, Debug)]
pub struct GeneWiseDispersionInput<'a> {
    /// Raw count matrix.
    pub counts: &'a CountMatrix,
    /// Design matrix generated by R or caller code.
    pub design: &'a DesignMatrix,
    /// Positive sample size factors.
    pub size_factors: &'a [f64],
    /// Optional gene/sample normalization factors, which preempt size factors.
    pub normalization_factors: Option<&'a RowMajorMatrix<f64>>,
    /// Counts normalized by size factors or normalization factors.
    pub normalized_counts: &'a RowMajorMatrix<f64>,
    /// Per-gene base means.
    pub base_mean: &'a [f64],
    /// Per-gene base variances.
    pub base_var: &'a [f64],
    /// Per-gene all-zero flags.
    pub all_zero: &'a [bool],
    /// Optional row-normalized observation weights.
    pub observation_weights: Option<&'a RowMajorMatrix<f64>>,
}

#[derive(Clone, Copy, Debug)]
struct GeneDispersionFitDiagnostics {
    estimate: f64,
    iterations: usize,
    converged: bool,
    initial_lp: f64,
    last_lp: f64,
}
