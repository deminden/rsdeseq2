use nalgebra::{DMatrix, DVector, Dyn, linalg::LU};
use rayon::prelude::*;
use rcompat_lbfgsb::{Bounds, OptimControl, optim_lbfgsb, optim_lbfgsb_with_gradient};
use std::cell::RefCell;

use crate::core::CountMatrix;
use crate::design::DesignMatrix;
use crate::errors::{DeseqError, invalid_dimensions};
use crate::glm::NbinomGlmFit;
use crate::glm::beta::fit_intercept_only_fixed_dispersion;
use crate::glm::fallback::{OptimFallbackRows, optim_fallback_rows, record_optim_fallback_row};
use crate::glm::nb::{
    FixedDispersionLikelihood, nbinom_log_likelihood_matrix, nbinom_log_likelihood_weighted,
};
use crate::glm::r_arithmetic::{nbinom_log, normal_log};
use crate::math::optim::BoundedOptimizationOutput;
use crate::matrix::RowMajorMatrix;

/// Internal retention policy; fitting and fallback arithmetic are identical.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum OutputRequirements {
    Inference,
    Means,
    Transform,
}
impl OutputRequirements {
    fn covariance(self) -> bool {
        matches!(self, Self::Inference)
    }
}

/// Fit a fixed-dispersion negative-binomial GLM with DESeq2-style dispatch.
///
/// Intercept-only designs with the default tiny ridge use DESeq2's closed-form
/// shortcut. Other designs use the general fixed-dispersion IRLS path.
pub fn fit_irls(
    counts: &CountMatrix,
    design: &DesignMatrix,
    size_factors: &[f64],
    dispersions: &[f64],
    options: IrlsOptions,
) -> Result<NbinomGlmFit, DeseqError> {
    if is_intercept_only_design(design)
        && options.uses_intercept_shortcut_for_coefficients(design.n_coefficients())?
    {
        fit_intercept_only_fixed_dispersion(counts, size_factors, dispersions)
    } else {
        fit_fixed_dispersion_irls(counts, design, size_factors, dispersions, options)
    }
}

/// Linear solver for the IRLS weighted least-squares update.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum IrlsSolver {
    /// Solve `(X' W X + ridge) beta = X' W z` directly.
    ///
    /// This preserves the initial Rust behavior and is useful for existing
    /// `useQR=FALSE` DESeq2 references.
    NormalEquations,
    /// Solve DESeq2's augmented QR problem `[sqrt(W) X; sqrt(ridge)] beta`.
    #[default]
    Qr,
}

/// Options for the initial fixed-dispersion IRLS implementation.
#[derive(Clone, Debug, PartialEq)]
pub struct IrlsOptions {
    /// Convergence tolerance matching DESeq2's `betaTol` criterion.
    pub beta_tol: f64,
    /// Maximum IRLS iterations.
    pub maxit: usize,
    /// Lower bound on fitted means during fitting.
    pub min_mu: f64,
    /// Mark rows as not converged if any beta exceeds this absolute value.
    pub max_beta_abs: f64,
    /// Natural-log-scale ridge value added to each coefficient.
    pub ridge_lambda: f64,
    /// Optional natural-log-scale ridge values, one per coefficient.
    ///
    /// When supplied, this vector overrides `ridge_lambda` and is used as the
    /// diagonal of DESeq2's `diag(lambda)` ridge matrix.
    pub ridge_lambda_by_coefficient: Option<Vec<f64>>,
    /// Weighted least-squares solver used inside each IRLS step.
    pub solver: IrlsSolver,
    /// Also refit rows that fail to converge within `maxit` using bounded optimization.
    pub use_optim: bool,
    /// Send every row through bounded optimization after IRLS.
    pub force_optim: bool,
    /// Maximum bounded-optimizer iterations for fallback rows.
    pub optim_maxit: usize,
    /// Projected-gradient tolerance for fallback optimization.
    pub optim_tol: f64,
    /// Match R `optim()` callback arithmetic and numerical gradients.
    pub r_optim_compat: bool,
}

impl Default for IrlsOptions {
    fn default() -> Self {
        let inv_ln2 = std::f64::consts::LOG2_E;
        Self {
            beta_tol: 1e-8,
            maxit: 100,
            min_mu: 0.5,
            max_beta_abs: 30.0,
            ridge_lambda: 1e-6 * inv_ln2 * inv_ln2,
            ridge_lambda_by_coefficient: None,
            solver: IrlsSolver::Qr,
            use_optim: true,
            force_optim: false,
            optim_maxit: 100,
            optim_tol: 0.0,
            r_optim_compat: true,
        }
    }
}

impl IrlsOptions {
    /// Set natural-log-scale ridge values, one per coefficient.
    pub fn ridge_lambda_by_coefficient(mut self, ridge_lambda: Vec<f64>) -> Self {
        self.ridge_lambda_by_coefficient = Some(ridge_lambda);
        self
    }

    fn ridge_lambdas_for_coefficients(&self, p: usize) -> Result<Vec<f64>, DeseqError> {
        let values = match &self.ridge_lambda_by_coefficient {
            Some(values) => {
                if values.len() != p {
                    return Err(invalid_dimensions(
                        "IRLS ridge lambda coefficients",
                        p,
                        values.len(),
                    ));
                }
                values.clone()
            }
            None => vec![self.ridge_lambda; p],
        };
        for (idx, value) in values.iter().copied().enumerate() {
            validate_nonnegative_finite(value, "ridge lambda", idx)?;
        }
        Ok(values)
    }

    pub(crate) fn uses_intercept_shortcut_for_coefficients(
        &self,
        p: usize,
    ) -> Result<bool, DeseqError> {
        let inv_ln2 = std::f64::consts::LOG2_E;
        let default_nat_log_lambda = 1e-6 * inv_ln2 * inv_ln2;
        Ok(self
            .ridge_lambdas_for_coefficients(p)?
            .into_iter()
            .all(|value| value <= default_nat_log_lambda))
    }
}

fn is_intercept_only_design(design: &DesignMatrix) -> bool {
    design.n_coefficients() == 1
        && design
            .matrix()
            .as_slice()
            .iter()
            .all(|value| (*value - 1.0).abs() <= f64::EPSILON)
}

/// Fit an unweighted fixed-dispersion NB GLM by IRLS.
///
/// This implements the standard-design-matrix branch of DESeq2's `fitBeta`
/// loop. The weighted least-squares update can use either the direct normal
/// equations branch or DESeq2's augmented QR branch. Observation weights are
/// supported, and fallback rows can be refit with bounded limited-memory
/// BFGS-style pure-Rust optimization. Explicit contrast testing is layered on top of the stored
/// beta covariance matrices.
pub fn fit_fixed_dispersion_irls(
    counts: &CountMatrix,
    design: &DesignMatrix,
    size_factors: &[f64],
    dispersions: &[f64],
    options: IrlsOptions,
) -> Result<NbinomGlmFit, DeseqError> {
    fit_fixed_dispersion_irls_with_weights(counts, design, size_factors, dispersions, None, options)
}

/// Fit an unweighted or observation-weighted fixed-dispersion NB GLM by IRLS.
///
/// Observation weights follow DESeq2's low-level `fitBeta` semantics: the
/// caller supplies non-negative gene/sample weights, and each row's working
/// weights and log likelihood are multiplied by those values.
pub fn fit_fixed_dispersion_irls_with_weights(
    counts: &CountMatrix,
    design: &DesignMatrix,
    size_factors: &[f64],
    dispersions: &[f64],
    weights: Option<&RowMajorMatrix<f64>>,
    options: IrlsOptions,
) -> Result<NbinomGlmFit, DeseqError> {
    let normalization_factors = normalization_factors_from_size_factors(counts, size_factors)?;
    fit_fixed_dispersion_irls_with_normalization_factors_and_weights(
        counts,
        design,
        &normalization_factors,
        dispersions,
        weights,
        options,
    )
}

/// Fit an unweighted fixed-dispersion NB GLM by IRLS with normalization factors.
pub fn fit_fixed_dispersion_irls_with_normalization_factors(
    counts: &CountMatrix,
    design: &DesignMatrix,
    normalization_factors: &RowMajorMatrix<f64>,
    dispersions: &[f64],
    options: IrlsOptions,
) -> Result<NbinomGlmFit, DeseqError> {
    fit_fixed_dispersion_irls_with_normalization_factors_and_weights(
        counts,
        design,
        normalization_factors,
        dispersions,
        None,
        options,
    )
}

/// Fit an unweighted or observation-weighted fixed-dispersion NB GLM by IRLS
/// with gene/sample normalization factors.
pub fn fit_fixed_dispersion_irls_with_normalization_factors_and_weights(
    counts: &CountMatrix,
    design: &DesignMatrix,
    normalization_factors: &RowMajorMatrix<f64>,
    dispersions: &[f64],
    weights: Option<&RowMajorMatrix<f64>>,
    options: IrlsOptions,
) -> Result<NbinomGlmFit, DeseqError> {
    fit_nf_irls_for_output(
        counts,
        design,
        normalization_factors,
        dispersions,
        weights,
        options,
        OutputRequirements::Inference,
    )
}

pub(crate) fn fit_sf_irls_for_output(
    counts: &CountMatrix,
    design: &DesignMatrix,
    size_factors: &[f64],
    dispersions: &[f64],
    weights: Option<&RowMajorMatrix<f64>>,
    options: IrlsOptions,
    requirements: OutputRequirements,
) -> Result<NbinomGlmFit, DeseqError> {
    let normalization_factors = normalization_factors_from_size_factors(counts, size_factors)?;
    fit_nf_irls_for_output(
        counts,
        design,
        &normalization_factors,
        dispersions,
        weights,
        options,
        requirements,
    )
}

pub(crate) fn fit_nf_irls_for_output(
    counts: &CountMatrix,
    design: &DesignMatrix,
    normalization_factors: &RowMajorMatrix<f64>,
    dispersions: &[f64],
    weights: Option<&RowMajorMatrix<f64>>,
    options: IrlsOptions,
    requirements: OutputRequirements,
) -> Result<NbinomGlmFit, DeseqError> {
    validate_nf_irls_inputs(
        counts,
        design,
        normalization_factors,
        dispersions,
        weights,
        &options,
    )?;

    let x = DMatrix::from_row_slice(
        design.n_samples(),
        design.n_coefficients(),
        design.matrix().as_slice(),
    );
    let initial_lu = (x.transpose() * &x).lu();
    let p = design.n_coefficients();
    let ridge_lambda = options.ridge_lambdas_for_coefficients(p)?;
    let n_beta = checked_output_elements(counts.n_genes(), p)?;
    let covariance_columns = checked_output_elements(p, p)?;
    let n_covariance = if requirements.covariance() {
        checked_output_elements(counts.n_genes(), covariance_columns)?
    } else {
        0
    };
    let n_means = checked_output_elements(counts.n_genes(), counts.n_samples())?;
    let mut beta_values = Vec::with_capacity(n_beta);
    let mut beta_var_values = Vec::with_capacity(n_beta);
    let mut beta_optim_start_values = vec![f64::NAN; n_beta];
    let mut beta_covariance_values = Vec::with_capacity(n_covariance);
    let mut mu_values = Vec::with_capacity(n_means);
    let mut hat_values = Vec::with_capacity(n_means);
    let mut beta_iter = Vec::with_capacity(counts.n_genes());
    let mut beta_converged = Vec::with_capacity(counts.n_genes());
    let mut beta_optim_iter = vec![f64::NAN; counts.n_genes()];
    let mut beta_optim_start_objective = vec![f64::NAN; counts.n_genes()];
    let mut beta_optim_objective = vec![f64::NAN; counts.n_genes()];
    let mut beta_optim_gradient_norm = vec![f64::NAN; counts.n_genes()];

    // Bound temporary output storage, retaining each gene's arithmetic and order.
    let row_bytes = covariance_columns
        .saturating_add(2 * p)
        .saturating_add(2 * counts.n_samples())
        .saturating_mul(std::mem::size_of::<f64>())
        .saturating_add(std::mem::size_of::<Result<IrlsRowOutput, DeseqError>>());
    let batch_rows = (8 * 1024 * 1024 / row_bytes.max(1)).clamp(1, 2048);
    let mut transient_fallback = OptimFallbackRows::default();
    for first in (0..counts.n_genes()).step_by(batch_rows) {
        let last = first.saturating_add(batch_rows).min(counts.n_genes());
        let row_outputs = (first..last)
            .into_par_iter()
            .map(|gene| {
                fit_fixed_dispersion_irls_row(IrlsRowInput {
                    gene,
                    counts,
                    x: &x,
                    initial_lu: &initial_lu,
                    normalization_factors,
                    dispersions,
                    weights,
                    ridge_lambda: &ridge_lambda,
                    options: &options,
                })
            })
            .collect::<Vec<_>>();
        // Indexed collection plus ordered extraction reports the lowest failed gene.
        for (offset, row) in row_outputs.into_iter().enumerate() {
            let row = row?;
            if !requirements.covariance() {
                record_optim_fallback_row(
                    &mut transient_fallback,
                    first + offset,
                    row.converged,
                    &row.beta,
                    &row.beta_covariance,
                    options.use_optim,
                    options.force_optim,
                );
            }
            beta_values.extend(row.beta);
            beta_var_values.extend(row.beta_se);
            if requirements.covariance() {
                beta_covariance_values.extend(row.beta_covariance);
            }
            mu_values.extend(row.mu);
            hat_values.extend(row.hat);
            beta_iter.push(row.iter);
            beta_converged.push(row.converged);
        }
    }

    let (mut beta_values, mut beta_covariance_values, fallback_rows) = if requirements.covariance()
    {
        let beta_for_routing = RowMajorMatrix::from_row_major(counts.n_genes(), p, beta_values)?;
        let covariance_for_routing =
            RowMajorMatrix::from_row_major(counts.n_genes(), p * p, beta_covariance_values)?;
        let rows = optim_fallback_rows(
            &beta_converged,
            &beta_for_routing,
            &covariance_for_routing,
            options.use_optim,
            options.force_optim,
        )?;
        (
            beta_for_routing.into_values(),
            covariance_for_routing.into_values(),
            rows,
        )
    } else {
        (beta_values, beta_covariance_values, transient_fallback)
    };
    let mut optim_log_like = vec![None; counts.n_genes()];
    if !fallback_rows.rows.is_empty() {
        refit_optim_fallback_rows(
            &fallback_rows.rows,
            &mut beta_values,
            &mut beta_var_values,
            &mut beta_optim_start_values,
            &mut beta_covariance_values,
            &mut mu_values,
            &mut hat_values,
            &mut beta_converged,
            &mut optim_log_like,
            &mut beta_optim_iter,
            &mut beta_optim_start_objective,
            &mut beta_optim_objective,
            &mut beta_optim_gradient_norm,
            OptimFallbackInput {
                counts,
                x: &x,
                normalization_factors,
                dispersions,
                weights,
                ridge_lambda: &ridge_lambda,
                options: &options,
            },
        )?;
    }
    let beta = RowMajorMatrix::from_row_major(counts.n_genes(), p, beta_values)?;
    let beta_se = RowMajorMatrix::from_row_major(counts.n_genes(), p, beta_var_values)?;
    let beta_optim_start =
        RowMajorMatrix::from_row_major(counts.n_genes(), p, beta_optim_start_values)?;
    let beta_covariance = if requirements.covariance() {
        Some(RowMajorMatrix::from_row_major(
            counts.n_genes(),
            p * p,
            beta_covariance_values,
        )?)
    } else {
        None
    };
    let mu = RowMajorMatrix::from_row_major(counts.n_genes(), counts.n_samples(), mu_values)?;
    let hat_diagonal =
        RowMajorMatrix::from_row_major(counts.n_genes(), counts.n_samples(), hat_values)?;
    let mut log_like = nbinom_log_likelihood_matrix(counts, &mu, dispersions, weights)?;
    for (gene, row_log_like) in optim_log_like.into_iter().enumerate() {
        if let Some(row_log_like) = row_log_like {
            log_like[gene] = row_log_like;
        }
    }

    Ok(NbinomGlmFit {
        log_like,
        beta_converged,
        beta,
        beta_se,
        beta_optim_start,
        beta_covariance,
        mu,
        beta_iter,
        beta_optim_iter,
        beta_optim_start_objective,
        beta_optim_objective,
        beta_optim_gradient_norm,
        model_matrix: design.clone(),
        n_terms: p,
        hat_diagonal,
    })
}

fn checked_output_elements(rows: usize, columns: usize) -> Result<usize, DeseqError> {
    rows.checked_mul(columns)
        .filter(|&length| length <= isize::MAX as usize / std::mem::size_of::<f64>())
        .ok_or_else(|| DeseqError::InvalidDimensions {
            context: "IRLS output allocation overflow".into(),
            expected: isize::MAX as usize / std::mem::size_of::<f64>(),
            actual: rows.saturating_mul(columns),
        })
}

struct OptimFallbackInput<'a> {
    counts: &'a CountMatrix,
    x: &'a DMatrix<f64>,
    normalization_factors: &'a RowMajorMatrix<f64>,
    dispersions: &'a [f64],
    weights: Option<&'a RowMajorMatrix<f64>>,
    ridge_lambda: &'a [f64],
    options: &'a IrlsOptions,
}

struct IrlsRowInput<'a> {
    gene: usize,
    counts: &'a CountMatrix,
    x: &'a DMatrix<f64>,
    initial_lu: &'a LU<f64, Dyn, Dyn>,
    normalization_factors: &'a RowMajorMatrix<f64>,
    dispersions: &'a [f64],
    weights: Option<&'a RowMajorMatrix<f64>>,
    ridge_lambda: &'a [f64],
    options: &'a IrlsOptions,
}

struct IrlsRowOutput {
    beta: Vec<f64>,
    beta_se: Vec<f64>,
    beta_covariance: Vec<f64>,
    mu: Vec<f64>,
    hat: Vec<f64>,
    iter: usize,
    converged: bool,
}

// Per-gene IRLS iterations and scratch reuse.
include!("irls/iteration.rs");

// Bounded-optimizer fallback, starting values and acceptance checks.
include!("irls/fallback.rs");

// Optimizer objectives, gradients and Hessians.
include!("irls/objective.rs");

// Projected-gradient convergence and checked arithmetic.
include!("irls/arithmetic.rs");

// Starting coefficients and stable fitted means.
include!("irls/means.rs");

// Working weights and responses.
include!("irls/working.rs");

// Weighted least squares, QR, covariance and design crossproducts.
include!("irls/linear_algebra.rs");

// Input validation, normalization factors and output scaling.
include!("irls/validation.rs");

// Private numerical regression tests.
include!("irls/tests.rs");
