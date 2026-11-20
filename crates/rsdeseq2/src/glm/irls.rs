use crate::numeric::checked::{self, CheckedExt};
use crate::numeric::math::ProfileMath;
use crate::numeric::validate::{
    nonnegative_finite as validate_nonnegative_finite, positive_finite as validate_positive_finite,
};
use nalgebra::{DMatrix, DVector, Dyn, linalg::LU};
use rayon::prelude::*;
use rcompat_lbfgsb::{Bounds, OptimControl, optim_lbfgsb, optim_lbfgsb_with_gradient};
use std::cell::RefCell;

use crate::design::DesignMatrix;
use crate::errors::{DeseqError, invalid_dimensions};
use crate::glm::NbinomGlmFit;
use crate::glm::beta::{GlmInput, fit_intercept_only_glm};
use crate::glm::fallback::requires_optim_fallback;
use crate::glm::nb::{
    FixedDispersionLikelihood, nbinom_log_likelihood_matrix, nbinom_log_likelihood_weighted,
};
use crate::glm::r_arithmetic::{nbinom_log, normal_log};
use crate::math::optim::BoundedOptimizationOutput;
use crate::matrix::CountMatrix;
use crate::matrix::RowMajorMatrix;
use crate::matrix::try_alloc_matrix;
use crate::normalization::Offsets;

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
        fit_intercept_only_glm(GlmInput::new(
            counts,
            Offsets::PerSample(size_factors),
            dispersions,
        ))
    } else {
        fit_fixed_dispersion_glm(
            GlmInput::new(counts, Offsets::PerSample(size_factors), dispersions),
            design,
            options,
        )
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
    /// The settings DESeq2 1.52 `rlogData` passes to `fitNbinomGLMs`:
    /// `betaTol = 1e-4`, `useOptim = FALSE`, `useQR = TRUE`; everything else
    /// as [`IrlsOptions::default`].
    pub fn deseq2_rlog() -> Self {
        Self {
            beta_tol: 1e-4,
            use_optim: false,
            solver: IrlsSolver::Qr,
            ..Self::default()
        }
    }

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
#[deprecated(
    since = "0.6.0",
    note = "use fit_fixed_dispersion_glm(GlmInput::new(counts, Offsets::PerSample(size_factors), dispersions), design, options)"
)]
pub fn fit_fixed_dispersion_irls(
    counts: &CountMatrix,
    design: &DesignMatrix,
    size_factors: &[f64],
    dispersions: &[f64],
    options: IrlsOptions,
) -> Result<NbinomGlmFit, DeseqError> {
    fit_fixed_dispersion_glm(
        GlmInput::new(counts, Offsets::PerSample(size_factors), dispersions),
        design,
        options,
    )
}

/// Fit an unweighted or observation-weighted fixed-dispersion NB GLM by IRLS.
///
/// Observation weights follow DESeq2's low-level `fitBeta` semantics: the
/// caller supplies non-negative gene/sample weights, and each row's working
/// weights and log likelihood are multiplied by those values.
#[deprecated(
    since = "0.6.0",
    note = "use fit_fixed_dispersion_glm with Offsets::PerSample and GlmInput::weights"
)]
pub fn fit_fixed_dispersion_irls_with_weights(
    counts: &CountMatrix,
    design: &DesignMatrix,
    size_factors: &[f64],
    dispersions: &[f64],
    weights: Option<&RowMajorMatrix<f64>>,
    options: IrlsOptions,
) -> Result<NbinomGlmFit, DeseqError> {
    fit_fixed_dispersion_glm(
        GlmInput {
            counts,
            offsets: Offsets::PerSample(size_factors),
            dispersions,
            weights,
        },
        design,
        options,
    )
}

/// Fit an unweighted fixed-dispersion NB GLM by IRLS with normalization factors.
#[deprecated(
    since = "0.6.0",
    note = "use fit_fixed_dispersion_glm(GlmInput::new(counts, Offsets::PerCell(normalization_factors), dispersions), design, options)"
)]
pub fn fit_fixed_dispersion_irls_with_normalization_factors(
    counts: &CountMatrix,
    design: &DesignMatrix,
    normalization_factors: &RowMajorMatrix<f64>,
    dispersions: &[f64],
    options: IrlsOptions,
) -> Result<NbinomGlmFit, DeseqError> {
    fit_fixed_dispersion_glm(
        GlmInput::new(counts, Offsets::PerCell(normalization_factors), dispersions),
        design,
        options,
    )
}

/// Fit an unweighted or observation-weighted fixed-dispersion NB GLM by IRLS
/// with gene/sample normalization factors.
#[deprecated(
    since = "0.6.0",
    note = "use fit_fixed_dispersion_glm with Offsets::PerCell and GlmInput::weights"
)]
pub fn fit_fixed_dispersion_irls_with_normalization_factors_and_weights(
    counts: &CountMatrix,
    design: &DesignMatrix,
    normalization_factors: &RowMajorMatrix<f64>,
    dispersions: &[f64],
    weights: Option<&RowMajorMatrix<f64>>,
    options: IrlsOptions,
) -> Result<NbinomGlmFit, DeseqError> {
    fit_fixed_dispersion_glm(
        GlmInput {
            counts,
            offsets: Offsets::PerCell(normalization_factors),
            dispersions,
            weights,
        },
        design,
        options,
    )
}

/// Fit a fixed-dispersion NB GLM by IRLS (DESeq2's `fitBeta` loop) with
/// size-factor or normalization-factor offsets and optional observation
/// weights.
pub fn fit_fixed_dispersion_glm(
    input: GlmInput<'_>,
    design: &DesignMatrix,
    options: IrlsOptions,
) -> Result<NbinomGlmFit, DeseqError> {
    match input.offsets {
        Offsets::PerSample(size_factors) => fit_sf_irls_for_output(
            input.counts,
            design,
            size_factors,
            input.dispersions,
            input.weights,
            options,
            OutputRequirements::Inference,
        ),
        Offsets::PerCell(normalization_factors) => fit_nf_irls_for_output(
            input.counts,
            design,
            normalization_factors,
            input.dispersions,
            input.weights,
            options,
            OutputRequirements::Inference,
        ),
    }
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
    validate_size_factor_offsets(counts, size_factors)?;
    fit_offsets_irls_for_output(
        counts,
        design,
        Offsets::PerSample(size_factors),
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
    fit_offsets_irls_for_output(
        counts,
        design,
        Offsets::PerCell(normalization_factors),
        dispersions,
        weights,
        options,
        requirements,
    )
}

pub(crate) fn fit_offsets_irls_for_output(
    counts: &CountMatrix,
    design: &DesignMatrix,
    normalization_factors: Offsets<'_>,
    dispersions: &[f64],
    weights: Option<&RowMajorMatrix<f64>>,
    options: IrlsOptions,
    requirements: OutputRequirements,
) -> Result<NbinomGlmFit, DeseqError> {
    let input = GlmInput {
        counts,
        offsets: normalization_factors,
        dispersions,
        weights,
    };
    fit_irls_rows(
        IrlsRows::<false> {
            input,
            selection: &[],
        },
        design,
        options,
        requirements,
    )
}

pub(crate) fn fit_selected_irls_means(
    input: GlmInput<'_>,
    design: &DesignMatrix,
    options: IrlsOptions,
    genes: &[usize],
) -> Result<RowMajorMatrix<f64>, DeseqError> {
    if let Offsets::PerSample(factors) = input.offsets {
        validate_size_factor_offsets(input.counts, factors)?;
    }
    fit_irls_rows(
        IrlsRows::<true> {
            input,
            selection: genes,
        },
        design,
        options,
        OutputRequirements::Means,
    )
    .map(|fit| fit.mu)
}

#[derive(Clone, Copy)]
struct IrlsRows<'a, const SELECTED: bool> {
    input: GlmInput<'a>,
    /// Source gene indices in fit order. Output and error positions are compact.
    selection: &'a [usize],
}

impl<'a, const SELECTED: bool> IrlsRows<'a, SELECTED> {
    fn len(self) -> usize {
        if SELECTED {
            self.selection.len()
        } else {
            self.input.counts.n_genes()
        }
    }

    fn source(self, gene: usize) -> usize {
        if SELECTED { self.selection[gene] } else { gene }
    }

    fn weights(self, source_gene: usize) -> Result<Option<&'a [f64]>, DeseqError> {
        self.input
            .weights
            .map(|matrix| matrix.row(source_gene))
            .transpose()
    }
}

fn fit_irls_rows<const SELECTED: bool>(
    rows: IrlsRows<'_, SELECTED>,
    design: &DesignMatrix,
    options: IrlsOptions,
    requirements: OutputRequirements,
) -> Result<NbinomGlmFit, DeseqError> {
    validate_nf_irls_inputs(rows, design, &options)?;
    let samples = rows.input.counts.n_samples();
    let genes = rows.len();

    let x = DMatrix::from_row_slice(
        design.n_samples(),
        design.n_coefficients(),
        design.matrix().as_slice(),
    );
    let initial_lu = (x.transpose() * &x).lu();
    let initial_beta = InitialBeta::for_design(design)?;
    let shared_log_offsets = match rows.input.offsets {
        Offsets::PerSample(factors) => {
            Some(factors.iter().copied().map(f64::ln_p).collect::<Vec<_>>())
        }
        Offsets::PerCell(_) => None,
    };
    let p = design.n_coefficients();
    let ridge_lambda = options.ridge_lambdas_for_coefficients(p)?;
    let n_beta = checked_output_elements(genes, p)?;
    let covariance_columns = checked_output_elements(p, p)?;
    let n_covariance = if requirements.covariance() {
        checked_output_elements(genes, covariance_columns)?
    } else {
        0
    };
    let n_means = checked_output_elements(genes, samples)?;
    // The sizes are checked above; these report a failed allocation as an
    // error instead of aborting.
    let mut beta_values = try_alloc_matrix(n_beta, 1, "IRLS output")?;
    beta_values.resize(n_beta, 0.0);
    let mut beta_var_values = try_alloc_matrix(n_beta, 1, "IRLS output")?;
    beta_var_values.resize(n_beta, f64::NAN);
    let mut beta_optim_start_values = try_alloc_matrix(n_beta, 1, "IRLS output")?;
    beta_optim_start_values.resize(n_beta, f64::NAN);
    let mut beta_covariance_values = try_alloc_matrix(n_covariance, 1, "IRLS output")?;
    beta_covariance_values.resize(n_covariance, f64::NAN);
    let mut mu_values = try_alloc_matrix(n_means, 1, "IRLS output")?;
    mu_values.resize(n_means, 0.0);
    let mut hat_values = try_alloc_matrix(n_means, 1, "IRLS output")?;
    hat_values.resize(n_means, f64::NAN);
    let mut beta_iter = Vec::with_capacity(genes);
    let mut beta_converged = Vec::with_capacity(genes);
    let mut beta_optim_iter = vec![f64::NAN; genes];
    let mut beta_optim_start_objective = vec![f64::NAN; genes];
    let mut beta_optim_objective = vec![f64::NAN; genes];
    let mut beta_optim_gradient_norm = vec![f64::NAN; genes];

    let covariance_rows = if requirements.covariance() {
        rayon::iter::Either::Left(
            beta_covariance_values
                .par_chunks_mut(covariance_columns)
                .map(Some),
        )
    } else {
        rayon::iter::Either::Right((0..genes).into_par_iter().map(|_| None))
    };
    let row_statuses = (
        beta_values.par_chunks_mut(p),
        beta_var_values.par_chunks_mut(p),
        covariance_rows,
        mu_values.par_chunks_mut(samples),
        hat_values.par_chunks_mut(samples),
    )
        .into_par_iter()
        .enumerate()
        .map_init(
            IrlsScratch::default,
            |scratch, (gene, (beta, beta_se, covariance, mu, hat))| {
                fit_fixed_dispersion_irls_row(
                    scratch,
                    IrlsRowInput {
                        gene,
                        shared_log_offsets: shared_log_offsets.as_deref(),
                        rows,
                        x: &x,
                        initial_lu: &initial_lu,
                        initial_beta,
                        ridge_lambda: &ridge_lambda,
                        options: &options,
                    },
                    IrlsRowOutput {
                        beta,
                        beta_se,
                        covariance,
                        mu,
                        hat,
                    },
                )
            },
        )
        .collect::<Vec<_>>();
    // Indexed collection and ordered extraction report the lowest failed gene.
    let mut fallback_rows = Vec::new();
    for (gene, status) in row_statuses.into_iter().enumerate() {
        let status = status?;
        beta_iter.push(status.iter);
        beta_converged.push(status.converged);
        if status.refit {
            fallback_rows.push(gene);
        }
    }

    let mut optim_log_like = vec![None; genes];
    if !fallback_rows.is_empty() {
        refit_optim_fallback_rows(
            &fallback_rows,
            &mut beta_values,
            &mut beta_var_values,
            &mut beta_optim_start_values,
            &mut beta_covariance_values,
            &mut mu_values,
            &mut beta_converged,
            &mut optim_log_like,
            &mut beta_optim_iter,
            &mut beta_optim_start_objective,
            &mut beta_optim_objective,
            &mut beta_optim_gradient_norm,
            OptimFallbackInput {
                rows,
                x: &x,
                ridge_lambda: &ridge_lambda,
                options: &options,
            },
        )?;
    }
    let beta = RowMajorMatrix::from_row_major(genes, p, beta_values)?;
    let beta_se = RowMajorMatrix::from_row_major(genes, p, beta_var_values)?;
    let beta_optim_start = RowMajorMatrix::from_row_major(genes, p, beta_optim_start_values)?;
    let beta_covariance = if requirements.covariance() {
        Some(RowMajorMatrix::from_row_major(
            genes,
            p * p,
            beta_covariance_values,
        )?)
    } else {
        None
    };
    let mu = RowMajorMatrix::from_row_major(genes, samples, mu_values)?;
    let hat_diagonal = RowMajorMatrix::from_row_major(genes, samples, hat_values)?;
    let mut log_like = if SELECTED {
        let mut result = Vec::with_capacity(genes);
        for (gene, &source_gene) in rows.selection.iter().enumerate() {
            result.push(nbinom_log_likelihood_weighted(
                rows.input.counts.row(source_gene)?,
                mu.row(gene)?,
                rows.input.dispersions[source_gene],
                rows.weights(source_gene)?,
            )?);
        }
        result
    } else {
        nbinom_log_likelihood_matrix(
            rows.input.counts,
            &mu,
            rows.input.dispersions,
            rows.input.weights,
        )?
    };
    for (likelihood, optimized) in log_like.iter_mut().zip(optim_log_like) {
        if let Some(value) = optimized {
            *likelihood = value;
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

struct OptimFallbackInput<'a, const SELECTED: bool> {
    rows: IrlsRows<'a, SELECTED>,
    x: &'a DMatrix<f64>,
    ridge_lambda: &'a [f64],
    options: &'a IrlsOptions,
}

struct IrlsRowInput<'a, const SELECTED: bool> {
    gene: usize,
    /// `ln` of the size factors when every gene shares them.
    shared_log_offsets: Option<&'a [f64]>,
    rows: IrlsRows<'a, SELECTED>,
    x: &'a DMatrix<f64>,
    initial_lu: &'a LU<f64, Dyn, Dyn>,
    initial_beta: InitialBeta,
    ridge_lambda: &'a [f64],
    options: &'a IrlsOptions,
}

struct IrlsRowOutput<'a> {
    beta: &'a mut [f64],
    beta_se: &'a mut [f64],
    covariance: Option<&'a mut [f64]>,
    mu: &'a mut [f64],
    hat: &'a mut [f64],
}

struct IrlsRowStatus {
    iter: usize,
    converged: bool,
    refit: bool,
}

// Per-gene IRLS iterations and scratch reuse.
#[path = "irls/iteration.rs"]
mod iteration;
use iteration::*;

// Bounded-optimizer fallback, starting values and acceptance checks.
#[path = "irls/fallback.rs"]
mod fallback;
use fallback::*;

// Optimizer objectives, gradients and Hessians.
#[path = "irls/objective.rs"]
mod objective;
use objective::*;

// Projected-gradient convergence and checked arithmetic.
#[path = "irls/arithmetic.rs"]
mod arithmetic;
use arithmetic::*;

// Starting coefficients and stable fitted means.
#[path = "irls/means.rs"]
mod means;
use means::*;

// Working weights and responses.
#[path = "irls/working.rs"]
mod working;
use working::*;

// Weighted least squares, QR, covariance and design crossproducts.
#[path = "irls/linear_algebra.rs"]
mod linear_algebra;
use linear_algebra::*;

// Input validation, normalization factors and output scaling.
#[path = "irls/validation.rs"]
mod validation;
pub(super) use validation::*;

// Private numerical regression tests.
#[cfg(test)]
#[path = "irls/tests.rs"]
mod tests;
