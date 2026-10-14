use nalgebra::{DMatrix, DVector};
use rayon::prelude::*;
use rcompat_lbfgsb::{Bounds, OptimControl, optim_lbfgsb, optim_lbfgsb_with_gradient};
use std::cell::RefCell;

use crate::core::CountMatrix;
use crate::design::DesignMatrix;
use crate::errors::{DeseqError, invalid_dimensions};
use crate::glm::NbinomGlmFit;
use crate::glm::beta::fit_intercept_only_fixed_dispersion;
use crate::glm::fallback::optim_fallback_rows;
use crate::glm::nb::{nbinom_log_likelihood_matrix, nbinom_log_likelihood_weighted};
use crate::glm::r_arithmetic::{nbinom_log, normal_log};
use crate::math::optim::BoundedOptimizationOutput;
use crate::matrix::RowMajorMatrix;

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
    let p = design.n_coefficients();
    let ridge_lambda = options.ridge_lambdas_for_coefficients(p)?;
    let mut beta_values = Vec::with_capacity(counts.n_genes() * p);
    let mut beta_var_values = Vec::with_capacity(counts.n_genes() * p);
    let mut beta_optim_start_values = vec![f64::NAN; counts.n_genes() * p];
    let mut beta_covariance_values = Vec::with_capacity(counts.n_genes() * p * p);
    let mut mu_values = Vec::with_capacity(counts.n_genes() * counts.n_samples());
    let mut hat_values = Vec::with_capacity(counts.n_genes() * counts.n_samples());
    let mut beta_iter = Vec::with_capacity(counts.n_genes());
    let mut beta_converged = Vec::with_capacity(counts.n_genes());
    let mut beta_optim_iter = vec![f64::NAN; counts.n_genes()];
    let mut beta_optim_start_objective = vec![f64::NAN; counts.n_genes()];
    let mut beta_optim_objective = vec![f64::NAN; counts.n_genes()];
    let mut beta_optim_gradient_norm = vec![f64::NAN; counts.n_genes()];

    let row_outputs = (0..counts.n_genes())
        .into_par_iter()
        .map(|gene| {
            fit_fixed_dispersion_irls_row(IrlsRowInput {
                gene,
                counts,
                x: &x,
                normalization_factors,
                dispersions,
                weights,
                ridge_lambda: &ridge_lambda,
                options: &options,
            })
        })
        .collect::<Result<Vec<_>, DeseqError>>()?;

    for row in row_outputs {
        beta_values.extend(row.beta);
        beta_var_values.extend(row.beta_se);
        beta_covariance_values.extend(row.beta_covariance);
        mu_values.extend(row.mu);
        hat_values.extend(row.hat);
        beta_iter.push(row.iter);
        beta_converged.push(row.converged);
    }

    let beta_for_routing =
        RowMajorMatrix::from_row_major(counts.n_genes(), p, beta_values.clone())?;
    let covariance_for_routing =
        RowMajorMatrix::from_row_major(counts.n_genes(), p * p, beta_covariance_values.clone())?;
    let fallback_rows = optim_fallback_rows(
        &beta_converged,
        &beta_for_routing,
        &covariance_for_routing,
        options.use_optim,
        options.force_optim,
    )?;
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
    let beta_covariance =
        RowMajorMatrix::from_row_major(counts.n_genes(), p * p, beta_covariance_values)?;
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
        beta_covariance: Some(beta_covariance),
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

fn fit_fixed_dispersion_irls_row(input: IrlsRowInput<'_>) -> Result<IrlsRowOutput, DeseqError> {
    let gene = input.gene;
    let p = input.x.ncols();
    if input.counts.is_all_zero_gene(gene)? {
        return Err(DeseqError::InvalidCounts {
            reason: format!("gene {gene} is all zero; DESeq2 GLM fitting excludes allZero rows"),
        });
    }
    let y = input
        .counts
        .row(gene)?
        .iter()
        .copied()
        .map(f64::from)
        .collect::<Vec<_>>();
    let nf = input.normalization_factors.row(gene)?;
    let weight_row = input.weights.map(|matrix| matrix.row(gene)).transpose()?;
    let dispersion = input.dispersions[gene];
    let mut beta = initial_beta(input.x, &y, nf)?;
    let mut mu = fitted_mu(input.x, &beta, nf, input.options.min_mu)?;
    let mut dev_old = 0.0;
    let mut dev = 0.0;
    let mut iter = 0_usize;
    let mut converged = false;

    for t in 0..input.options.maxit {
        iter += 1;
        let w = working_weights(&mu, dispersion, weight_row)?;
        let z = working_response(&mu, nf, &y)?;
        let Some(next_beta) =
            solve_weighted_least_squares(input.x, &w, &z, input.ridge_lambda, input.options.solver)
        else {
            iter = input.options.maxit;
            break;
        };
        beta = next_beta;
        if beta
            .iter()
            .any(|value| !value.is_finite() || value.abs() > input.options.max_beta_abs)
        {
            iter = input.options.maxit;
            break;
        }
        mu = fitted_mu(input.x, &beta, nf, input.options.min_mu)?;
        dev = -2.0
            * nbinom_log_likelihood_weighted(input.counts.row(gene)?, &mu, dispersion, weight_row)?;
        let Some(conv_test) = irls_deviance_convergence_stat(dev, dev_old) else {
            iter = input.options.maxit;
            break;
        };
        if t > 0 && conv_test < input.options.beta_tol {
            converged = true;
            break;
        }
        dev_old = dev;
    }

    let mut row = IrlsRowOutput {
        beta: Vec::with_capacity(p),
        beta_se: Vec::with_capacity(p),
        beta_covariance: Vec::with_capacity(p * p),
        mu: Vec::with_capacity(input.counts.n_samples()),
        hat: Vec::with_capacity(input.counts.n_samples()),
        iter,
        converged: false,
    };

    let w = working_weights(&mu, dispersion, weight_row)?;
    let Some((beta_covariance, hat_diag)) =
        covariance_and_hat_diagonal(input.x, &w, input.ridge_lambda)
    else {
        row.iter = input.options.maxit;
        row.beta_se.extend(std::iter::repeat_n(f64::NAN, p));
        row.beta_covariance
            .extend(std::iter::repeat_n(f64::NAN, p * p));
        row.hat
            .extend(std::iter::repeat_n(f64::NAN, input.counts.n_samples()));
        for (col, value) in beta.iter().copied().enumerate() {
            row.beta
                .push(checked_log2_scale(value, col, "IRLS beta log2 scale")?);
        }
        row.mu
            .extend(fitted_mu_unfloored(input.x, &beta, nf)?.iter().copied());
        return Ok(row);
    };

    for (col, value) in beta.iter().copied().enumerate() {
        row.beta
            .push(checked_log2_scale(value, col, "IRLS beta log2 scale")?);
    }
    for diagonal in 0..p {
        let value = beta_covariance[diagonal * p + diagonal];
        row.beta_se.push(checked_log2_standard_error(
            value,
            diagonal,
            "IRLS beta standard error",
        )?);
    }
    for (idx, value) in beta_covariance.into_iter().enumerate() {
        row.beta_covariance.push(checked_log2_covariance(
            value,
            idx,
            "IRLS beta covariance log2 scale",
        )?);
    }
    row.mu
        .extend(fitted_mu_unfloored(input.x, &beta, nf)?.iter().copied());
    row.hat.extend(hat_diag);
    row.converged = converged && row.iter < input.options.maxit;
    let _ = dev;
    Ok(row)
}

#[allow(clippy::too_many_arguments)]
fn refit_optim_fallback_rows(
    rows: &[usize],
    beta_values: &mut [f64],
    beta_var_values: &mut [f64],
    beta_optim_start_values: &mut [f64],
    beta_covariance_values: &mut [f64],
    mu_values: &mut [f64],
    _hat_values: &mut [f64],
    beta_converged: &mut [bool],
    optim_log_like: &mut [Option<f64>],
    beta_optim_iter: &mut [f64],
    beta_optim_start_objective: &mut [f64],
    beta_optim_objective: &mut [f64],
    beta_optim_gradient_norm: &mut [f64],
    input: OptimFallbackInput<'_>,
) -> Result<(), DeseqError> {
    let p = input.x.ncols();
    let n = input.x.nrows();
    for gene in rows.iter().copied() {
        let counts_row = input.counts.row(gene)?;
        let nf = input.normalization_factors.row(gene)?;
        let weight_row = input.weights.map(|matrix| matrix.row(gene)).transpose()?;
        let dispersion = input.dispersions[gene];
        let beta_start = optim_start_beta_log2(
            &beta_values[gene * p..(gene + 1) * p],
            input.x,
            counts_row,
            nf,
            input.options.max_beta_abs,
        )?;
        let beta_input = BetaOptimInput {
            x: input.x,
            counts: counts_row,
            nf,
            dispersion,
            weights: weight_row,
            ridge_lambda: input.ridge_lambda,
        };
        let start_objective = if input.options.r_optim_compat {
            beta_log2_r_objective(&beta_input, &beta_start)?
        } else {
            beta_log2_objective(&beta_input, &beta_start)?
        };
        let output = optimize_beta_log2(beta_input, &beta_start, input.options)?;
        let (_, final_gradient, _) =
            beta_log2_objective_gradient_hessian(&beta_input, &output.parameters)?;
        let final_gradient_norm = projected_gradient_norm(
            &output.parameters,
            &final_gradient,
            -input.options.max_beta_abs,
            input.options.max_beta_abs,
        )
        .unwrap_or(f64::NAN);

        let mu_unfloored = fitted_mu_log2_unfloored(input.x, &output.parameters, nf)?;
        let mu_for_inference = mu_unfloored
            .iter()
            .copied()
            .map(|value| value.max(input.options.min_mu))
            .collect::<Vec<_>>();
        let w = working_weights(&mu_for_inference, dispersion, weight_row)?;
        let Some((beta_covariance, _hat_diag)) =
            covariance_and_hat_diagonal(input.x, &w, input.ridge_lambda)
        else {
            return Err(DeseqError::UnsupportedFeature {
                feature: "optim fallback covariance is singular".to_string(),
            });
        };

        for (col, value) in output.parameters.iter().copied().enumerate() {
            beta_optim_start_values[gene * p + col] = beta_start[col];
            beta_values[gene * p + col] = value;
            let covariance_value = beta_covariance[col * p + col];
            beta_var_values[gene * p + col] = checked_log2_standard_error(
                covariance_value,
                gene * p + col,
                "optim fallback beta standard error",
            )?;
        }
        for (idx, value) in beta_covariance.into_iter().enumerate() {
            beta_covariance_values[gene * p * p + idx] = checked_log2_covariance(
                value,
                gene * p * p + idx,
                "optim fallback beta covariance log2 scale",
            )?;
        }
        for (sample, value) in mu_unfloored.iter().copied().enumerate() {
            mu_values[gene * n + sample] = value;
        }
        // DESeq2 updates beta, SE, mu, and log-likelihood after L-BFGS-B
        // fallback, but preserves the original IRLS hat diagonals returned by
        // fitBetaWrapper. Cook's replacement relies on that pre-optim leverage.
        beta_converged[gene] = output.converged;
        beta_optim_iter[gene] = output.iterations as f64;
        beta_optim_start_objective[gene] = start_objective;
        beta_optim_objective[gene] = output.value;
        beta_optim_gradient_norm[gene] = final_gradient_norm;
        optim_log_like[gene] = Some(nbinom_log_likelihood_weighted(
            counts_row,
            &mu_for_inference,
            dispersion,
            weight_row,
        )?);
    }
    Ok(())
}

fn optim_start_beta_log2(
    current_beta_log2: &[f64],
    x: &DMatrix<f64>,
    counts: &[u32],
    nf: &[f64],
    bound: f64,
) -> Result<Vec<f64>, DeseqError> {
    if current_beta_log2
        .iter()
        .copied()
        .all(|value| value.is_finite() && value.abs() <= bound)
    {
        return Ok(current_beta_log2
            .iter()
            .copied()
            .map(|value| value.clamp(-bound, bound))
            .collect());
    }

    let y = counts.iter().copied().map(f64::from).collect::<Vec<_>>();
    Ok(initial_beta(x, &y, nf)?
        .iter()
        .copied()
        .map(|value| value.clamp(-bound, bound))
        .collect())
}

#[derive(Clone, Copy)]
struct BetaOptimInput<'a> {
    x: &'a DMatrix<f64>,
    counts: &'a [u32],
    nf: &'a [f64],
    dispersion: f64,
    weights: Option<&'a [f64]>,
    ridge_lambda: &'a [f64],
}

fn optimize_beta_log2(
    input: BetaOptimInput<'_>,
    start: &[f64],
    options: &IrlsOptions,
) -> Result<BoundedOptimizationOutput, DeseqError> {
    minimize_beta_log2_lbfgsb(
        input,
        start,
        -options.max_beta_abs,
        options.max_beta_abs,
        options.optim_maxit,
        options.optim_tol,
        options.r_optim_compat,
    )
}

fn minimize_beta_log2_lbfgsb(
    input: BetaOptimInput<'_>,
    start: &[f64],
    lower: f64,
    upper: f64,
    maxit: usize,
    gradient_tol: f64,
    r_optim_compat: bool,
) -> Result<BoundedOptimizationOutput, DeseqError> {
    let parameters = start
        .iter()
        .copied()
        .map(|value| value.clamp(lower, upper))
        .collect::<Vec<_>>();
    let start_parameters = parameters.clone();
    let deferred_error = RefCell::new(None);
    let mut objective = |beta: &[f64]| match if r_optim_compat {
        beta_log2_r_objective(&input, beta)
    } else {
        beta_log2_objective(&input, beta)
    } {
        Ok(value) => value,
        Err(error) => {
            deferred_error.borrow_mut().get_or_insert(error);
            1.0e300
        }
    };
    let start_value = if r_optim_compat {
        beta_log2_r_objective(&input, &parameters)?
    } else {
        beta_log2_objective(&input, &parameters)?
    };
    let bounds = Bounds::new(vec![lower; parameters.len()], vec![upper; parameters.len()])
        .map_err(|error| DeseqError::UnsupportedFeature {
            feature: format!("optim fallback L-BFGS-B bounds failed: {error}"),
        })?;
    let control = optim_fallback_control(parameters.len(), maxit, gradient_tol);
    let solution = if r_optim_compat {
        let numeric = match optim_lbfgsb(
            parameters.clone(),
            bounds.clone(),
            &mut objective,
            control.clone(),
        ) {
            Ok(solution) => solution,
            Err(error) => {
                return Err(DeseqError::UnsupportedFeature {
                    feature: format!("optim fallback L-BFGS-B failed: {error}"),
                });
            }
        };
        let mut chosen = numeric;
        let mut analytic_objective =
            |beta: &[f64]| beta_log2_objective(&input, beta).unwrap_or(1.0e300);
        let mut analytic_gradient = |beta: &[f64]| {
            beta_log2_value_gradient(&input, beta)
                .map(|(_, gradient)| gradient)
                .unwrap_or_else(|_| vec![0.0; beta.len()])
        };
        if let Ok(mut analytic) = optim_lbfgsb_with_gradient(
            parameters,
            bounds,
            &mut analytic_objective,
            &mut analytic_gradient,
            control.clone(),
        ) {
            let numeric_reduced = beta_log2_objective(&input, &chosen.par)?;
            let analytic_reduced = beta_log2_objective(&input, &analytic.par)?;
            if analytic.is_success()
                && !r_optim_candidate_within_analytic_tolerance(
                    numeric_reduced,
                    analytic_reduced,
                    control.factr,
                )
                && r_optim_candidates_materially_disagree(
                    &chosen.par,
                    &analytic.par,
                    &control.ndeps,
                    &control.parscale,
                )
            {
                analytic.value = beta_log2_r_objective(&input, &analytic.par)?;
                chosen = analytic;
            }
        }
        Ok(chosen)
    } else {
        let mut gradient = |beta: &[f64]| match beta_log2_value_gradient(&input, beta) {
            Ok((_, gradient)) => gradient,
            Err(error) => {
                deferred_error.borrow_mut().get_or_insert(error);
                vec![0.0; beta.len()]
            }
        };
        optim_lbfgsb_with_gradient(parameters, bounds, &mut objective, &mut gradient, control)
    }
    .map_err(|error| DeseqError::UnsupportedFeature {
        feature: format!("optim fallback L-BFGS-B failed: {error}"),
    })?;
    if let Some(error) = deferred_error.into_inner() {
        return Err(error);
    }
    let (_, final_gradient) = beta_log2_value_gradient(&input, &solution.par)?;
    let final_gradient_norm = projected_gradient_norm(&solution.par, &final_gradient, lower, upper);
    let converged_by_gradient = final_gradient_norm
        .is_some_and(|norm| norm <= gradient_tol.max(1.0e-4) && solution.value.is_finite());
    let converged = solution.is_success() || converged_by_gradient;
    let iterations = solution.gradient_count();
    if !converged && final_gradient_norm.is_some_and(|norm| norm > 0.1) {
        return Ok(BoundedOptimizationOutput {
            parameters: start_parameters,
            value: start_value,
            converged: false,
            iterations,
        });
    }
    Ok(BoundedOptimizationOutput {
        parameters: solution.par,
        value: solution.value,
        converged,
        iterations,
    })
}

fn optim_fallback_control(dimension: usize, maxit: usize, gradient_tol: f64) -> OptimControl {
    let mut control = OptimControl::default_for_dimension(dimension);
    control.maxit = maxit;
    control.factr = 1.0e7;
    control.pgtol = gradient_tol;
    control.lmm = 5;
    control
}

fn r_optim_candidate_within_analytic_tolerance(
    compatible_objective: f64,
    analytic_objective: f64,
    factr: f64,
) -> bool {
    if !compatible_objective.is_finite()
        || !analytic_objective.is_finite()
        || !factr.is_finite()
        || factr < 0.0
    {
        return false;
    }
    // `factr * epsilon` is the L-BFGS-B relative-reduction tolerance. Allow a
    // tenfold band for the deliberately finite-difference-compatible endpoint;
    // beyond it the callback path and analytic solution disagree materially.
    let scale = compatible_objective
        .abs()
        .max(analytic_objective.abs())
        .max(1.0);
    let tolerance = 10.0 * factr * f64::EPSILON * scale;
    compatible_objective <= analytic_objective + tolerance
}

fn r_optim_candidates_materially_disagree(
    compatible: &[f64],
    analytic: &[f64],
    ndeps: &[f64],
    parscale: &[f64],
) -> bool {
    if compatible.len() != analytic.len()
        || compatible.len() != ndeps.len()
        || compatible.len() != parscale.len()
    {
        return true;
    }
    compatible
        .iter()
        .zip(analytic)
        .zip(ndeps.iter().zip(parscale))
        .any(|((&r_value, &analytic_value), (&step, &scale))| {
            let resolution = step * scale;
            !r_value.is_finite()
                || !analytic_value.is_finite()
                || !resolution.is_finite()
                || resolution <= 0.0
                || (r_value - analytic_value).abs() > 10.0 * resolution
        })
}

fn beta_log2_value_gradient(
    input: &BetaOptimInput<'_>,
    beta: &[f64],
) -> Result<(f64, Vec<f64>), DeseqError> {
    let (value, gradient, _) = beta_log2_objective_gradient_hessian(input, beta)?;
    Ok((value, gradient))
}

fn beta_log2_objective_gradient_hessian(
    input: &BetaOptimInput<'_>,
    beta: &[f64],
) -> Result<(f64, Vec<f64>, DMatrix<f64>), DeseqError> {
    let p = input.x.ncols();
    if beta.len() != p {
        return Err(invalid_dimensions("optim beta coefficients", p, beta.len()));
    }
    if input.ridge_lambda.len() != p {
        return Err(invalid_dimensions(
            "optim ridge coefficients",
            p,
            input.ridge_lambda.len(),
        ));
    }
    if input.counts.len() != input.x.nrows() || input.nf.len() != input.x.nrows() {
        return Err(invalid_dimensions(
            "optim samples",
            input.x.nrows(),
            input.counts.len().min(input.nf.len()),
        ));
    }
    if let Some(weights) = input.weights
        && weights.len() != input.x.nrows()
    {
        return Err(invalid_dimensions(
            "optim weights",
            input.x.nrows(),
            weights.len(),
        ));
    }
    validate_positive_finite(input.dispersion, "dispersion", 0)?;

    let mut objective = 0.0;
    let mut gradient = vec![0.0; p];
    let mut hessian = DMatrix::zeros(p, p);
    let ln2 = std::f64::consts::LN_2;
    let ln2_squared = ln2 * ln2;
    for sample in 0..input.x.nrows() {
        validate_positive_finite(input.nf[sample], "normalization factor", sample)?;
        let weight = input.weights.map_or(1.0, |weights| weights[sample]);
        validate_nonnegative_finite(weight, "weight", sample)?;
        let mut eta = 0.0;
        for (col, beta_value) in beta.iter().copied().enumerate().take(p) {
            let next_eta = input.x[(sample, col)].mul_add(beta_value, eta);
            if !next_eta.is_finite() {
                return Ok(beta_optim_penalty(p));
            }
            eta = next_eta;
        }
        let mu = input.nf[sample] * 2.0_f64.powf(eta);
        if !mu.is_finite() || mu <= 0.0 {
            return Ok(beta_optim_penalty(p));
        }
        let objective_term =
            match beta_log2_nb_objective_term(input.counts[sample], mu, input.dispersion) {
                Some(value) => value,
                None => return Ok(beta_optim_penalty(p)),
            };
        let Some(weighted_objective_term) = checked_product2(weight, objective_term) else {
            return Ok(beta_optim_penalty(p));
        };
        let Some(next_objective) = checked_sum2(objective, weighted_objective_term) else {
            return Ok(beta_optim_penalty(p));
        };
        objective = next_objective;
        let Some(disp_mu) = checked_product2(input.dispersion, mu) else {
            return Ok(beta_optim_penalty(p));
        };
        let Some(one_plus_disp_mu) = checked_sum2(1.0, disp_mu) else {
            return Ok(beta_optim_penalty(p));
        };
        if !one_plus_disp_mu.is_finite() || one_plus_disp_mu <= 0.0 {
            return Ok(beta_optim_penalty(p));
        }
        let inv_one_plus_disp_mu = one_plus_disp_mu.recip();
        let count_residual = f64::from(input.counts[sample]) - mu;
        if !count_residual.is_finite() {
            return Ok(beta_optim_penalty(p));
        }
        let Some(sample_score) =
            checked_product4(weight, ln2, count_residual, inv_one_plus_disp_mu)
        else {
            return Ok(beta_optim_penalty(p));
        };
        for (col, gradient_value) in gradient.iter_mut().enumerate().take(p) {
            let Some(term) = checked_product2(-input.x[(sample, col)], sample_score) else {
                return Ok(beta_optim_penalty(p));
            };
            let Some(next_gradient) = checked_sum2(*gradient_value, term) else {
                return Ok(beta_optim_penalty(p));
            };
            *gradient_value = next_gradient;
        }
        let Some(disp_count) = checked_product2(input.dispersion, f64::from(input.counts[sample]))
        else {
            return Ok(beta_optim_penalty(p));
        };
        let Some(one_plus_disp_count) = checked_sum2(1.0, disp_count) else {
            return Ok(beta_optim_penalty(p));
        };
        let Some(sample_hessian_weight) = checked_product6(
            weight,
            ln2_squared,
            mu,
            one_plus_disp_count,
            inv_one_plus_disp_mu,
            inv_one_plus_disp_mu,
        ) else {
            return Ok(beta_optim_penalty(p));
        };
        for row in 0..p {
            for col in 0..p {
                let Some(term) = checked_product3(
                    input.x[(sample, row)],
                    sample_hessian_weight,
                    input.x[(sample, col)],
                ) else {
                    return Ok(beta_optim_penalty(p));
                };
                let Some(next_hessian) = checked_sum2(hessian[(row, col)], term) else {
                    return Ok(beta_optim_penalty(p));
                };
                hessian[(row, col)] = next_hessian;
            }
        }
    }

    if !objective.is_finite() {
        return Ok(beta_optim_penalty(p));
    }
    for col in 0..p {
        validate_nonnegative_finite(input.ridge_lambda[col], "ridge lambda", col)?;
        let Some(ridge_log2) = checked_product2(input.ridge_lambda[col], ln2_squared) else {
            return Ok(beta_optim_penalty(p));
        };
        let Some(objective_term) = checked_product4(0.5, ridge_log2, beta[col], beta[col]) else {
            return Ok(beta_optim_penalty(p));
        };
        let Some(gradient_term) = checked_product2(ridge_log2, beta[col]) else {
            return Ok(beta_optim_penalty(p));
        };
        let Some(next_objective) = checked_sum2(objective, objective_term) else {
            return Ok(beta_optim_penalty(p));
        };
        let Some(next_gradient) = checked_sum2(gradient[col], gradient_term) else {
            return Ok(beta_optim_penalty(p));
        };
        let Some(next_hessian) = checked_sum2(hessian[(col, col)], ridge_log2) else {
            return Ok(beta_optim_penalty(p));
        };
        objective = next_objective;
        gradient[col] = next_gradient;
        hessian[(col, col)] = next_hessian;
    }
    Ok((objective, gradient, hessian))
}

fn beta_log2_nb_objective_term(count: u32, mu: f64, dispersion: f64) -> Option<f64> {
    validate_positive_finite(mu, "optim fallback mean", 0).ok()?;
    validate_positive_finite(dispersion, "dispersion", 0).ok()?;
    let count = f64::from(count);
    let size = dispersion.recip();
    if !size.is_finite() || size <= 0.0 {
        return None;
    }
    let size_plus_mu = checked_sum2(size, mu)?;
    let count_plus_size = checked_sum2(count, size)?;
    let first = checked_product2(count_plus_size, size_plus_mu.ln())?;
    let second = checked_product2(size, size.ln())?;
    let third = checked_product2(count, mu.ln())?;
    checked_sum2(checked_sum2(first, -second)?, -third)
}

fn beta_optim_penalty(p: usize) -> (f64, Vec<f64>, DMatrix<f64>) {
    (1.0e300, vec![0.0; p], DMatrix::identity(p, p))
}

fn beta_log2_objective(input: &BetaOptimInput<'_>, beta: &[f64]) -> Result<f64, DeseqError> {
    beta_log2_objective_gradient_hessian(input, beta).map(|(objective, _, _)| objective)
}

fn beta_log2_r_objective(input: &BetaOptimInput<'_>, beta: &[f64]) -> Result<f64, DeseqError> {
    let mut ll = 0.0;
    for sample in 0..input.x.nrows() {
        let mut eta = 0.0;
        for (col, b) in beta.iter().copied().enumerate() {
            eta = input.x[(sample, col)].mul_add(b, eta)
        }
        let mu = input.nf[sample] * 2.0_f64.powf(eta);
        if !mu.is_finite() || mu <= 0.0 {
            return Ok(1e300);
        }
        ll += input.weights.map_or(1.0, |w| w[sample])
            * nbinom_log(input.counts[sample], mu, input.dispersion)
    }
    let l2 = std::f64::consts::LN_2.powi(2);
    let mut lp = 0.0;
    for (col, b) in beta.iter().copied().enumerate() {
        let precision = input.ridge_lambda[col] * l2;
        if precision > 0.0 {
            lp += normal_log(b, precision);
        }
    }
    Ok(-(ll + lp))
}

#[cfg(test)]
fn beta_log2_numeric_gradient(
    input: &BetaOptimInput<'_>,
    beta: &[f64],
    lower: f64,
    upper: f64,
) -> Result<Vec<f64>, DeseqError> {
    const R_OPTIM_NDEPS: f64 = 1.0e-3;
    let mut gradient = Vec::with_capacity(beta.len());
    for col in 0..beta.len() {
        let forward = (beta[col] + R_OPTIM_NDEPS).min(upper);
        let forward_eps = forward - beta[col];
        let backward = (beta[col] - R_OPTIM_NDEPS).max(lower);
        let backward_eps = beta[col] - backward;
        let denominator = forward_eps + backward_eps;
        if !denominator.is_finite() || denominator <= 0.0 {
            return Err(DeseqError::UnsupportedFeature {
                feature: "optim fallback finite-difference gradient at degenerate bounds"
                    .to_string(),
            });
        }
        let mut forward_beta = beta.to_vec();
        forward_beta[col] = forward;
        let mut backward_beta = beta.to_vec();
        backward_beta[col] = backward;
        let forward_value = beta_log2_objective(input, &forward_beta)?;
        let backward_value = beta_log2_objective(input, &backward_beta)?;
        let Some(difference) = checked_sum2(forward_value, -backward_value) else {
            return Ok(vec![0.0; beta.len()]);
        };
        let Some(value) = checked_div2(difference, denominator) else {
            return Ok(vec![0.0; beta.len()]);
        };
        gradient.push(value);
    }
    Ok(gradient)
}

fn projected_gradient_norm(
    parameters: &[f64],
    gradient: &[f64],
    lower: f64,
    upper: f64,
) -> Option<f64> {
    let scale = parameters
        .iter()
        .copied()
        .zip(gradient.iter().copied())
        .map(|(parameter, gradient)| {
            if (parameter <= lower && gradient > 0.0) || (parameter >= upper && gradient < 0.0) {
                0.0
            } else {
                gradient.abs()
            }
        })
        .try_fold(0.0_f64, |scale, value| {
            value.is_finite().then_some(scale.max(value))
        })?;
    if scale == 0.0 {
        return Some(0.0);
    }
    let mut sum = 0.0;
    for (parameter, gradient) in parameters.iter().copied().zip(gradient.iter().copied()) {
        let value =
            if (parameter <= lower && gradient > 0.0) || (parameter >= upper && gradient < 0.0) {
                0.0
            } else {
                gradient
            };
        let scaled = checked_div2(value, scale)?;
        let term = checked_product2(scaled, scaled)?;
        let next = sum + term;
        if !term.is_finite() || !next.is_finite() {
            return None;
        }
        sum = next;
    }
    let norm = scale * sum.sqrt();
    norm.is_finite().then_some(norm)
}

fn checked_sum2(left: f64, right: f64) -> Option<f64> {
    let sum = left + right;
    (left.is_finite() && right.is_finite() && sum.is_finite()).then_some(sum)
}

fn checked_product2(left: f64, right: f64) -> Option<f64> {
    let product = left * right;
    (left.is_finite() && right.is_finite() && product.is_finite()).then_some(product)
}

fn checked_div2(left: f64, right: f64) -> Option<f64> {
    let quotient = left / right;
    (left.is_finite() && right.is_finite() && right != 0.0 && quotient.is_finite())
        .then_some(quotient)
}

fn checked_product3(left: f64, middle: f64, right: f64) -> Option<f64> {
    checked_product2(checked_product2(left, middle)?, right)
}

fn checked_product4(first: f64, second: f64, third: f64, fourth: f64) -> Option<f64> {
    checked_product2(
        checked_product2(first, second)?,
        checked_product2(third, fourth)?,
    )
}

fn checked_product6(
    first: f64,
    second: f64,
    third: f64,
    fourth: f64,
    fifth: f64,
    sixth: f64,
) -> Option<f64> {
    checked_product2(
        checked_product3(first, second, third)?,
        checked_product3(fourth, fifth, sixth)?,
    )
}

fn irls_deviance_convergence_stat(dev: f64, dev_old: f64) -> Option<f64> {
    let delta = checked_sum2(dev, -dev_old)?.abs();
    let scale = checked_sum2(dev.abs(), 0.1)?;
    checked_div2(delta, scale)
}

fn initial_beta(x: &DMatrix<f64>, y: &[f64], nf: &[f64]) -> Result<DVector<f64>, DeseqError> {
    let response = y
        .iter()
        .copied()
        .zip(nf.iter().copied())
        .map(|(count, factor)| {
            validate_positive_finite(factor, "normalization factor", 0)?;
            Ok((count / factor + 0.1).ln())
        })
        .collect::<Result<Vec<_>, DeseqError>>()?;
    let xtx = x.transpose() * x;
    let xty = x.transpose() * DVector::from_vec(response);
    if let Some(beta) = xtx.lu().solve(&xty) {
        Ok(beta)
    } else {
        Ok(DVector::from_element(x.ncols(), 0.0))
    }
}

fn fitted_mu(
    x: &DMatrix<f64>,
    beta: &DVector<f64>,
    nf: &[f64],
    min_mu: f64,
) -> Result<Vec<f64>, DeseqError> {
    fitted_mu_impl(x, beta, nf, Some(min_mu))
}

fn fitted_mu_unfloored(
    x: &DMatrix<f64>,
    beta: &DVector<f64>,
    nf: &[f64],
) -> Result<Vec<f64>, DeseqError> {
    fitted_mu_impl(x, beta, nf, None)
}

fn fitted_mu_log2_unfloored(
    x: &DMatrix<f64>,
    beta: &[f64],
    nf: &[f64],
) -> Result<Vec<f64>, DeseqError> {
    if beta.len() != x.ncols() {
        return Err(invalid_dimensions(
            "log2 beta coefficients",
            x.ncols(),
            beta.len(),
        ));
    }
    (0..x.nrows())
        .map(|sample| {
            validate_positive_finite(nf[sample], "normalization factor", sample)?;
            let eta = checked_row_dot_slice(x, sample, beta).ok_or(DeseqError::NonFiniteValue {
                context: "optim fallback linear predictor".to_string(),
                index: Some(sample),
                value: f64::NAN,
            })?;
            let mu = nf[sample] * 2.0_f64.powf(eta);
            if !mu.is_finite() || mu <= 0.0 {
                return Err(DeseqError::NonFiniteValue {
                    context: "optim fallback fitted mean".to_string(),
                    index: Some(sample),
                    value: mu,
                });
            }
            Ok(mu)
        })
        .collect()
}

fn fitted_mu_impl(
    x: &DMatrix<f64>,
    beta: &DVector<f64>,
    nf: &[f64],
    min_mu: Option<f64>,
) -> Result<Vec<f64>, DeseqError> {
    (0..x.nrows())
        .zip(nf.iter().copied())
        .enumerate()
        .map(|(sample, (row, factor))| {
            validate_positive_finite(factor, "normalization factor", sample)?;
            let eta = checked_row_dot_vector(x, row, beta).ok_or(DeseqError::NonFiniteValue {
                context: "IRLS linear predictor".to_string(),
                index: Some(sample),
                value: f64::NAN,
            })?;
            let mu = finite_scaled_exp(factor, eta).ok_or(DeseqError::NonFiniteValue {
                context: "IRLS fitted mean".to_string(),
                index: Some(sample),
                value: f64::NAN,
            })?;
            if mu <= 0.0 {
                return Err(DeseqError::NonFiniteValue {
                    context: "IRLS fitted mean".to_string(),
                    index: Some(sample),
                    value: mu,
                });
            }
            Ok(min_mu.map_or(mu, |min_mu| mu.max(min_mu)))
        })
        .collect()
}

fn checked_row_dot_slice(x: &DMatrix<f64>, row: usize, beta: &[f64]) -> Option<f64> {
    let mut sum = 0.0;
    for col in 0..x.ncols() {
        sum = x[(row, col)].mul_add(beta[col], sum);
        if !sum.is_finite() {
            return None;
        }
    }
    Some(sum)
}

fn checked_row_dot_vector(x: &DMatrix<f64>, row: usize, beta: &DVector<f64>) -> Option<f64> {
    let mut sum = ScaledSum::default();
    for col in 0..x.ncols() {
        sum.add(checked_product2(x[(row, col)], beta[col])?)?;
    }
    sum.finish()
}

#[derive(Default)]
struct ScaledSum {
    scale: f64,
    normalized_sum: f64,
}

impl ScaledSum {
    fn add(&mut self, value: f64) -> Option<()> {
        if !value.is_finite() {
            return None;
        }
        let abs = value.abs();
        if abs == 0.0 {
            return Some(());
        }
        if self.scale == 0.0 {
            self.scale = abs;
            self.normalized_sum = value / abs;
            return Some(());
        }
        if abs > self.scale {
            let rescaled = self.normalized_sum * (self.scale / abs);
            self.normalized_sum = checked_sum2(rescaled, value / abs)?;
            self.scale = abs;
        } else {
            self.normalized_sum = checked_sum2(self.normalized_sum, value / self.scale)?;
        }
        Some(())
    }

    fn finish(self) -> Option<f64> {
        if self.scale == 0.0 {
            return Some(0.0);
        }
        checked_product2(self.normalized_sum, self.scale)
    }
}

fn finite_scaled_exp(factor: f64, eta: f64) -> Option<f64> {
    if !factor.is_finite() || factor <= 0.0 || !eta.is_finite() {
        return None;
    }
    let log_mu = factor.ln() + eta;
    if log_mu >= f64::MAX.ln() {
        return Some(f64::MAX);
    }
    if log_mu <= f64::MIN_POSITIVE.ln() {
        return Some(f64::MIN_POSITIVE);
    }
    let mu = log_mu.exp();
    (mu.is_finite() && mu > 0.0).then_some(mu)
}

fn working_weights(
    mu: &[f64],
    dispersion: f64,
    weights: Option<&[f64]>,
) -> Result<Vec<f64>, DeseqError> {
    validate_positive_finite(dispersion, "dispersion", 0)?;
    mu.iter()
        .copied()
        .enumerate()
        .map(|(sample, value)| {
            validate_positive_finite(value, "mu", sample)?;
            let working_weight =
                stable_working_weight(value, dispersion).ok_or(DeseqError::NonFiniteValue {
                    context: "IRLS working weight".to_string(),
                    index: Some(sample),
                    value: f64::NAN,
                })?;
            Ok(match weights {
                Some(weights) => {
                    let weight = weights[sample];
                    validate_nonnegative_finite(weight, "weight", sample)?;
                    checked_product2(weight, working_weight).ok_or(DeseqError::NonFiniteValue {
                        context: "IRLS weighted working weight".to_string(),
                        index: Some(sample),
                        value: f64::NAN,
                    })?
                }
                None => working_weight,
            })
        })
        .collect()
}

fn stable_working_weight(mu: f64, dispersion: f64) -> Option<f64> {
    validate_positive_finite(mu, "mu", 0).ok()?;
    validate_positive_finite(dispersion, "dispersion", 0).ok()?;
    let disp_mu = dispersion * mu;
    if disp_mu.is_infinite() && disp_mu.is_sign_positive() {
        return Some(dispersion.recip());
    }
    let denominator = checked_sum2(1.0, disp_mu)?;
    checked_div2(mu, denominator)
}

fn working_response(mu: &[f64], nf: &[f64], y: &[f64]) -> Result<Vec<f64>, DeseqError> {
    mu.iter()
        .copied()
        .zip(nf.iter().copied())
        .zip(y.iter().copied())
        .enumerate()
        .map(|(sample, ((mu, factor), count))| {
            validate_positive_finite(mu, "mu", sample)?;
            validate_positive_finite(factor, "normalization factor", sample)?;
            let log_normalized_mu =
                stable_log_normalized_mean(mu, factor).ok_or(DeseqError::NonFiniteValue {
                    context: "IRLS working response normalized mean".to_string(),
                    index: Some(sample),
                    value: f64::NAN,
                })?;
            let residual = checked_sum2(count, -mu).ok_or(DeseqError::NonFiniteValue {
                context: "IRLS working response residual".to_string(),
                index: Some(sample),
                value: f64::NAN,
            })?;
            let scaled_residual = checked_div2(residual, mu).ok_or(DeseqError::NonFiniteValue {
                context: "IRLS working response scaled residual".to_string(),
                index: Some(sample),
                value: f64::NAN,
            })?;
            checked_sum2(log_normalized_mu, scaled_residual).ok_or(DeseqError::NonFiniteValue {
                context: "IRLS working response".to_string(),
                index: Some(sample),
                value: f64::NAN,
            })
        })
        .collect()
}

fn stable_log_normalized_mean(mu: f64, factor: f64) -> Option<f64> {
    validate_positive_finite(mu, "mu", 0).ok()?;
    validate_positive_finite(factor, "normalization factor", 0).ok()?;
    checked_sum2(mu.ln(), -factor.ln())
}

fn solve_weighted_least_squares(
    x: &DMatrix<f64>,
    w: &[f64],
    z: &[f64],
    ridge_lambda: &[f64],
    solver: IrlsSolver,
) -> Option<DVector<f64>> {
    match solver {
        IrlsSolver::NormalEquations => {
            let (xtwx, xtwz) = xtwx_and_xtwz(x, w, z, ridge_lambda)?;
            xtwx.lu().solve(&xtwz)
        }
        IrlsSolver::Qr => solve_weighted_least_squares_qr(x, w, z, ridge_lambda),
    }
}

fn solve_weighted_least_squares_qr(
    x: &DMatrix<f64>,
    w: &[f64],
    z: &[f64],
    ridge_lambda: &[f64],
) -> Option<DVector<f64>> {
    let n = x.nrows();
    let p = x.ncols();
    let mut augmented_x = DMatrix::zeros(n + p, p);
    let mut augmented_z = DVector::zeros(n + p);
    for row in 0..n {
        let sqrt_weight = w[row].sqrt();
        if !sqrt_weight.is_finite() {
            return None;
        }
        augmented_z[row] = checked_product2(z[row], sqrt_weight)?;
        for col in 0..p {
            augmented_x[(row, col)] = checked_product2(x[(row, col)], sqrt_weight)?;
        }
    }
    for col in 0..p {
        let sqrt_ridge = ridge_lambda[col].sqrt();
        if !sqrt_ridge.is_finite() {
            return None;
        }
        augmented_x[(n + col, col)] = sqrt_ridge;
    }
    let (q, r) = augmented_x.qr().unpack();
    let rhs = q.transpose() * augmented_z;
    let r_econ = r.view((0, 0), (p, p)).into_owned();
    let rhs_econ = rhs.rows(0, p).into_owned();
    r_econ.lu().solve(&rhs_econ)
}

fn covariance_and_hat_diagonal(
    x: &DMatrix<f64>,
    w: &[f64],
    ridge_lambda: &[f64],
) -> Option<(Vec<f64>, Vec<f64>)> {
    let zeros = vec![0.0; x.nrows()];
    let (xtwx_ridge, _) = xtwx_and_xtwz(x, w, &zeros, ridge_lambda)?;
    let xtwx = xtwx_without_ridge(x, w)?;
    let inverse = xtwx_ridge.try_inverse()?;
    let sigma = &inverse * xtwx * &inverse;
    let mut beta_covariance = Vec::with_capacity(x.ncols() * x.ncols());
    for row in 0..x.ncols() {
        for col in 0..x.ncols() {
            let value = sigma[(row, col)];
            if !value.is_finite() {
                return None;
            }
            beta_covariance.push(value);
        }
    }
    let mut hat = Vec::with_capacity(x.nrows());
    for (sample, &weight) in w.iter().enumerate().take(x.nrows()) {
        let mut value = 0.0;
        if !weight.is_finite() {
            return None;
        }
        let active = active_design_columns(x, sample);
        for (left, left_value) in active.iter().copied() {
            for (right, right_value) in active.iter().copied() {
                value = checked_sum2(
                    value,
                    checked_product4(left_value, weight, right_value, inverse[(right, left)])?,
                )?;
            }
        }
        hat.push(value);
    }
    Some((beta_covariance, hat))
}

fn xtwx_and_xtwz(
    x: &DMatrix<f64>,
    w: &[f64],
    z: &[f64],
    ridge_lambda: &[f64],
) -> Option<(DMatrix<f64>, DVector<f64>)> {
    let mut xtwx = DMatrix::zeros(x.ncols(), x.ncols());
    let mut xtwz = DVector::zeros(x.ncols());
    for (sample, &weight) in w.iter().enumerate().take(x.nrows()) {
        let active = active_design_columns(x, sample);
        for (col, col_value) in active.iter().copied() {
            xtwz[col] = checked_sum2(xtwz[col], checked_product3(col_value, weight, z[sample])?)?;
        }
        for (left_pos, (row, row_value)) in active.iter().copied().enumerate() {
            for (col, col_value) in active.iter().copied().skip(left_pos) {
                xtwx[(row, col)] = checked_sum2(
                    xtwx[(row, col)],
                    checked_product3(row_value, weight, col_value)?,
                )?;
                if row != col {
                    xtwx[(col, row)] = xtwx[(row, col)];
                }
            }
        }
    }
    for diagonal in 0..x.ncols() {
        xtwx[(diagonal, diagonal)] =
            checked_sum2(xtwx[(diagonal, diagonal)], ridge_lambda[diagonal])?;
    }
    Some((xtwx, xtwz))
}

fn xtwx_without_ridge(x: &DMatrix<f64>, w: &[f64]) -> Option<DMatrix<f64>> {
    let mut xtwx = DMatrix::zeros(x.ncols(), x.ncols());
    for (sample, &weight) in w.iter().enumerate().take(x.nrows()) {
        let active = active_design_columns(x, sample);
        for (left_pos, (row, row_value)) in active.iter().copied().enumerate() {
            for (col, col_value) in active.iter().copied().skip(left_pos) {
                xtwx[(row, col)] = checked_sum2(
                    xtwx[(row, col)],
                    checked_product3(row_value, weight, col_value)?,
                )?;
                if row != col {
                    xtwx[(col, row)] = xtwx[(row, col)];
                }
            }
        }
    }
    Some(xtwx)
}

fn active_design_columns(x: &DMatrix<f64>, row: usize) -> Vec<(usize, f64)> {
    (0..x.ncols())
        .filter_map(|col| {
            let value = x[(row, col)];
            (value != 0.0).then_some((col, value))
        })
        .collect()
}

fn normalization_factors_from_size_factors(
    counts: &CountMatrix,
    size_factors: &[f64],
) -> Result<RowMajorMatrix<f64>, DeseqError> {
    if size_factors.len() != counts.n_samples() {
        return Err(invalid_dimensions(
            "size factors",
            counts.n_samples(),
            size_factors.len(),
        ));
    }
    let mut values = Vec::with_capacity(counts.n_genes() * counts.n_samples());
    for _ in 0..counts.n_genes() {
        for (idx, factor) in size_factors.iter().copied().enumerate() {
            validate_positive_finite(factor, "size factor", idx)?;
            values.push(factor);
        }
    }
    RowMajorMatrix::from_row_major(counts.n_genes(), counts.n_samples(), values)
}

fn validate_nf_irls_inputs(
    counts: &CountMatrix,
    design: &DesignMatrix,
    normalization_factors: &RowMajorMatrix<f64>,
    dispersions: &[f64],
    weights: Option<&RowMajorMatrix<f64>>,
    options: &IrlsOptions,
) -> Result<(), DeseqError> {
    if design.n_samples() != counts.n_samples() {
        return Err(invalid_dimensions(
            "design rows",
            counts.n_samples(),
            design.n_samples(),
        ));
    }
    if normalization_factors.n_rows() != counts.n_genes()
        || normalization_factors.n_cols() != counts.n_samples()
    {
        return Err(DeseqError::InvalidDimensions {
            context: "normalization factors".to_string(),
            expected: counts.n_genes() * counts.n_samples(),
            actual: normalization_factors.len(),
        });
    }
    if dispersions.len() != counts.n_genes() {
        return Err(invalid_dimensions(
            "dispersions",
            counts.n_genes(),
            dispersions.len(),
        ));
    }
    for (idx, dispersion) in dispersions.iter().copied().enumerate() {
        validate_positive_finite(dispersion, "dispersion", idx)?;
    }
    if let Some(weights) = weights {
        if weights.n_rows() != counts.n_genes() || weights.n_cols() != counts.n_samples() {
            return Err(DeseqError::InvalidDimensions {
                context: "weights".to_string(),
                expected: counts.n_genes() * counts.n_samples(),
                actual: weights.len(),
            });
        }
        for (idx, weight) in weights.as_slice().iter().copied().enumerate() {
            validate_nonnegative_finite(weight, "weight", idx)?;
        }
    }
    if !options.beta_tol.is_finite()
        || options.beta_tol <= 0.0
        || options.maxit == 0
        || options.optim_maxit == 0
        || !options.optim_tol.is_finite()
        || options.optim_tol < 0.0
        || !options.min_mu.is_finite()
        || options.min_mu <= 0.0
        || !options.max_beta_abs.is_finite()
        || options.max_beta_abs <= 0.0
        || !options.ridge_lambda.is_finite()
        || options.ridge_lambda < 0.0
    {
        return Err(DeseqError::UnsupportedFeature {
            feature: "invalid IRLS options".to_string(),
        });
    }
    normalization_factors.validate_finite("normalization factors")?;
    Ok(())
}

fn validate_positive_finite(value: f64, context: &str, index: usize) -> Result<(), DeseqError> {
    if !value.is_finite() || value <= 0.0 {
        return Err(DeseqError::NonFiniteValue {
            context: context.to_string(),
            index: Some(index),
            value,
        });
    }
    Ok(())
}

fn validate_nonnegative_finite(value: f64, context: &str, index: usize) -> Result<(), DeseqError> {
    if !value.is_finite() || value < 0.0 {
        return Err(DeseqError::NonFiniteValue {
            context: context.to_string(),
            index: Some(index),
            value,
        });
    }
    Ok(())
}

fn checked_log2_scale(value: f64, index: usize, context: &str) -> Result<f64, DeseqError> {
    let scaled = std::f64::consts::LOG2_E * value;
    if scaled.is_finite() {
        Ok(scaled)
    } else {
        Err(DeseqError::NonFiniteValue {
            context: context.to_string(),
            index: Some(index),
            value: scaled,
        })
    }
}

fn checked_log2_standard_error(
    covariance: f64,
    index: usize,
    context: &str,
) -> Result<f64, DeseqError> {
    let se = covariance.max(0.0).sqrt();
    checked_log2_scale(se, index, context)
}

fn checked_log2_covariance(value: f64, index: usize, context: &str) -> Result<f64, DeseqError> {
    let log2_e = std::f64::consts::LOG2_E;
    let log2_e_squared = checked_log2_scale(log2_e, index, context)?;
    let scaled = log2_e_squared * value;
    if scaled.is_finite() {
        Ok(scaled)
    } else {
        Err(DeseqError::NonFiniteValue {
            context: context.to_string(),
            index: Some(index),
            value: scaled,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use serde::Deserialize;
    use std::{env, fs, path::PathBuf};

    #[test]
    fn log2_output_scaling_rejects_nonfinite_values() {
        assert!(matches!(
            checked_log2_scale(f64::MAX, 0, "test beta scale"),
            Err(DeseqError::NonFiniteValue { context, index, .. })
                if context == "test beta scale" && index == Some(0)
        ));
        assert!(matches!(
            checked_log2_covariance(f64::MAX, 1, "test covariance scale"),
            Err(DeseqError::NonFiniteValue { context, index, .. })
                if context == "test covariance scale" && index == Some(1)
        ));
    }

    #[test]
    fn optim_start_keeps_finite_beta_on_bounds() {
        let x = DMatrix::from_row_slice(2, 2, &[1.0, 0.0, 0.0, 1.0]);
        let counts = [1_u32, 2_u32];
        let nf = [1.0, 1.0];

        let start = optim_start_beta_log2(&[30.0, -30.0], &x, &counts, &nf, 30.0).unwrap();

        assert_eq!(start, vec![30.0, -30.0]);
    }

    #[test]
    fn optim_start_preserves_backup_qr_natural_log_values() {
        let x = DMatrix::from_element(2, 1, 1.0);
        let counts = [9, 9];
        let nf = [1.0, 1.0];
        let start = optim_start_beta_log2(&[f64::NAN], &x, &counts, &nf, 30.0).unwrap();

        assert_relative_eq!(start[0], 9.1_f64.ln(), epsilon = 1e-14);
    }

    #[test]
    fn optim_fallback_control_matches_r_optim_lbfgsb_defaults() {
        let options = IrlsOptions::default();
        assert_eq!(options.max_beta_abs, 30.0);
        let control = optim_fallback_control(3, 100, 0.0);

        assert_eq!(control.maxit, 100);
        assert_eq!(control.factr, 1.0e7);
        assert_eq!(control.pgtol, 0.0);
        assert_eq!(control.lmm, 5);
        assert_eq!(control.fnscale, 1.0);
        assert_eq!(control.parscale, vec![1.0; 3]);
        assert_eq!(control.ndeps, vec![1.0e-3; 3]);
    }

    #[test]
    fn r_optim_candidate_accepts_objective_gap_within_factr_band() {
        let analytic = 427.391043312651;
        let tolerance = 10.0 * 1.0e7 * f64::EPSILON * analytic;

        assert!(r_optim_candidate_within_analytic_tolerance(
            analytic + tolerance,
            analytic,
            1.0e7,
        ));
    }

    #[test]
    fn r_optim_candidate_rejects_material_objective_gap() {
        let analytic = 45.11460017889003;
        let tolerance = 10.0 * 1.0e7 * f64::EPSILON * analytic;

        assert!(!r_optim_candidate_within_analytic_tolerance(
            analytic + 1.01 * tolerance,
            analytic,
            1.0e7,
        ));
        assert!(!r_optim_candidate_within_analytic_tolerance(
            f64::INFINITY,
            analytic,
            1.0e7,
        ));
    }

    #[test]
    fn r_optim_candidate_disagreement_is_scaled_by_numeric_resolution() {
        let ndeps = [1.0e-3, 1.0e-3];
        let parscale = [1.0, 2.0];

        assert!(!r_optim_candidates_materially_disagree(
            &[1.0, 2.0],
            &[1.009, 2.019],
            &ndeps,
            &parscale,
        ));
        assert!(r_optim_candidates_materially_disagree(
            &[1.0, 2.0],
            &[1.010_001, 2.0],
            &ndeps,
            &parscale,
        ));
    }

    #[test]
    fn beta_objective_uses_log2_mu_and_nb_kernel_without_gamma_constants() {
        let x = DMatrix::from_row_slice(1, 1, &[1.0]);
        let counts = [3_u32];
        let nf = [2.0];
        let ridge_lambda = [0.0];
        let input = BetaOptimInput {
            x: &x,
            counts: &counts,
            nf: &nf,
            dispersion: 0.5,
            weights: None,
            ridge_lambda: &ridge_lambda,
        };

        let (objective, gradient, hessian) =
            beta_log2_objective_gradient_hessian(&input, &[1.0]).unwrap();
        let size = 2.0_f64;
        let mu = 4.0_f64;
        let expected = (3.0 + size) * (size + mu).ln() - size * size.ln() - 3.0 * mu.ln();
        let expected_gradient = std::f64::consts::LN_2 * (mu - 3.0) / (1.0 + 0.5 * mu);
        let expected_hessian =
            std::f64::consts::LN_2.powi(2) * mu * (1.0 + 0.5 * 3.0) / (1.0 + 0.5 * mu).powi(2);

        assert_relative_eq!(objective, expected, epsilon = 1e-12);
        assert_relative_eq!(gradient[0], expected_gradient, epsilon = 1e-12);
        assert_relative_eq!(hessian[(0, 0)], expected_hessian, epsilon = 1e-12);
    }

    #[test]
    fn beta_objective_converts_natural_log_ridge_to_log2_penalty() {
        let x = DMatrix::from_row_slice(1, 1, &[0.0]);
        let counts = [0_u32];
        let nf = [1.0];
        let natural_log_ridge = 2.5;
        let ridge_lambda = [natural_log_ridge];
        let input = BetaOptimInput {
            x: &x,
            counts: &counts,
            nf: &nf,
            dispersion: 0.5,
            weights: None,
            ridge_lambda: &ridge_lambda,
        };

        let beta = [3.0];
        let (objective, gradient, hessian) =
            beta_log2_objective_gradient_hessian(&input, &beta).unwrap();
        let base = beta_log2_nb_objective_term(0, 1.0, 0.5).unwrap();
        let ridge_log2 = natural_log_ridge * std::f64::consts::LN_2.powi(2);

        assert_relative_eq!(
            objective,
            base + 0.5 * ridge_log2 * beta[0] * beta[0],
            epsilon = 1e-12
        );
        assert_relative_eq!(gradient[0], ridge_log2 * beta[0], epsilon = 1e-12);
        assert_relative_eq!(hessian[(0, 0)], ridge_log2, epsilon = 1e-12);
    }

    #[test]
    fn beta_numeric_gradient_uses_r_optim_bounded_difference() {
        let x = DMatrix::from_row_slice(2, 1, &[1.0, 1.0]);
        let counts = [1_u32, 8_u32];
        let nf = [1.0, 1.5];
        let ridge_lambda = [0.25];
        let input = BetaOptimInput {
            x: &x,
            counts: &counts,
            nf: &nf,
            dispersion: 0.3,
            weights: None,
            ridge_lambda: &ridge_lambda,
        };
        let beta = [30.0];
        let gradient = beta_log2_numeric_gradient(&input, &beta, -30.0, 30.0).unwrap();
        let forward_value = beta_log2_objective(&input, &[30.0]).unwrap();
        let backward_value = beta_log2_objective(&input, &[29.999]).unwrap();

        assert_relative_eq!(
            gradient[0],
            (forward_value - backward_value) / 0.001,
            epsilon = 1e-5
        );
    }

    #[test]
    fn optim_fallback_replays_objective_only_r_optim_hard_real_case_when_available() {
        let Some(fixture) = optional_hard_real_r_optim_fixture() else {
            return;
        };
        let contrast = fixture
            .contrasts
            .iter()
            .find(|contrast| contrast.contrast == "heart_blocked_permutation_rep01")
            .unwrap_or_else(|| fixture.contrasts.first().expect("hard fixture contrast"));
        let case = contrast
            .cases
            .iter()
            .find(|case| case.case_kind == "actual_or_rough_optimizer_row")
            .unwrap_or_else(|| contrast.cases.first().expect("hard fixture case"));
        assert_eq!(case.lower, vec![-30.0; contrast.coefficients]);
        assert_eq!(case.upper, vec![30.0; contrast.coefficients]);
        assert_eq!(case.control.maxit, 100);
        assert_eq!(case.control.factr, 1.0e7);
        assert_eq!(case.control.pgtol, 0.0);
        assert_eq!(case.control.lmm, 5);
        assert_eq!(case.control.ndeps, vec![1.0e-3; contrast.coefficients]);

        let design_values = contrast
            .design
            .iter()
            .flat_map(|row| row.iter().copied())
            .collect::<Vec<_>>();
        let x = DMatrix::from_row_slice(contrast.samples, contrast.coefficients, &design_values);
        let counts = case
            .counts
            .iter()
            .copied()
            .map(|value| {
                assert!(value.fract() == 0.0 && value >= 0.0 && value <= f64::from(u32::MAX));
                value as u32
            })
            .collect::<Vec<_>>();
        let ridge_lambda = vec![0.0; contrast.coefficients];
        let input = BetaOptimInput {
            x: &x,
            counts: &counts,
            nf: &contrast.size_factors,
            dispersion: case.dispersion,
            weights: None,
            ridge_lambda: &ridge_lambda,
        };

        let output = minimize_beta_log2_lbfgsb(
            input,
            &case.initial_par,
            case.lower[0],
            case.upper[0],
            case.control.maxit,
            case.control.pgtol,
            false,
        )
        .unwrap();

        assert!(output.converged, "{}", case.fixture);
        assert!(
            output.iterations.abs_diff(case.result.counts.gradient) <= 4,
            "{} optimizer count drift: actual={} expected={}",
            case.fixture,
            output.iterations,
            case.result.counts.gradient
        );
        assert_relative_eq!(output.value, case.result.value, epsilon = 5e-6);
        for (actual, expected) in output.parameters.iter().zip(case.result.par.iter()) {
            assert_relative_eq!(*actual, *expected, epsilon = 2e-3);
        }
    }

    #[test]
    fn bounded_beta_numeric_helpers_reject_nonfinite_accumulation() {
        assert_relative_eq!(
            projected_gradient_norm(&[0.0, 0.0], &[f64::MAX / 2.0, f64::MAX / 2.0], -30.0, 30.0)
                .unwrap(),
            f64::MAX / 2.0 * 2.0_f64.sqrt(),
            epsilon = 1e292
        );
        assert_eq!(
            projected_gradient_norm(&[0.0, 0.0], &[f64::MAX, f64::MAX], -30.0, 30.0),
            None
        );
        assert_eq!(checked_div2(1.0, 0.0), None);
        assert_eq!(checked_div2(f64::NAN, 1.0), None);
    }

    #[test]
    fn irls_deviance_convergence_stat_rejects_nonfinite_arithmetic() {
        assert_eq!(irls_deviance_convergence_stat(f64::MAX, -f64::MAX), None);
        assert_eq!(irls_deviance_convergence_stat(f64::INFINITY, 1.0), None);
        assert_eq!(irls_deviance_convergence_stat(0.0, 0.0), Some(0.0));
    }

    #[test]
    fn weighted_least_squares_helpers_reject_nonfinite_accumulation() {
        let x = DMatrix::from_row_slice(1, 1, &[f64::MAX]);
        assert!(xtwx_and_xtwz(&x, &[2.0], &[1.0], &[0.0]).is_none());
        assert!(xtwx_without_ridge(&x, &[2.0]).is_none());

        let unit_x = DMatrix::from_row_slice(1, 1, &[1.0]);
        assert!(solve_weighted_least_squares_qr(&unit_x, &[4.0], &[f64::MAX], &[0.0]).is_none());
    }

    #[test]
    fn beta_objective_returns_penalty_for_nonfinite_accumulation() {
        let x = DMatrix::from_row_slice(1, 1, &[f64::MAX]);
        let counts = [1_u32];
        let nf = [1.0];
        let ridge_lambda = [0.0];
        let input = BetaOptimInput {
            x: &x,
            counts: &counts,
            nf: &nf,
            dispersion: 0.2,
            weights: None,
            ridge_lambda: &ridge_lambda,
        };

        let (objective, gradient, hessian) =
            beta_log2_objective_gradient_hessian(&input, &[20.0]).unwrap();

        assert_eq!(objective, 1.0e300);
        assert_eq!(gradient, vec![0.0]);
        assert_relative_eq!(hessian[(0, 0)], 1.0, epsilon = 1e-12);
    }

    #[test]
    fn beta_objective_returns_penalty_for_overflowed_ridge_term() {
        let x = DMatrix::from_row_slice(1, 1, &[1.0]);
        let counts = [1_u32];
        let nf = [1.0];
        let ridge_lambda = [f64::MAX];
        let input = BetaOptimInput {
            x: &x,
            counts: &counts,
            nf: &nf,
            dispersion: 0.2,
            weights: None,
            ridge_lambda: &ridge_lambda,
        };

        let (objective, gradient, hessian) =
            beta_log2_objective_gradient_hessian(&input, &[20.0]).unwrap();

        assert_eq!(objective, 1.0e300);
        assert_eq!(gradient, vec![0.0]);
        assert_relative_eq!(hessian[(0, 0)], 1.0, epsilon = 1e-12);
    }

    #[test]
    fn fitted_mu_rejects_overflowed_linear_predictor() {
        let x = DMatrix::from_row_slice(1, 1, &[f64::MAX]);
        let beta = DVector::from_vec(vec![2.0]);
        let nf = [1.0];

        let err = fitted_mu(&x, &beta, &nf, 0.5).unwrap_err();

        match err {
            DeseqError::NonFiniteValue { context, index, .. } => {
                assert_eq!(context, "IRLS linear predictor");
                assert_eq!(index, Some(0));
            }
            other => panic!("unexpected error: {other:?}"),
        }
    }

    #[test]
    fn fitted_mu_saturates_overflowed_exponential_mean() {
        let x = DMatrix::from_row_slice(1, 1, &[1000.0]);
        let beta = DVector::from_vec(vec![1.0]);
        let nf = [1.0];

        let mu = fitted_mu(&x, &beta, &nf, 0.5).unwrap();

        assert_eq!(mu, vec![f64::MAX]);
    }

    #[test]
    fn fitted_mu_avoids_intermediate_exponential_overflow_with_small_factor() {
        let x = DMatrix::from_row_slice(1, 1, &[1000.0]);
        let beta = DVector::from_vec(vec![1.0]);
        let nf = [f64::MIN_POSITIVE];

        let mu = fitted_mu(&x, &beta, &nf, 0.5).unwrap();

        assert!(mu[0].is_finite());
        assert!(mu[0] > 1.0e126);
    }

    #[test]
    fn fitted_mu_saturates_final_overflowed_mean_with_small_factor() {
        let x = DMatrix::from_row_slice(1, 1, &[2000.0]);
        let beta = DVector::from_vec(vec![1.0]);
        let nf = [f64::MIN_POSITIVE];

        let mu = fitted_mu(&x, &beta, &nf, 0.5).unwrap();

        assert_eq!(mu, vec![f64::MAX]);
    }

    #[test]
    fn fitted_mu_log2_rejects_overflowed_linear_predictor() {
        let x = DMatrix::from_row_slice(1, 1, &[f64::MAX]);
        let nf = [1.0];

        let err = fitted_mu_log2_unfloored(&x, &[2.0], &nf).unwrap_err();

        match err {
            DeseqError::NonFiniteValue { context, index, .. } => {
                assert_eq!(context, "optim fallback linear predictor");
                assert_eq!(index, Some(0));
            }
            other => panic!("unexpected error: {other:?}"),
        }
    }

    #[test]
    fn working_weights_reject_nonfinite_arithmetic() {
        assert_eq!(working_weights(&[f64::MAX], 2.0, None).unwrap(), vec![0.5]);

        let err = working_weights(&[f64::MAX], f64::MIN_POSITIVE, Some(&[10.0])).unwrap_err();
        assert!(matches!(
            err,
            DeseqError::NonFiniteValue { context, index, .. }
                if context == "IRLS weighted working weight" && index == Some(0)
        ));
    }

    #[test]
    fn working_response_avoids_normalized_mean_overflow() {
        let z = working_response(&[f64::MAX], &[f64::MIN_POSITIVE], &[1.0]).unwrap();

        assert!(z[0].is_finite());
        assert!(z[0] > 1400.0);
    }

    #[test]
    fn working_response_rejects_nonfinite_arithmetic() {
        let err = working_response(&[f64::MIN_POSITIVE], &[1.0], &[f64::MAX]).unwrap_err();
        assert!(matches!(
            err,
            DeseqError::NonFiniteValue { context, index, .. }
                if context == "IRLS working response scaled residual" && index == Some(0)
        ));
    }

    #[derive(Debug, Deserialize)]
    struct HardROptimFixture {
        contrasts: Vec<HardROptimContrast>,
    }

    #[derive(Debug, Deserialize)]
    struct HardROptimContrast {
        contrast: String,
        samples: usize,
        coefficients: usize,
        design: Vec<Vec<f64>>,
        size_factors: Vec<f64>,
        cases: Vec<HardROptimCase>,
    }

    #[derive(Debug, Deserialize)]
    struct HardROptimCase {
        fixture: String,
        case_kind: String,
        dispersion: f64,
        counts: Vec<f64>,
        initial_par: Vec<f64>,
        lower: Vec<f64>,
        upper: Vec<f64>,
        control: HardROptimControl,
        result: HardROptimResult,
    }

    #[derive(Debug, Deserialize)]
    struct HardROptimControl {
        maxit: usize,
        ndeps: Vec<f64>,
        factr: f64,
        pgtol: f64,
        lmm: usize,
    }

    #[derive(Debug, Deserialize)]
    struct HardROptimResult {
        par: Vec<f64>,
        value: f64,
        counts: HardROptimCounts,
    }

    #[derive(Debug, Deserialize)]
    struct HardROptimCounts {
        gradient: usize,
    }

    fn optional_hard_real_r_optim_fixture() -> Option<HardROptimFixture> {
        let mut candidates = Vec::new();
        if let Some(path) = env::var_os("RSDESEQ2_RCOMPAT_HARD_OPTIM_FIXTURE") {
            candidates.push(PathBuf::from(path));
        }
        candidates.push(
            PathBuf::from(env!("CARGO_MANIFEST_DIR"))
                .join("../../../rcompat-lbfgsb/fixtures/deseq_hard_real_subset/optim_cases.json"),
        );

        for path in candidates {
            if path.exists() {
                let text = fs::read_to_string(&path)
                    .unwrap_or_else(|error| panic!("failed to read {}: {error}", path.display()));
                return Some(serde_json::from_str(&text).unwrap_or_else(|error| {
                    panic!("failed to parse {}: {error}", path.display())
                }));
            }
        }
        eprintln!(
            "skipping hard real L-BFGS-B replay; set RSDESEQ2_RCOMPAT_HARD_OPTIM_FIXTURE or keep the rcompat-lbfgsb fixture bundle available"
        );
        None
    }
}
