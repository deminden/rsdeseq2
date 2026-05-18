use nalgebra::{DMatrix, DVector};

use crate::core::CountMatrix;
use crate::design::DesignMatrix;
use crate::errors::{invalid_dimensions, DeseqError};
use crate::glm::nb::{nbinom_log_likelihood_matrix, nbinom_log_likelihood_weighted};
use crate::glm::NbinomGlmFit;
use crate::matrix::RowMajorMatrix;

/// Placeholder for future IRLS beta fitting.
pub fn fit_irls() -> Result<(), DeseqError> {
    Err(DeseqError::UnsupportedFeature {
        feature: "full beta-prior, contrast-output, and optim-fallback IRLS variants".to_string(),
    })
}

/// Linear solver for the IRLS weighted least-squares update.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum IrlsSolver {
    /// Solve `(X' W X + ridge) beta = X' W z` directly.
    ///
    /// This preserves the initial Rust behavior and is useful for existing
    /// `useQR=FALSE` DESeq2 references.
    #[default]
    NormalEquations,
    /// Solve DESeq2's augmented QR problem `[sqrt(W) X; sqrt(ridge)] beta`.
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
}

impl Default for IrlsOptions {
    fn default() -> Self {
        Self {
            beta_tol: 1e-8,
            maxit: 100,
            min_mu: 0.5,
            max_beta_abs: 30.0,
            ridge_lambda: 1e-6 / std::f64::consts::LN_2.powi(2),
            ridge_lambda_by_coefficient: None,
            solver: IrlsSolver::NormalEquations,
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
        let default_nat_log_lambda = 1e-6 / std::f64::consts::LN_2.powi(2);
        Ok(self
            .ridge_lambdas_for_coefficients(p)?
            .into_iter()
            .all(|value| value <= default_nat_log_lambda))
    }
}

/// Fit an unweighted fixed-dispersion NB GLM by IRLS.
///
/// This implements the standard-design-matrix branch of DESeq2's `fitBeta`
/// loop. The weighted least-squares update can use either the direct normal
/// equations branch or DESeq2's augmented QR branch. Observation weights,
/// contrast output, and optim fallback remain future work.
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
    let mut beta_covariance_values = Vec::with_capacity(counts.n_genes() * p * p);
    let mut mu_values = Vec::with_capacity(counts.n_genes() * counts.n_samples());
    let mut hat_values = Vec::with_capacity(counts.n_genes() * counts.n_samples());
    let mut beta_iter = Vec::with_capacity(counts.n_genes());
    let mut beta_converged = Vec::with_capacity(counts.n_genes());

    for (gene, dispersion) in dispersions.iter().copied().enumerate() {
        if counts.is_all_zero_gene(gene)? {
            return Err(DeseqError::InvalidCounts {
                reason: format!(
                    "gene {gene} is all zero; DESeq2 GLM fitting excludes allZero rows"
                ),
            });
        }
        let y = counts
            .row(gene)?
            .iter()
            .copied()
            .map(f64::from)
            .collect::<Vec<_>>();
        let nf = normalization_factors.row(gene)?;
        let weight_row = weights.map(|matrix| matrix.row(gene)).transpose()?;
        let mut beta = initial_beta(&x, &y, nf)?;
        let mut mu = fitted_mu(&x, &beta, nf, options.min_mu)?;
        let mut dev_old = 0.0;
        let mut dev = 0.0;
        let mut iter = 0_usize;
        let mut converged = false;

        for t in 0..options.maxit {
            iter += 1;
            let w = working_weights(&mu, dispersion, weight_row)?;
            let z = working_response(&mu, nf, &y)?;
            let Some(next_beta) =
                solve_weighted_least_squares(&x, &w, &z, &ridge_lambda, options.solver)
            else {
                iter = options.maxit;
                break;
            };
            beta = next_beta;
            if beta
                .iter()
                .any(|value| !value.is_finite() || value.abs() > options.max_beta_abs)
            {
                iter = options.maxit;
                break;
            }
            mu = fitted_mu(&x, &beta, nf, options.min_mu)?;
            dev = -2.0
                * nbinom_log_likelihood_weighted(counts.row(gene)?, &mu, dispersion, weight_row)?;
            let conv_test = (dev - dev_old).abs() / (dev.abs() + 0.1);
            if !conv_test.is_finite() {
                iter = options.maxit;
                break;
            }
            if t > 0 && conv_test < options.beta_tol {
                converged = true;
                break;
            }
            dev_old = dev;
        }

        let w = working_weights(&mu, dispersion, weight_row)?;
        let Some((beta_covariance, hat_diag)) = covariance_and_hat_diagonal(&x, &w, &ridge_lambda)
        else {
            iter = options.maxit;
            (0..p).for_each(|_| beta_var_values.push(f64::NAN));
            (0..p * p).for_each(|_| beta_covariance_values.push(f64::NAN));
            hat_values.extend(vec![f64::NAN; counts.n_samples()]);
            beta_values.extend(beta.iter().map(|value| std::f64::consts::LOG2_E * value));
            let output_mu = fitted_mu_unfloored(&x, &beta, nf)?;
            mu_values.extend(output_mu.iter().copied());
            beta_iter.push(iter);
            beta_converged.push(false);
            continue;
        };

        beta_values.extend(beta.iter().map(|value| std::f64::consts::LOG2_E * value));
        for diagonal in 0..p {
            let value = beta_covariance[diagonal * p + diagonal];
            beta_var_values.push(std::f64::consts::LOG2_E * value.max(0.0).sqrt());
        }
        beta_covariance_values.extend(
            beta_covariance
                .into_iter()
                .map(|value| std::f64::consts::LOG2_E.powi(2) * value),
        );
        let output_mu = fitted_mu_unfloored(&x, &beta, nf)?;
        mu_values.extend(output_mu.iter().copied());
        hat_values.extend(hat_diag);
        beta_iter.push(iter);
        beta_converged.push(converged && iter < options.maxit);
        let _ = dev;
    }

    let beta = RowMajorMatrix::from_row_major(counts.n_genes(), p, beta_values)?;
    let beta_se = RowMajorMatrix::from_row_major(counts.n_genes(), p, beta_var_values)?;
    let beta_covariance =
        RowMajorMatrix::from_row_major(counts.n_genes(), p * p, beta_covariance_values)?;
    let mu = RowMajorMatrix::from_row_major(counts.n_genes(), counts.n_samples(), mu_values)?;
    let hat_diagonal =
        RowMajorMatrix::from_row_major(counts.n_genes(), counts.n_samples(), hat_values)?;
    let log_like = nbinom_log_likelihood_matrix(counts, &mu, dispersions, weights)?;

    Ok(NbinomGlmFit {
        log_like,
        beta_converged,
        beta,
        beta_se,
        beta_covariance: Some(beta_covariance),
        mu,
        beta_iter,
        model_matrix: design.clone(),
        n_terms: p,
        hat_diagonal,
    })
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

fn fitted_mu_impl(
    x: &DMatrix<f64>,
    beta: &DVector<f64>,
    nf: &[f64],
    min_mu: Option<f64>,
) -> Result<Vec<f64>, DeseqError> {
    let eta = x * beta;
    eta.iter()
        .copied()
        .zip(nf.iter().copied())
        .enumerate()
        .map(|(sample, (eta, factor))| {
            validate_positive_finite(factor, "normalization factor", sample)?;
            let mu = factor * eta.exp();
            if !mu.is_finite() {
                return Err(DeseqError::NonFiniteValue {
                    context: "IRLS fitted mean".to_string(),
                    index: Some(sample),
                    value: mu,
                });
            }
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
            let working_weight = value / (1.0 + dispersion * value);
            Ok(match weights {
                Some(weights) => {
                    let weight = weights[sample];
                    validate_nonnegative_finite(weight, "weight", sample)?;
                    weight * working_weight
                }
                None => working_weight,
            })
        })
        .collect()
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
            Ok((mu / factor).ln() + (count - mu) / mu)
        })
        .collect()
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
            let (xtwx, xtwz) = xtwx_and_xtwz(x, w, z, ridge_lambda);
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
        augmented_z[row] = z[row] * sqrt_weight;
        for col in 0..p {
            augmented_x[(row, col)] = x[(row, col)] * sqrt_weight;
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
    r.lu().solve(&rhs)
}

fn covariance_and_hat_diagonal(
    x: &DMatrix<f64>,
    w: &[f64],
    ridge_lambda: &[f64],
) -> Option<(Vec<f64>, Vec<f64>)> {
    let zeros = vec![0.0; x.nrows()];
    let (xtwx_ridge, _) = xtwx_and_xtwz(x, w, &zeros, ridge_lambda);
    let xtwx = xtwx_without_ridge(x, w);
    let inverse = xtwx_ridge.try_inverse()?;
    let sigma = &inverse * xtwx * &inverse;
    let mut beta_covariance = Vec::with_capacity(x.ncols() * x.ncols());
    for row in 0..x.ncols() {
        for col in 0..x.ncols() {
            beta_covariance.push(sigma[(row, col)]);
        }
    }
    let mut hat = Vec::with_capacity(x.nrows());
    for sample in 0..x.nrows() {
        let mut value = 0.0;
        for left in 0..x.ncols() {
            for right in 0..x.ncols() {
                value += x[(sample, left)]
                    * w[sample].sqrt()
                    * x[(sample, right)]
                    * w[sample].sqrt()
                    * inverse[(right, left)];
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
) -> (DMatrix<f64>, DVector<f64>) {
    let mut xtwx = DMatrix::zeros(x.ncols(), x.ncols());
    let mut xtwz = DVector::zeros(x.ncols());
    for sample in 0..x.nrows() {
        for col in 0..x.ncols() {
            xtwz[col] += x[(sample, col)] * w[sample] * z[sample];
            for row in 0..x.ncols() {
                xtwx[(row, col)] += x[(sample, row)] * w[sample] * x[(sample, col)];
            }
        }
    }
    for diagonal in 0..x.ncols() {
        xtwx[(diagonal, diagonal)] += ridge_lambda[diagonal];
    }
    (xtwx, xtwz)
}

fn xtwx_without_ridge(x: &DMatrix<f64>, w: &[f64]) -> DMatrix<f64> {
    let mut xtwx = DMatrix::zeros(x.ncols(), x.ncols());
    for sample in 0..x.nrows() {
        for col in 0..x.ncols() {
            for row in 0..x.ncols() {
                xtwx[(row, col)] += x[(sample, row)] * w[sample] * x[(sample, col)];
            }
        }
    }
    xtwx
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
