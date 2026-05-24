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

/// Estimate gene-wise dispersions using DESeq2's linear-mu branch shape.
///
/// This implements a clean Rust subset of `estimateDispersionsGeneEst`:
/// base normalized counts are projected through the supplied design matrix,
/// raw means are reconstructed from size factors or gene/sample normalization
/// factors, rough/moments starts are bounded, and each gene's dispersion is
/// optimized on a two-pass log-alpha grid with optional Cox-Reid correction and
/// without priors. General iterative GLM mean refitting remains future work.
pub fn estimate_gene_wise_dispersions_linear_mu(
    input: GeneWiseDispersionInput<'_>,
    options: GeneWiseDispersionOptions,
) -> Result<GeneWiseDispersionOutput, DeseqError> {
    validate_gene_est_inputs(&input, options)?;
    let max_disp = max_dispersion(options, input.counts.n_samples());
    let normalized_mu = linear_model_mu(input.normalized_counts, input.design)?;
    let rough_disp = rough_dispersion_estimates(input.normalized_counts, input.design)?;
    let moments_disp = match input.normalization_factors {
        Some(normalization_factors) => moments_dispersion_estimates_with_normalization_factors(
            input.base_mean,
            input.base_var,
            normalization_factors,
            Some(input.all_zero),
        )?,
        None => moments_dispersion_estimates(input.base_mean, input.base_var, input.size_factors)?,
    };
    let initial_disp =
        initial_dispersion_estimates(&rough_disp, &moments_disp, options.min_disp, max_disp)?;

    let mut mu_values = vec![f64::NAN; input.counts.n_genes() * input.counts.n_samples()];
    let mut disp_gene_est = vec![f64::NAN; input.counts.n_genes()];
    let mut disp_iter = vec![0; input.counts.n_genes()];
    let mut converged = vec![false; input.counts.n_genes()];

    for gene in 0..input.counts.n_genes() {
        if input.all_zero[gene] {
            continue;
        }
        let mu_start = gene * input.counts.n_samples();
        let normalization_factor_row = input
            .normalization_factors
            .map(|normalization_factors| normalization_factors.row(gene))
            .transpose()?;
        for sample in 0..input.counts.n_samples() {
            let factor = normalization_factor_row
                .map(|row| row[sample])
                .unwrap_or(input.size_factors[sample]);
            let value = normalized_mu.row(gene)?[sample] * factor;
            mu_values[mu_start + sample] = value.max(options.min_mu);
        }
        let row_mu = &mu_values[mu_start..mu_start + input.counts.n_samples()];
        let (estimate, iterations, is_converged) = fit_dispersion_for_gene(
            input.counts.row(gene)?,
            row_mu,
            input.design,
            initial_disp[gene],
            options,
            input.counts.n_samples(),
        )?;
        disp_gene_est[gene] = estimate.clamp(options.min_disp, max_disp);
        disp_iter[gene] = iterations;
        converged[gene] = is_converged;
    }

    Ok(GeneWiseDispersionOutput {
        disp_gene_est,
        disp_iter,
        rough_disp,
        moments_disp,
        initial_disp,
        mu: RowMajorMatrix::from_row_major(
            input.counts.n_genes(),
            input.counts.n_samples(),
            mu_values,
        )?,
        converged,
    })
}

/// Estimate gene-wise dispersions using one or more GLM mean-refit iterations.
///
/// This follows the non-`linearMu` branch shape of DESeq2's
/// `estimateDispersionsGeneEst`: rough/moments estimates initialize
/// `alpha_hat`, non-all-zero rows alternate between fixed-dispersion NB GLM
/// mean fitting and fixed-mean dispersion optimization, and rows stop
/// refitting when the log-dispersion move is at most `0.05`. When
/// row-normalized observation weights are supplied, they are used in the
/// fixed-dispersion IRLS mean fit and the fixed-mean likelihood objective;
/// Cox-Reid terms use DESeq2's thresholded weighted sample subset.
/// glmGamPoi fitting remains a future high-level branch.
pub fn estimate_gene_wise_dispersions_glm_mu(
    input: GeneWiseDispersionInput<'_>,
    options: GeneWiseDispersionOptions,
    irls_options: IrlsOptions,
) -> Result<GeneWiseDispersionOutput, DeseqError> {
    validate_gene_est_inputs(&input, options)?;
    let max_disp = max_dispersion(options, input.counts.n_samples());
    let rough_disp = rough_dispersion_estimates(input.normalized_counts, input.design)?;
    let moments_disp = match input.normalization_factors {
        Some(normalization_factors) => moments_dispersion_estimates_with_normalization_factors(
            input.base_mean,
            input.base_var,
            normalization_factors,
            Some(input.all_zero),
        )?,
        None => moments_dispersion_estimates(input.base_mean, input.base_var, input.size_factors)?,
    };
    let initial_disp =
        initial_dispersion_estimates(&rough_disp, &moments_disp, options.min_disp, max_disp)?;

    let mut alpha_hat = initial_disp.clone();
    let mut alpha_hat_new = initial_disp.clone();
    let alpha_init = initial_disp.clone();
    let fitting_gene_order = input
        .all_zero
        .iter()
        .copied()
        .enumerate()
        .filter_map(|(gene, all_zero)| (!all_zero).then_some(gene))
        .collect::<Vec<_>>();
    let mut fitidx = input
        .all_zero
        .iter()
        .map(|all_zero| !all_zero)
        .collect::<Vec<_>>();
    let mut mu_values = vec![f64::NAN; input.counts.n_genes() * input.counts.n_samples()];
    let mut disp_iter = vec![0; input.counts.n_genes()];
    let mut initial_lp = vec![f64::NAN; input.counts.n_genes()];
    let mut last_lp = vec![f64::NAN; input.counts.n_genes()];

    let mut mean_options = irls_options;
    mean_options.min_mu = options.min_mu;

    for _ in 0..options.niter {
        let fit_genes = fitidx
            .iter()
            .copied()
            .enumerate()
            .filter_map(|(gene, should_fit)| should_fit.then_some(gene))
            .collect::<Vec<_>>();
        if fit_genes.is_empty() {
            break;
        }

        let compact_counts = compact_counts_rows(input.counts, &fit_genes)?;
        let compact_disp = compact_gene_values(&alpha_hat, &fit_genes)?;
        let compact_weights = input
            .observation_weights
            .map(|weights| compact_matrix_rows(weights, &fit_genes))
            .transpose()?;
        let fit = match input.normalization_factors {
            Some(normalization_factors) => {
                let compact_factors = compact_matrix_rows(normalization_factors, &fit_genes)?;
                fit_fixed_dispersion_irls_with_normalization_factors_and_weights(
                    &compact_counts,
                    input.design,
                    &compact_factors,
                    &compact_disp,
                    compact_weights.as_ref(),
                    mean_options.clone(),
                )?
            }
            None => fit_fixed_dispersion_irls_with_weights(
                &compact_counts,
                input.design,
                input.size_factors,
                &compact_disp,
                compact_weights.as_ref(),
                mean_options.clone(),
            )?,
        };

        for (compact_row, gene) in fit_genes.iter().copied().enumerate() {
            let fit_mu_raw = fit.mu.row(compact_row)?;
            let fit_mu = fit_mu_raw
                .iter()
                .copied()
                .map(|value| value.max(options.min_mu))
                .collect::<Vec<_>>();
            let start = gene * input.counts.n_samples();
            mu_values[start..start + input.counts.n_samples()].copy_from_slice(&fit_mu);
            // DESeq2 passes the full non-all-zero weight matrix into fitDisp
            // even when counts/mu are subset by fitidx; the C++ then indexes
            // weights by compact row position.
            let weight_row = input
                .observation_weights
                .map(|weights| weights.row(fitting_gene_order[compact_row]))
                .transpose()?;
            let diagnostics = fit_dispersion_for_gene_detailed_with_weights(
                input.counts.row(gene)?,
                &fit_mu,
                input.design,
                alpha_hat[gene],
                options,
                input.counts.n_samples(),
                weight_row,
            )?;
            alpha_hat_new[gene] = diagnostics.estimate.min(max_disp);
            disp_iter[gene] = diagnostics.iterations;
            initial_lp[gene] = diagnostics.initial_lp;
            last_lp[gene] = diagnostics.last_lp;
        }

        fitidx = input
            .all_zero
            .iter()
            .copied()
            .enumerate()
            .map(|(gene, all_zero)| {
                if all_zero {
                    return false;
                }
                let move_size = (alpha_hat_new[gene].ln() - alpha_hat[gene].ln()).abs();
                move_size.is_finite() && move_size > 0.05
            })
            .collect();
        alpha_hat.clone_from(&alpha_hat_new);
        if !fitidx.iter().any(|should_fit| *should_fit) {
            break;
        }
    }

    let mut disp_gene_est = alpha_hat;
    if options.niter == 1 {
        for gene in 0..input.counts.n_genes() {
            if input.all_zero[gene] || !initial_lp[gene].is_finite() || !last_lp[gene].is_finite() {
                continue;
            }
            if last_lp[gene] < initial_lp[gene] + initial_lp[gene].abs() / 1.0e6 {
                disp_gene_est[gene] = alpha_init[gene];
            }
        }
    }

    let mut converged = vec![false; input.counts.n_genes()];
    for gene in 0..input.counts.n_genes() {
        if input.all_zero[gene] {
            disp_gene_est[gene] = f64::NAN;
            continue;
        }
        converged[gene] = disp_iter[gene] < options.maxit && disp_iter[gene] != 1;
        if !converged[gene] && disp_gene_est[gene] > options.min_disp * 10.0 {
            let mu = &mu_values[gene * input.counts.n_samples()
                ..gene * input.counts.n_samples() + input.counts.n_samples()];
            let weight_row = input
                .observation_weights
                .map(|weights| weights.row(gene))
                .transpose()?;
            disp_gene_est[gene] = fit_dispersion_grid_inner(DispersionOptimizerInput {
                counts: input.counts.row(gene)?,
                mu,
                design: Some(input.design),
                initial_dispersion: disp_gene_est[gene],
                options,
                n_samples: input.counts.n_samples(),
                prior: None,
                weights: weight_row,
            })?
            .0;
        }
        disp_gene_est[gene] = disp_gene_est[gene].clamp(options.min_disp, max_disp);
    }

    Ok(GeneWiseDispersionOutput {
        disp_gene_est,
        disp_iter,
        rough_disp,
        moments_disp,
        initial_disp,
        mu: RowMajorMatrix::from_row_major(
            input.counts.n_genes(),
            input.counts.n_samples(),
            mu_values,
        )?,
        converged,
    })
}

/// Project normalized counts onto the supplied design matrix.
///
/// This is the Rust analogue of DESeq2's `linearModelMu` helper for row-wise
/// fitted values, using `Y X (X'X)^-1 X'`.
pub fn linear_model_mu(
    normalized_counts: &RowMajorMatrix<f64>,
    design: &DesignMatrix,
) -> Result<RowMajorMatrix<f64>, DeseqError> {
    if normalized_counts.n_cols() != design.n_samples() {
        return Err(invalid_dimensions(
            "linear mu samples",
            design.n_samples(),
            normalized_counts.n_cols(),
        ));
    }
    let y = DMatrix::from_row_slice(
        normalized_counts.n_rows(),
        normalized_counts.n_cols(),
        normalized_counts.as_slice(),
    );
    let x = DMatrix::from_row_slice(
        design.n_samples(),
        design.n_coefficients(),
        design.matrix().as_slice(),
    );
    let xtx = x.transpose() * &x;
    let Some(xtx_inverse) = xtx.try_inverse() else {
        return Err(DeseqError::InvalidDimensions {
            context: "linear mu design rank".to_string(),
            expected: design.n_coefficients(),
            actual: 0,
        });
    };
    let hat = &x * xtx_inverse * x.transpose();
    let fitted = y * hat;
    let mut values = Vec::with_capacity(normalized_counts.n_rows() * normalized_counts.n_cols());
    for row in 0..normalized_counts.n_rows() {
        for col in 0..normalized_counts.n_cols() {
            values.push(fitted[(row, col)]);
        }
    }
    RowMajorMatrix::from_row_major(
        normalized_counts.n_rows(),
        normalized_counts.n_cols(),
        values,
    )
}

/// DESeq2-style rough dispersion estimates from normalized counts.
pub fn rough_dispersion_estimates(
    normalized_counts: &RowMajorMatrix<f64>,
    design: &DesignMatrix,
) -> Result<Vec<f64>, DeseqError> {
    if design.n_samples() <= design.n_coefficients() {
        return Err(DeseqError::InvalidDimensions {
            context: "dispersion residual degrees of freedom".to_string(),
            expected: design.n_coefficients() + 1,
            actual: design.n_samples(),
        });
    }
    let mu = linear_model_mu(normalized_counts, design)?;
    let residual_df = (design.n_samples() - design.n_coefficients()) as f64;
    let mut estimates = Vec::with_capacity(normalized_counts.n_rows());
    for gene in 0..normalized_counts.n_rows() {
        let y = normalized_counts.row(gene)?;
        let mu = mu.row(gene)?;
        let sum = y
            .iter()
            .copied()
            .zip(mu.iter().copied())
            .map(|(count, fitted)| {
                let fitted = fitted.max(1.0);
                ((count - fitted).powi(2) - fitted) / fitted.powi(2)
            })
            .sum::<f64>();
        estimates.push((sum / residual_df).max(0.0));
    }
    Ok(estimates)
}

/// DESeq2-style moments dispersion estimates.
pub fn moments_dispersion_estimates(
    base_mean: &[f64],
    base_var: &[f64],
    size_factors: &[f64],
) -> Result<Vec<f64>, DeseqError> {
    if base_mean.len() != base_var.len() {
        return Err(invalid_dimensions(
            "moments dispersion base statistics",
            base_mean.len(),
            base_var.len(),
        ));
    }
    validate_size_factors(size_factors)?;
    let xim =
        size_factors.iter().map(|value| value.recip()).sum::<f64>() / size_factors.len() as f64;
    Ok(base_mean
        .iter()
        .copied()
        .zip(base_var.iter().copied())
        .map(|(mean, variance)| {
            if mean > 0.0 {
                (variance - xim * mean) / mean.powi(2)
            } else {
                f64::NAN
            }
        })
        .collect())
}

/// DESeq2-style moments dispersion estimates with gene/sample normalization factors.
///
/// This follows `momentsDispEstimate`: when normalization factors are present,
/// `xim = mean(1 / colMeans(normalizationFactors))`. If `all_zero` is supplied,
/// all-zero rows are excluded from the column means, matching the fact that
/// DESeq2 calls the helper on `objectNZ`.
pub fn moments_dispersion_estimates_with_normalization_factors(
    base_mean: &[f64],
    base_var: &[f64],
    normalization_factors: &RowMajorMatrix<f64>,
    all_zero: Option<&[bool]>,
) -> Result<Vec<f64>, DeseqError> {
    if base_mean.len() != base_var.len() {
        return Err(invalid_dimensions(
            "moments dispersion base statistics",
            base_mean.len(),
            base_var.len(),
        ));
    }
    if normalization_factors.n_rows() != base_mean.len() {
        return Err(invalid_dimensions(
            "moments dispersion normalization-factor rows",
            base_mean.len(),
            normalization_factors.n_rows(),
        ));
    }
    if let Some(all_zero) = all_zero {
        if all_zero.len() != base_mean.len() {
            return Err(invalid_dimensions(
                "moments dispersion allZero",
                base_mean.len(),
                all_zero.len(),
            ));
        }
    }
    let xim = normalization_factor_moments_xim(normalization_factors, all_zero)?;
    moments_dispersion_estimates_with_xim(base_mean, base_var, xim)
}

fn moments_dispersion_estimates_with_xim(
    base_mean: &[f64],
    base_var: &[f64],
    xim: f64,
) -> Result<Vec<f64>, DeseqError> {
    if base_mean.len() != base_var.len() {
        return Err(invalid_dimensions(
            "moments dispersion base statistics",
            base_mean.len(),
            base_var.len(),
        ));
    }
    if !xim.is_finite() || xim <= 0.0 {
        return Err(DeseqError::InvalidSizeFactors {
            reason: "moments dispersion normalization factor summary must be finite and positive"
                .to_string(),
        });
    }
    Ok(base_mean
        .iter()
        .copied()
        .zip(base_var.iter().copied())
        .map(|(mean, variance)| {
            if mean > 0.0 {
                (variance - xim * mean) / mean.powi(2)
            } else {
                f64::NAN
            }
        })
        .collect())
}

fn normalization_factor_moments_xim(
    normalization_factors: &RowMajorMatrix<f64>,
    all_zero: Option<&[bool]>,
) -> Result<f64, DeseqError> {
    let mut col_sums = vec![0.0; normalization_factors.n_cols()];
    let mut n_rows_used = 0_usize;
    for row in 0..normalization_factors.n_rows() {
        if all_zero.is_some_and(|flags| flags[row]) {
            continue;
        }
        for (sample, value) in normalization_factors.row(row)?.iter().copied().enumerate() {
            validate_normalization_factor(value, sample)?;
            col_sums[sample] += value;
        }
        n_rows_used += 1;
    }
    if n_rows_used == 0 {
        return Err(DeseqError::InvalidCounts {
            reason: "no non-all-zero rows available for normalization-factor moments estimate"
                .to_string(),
        });
    }
    let inverse_col_mean_sum = col_sums
        .iter()
        .copied()
        .enumerate()
        .map(|(sample, sum)| {
            let col_mean = sum / n_rows_used as f64;
            if !col_mean.is_finite() || col_mean <= 0.0 {
                Err(DeseqError::InvalidSizeFactors {
                    reason: format!(
                        "normalization-factor column mean at sample {sample} must be finite and positive"
                    ),
                })
            } else {
                Ok(col_mean.recip())
            }
        })
        .collect::<Result<Vec<_>, DeseqError>>()?
        .into_iter()
        .sum::<f64>();
    Ok(inverse_col_mean_sum / normalization_factors.n_cols() as f64)
}

/// Combine rough and moments estimates using DESeq2's bounded start shape.
pub fn initial_dispersion_estimates(
    rough_disp: &[f64],
    moments_disp: &[f64],
    min_disp: f64,
    max_disp: f64,
) -> Result<Vec<f64>, DeseqError> {
    if rough_disp.len() != moments_disp.len() {
        return Err(invalid_dimensions(
            "initial dispersion starts",
            rough_disp.len(),
            moments_disp.len(),
        ));
    }
    validate_dispersion_bounds(min_disp, max_disp)?;
    Ok(rough_disp
        .iter()
        .copied()
        .zip(moments_disp.iter().copied())
        .map(|(rough, moments)| {
            if !rough.is_finite() || !moments.is_finite() {
                f64::NAN
            } else {
                rough.min(moments).clamp(min_disp, max_disp)
            }
        })
        .collect())
}

fn fit_dispersion_for_gene(
    counts: &[u32],
    mu: &[f64],
    design: &DesignMatrix,
    initial_dispersion: f64,
    options: GeneWiseDispersionOptions,
    n_samples: usize,
) -> Result<(f64, usize, bool), DeseqError> {
    let diagnostics = fit_dispersion_for_gene_detailed(
        counts,
        mu,
        design,
        initial_dispersion,
        options,
        n_samples,
    )?;
    Ok((
        diagnostics.estimate,
        diagnostics.iterations,
        diagnostics.converged,
    ))
}

fn fit_dispersion_for_gene_detailed(
    counts: &[u32],
    mu: &[f64],
    design: &DesignMatrix,
    initial_dispersion: f64,
    options: GeneWiseDispersionOptions,
    n_samples: usize,
) -> Result<GeneDispersionFitDiagnostics, DeseqError> {
    fit_dispersion_for_gene_detailed_with_weights(
        counts,
        mu,
        design,
        initial_dispersion,
        options,
        n_samples,
        None,
    )
}

fn fit_dispersion_for_gene_detailed_with_weights(
    counts: &[u32],
    mu: &[f64],
    design: &DesignMatrix,
    initial_dispersion: f64,
    options: GeneWiseDispersionOptions,
    n_samples: usize,
    weights: Option<&[f64]>,
) -> Result<GeneDispersionFitDiagnostics, DeseqError> {
    match options.fit_method {
        GeneWiseDispersionFitMethod::Grid => {
            let (dispersion, evaluations) = fit_dispersion_grid_inner(DispersionOptimizerInput {
                counts,
                mu,
                design: Some(design),
                initial_dispersion,
                options,
                n_samples,
                prior: None,
                weights,
            })?;
            let last_lp = dispersion_log_posterior_objective(
                DispersionObjectiveInput {
                    counts,
                    mu,
                    design: Some(design),
                    use_cox_reid: options.use_cox_reid,
                    prior: None,
                    weights,
                    weight_threshold: options.weight_threshold,
                },
                dispersion.ln(),
            )?;
            Ok(GeneDispersionFitDiagnostics {
                estimate: dispersion,
                iterations: evaluations,
                converged: true,
                initial_lp: last_lp,
                last_lp,
            })
        }
        GeneWiseDispersionFitMethod::LineSearch => {
            let line_search = fit_dispersion_line_search_inner(DispersionOptimizerInput {
                counts,
                mu,
                design: Some(design),
                initial_dispersion,
                options,
                n_samples,
                prior: None,
                weights,
            })?;
            let mut dispersion = line_search.dispersion;
            if !line_search.converged && dispersion > options.min_disp * 10.0 {
                dispersion = fit_dispersion_grid_inner(DispersionOptimizerInput {
                    counts,
                    mu,
                    design: Some(design),
                    initial_dispersion: dispersion,
                    options,
                    n_samples,
                    prior: None,
                    weights,
                })?
                .0;
            }
            Ok(GeneDispersionFitDiagnostics {
                estimate: dispersion,
                iterations: line_search.iter,
                converged: line_search.converged,
                initial_lp: line_search.initial_lp,
                last_lp: line_search.last_lp,
            })
        }
    }
}

/// Fit one dispersion by DESeq2's Armijo line-search shape.
pub fn fit_dispersion_line_search(
    counts: &[u32],
    mu: &[f64],
    design: &DesignMatrix,
    initial_dispersion: f64,
    options: GeneWiseDispersionOptions,
    n_samples: usize,
) -> Result<DispersionLineSearchOutput, DeseqError> {
    fit_dispersion_line_search_inner(DispersionOptimizerInput {
        counts,
        mu,
        design: Some(design),
        initial_dispersion,
        options,
        n_samples,
        prior: None,
        weights: None,
    })
}

/// Fit one dispersion by DESeq2's Armijo line-search shape with a log-alpha prior.
pub fn fit_dispersion_line_search_with_prior(
    counts: &[u32],
    mu: &[f64],
    design: &DesignMatrix,
    initial_dispersion: f64,
    options: GeneWiseDispersionOptions,
    n_samples: usize,
    prior: DispersionPrior,
) -> Result<DispersionLineSearchOutput, DeseqError> {
    fit_dispersion_line_search_inner(DispersionOptimizerInput {
        counts,
        mu,
        design: Some(design),
        initial_dispersion,
        options,
        n_samples,
        prior: Some(prior),
        weights: None,
    })
}

/// Fit one dispersion by DESeq2's weighted Armijo line-search shape with a log-alpha prior.
pub fn fit_dispersion_line_search_with_prior_and_weights(
    input: WeightedDispersionFitInput<'_>,
) -> Result<DispersionLineSearchOutput, DeseqError> {
    fit_dispersion_line_search_inner(DispersionOptimizerInput {
        counts: input.counts,
        mu: input.mu,
        design: Some(input.design),
        initial_dispersion: input.initial_dispersion,
        options: input.options,
        n_samples: input.n_samples,
        prior: Some(input.prior),
        weights: Some(input.weights),
    })
}

/// Fit one dispersion by line search without Cox-Reid correction.
pub fn fit_dispersion_line_search_no_cr(
    counts: &[u32],
    mu: &[f64],
    initial_dispersion: f64,
    options: GeneWiseDispersionOptions,
    n_samples: usize,
) -> Result<DispersionLineSearchOutput, DeseqError> {
    fit_dispersion_line_search_inner(DispersionOptimizerInput {
        counts,
        mu,
        design: None,
        initial_dispersion,
        options,
        n_samples,
        prior: None,
        weights: None,
    })
}

/// Fit one dispersion by line search with a log-alpha prior and without Cox-Reid correction.
pub fn fit_dispersion_line_search_no_cr_with_prior(
    counts: &[u32],
    mu: &[f64],
    initial_dispersion: f64,
    options: GeneWiseDispersionOptions,
    n_samples: usize,
    prior: DispersionPrior,
) -> Result<DispersionLineSearchOutput, DeseqError> {
    fit_dispersion_line_search_inner(DispersionOptimizerInput {
        counts,
        mu,
        design: None,
        initial_dispersion,
        options,
        n_samples,
        prior: Some(prior),
        weights: None,
    })
}

/// Fit one dispersion by weighted line search with a log-alpha prior and without Cox-Reid correction.
pub fn fit_dispersion_line_search_no_cr_with_prior_and_weights(
    counts: &[u32],
    mu: &[f64],
    initial_dispersion: f64,
    options: GeneWiseDispersionOptions,
    n_samples: usize,
    prior: DispersionPrior,
    weights: &[f64],
) -> Result<DispersionLineSearchOutput, DeseqError> {
    fit_dispersion_line_search_inner(DispersionOptimizerInput {
        counts,
        mu,
        design: None,
        initial_dispersion,
        options,
        n_samples,
        prior: Some(prior),
        weights: Some(weights),
    })
}

#[derive(Clone, Copy)]
struct DispersionOptimizerInput<'a> {
    counts: &'a [u32],
    mu: &'a [f64],
    design: Option<&'a DesignMatrix>,
    initial_dispersion: f64,
    options: GeneWiseDispersionOptions,
    n_samples: usize,
    prior: Option<DispersionPrior>,
    weights: Option<&'a [f64]>,
}

#[derive(Clone, Copy)]
struct DispersionObjectiveInput<'a> {
    counts: &'a [u32],
    mu: &'a [f64],
    design: Option<&'a DesignMatrix>,
    use_cox_reid: bool,
    prior: Option<DispersionPrior>,
    weights: Option<&'a [f64]>,
    weight_threshold: f64,
}

fn fit_dispersion_line_search_inner(
    input: DispersionOptimizerInput<'_>,
) -> Result<DispersionLineSearchOutput, DeseqError> {
    let counts = input.counts;
    let mu = input.mu;
    let design = input.design;
    let initial_dispersion = input.initial_dispersion;
    let mut options = input.options;
    let n_samples = input.n_samples;
    let prior = input.prior;
    let weights = input.weights;
    if counts.len() != mu.len() {
        return Err(invalid_dimensions(
            "dispersion line-search mu",
            counts.len(),
            mu.len(),
        ));
    }
    if design.is_none() {
        options.use_cox_reid = false;
    }
    validate_gene_est_options(options)?;
    validate_dispersion_prior(prior)?;
    let max_disp = max_dispersion(options, n_samples);
    validate_dispersion_bounds(options.min_disp, max_disp)?;
    if !initial_dispersion.is_finite() || initial_dispersion <= 0.0 {
        return Err(DeseqError::InvalidDispersion {
            reason: "initial dispersion must be finite and positive".to_string(),
        });
    }

    let min_log_alpha = (options.min_disp / 10.0).ln();
    let mut log_alpha = initial_dispersion.clamp(options.min_disp, max_disp).ln();
    let objective = DispersionObjectiveInput {
        counts,
        mu,
        design,
        use_cox_reid: options.use_cox_reid,
        prior,
        weights,
        weight_threshold: options.weight_threshold,
    };
    let mut lp = dispersion_log_posterior_objective(objective, log_alpha)?;
    let mut dlp = dispersion_log_posterior_derivative_objective(objective, log_alpha)?;
    let initial_lp = lp;
    let initial_dlp = dlp;
    let mut kappa = options.kappa_0;
    let mut iter = 0_usize;
    let mut iter_accept = 0_usize;
    let mut last_change = -1.0;
    let epsilon = 1.0e-4;

    for _ in 0..options.maxit {
        iter += 1;
        if !dlp.is_finite() || dlp.abs() <= f64::EPSILON || !kappa.is_finite() || kappa <= 0.0 {
            break;
        }

        let proposal = log_alpha + kappa * dlp;
        if proposal < -30.0 {
            kappa = (-30.0 - log_alpha) / dlp;
        }
        if proposal > 10.0 {
            kappa = (10.0 - log_alpha) / dlp;
        }
        if !kappa.is_finite() || kappa <= 0.0 {
            break;
        }

        let proposed_log_alpha = log_alpha + kappa * dlp;
        let theta_kappa = -dispersion_log_posterior_objective(objective, proposed_log_alpha)?;
        let theta_hat_kappa = -lp - kappa * epsilon * dlp.powi(2);
        if theta_kappa <= theta_hat_kappa {
            iter_accept += 1;
            log_alpha = proposed_log_alpha;
            let lp_new = dispersion_log_posterior_objective(objective, log_alpha)?;
            last_change = lp_new - lp;
            lp = lp_new;
            if last_change < options.disp_tol {
                break;
            }
            if log_alpha < min_log_alpha {
                break;
            }
            dlp = dispersion_log_posterior_derivative_objective(objective, log_alpha)?;
            kappa = (kappa * 1.1).min(options.kappa_0);
            if iter_accept % 5 == 0 {
                kappa /= 2.0;
            }
        } else {
            kappa /= 2.0;
        }
    }

    let dispersion = log_alpha.exp().clamp(options.min_disp, max_disp);
    let last_dlp = dispersion_log_posterior_derivative_objective(objective, log_alpha)?;
    let last_d2lp = dispersion_log_posterior_second_derivative_objective(objective, log_alpha)?;
    Ok(DispersionLineSearchOutput {
        dispersion,
        log_alpha,
        iter,
        iter_accept,
        initial_lp,
        initial_dlp,
        last_lp: lp,
        last_dlp,
        last_d2lp,
        last_change,
        converged: iter < options.maxit && iter != 1,
    })
}

/// Fit a dispersion for one gene by DESeq2-style two-pass log-alpha grid search.
pub fn fit_dispersion_grid(
    counts: &[u32],
    mu: &[f64],
    design: &DesignMatrix,
    initial_dispersion: f64,
    options: GeneWiseDispersionOptions,
    n_samples: usize,
) -> Result<(f64, usize), DeseqError> {
    fit_dispersion_grid_inner(DispersionOptimizerInput {
        counts,
        mu,
        design: Some(design),
        initial_dispersion,
        options,
        n_samples,
        prior: None,
        weights: None,
    })
}

/// Fit a dispersion by two-pass log-alpha grid search with a log-alpha prior.
pub fn fit_dispersion_grid_with_prior(
    counts: &[u32],
    mu: &[f64],
    design: &DesignMatrix,
    initial_dispersion: f64,
    options: GeneWiseDispersionOptions,
    n_samples: usize,
    prior: DispersionPrior,
) -> Result<(f64, usize), DeseqError> {
    fit_dispersion_grid_inner(DispersionOptimizerInput {
        counts,
        mu,
        design: Some(design),
        initial_dispersion,
        options,
        n_samples,
        prior: Some(prior),
        weights: None,
    })
}

/// Fit a dispersion by weighted two-pass log-alpha grid search with a log-alpha prior.
pub fn fit_dispersion_grid_with_prior_and_weights(
    input: WeightedDispersionFitInput<'_>,
) -> Result<(f64, usize), DeseqError> {
    fit_dispersion_grid_inner(DispersionOptimizerInput {
        counts: input.counts,
        mu: input.mu,
        design: Some(input.design),
        initial_dispersion: input.initial_dispersion,
        options: input.options,
        n_samples: input.n_samples,
        prior: Some(input.prior),
        weights: Some(input.weights),
    })
}

/// Fit a dispersion for one gene without Cox-Reid correction.
pub fn fit_dispersion_grid_no_cr(
    counts: &[u32],
    mu: &[f64],
    initial_dispersion: f64,
    options: GeneWiseDispersionOptions,
    n_samples: usize,
) -> Result<(f64, usize), DeseqError> {
    fit_dispersion_grid_inner(DispersionOptimizerInput {
        counts,
        mu,
        design: None,
        initial_dispersion,
        options,
        n_samples,
        prior: None,
        weights: None,
    })
}

/// Fit a dispersion by two-pass log-alpha grid search with a prior and without Cox-Reid correction.
pub fn fit_dispersion_grid_no_cr_with_prior(
    counts: &[u32],
    mu: &[f64],
    initial_dispersion: f64,
    options: GeneWiseDispersionOptions,
    n_samples: usize,
    prior: DispersionPrior,
) -> Result<(f64, usize), DeseqError> {
    fit_dispersion_grid_inner(DispersionOptimizerInput {
        counts,
        mu,
        design: None,
        initial_dispersion,
        options,
        n_samples,
        prior: Some(prior),
        weights: None,
    })
}

/// Fit a dispersion by weighted two-pass log-alpha grid search with a prior and without Cox-Reid correction.
pub fn fit_dispersion_grid_no_cr_with_prior_and_weights(
    counts: &[u32],
    mu: &[f64],
    initial_dispersion: f64,
    options: GeneWiseDispersionOptions,
    n_samples: usize,
    prior: DispersionPrior,
    weights: &[f64],
) -> Result<(f64, usize), DeseqError> {
    fit_dispersion_grid_inner(DispersionOptimizerInput {
        counts,
        mu,
        design: None,
        initial_dispersion,
        options,
        n_samples,
        prior: Some(prior),
        weights: Some(weights),
    })
}

fn fit_dispersion_grid_inner(
    input: DispersionOptimizerInput<'_>,
) -> Result<(f64, usize), DeseqError> {
    let counts = input.counts;
    let mu = input.mu;
    let design = input.design;
    let mut options = input.options;
    let n_samples = input.n_samples;
    let prior = input.prior;
    let weights = input.weights;
    if counts.len() != mu.len() {
        return Err(invalid_dimensions(
            "dispersion grid mu",
            counts.len(),
            mu.len(),
        ));
    }
    if design.is_none() {
        options.use_cox_reid = false;
    }
    validate_gene_est_options(options)?;
    validate_dispersion_prior(prior)?;
    if options.use_cox_reid && design.is_none() {
        return Err(DeseqError::UnsupportedFeature {
            feature: "Cox-Reid dispersion fitting requires a design matrix".to_string(),
        });
    }
    let max_disp = max_dispersion(options, n_samples);
    validate_dispersion_bounds(options.min_disp, max_disp)?;
    let objective = DispersionObjectiveInput {
        counts,
        mu,
        design,
        use_cox_reid: options.use_cox_reid,
        prior,
        weights,
        weight_threshold: options.weight_threshold,
    };
    let min_log = options.min_disp.ln();
    let max_log = max_disp.ln();
    let coarse = linspace(min_log, max_log, options.grid_points);
    let (best_log, _) = best_log_alpha(objective, &coarse)?;
    let delta = coarse[1] - coarse[0];
    let fine = linspace(best_log - delta, best_log + delta, options.grid_points);
    let (best_fine_log, _) = best_log_alpha(objective, &fine)?;
    Ok((
        best_fine_log.exp().clamp(options.min_disp, max_disp),
        options.grid_points * 2,
    ))
}

/// DESeq2's alpha-dependent NB log-likelihood kernel.
///
/// Terms independent of alpha are omitted, matching the objective used inside
/// DESeq2's dispersion optimizer.
pub fn dispersion_nb_log_likelihood_kernel(
    counts: &[u32],
    mu: &[f64],
    log_alpha: f64,
) -> Result<f64, DeseqError> {
    dispersion_nb_log_likelihood_kernel_weighted(counts, mu, log_alpha, None)
}

/// DESeq2's alpha-dependent NB log-likelihood kernel with optional observation weights.
///
/// Terms independent of alpha are omitted, matching the objective used inside
/// DESeq2's dispersion optimizer. When supplied, observation weights multiply
/// the per-sample terms.
pub fn dispersion_nb_log_likelihood_kernel_weighted(
    counts: &[u32],
    mu: &[f64],
    log_alpha: f64,
    weights: Option<&[f64]>,
) -> Result<f64, DeseqError> {
    if counts.len() != mu.len() {
        return Err(invalid_dimensions(
            "dispersion objective mu",
            counts.len(),
            mu.len(),
        ));
    }
    validate_observation_weight_slice(weights, counts.len(), "dispersion objective weights")?;
    if !log_alpha.is_finite() {
        return Err(DeseqError::InvalidDispersion {
            reason: "log dispersion must be finite".to_string(),
        });
    }
    let alpha = log_alpha.exp();
    if !alpha.is_finite() || alpha <= 0.0 {
        return Err(DeseqError::InvalidDispersion {
            reason: "dispersion must be finite and positive".to_string(),
        });
    }
    let inv_alpha = alpha.recip();
    let mut total = 0.0;
    for (sample, (count, mu)) in counts.iter().copied().zip(mu.iter().copied()).enumerate() {
        validate_positive_mu(mu, sample)?;
        let observation_weight = weights.map(|values| values[sample]).unwrap_or(1.0);
        let y = f64::from(count);
        total += observation_weight
            * (ln_gamma(y + inv_alpha)
                - ln_gamma(inv_alpha)
                - y * (mu + inv_alpha).ln()
                - inv_alpha * (1.0 + mu * alpha).ln());
    }
    Ok(total)
}

/// Derivative of DESeq2's alpha-dependent NB likelihood kernel with respect to log alpha.
pub fn dispersion_nb_log_likelihood_kernel_derivative(
    counts: &[u32],
    mu: &[f64],
    log_alpha: f64,
) -> Result<f64, DeseqError> {
    dispersion_nb_log_likelihood_kernel_derivative_weighted(counts, mu, log_alpha, None)
}

/// Derivative of the weighted NB likelihood kernel with respect to log alpha.
pub fn dispersion_nb_log_likelihood_kernel_derivative_weighted(
    counts: &[u32],
    mu: &[f64],
    log_alpha: f64,
    weights: Option<&[f64]>,
) -> Result<f64, DeseqError> {
    if counts.len() != mu.len() {
        return Err(invalid_dimensions(
            "dispersion objective derivative mu",
            counts.len(),
            mu.len(),
        ));
    }
    validate_observation_weight_slice(
        weights,
        counts.len(),
        "dispersion objective derivative weights",
    )?;
    if !log_alpha.is_finite() {
        return Err(DeseqError::InvalidDispersion {
            reason: "log dispersion must be finite".to_string(),
        });
    }
    let alpha = log_alpha.exp();
    if !alpha.is_finite() || alpha <= 0.0 {
        return Err(DeseqError::InvalidDispersion {
            reason: "dispersion must be finite and positive".to_string(),
        });
    }
    let inv_alpha = alpha.recip();
    let inv_alpha_squared = inv_alpha.powi(2);
    let mut derivative_alpha = 0.0;
    for (sample, (count, mu)) in counts.iter().copied().zip(mu.iter().copied()).enumerate() {
        validate_positive_mu(mu, sample)?;
        let observation_weight = weights.map(|values| values[sample]).unwrap_or(1.0);
        let y = f64::from(count);
        derivative_alpha += observation_weight
            * (digamma(inv_alpha) + (1.0 + mu * alpha).ln()
                - mu * alpha / (1.0 + mu * alpha)
                - digamma(y + inv_alpha)
                + y / (mu + inv_alpha));
    }
    Ok(inv_alpha_squared * derivative_alpha * alpha)
}

/// Second derivative of DESeq2's NB likelihood kernel with respect to log alpha.
pub fn dispersion_nb_log_likelihood_kernel_second_derivative(
    counts: &[u32],
    mu: &[f64],
    log_alpha: f64,
) -> Result<f64, DeseqError> {
    dispersion_nb_log_likelihood_kernel_second_derivative_weighted(counts, mu, log_alpha, None)
}

/// Second derivative of the weighted NB likelihood kernel with respect to log alpha.
pub fn dispersion_nb_log_likelihood_kernel_second_derivative_weighted(
    counts: &[u32],
    mu: &[f64],
    log_alpha: f64,
    weights: Option<&[f64]>,
) -> Result<f64, DeseqError> {
    if counts.len() != mu.len() {
        return Err(invalid_dimensions(
            "dispersion objective second derivative mu",
            counts.len(),
            mu.len(),
        ));
    }
    validate_observation_weight_slice(
        weights,
        counts.len(),
        "dispersion objective second derivative weights",
    )?;
    if !log_alpha.is_finite() {
        return Err(DeseqError::InvalidDispersion {
            reason: "log dispersion must be finite".to_string(),
        });
    }
    let alpha = log_alpha.exp();
    if !alpha.is_finite() || alpha <= 0.0 {
        return Err(DeseqError::InvalidDispersion {
            reason: "dispersion must be finite and positive".to_string(),
        });
    }
    let inv_alpha = alpha.recip();
    let inv_alpha_squared = inv_alpha.powi(2);
    let mut first_alpha_sum = 0.0;
    let mut second_alpha_sum = 0.0;
    for (sample, (count, mu)) in counts.iter().copied().zip(mu.iter().copied()).enumerate() {
        validate_positive_mu(mu, sample)?;
        let observation_weight = weights.map(|values| values[sample]).unwrap_or(1.0);
        let y = f64::from(count);
        let first_term = digamma(inv_alpha) + (1.0 + mu * alpha).ln()
            - mu * alpha / (1.0 + mu * alpha)
            - digamma(y + inv_alpha)
            + y / (mu + inv_alpha);
        let second_term = -inv_alpha_squared * trigamma(inv_alpha)?
            + mu.powi(2) * alpha * (1.0 + mu * alpha).powi(-2)
            + inv_alpha_squared * trigamma(y + inv_alpha)?
            + inv_alpha_squared * y * (mu + inv_alpha).powi(-2);
        first_alpha_sum += observation_weight * first_term;
        second_alpha_sum += observation_weight * second_term;
    }
    let second_alpha =
        -2.0 * inv_alpha.powi(3) * first_alpha_sum + inv_alpha_squared * second_alpha_sum;
    let first_log_alpha =
        dispersion_nb_log_likelihood_kernel_derivative_weighted(counts, mu, log_alpha, weights)?;
    Ok(second_alpha * alpha.powi(2) + first_log_alpha)
}

/// Cox-Reid adjustment term for one gene and design matrix.
pub fn cox_reid_adjustment(
    design: &DesignMatrix,
    mu: &[f64],
    log_alpha: f64,
) -> Result<f64, DeseqError> {
    cox_reid_adjustment_weighted(design, mu, log_alpha, None)
}

/// Cox-Reid adjustment term with optional DESeq2-style weighted sample subset.
pub fn cox_reid_adjustment_weighted(
    design: &DesignMatrix,
    mu: &[f64],
    log_alpha: f64,
    weights: Option<&[f64]>,
) -> Result<f64, DeseqError> {
    cox_reid_adjustment_weighted_with_threshold(
        design,
        mu,
        log_alpha,
        weights,
        GeneWiseDispersionOptions::default().weight_threshold,
    )
}

fn cox_reid_adjustment_weighted_with_threshold(
    design: &DesignMatrix,
    mu: &[f64],
    log_alpha: f64,
    weights: Option<&[f64]>,
    weight_threshold: f64,
) -> Result<f64, DeseqError> {
    if design.n_samples() != mu.len() {
        return Err(invalid_dimensions(
            "Cox-Reid design samples",
            mu.len(),
            design.n_samples(),
        ));
    }
    validate_observation_weight_slice(weights, mu.len(), "Cox-Reid weights")?;
    if !log_alpha.is_finite() {
        return Err(DeseqError::InvalidDispersion {
            reason: "log dispersion must be finite".to_string(),
        });
    }
    let alpha = log_alpha.exp();
    if !alpha.is_finite() || alpha <= 0.0 {
        return Err(DeseqError::InvalidDispersion {
            reason: "dispersion must be finite and positive".to_string(),
        });
    }
    let matrices = cox_reid_weighted_design_matrices_with_threshold(
        design,
        mu,
        alpha,
        weights,
        weight_threshold,
    )?;
    let determinant = matrices.xtwx.determinant();
    if !determinant.is_finite() || determinant <= 0.0 {
        return Err(DeseqError::InvalidDimensions {
            context: "Cox-Reid weighted design determinant".to_string(),
            expected: design.n_coefficients(),
            actual: 0,
        });
    }
    Ok(-0.5 * determinant.ln())
}

struct CoxReidDesignMatrices {
    xtwx: DMatrix<f64>,
    d_xtwx: DMatrix<f64>,
    d2_xtwx: DMatrix<f64>,
}

fn cox_reid_weighted_design_matrices_with_threshold(
    design: &DesignMatrix,
    mu: &[f64],
    alpha: f64,
    weights: Option<&[f64]>,
    weight_threshold: f64,
) -> Result<CoxReidDesignMatrices, DeseqError> {
    if design.n_samples() != mu.len() {
        return Err(invalid_dimensions(
            "Cox-Reid design samples",
            mu.len(),
            design.n_samples(),
        ));
    }
    if !alpha.is_finite() || alpha <= 0.0 {
        return Err(DeseqError::InvalidDispersion {
            reason: "dispersion must be finite and positive".to_string(),
        });
    }
    validate_observation_weight_slice(weights, mu.len(), "Cox-Reid weights")?;
    validate_weight_threshold(weight_threshold, "Cox-Reid weight threshold")?;
    let selected_samples = cox_reid_sample_indices(weights, mu.len(), weight_threshold);
    let selected_columns = match weights {
        Some(_) => cox_reid_column_indices(design, &selected_samples)?,
        None => (0..design.n_coefficients()).collect(),
    };
    if selected_samples.is_empty() || selected_columns.is_empty() {
        return Err(DeseqError::InvalidDimensions {
            context: "Cox-Reid weighted design subset".to_string(),
            expected: design.n_coefficients(),
            actual: 0,
        });
    }
    let p = selected_columns.len();
    let mut xtwx = DMatrix::<f64>::zeros(p, p);
    let mut d_xtwx = DMatrix::<f64>::zeros(p, p);
    let mut d2_xtwx = DMatrix::<f64>::zeros(p, p);
    for sample in selected_samples {
        let mu = mu[sample];
        validate_positive_mu(mu, sample)?;
        let inverse_variance = mu.recip() + alpha;
        let weight = inverse_variance.recip();
        let d_weight = -inverse_variance.powi(-2);
        let d2_weight = 2.0 * inverse_variance.powi(-3);
        let row = design.matrix().row(sample)?;
        for (left_idx, left_col) in selected_columns.iter().copied().enumerate() {
            for (right_idx, right_col) in selected_columns.iter().copied().enumerate() {
                let x_product = row[left_col] * row[right_col];
                xtwx[(left_idx, right_idx)] += x_product * weight;
                d_xtwx[(left_idx, right_idx)] += x_product * d_weight;
                d2_xtwx[(left_idx, right_idx)] += x_product * d2_weight;
            }
        }
    }
    Ok(CoxReidDesignMatrices {
        xtwx,
        d_xtwx,
        d2_xtwx,
    })
}

fn cox_reid_sample_indices(
    weights: Option<&[f64]>,
    n_samples: usize,
    weight_threshold: f64,
) -> Vec<usize> {
    match weights {
        Some(weights) => weights
            .iter()
            .copied()
            .enumerate()
            .filter_map(|(sample, weight)| (weight > weight_threshold).then_some(sample))
            .collect(),
        None => (0..n_samples).collect(),
    }
}

fn cox_reid_column_indices(
    design: &DesignMatrix,
    selected_samples: &[usize],
) -> Result<Vec<usize>, DeseqError> {
    let mut selected = Vec::with_capacity(design.n_coefficients());
    for column in 0..design.n_coefficients() {
        let mut sum_abs = 0.0;
        for sample in selected_samples {
            sum_abs += design.matrix().row(*sample)?[column].abs();
        }
        if sum_abs > 0.0 {
            selected.push(column);
        }
    }
    Ok(selected)
}

fn trace_product(left: &DMatrix<f64>, right: &DMatrix<f64>) -> f64 {
    let product = left * right;
    product.diagonal().iter().sum()
}

/// Derivative of the Cox-Reid adjustment with respect to log alpha.
pub fn cox_reid_adjustment_derivative(
    design: &DesignMatrix,
    mu: &[f64],
    log_alpha: f64,
) -> Result<f64, DeseqError> {
    cox_reid_adjustment_derivative_weighted(design, mu, log_alpha, None)
}

/// Derivative of the weighted Cox-Reid adjustment with respect to log alpha.
///
/// Observation weights define the DESeq2 `weightThreshold` sample subset for
/// the determinant; they do not multiply the Cox-Reid working weights.
pub fn cox_reid_adjustment_derivative_weighted(
    design: &DesignMatrix,
    mu: &[f64],
    log_alpha: f64,
    weights: Option<&[f64]>,
) -> Result<f64, DeseqError> {
    cox_reid_adjustment_derivative_weighted_with_threshold(
        design,
        mu,
        log_alpha,
        weights,
        GeneWiseDispersionOptions::default().weight_threshold,
    )
}

fn cox_reid_adjustment_derivative_weighted_with_threshold(
    design: &DesignMatrix,
    mu: &[f64],
    log_alpha: f64,
    weights: Option<&[f64]>,
    weight_threshold: f64,
) -> Result<f64, DeseqError> {
    let alpha = log_alpha.exp();
    let matrices = cox_reid_weighted_design_matrices_with_threshold(
        design,
        mu,
        alpha,
        weights,
        weight_threshold,
    )?;
    let Some(inverse) = matrices.xtwx.try_inverse() else {
        return Err(DeseqError::InvalidDimensions {
            context: "Cox-Reid weighted design inverse".to_string(),
            expected: design.n_coefficients(),
            actual: 0,
        });
    };
    Ok(-0.5 * trace_product(&inverse, &matrices.d_xtwx) * alpha)
}

/// Second derivative of the Cox-Reid adjustment with respect to log alpha.
pub fn cox_reid_adjustment_second_derivative(
    design: &DesignMatrix,
    mu: &[f64],
    log_alpha: f64,
) -> Result<f64, DeseqError> {
    cox_reid_adjustment_second_derivative_weighted(design, mu, log_alpha, None)
}

/// Second derivative of the weighted Cox-Reid adjustment with respect to log alpha.
///
/// Observation weights define the DESeq2 `weightThreshold` sample subset for
/// the determinant; they do not multiply the Cox-Reid working weights.
pub fn cox_reid_adjustment_second_derivative_weighted(
    design: &DesignMatrix,
    mu: &[f64],
    log_alpha: f64,
    weights: Option<&[f64]>,
) -> Result<f64, DeseqError> {
    cox_reid_adjustment_second_derivative_weighted_with_threshold(
        design,
        mu,
        log_alpha,
        weights,
        GeneWiseDispersionOptions::default().weight_threshold,
    )
}

fn cox_reid_adjustment_second_derivative_weighted_with_threshold(
    design: &DesignMatrix,
    mu: &[f64],
    log_alpha: f64,
    weights: Option<&[f64]>,
    weight_threshold: f64,
) -> Result<f64, DeseqError> {
    let alpha = log_alpha.exp();
    let matrices = cox_reid_weighted_design_matrices_with_threshold(
        design,
        mu,
        alpha,
        weights,
        weight_threshold,
    )?;
    let Some(inverse) = matrices.xtwx.try_inverse() else {
        return Err(DeseqError::InvalidDimensions {
            context: "Cox-Reid weighted design inverse".to_string(),
            expected: design.n_coefficients(),
            actual: 0,
        });
    };
    let trace_bi_db = trace_product(&inverse, &matrices.d_xtwx);
    let second_trace = (&inverse * &matrices.d_xtwx * &inverse * &matrices.d_xtwx)
        .diagonal()
        .iter()
        .sum::<f64>();
    let trace_bi_d2b = trace_product(&inverse, &matrices.d2_xtwx);
    let second_alpha =
        0.5 * trace_bi_db.powi(2) - 0.5 * (trace_bi_db.powi(2) - second_trace + trace_bi_d2b);
    let first_log_alpha = cox_reid_adjustment_derivative_weighted_with_threshold(
        design,
        mu,
        log_alpha,
        weights,
        weight_threshold,
    )?;
    Ok(second_alpha * alpha.powi(2) + first_log_alpha)
}

/// DESeq2's log-dispersion prior kernel, omitting additive constants.
pub fn dispersion_prior_log_density(
    log_alpha: f64,
    prior: DispersionPrior,
) -> Result<f64, DeseqError> {
    if !log_alpha.is_finite() {
        return Err(DeseqError::InvalidDispersion {
            reason: "log dispersion must be finite".to_string(),
        });
    }
    validate_dispersion_prior(Some(prior))?;
    Ok(-0.5 * (log_alpha - prior.log_mean).powi(2) / prior.variance)
}

/// Derivative of the log-dispersion prior kernel with respect to log alpha.
pub fn dispersion_prior_derivative(
    log_alpha: f64,
    prior: DispersionPrior,
) -> Result<f64, DeseqError> {
    if !log_alpha.is_finite() {
        return Err(DeseqError::InvalidDispersion {
            reason: "log dispersion must be finite".to_string(),
        });
    }
    validate_dispersion_prior(Some(prior))?;
    Ok(-(log_alpha - prior.log_mean) / prior.variance)
}

/// Second derivative of the log-dispersion prior kernel with respect to log alpha.
pub fn dispersion_prior_second_derivative(
    log_alpha: f64,
    prior: DispersionPrior,
) -> Result<f64, DeseqError> {
    if !log_alpha.is_finite() {
        return Err(DeseqError::InvalidDispersion {
            reason: "log dispersion must be finite".to_string(),
        });
    }
    validate_dispersion_prior(Some(prior))?;
    Ok(-prior.variance.recip())
}

/// Dispersion log posterior without prior and with optional Cox-Reid correction.
pub fn dispersion_log_posterior(
    counts: &[u32],
    mu: &[f64],
    design: Option<&DesignMatrix>,
    log_alpha: f64,
    use_cox_reid: bool,
) -> Result<f64, DeseqError> {
    dispersion_log_posterior_with_prior(counts, mu, design, log_alpha, use_cox_reid, None)
}

/// Dispersion log posterior with optional Cox-Reid correction and log-alpha prior.
pub fn dispersion_log_posterior_with_prior(
    counts: &[u32],
    mu: &[f64],
    design: Option<&DesignMatrix>,
    log_alpha: f64,
    use_cox_reid: bool,
    prior: Option<DispersionPrior>,
) -> Result<f64, DeseqError> {
    dispersion_log_posterior_with_prior_and_weights(
        counts,
        mu,
        design,
        log_alpha,
        use_cox_reid,
        prior,
        None,
    )
}

/// Dispersion log posterior with optional Cox-Reid correction, log-alpha prior, and weights.
pub fn dispersion_log_posterior_with_prior_and_weights(
    counts: &[u32],
    mu: &[f64],
    design: Option<&DesignMatrix>,
    log_alpha: f64,
    use_cox_reid: bool,
    prior: Option<DispersionPrior>,
    weights: Option<&[f64]>,
) -> Result<f64, DeseqError> {
    dispersion_log_posterior_objective(
        DispersionObjectiveInput {
            counts,
            mu,
            design,
            use_cox_reid,
            prior,
            weights,
            weight_threshold: GeneWiseDispersionOptions::default().weight_threshold,
        },
        log_alpha,
    )
}

fn dispersion_log_posterior_objective(
    input: DispersionObjectiveInput<'_>,
    log_alpha: f64,
) -> Result<f64, DeseqError> {
    let likelihood = dispersion_nb_log_likelihood_kernel_weighted(
        input.counts,
        input.mu,
        log_alpha,
        input.weights,
    )?;
    let posterior = if input.use_cox_reid {
        let Some(design) = input.design else {
            return Err(DeseqError::UnsupportedFeature {
                feature: "Cox-Reid dispersion objective requires a design matrix".to_string(),
            });
        };
        likelihood
            + cox_reid_adjustment_weighted_with_threshold(
                design,
                input.mu,
                log_alpha,
                input.weights,
                input.weight_threshold,
            )?
    } else {
        likelihood
    };
    if let Some(prior) = input.prior {
        Ok(posterior + dispersion_prior_log_density(log_alpha, prior)?)
    } else {
        Ok(posterior)
    }
}

/// Derivative of the dispersion log posterior with respect to log alpha.
pub fn dispersion_log_posterior_derivative(
    counts: &[u32],
    mu: &[f64],
    design: Option<&DesignMatrix>,
    log_alpha: f64,
    use_cox_reid: bool,
) -> Result<f64, DeseqError> {
    dispersion_log_posterior_derivative_with_prior(
        counts,
        mu,
        design,
        log_alpha,
        use_cox_reid,
        None,
    )
}

/// Derivative of the dispersion log posterior with an optional log-alpha prior.
pub fn dispersion_log_posterior_derivative_with_prior(
    counts: &[u32],
    mu: &[f64],
    design: Option<&DesignMatrix>,
    log_alpha: f64,
    use_cox_reid: bool,
    prior: Option<DispersionPrior>,
) -> Result<f64, DeseqError> {
    dispersion_log_posterior_derivative_with_prior_and_weights(
        counts,
        mu,
        design,
        log_alpha,
        use_cox_reid,
        prior,
        None,
    )
}

/// Derivative of the dispersion log posterior with optional prior and weights.
pub fn dispersion_log_posterior_derivative_with_prior_and_weights(
    counts: &[u32],
    mu: &[f64],
    design: Option<&DesignMatrix>,
    log_alpha: f64,
    use_cox_reid: bool,
    prior: Option<DispersionPrior>,
    weights: Option<&[f64]>,
) -> Result<f64, DeseqError> {
    dispersion_log_posterior_derivative_objective(
        DispersionObjectiveInput {
            counts,
            mu,
            design,
            use_cox_reid,
            prior,
            weights,
            weight_threshold: GeneWiseDispersionOptions::default().weight_threshold,
        },
        log_alpha,
    )
}

fn dispersion_log_posterior_derivative_objective(
    input: DispersionObjectiveInput<'_>,
    log_alpha: f64,
) -> Result<f64, DeseqError> {
    let likelihood = dispersion_nb_log_likelihood_kernel_derivative_weighted(
        input.counts,
        input.mu,
        log_alpha,
        input.weights,
    )?;
    let derivative = if input.use_cox_reid {
        let Some(design) = input.design else {
            return Err(DeseqError::UnsupportedFeature {
                feature: "Cox-Reid dispersion derivative requires a design matrix".to_string(),
            });
        };
        likelihood
            + cox_reid_adjustment_derivative_weighted_with_threshold(
                design,
                input.mu,
                log_alpha,
                input.weights,
                input.weight_threshold,
            )?
    } else {
        likelihood
    };
    if let Some(prior) = input.prior {
        Ok(derivative + dispersion_prior_derivative(log_alpha, prior)?)
    } else {
        Ok(derivative)
    }
}

/// Second derivative of the dispersion log posterior with respect to log alpha.
pub fn dispersion_log_posterior_second_derivative(
    counts: &[u32],
    mu: &[f64],
    design: Option<&DesignMatrix>,
    log_alpha: f64,
    use_cox_reid: bool,
) -> Result<f64, DeseqError> {
    dispersion_log_posterior_second_derivative_with_prior(
        counts,
        mu,
        design,
        log_alpha,
        use_cox_reid,
        None,
    )
}

/// Second derivative of the dispersion log posterior with an optional log-alpha prior.
pub fn dispersion_log_posterior_second_derivative_with_prior(
    counts: &[u32],
    mu: &[f64],
    design: Option<&DesignMatrix>,
    log_alpha: f64,
    use_cox_reid: bool,
    prior: Option<DispersionPrior>,
) -> Result<f64, DeseqError> {
    dispersion_log_posterior_second_derivative_with_prior_and_weights(
        counts,
        mu,
        design,
        log_alpha,
        use_cox_reid,
        prior,
        None,
    )
}

/// Second derivative of the dispersion log posterior with optional prior and weights.
pub fn dispersion_log_posterior_second_derivative_with_prior_and_weights(
    counts: &[u32],
    mu: &[f64],
    design: Option<&DesignMatrix>,
    log_alpha: f64,
    use_cox_reid: bool,
    prior: Option<DispersionPrior>,
    weights: Option<&[f64]>,
) -> Result<f64, DeseqError> {
    dispersion_log_posterior_second_derivative_objective(
        DispersionObjectiveInput {
            counts,
            mu,
            design,
            use_cox_reid,
            prior,
            weights,
            weight_threshold: GeneWiseDispersionOptions::default().weight_threshold,
        },
        log_alpha,
    )
}

fn dispersion_log_posterior_second_derivative_objective(
    input: DispersionObjectiveInput<'_>,
    log_alpha: f64,
) -> Result<f64, DeseqError> {
    let likelihood = dispersion_nb_log_likelihood_kernel_second_derivative_weighted(
        input.counts,
        input.mu,
        log_alpha,
        input.weights,
    )?;
    let second_derivative = if input.use_cox_reid {
        let Some(design) = input.design else {
            return Err(DeseqError::UnsupportedFeature {
                feature: "Cox-Reid dispersion second derivative requires a design matrix"
                    .to_string(),
            });
        };
        likelihood
            + cox_reid_adjustment_second_derivative_weighted_with_threshold(
                design,
                input.mu,
                log_alpha,
                input.weights,
                input.weight_threshold,
            )?
    } else {
        likelihood
    };
    if let Some(prior) = input.prior {
        Ok(second_derivative + dispersion_prior_second_derivative(log_alpha, prior)?)
    } else {
        Ok(second_derivative)
    }
}

fn best_log_alpha(
    objective: DispersionObjectiveInput<'_>,
    grid: &[f64],
) -> Result<(f64, f64), DeseqError> {
    let mut best_log = grid[0];
    let mut best_score = dispersion_log_posterior_objective(objective, best_log)?;
    for log_alpha in grid.iter().copied().skip(1) {
        let score = dispersion_log_posterior_objective(objective, log_alpha)?;
        if score > best_score {
            best_log = log_alpha;
            best_score = score;
        }
    }
    Ok((best_log, best_score))
}

fn linspace(start: f64, end: f64, len: usize) -> Vec<f64> {
    let step = (end - start) / (len as f64 - 1.0);
    (0..len).map(|idx| start + step * idx as f64).collect()
}

fn compact_counts_rows(
    counts: &CountMatrix,
    row_indices: &[usize],
) -> Result<CountMatrix, DeseqError> {
    let mut values = Vec::with_capacity(row_indices.len() * counts.n_samples());
    for row in row_indices {
        values.extend_from_slice(counts.row(*row)?);
    }
    let gene_names = counts.gene_names().map(|names| {
        row_indices
            .iter()
            .map(|row| names[*row].clone())
            .collect::<Vec<_>>()
    });
    let sample_names = counts.sample_names().map(<[String]>::to_vec);
    CountMatrix::from_row_major_u32_with_names(
        row_indices.len(),
        counts.n_samples(),
        values,
        gene_names,
        sample_names,
    )
}

fn compact_matrix_rows(
    matrix: &RowMajorMatrix<f64>,
    row_indices: &[usize],
) -> Result<RowMajorMatrix<f64>, DeseqError> {
    let mut values = Vec::with_capacity(row_indices.len() * matrix.n_cols());
    for row in row_indices {
        values.extend_from_slice(matrix.row(*row)?);
    }
    RowMajorMatrix::from_row_major(row_indices.len(), matrix.n_cols(), values)
}

fn compact_gene_values(values: &[f64], row_indices: &[usize]) -> Result<Vec<f64>, DeseqError> {
    let mut compact = Vec::with_capacity(row_indices.len());
    for row in row_indices {
        let Some(value) = values.get(*row).copied() else {
            return Err(invalid_dimensions(
                "compact gene values",
                row + 1,
                values.len(),
            ));
        };
        compact.push(value);
    }
    Ok(compact)
}

fn validate_gene_est_inputs(
    input: &GeneWiseDispersionInput<'_>,
    options: GeneWiseDispersionOptions,
) -> Result<(), DeseqError> {
    if input.design.n_samples() != input.counts.n_samples() {
        return Err(invalid_dimensions(
            "dispersion design samples",
            input.counts.n_samples(),
            input.design.n_samples(),
        ));
    }
    if input.normalized_counts.n_rows() != input.counts.n_genes()
        || input.normalized_counts.n_cols() != input.counts.n_samples()
    {
        return Err(DeseqError::InvalidDimensions {
            context: "dispersion normalized counts".to_string(),
            expected: input.counts.n_genes() * input.counts.n_samples(),
            actual: input.normalized_counts.len(),
        });
    }
    if input.base_mean.len() != input.counts.n_genes() {
        return Err(invalid_dimensions(
            "dispersion baseMean",
            input.counts.n_genes(),
            input.base_mean.len(),
        ));
    }
    if input.base_var.len() != input.counts.n_genes() {
        return Err(invalid_dimensions(
            "dispersion baseVar",
            input.counts.n_genes(),
            input.base_var.len(),
        ));
    }
    if input.all_zero.len() != input.counts.n_genes() {
        return Err(invalid_dimensions(
            "dispersion allZero",
            input.counts.n_genes(),
            input.all_zero.len(),
        ));
    }
    if input.size_factors.len() != input.counts.n_samples() {
        return Err(invalid_dimensions(
            "dispersion size factors",
            input.counts.n_samples(),
            input.size_factors.len(),
        ));
    }
    validate_size_factors(input.size_factors)?;
    if let Some(normalization_factors) = input.normalization_factors {
        validate_normalization_factors(
            normalization_factors,
            input.counts.n_genes(),
            input.counts.n_samples(),
        )?;
    }
    if let Some(observation_weights) = input.observation_weights {
        if observation_weights.n_rows() != input.counts.n_genes()
            || observation_weights.n_cols() != input.counts.n_samples()
        {
            return Err(DeseqError::InvalidDimensions {
                context: "dispersion observation weights".to_string(),
                expected: input.counts.n_genes() * input.counts.n_samples(),
                actual: observation_weights.len(),
            });
        }
        for (idx, weight) in observation_weights.as_slice().iter().copied().enumerate() {
            if !weight.is_finite() || weight < 0.0 {
                return Err(DeseqError::NonFiniteValue {
                    context: "dispersion observation weight".to_string(),
                    index: Some(idx),
                    value: weight,
                });
            }
        }
    }
    validate_gene_est_options(options)?;
    validate_dispersion_bounds(
        options.min_disp,
        max_dispersion(options, input.counts.n_samples()),
    )?;
    Ok(())
}

fn validate_gene_est_options(options: GeneWiseDispersionOptions) -> Result<(), DeseqError> {
    if !options.min_disp.is_finite() || options.min_disp <= 0.0 {
        return Err(DeseqError::InvalidDispersion {
            reason: "min_disp must be finite and positive".to_string(),
        });
    }
    if (options.min_disp / 10.0).ln() <= -30.0 {
        return Err(DeseqError::InvalidDispersion {
            reason: "log(min_disp / 10) must be above -30 for numerical stability".to_string(),
        });
    }
    if !options.min_mu.is_finite() || options.min_mu <= 0.0 {
        return Err(DeseqError::NonFiniteValue {
            context: "dispersion min_mu".to_string(),
            index: None,
            value: options.min_mu,
        });
    }
    if options.grid_points < 3 {
        return Err(DeseqError::InvalidDimensions {
            context: "dispersion grid points".to_string(),
            expected: 3,
            actual: options.grid_points,
        });
    }
    if !options.kappa_0.is_finite() || options.kappa_0 <= 0.0 {
        return Err(DeseqError::InvalidDispersion {
            reason: "kappa_0 must be finite and positive".to_string(),
        });
    }
    if !options.disp_tol.is_finite() || options.disp_tol <= 0.0 {
        return Err(DeseqError::InvalidDispersion {
            reason: "disp_tol must be finite and positive".to_string(),
        });
    }
    if options.maxit == 0 {
        return Err(DeseqError::InvalidDimensions {
            context: "dispersion maxit".to_string(),
            expected: 1,
            actual: 0,
        });
    }
    if options.niter == 0 {
        return Err(DeseqError::InvalidDimensions {
            context: "dispersion niter".to_string(),
            expected: 1,
            actual: 0,
        });
    }
    validate_weight_threshold(options.weight_threshold, "dispersion weight_threshold")?;
    Ok(())
}

fn validate_weight_threshold(value: f64, context: &str) -> Result<(), DeseqError> {
    if !value.is_finite() || value < 0.0 {
        return Err(DeseqError::InvalidDispersion {
            reason: format!("{context} must be finite and non-negative"),
        });
    }
    Ok(())
}

fn validate_dispersion_prior(prior: Option<DispersionPrior>) -> Result<(), DeseqError> {
    if let Some(prior) = prior {
        if !prior.log_mean.is_finite() {
            return Err(DeseqError::InvalidDispersion {
                reason: "dispersion prior log_mean must be finite".to_string(),
            });
        }
        if !prior.variance.is_finite() || prior.variance <= 0.0 {
            return Err(DeseqError::InvalidDispersion {
                reason: "dispersion prior variance must be finite and positive".to_string(),
            });
        }
    }
    Ok(())
}

fn validate_size_factors(size_factors: &[f64]) -> Result<(), DeseqError> {
    for (idx, value) in size_factors.iter().copied().enumerate() {
        if !value.is_finite() || value <= 0.0 {
            return Err(DeseqError::InvalidSizeFactors {
                reason: format!("size factor at sample {idx} must be finite and positive"),
            });
        }
    }
    Ok(())
}

fn validate_normalization_factors(
    normalization_factors: &RowMajorMatrix<f64>,
    n_genes: usize,
    n_samples: usize,
) -> Result<(), DeseqError> {
    if normalization_factors.n_rows() != n_genes || normalization_factors.n_cols() != n_samples {
        return Err(DeseqError::InvalidDimensions {
            context: "dispersion normalization factors".to_string(),
            expected: n_genes * n_samples,
            actual: normalization_factors.len(),
        });
    }
    for (idx, value) in normalization_factors.as_slice().iter().copied().enumerate() {
        validate_normalization_factor(value, idx)?;
    }
    Ok(())
}

fn validate_normalization_factor(value: f64, index: usize) -> Result<(), DeseqError> {
    if !value.is_finite() || value <= 0.0 {
        return Err(DeseqError::InvalidSizeFactors {
            reason: format!("normalization factor at index {index} must be finite and positive"),
        });
    }
    Ok(())
}

fn validate_positive_mu(mu: f64, sample: usize) -> Result<(), DeseqError> {
    if !mu.is_finite() || mu <= 0.0 {
        return Err(DeseqError::NonFiniteValue {
            context: "dispersion mean".to_string(),
            index: Some(sample),
            value: mu,
        });
    }
    Ok(())
}

fn validate_observation_weight_slice(
    weights: Option<&[f64]>,
    expected_len: usize,
    context: &str,
) -> Result<(), DeseqError> {
    let Some(weights) = weights else {
        return Ok(());
    };
    if weights.len() != expected_len {
        return Err(invalid_dimensions(context, expected_len, weights.len()));
    }
    for (idx, weight) in weights.iter().copied().enumerate() {
        if !weight.is_finite() || weight < 0.0 {
            return Err(DeseqError::NonFiniteValue {
                context: context.to_string(),
                index: Some(idx),
                value: weight,
            });
        }
    }
    Ok(())
}

fn validate_dispersion_bounds(min_disp: f64, max_disp: f64) -> Result<(), DeseqError> {
    if !max_disp.is_finite() || max_disp <= min_disp {
        return Err(DeseqError::InvalidDispersion {
            reason: "max dispersion must be finite and greater than min dispersion".to_string(),
        });
    }
    Ok(())
}

fn max_dispersion(options: GeneWiseDispersionOptions, n_samples: usize) -> f64 {
    options
        .max_disp
        .unwrap_or_else(|| 10.0_f64.max(n_samples as f64))
}
