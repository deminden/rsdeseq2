//! Batch orchestration; likelihood, prior and optimizers have no R dependency.
use super::{
    analysis::{AnalysisRow, GridPosteriorRow, PosteriorAnalysis},
    beta_binomial,
    negative_binomial::Row,
    prior::{CauchyCoefficientPrior, ERROR_SCOPE, estimate_cauchy_prior_variance},
    regression::{self, GeneralRow, RegressionMethod, RegressionPosteriorOptions},
};
use crate::errors::{PosteriorArray, PosteriorArrayAxis, scoped_invalid_options};
use crate::numeric::math::ProfileMath;
use crate::{
    Error,
    matrix::{MatrixView, RowMajorMatrix},
    probability::normal_upper_tail,
    shrinkage::svalues,
};
use nalgebra::{DMatrix, Matrix4, linalg::try_invert_to};
use rcompat_lbfgsb::OptimResult;

/// MLE scale supplied for empirical prior adaptation.
#[derive(Clone, Copy, Debug)]
pub enum CoefficientScale {
    /// Natural-log units.
    Natural,
    /// Log2 units (DESeq2's fold changes).
    Log2,
}
impl CoefficientScale {
    /// Convert a value on this scale to natural-log units.
    pub fn to_natural(self, value: f64) -> f64 {
        match self {
            Self::Natural => value,
            Self::Log2 => value * std::f64::consts::LN_2,
        }
    }
}
/// Coefficient estimates and SEs for empirical prior adaptation.
#[derive(Clone, Copy, Debug)]
pub struct CoefficientMle<'a> {
    /// MLE estimates, one per row.
    pub estimates: &'a [f64],
    /// Their standard errors.
    pub standard_errors: &'a [f64],
    /// The scale of both.
    pub scale: CoefficientScale,
}
/// Log offsets or normalization values; scale conversion stays in the core.
#[derive(Clone, Copy, Debug, Default)]
pub enum PosteriorOffsets<'a> {
    #[default]
    /// No offsets.
    None,
    /// Row-by-sample offsets on the natural-log scale.
    Log(MatrixView<'a, f64>),
    /// Row-by-sample normalization factors (log taken internally).
    Factors(MatrixView<'a, f64>),
    /// Per-sample size factors.
    SizeFactors(&'a [f64]),
}
/// Borrowed primitive data for a negative-binomial posterior batch.
#[derive(Clone, Copy, Debug)]
pub struct NegativeBinomialPosteriorInput<'a> {
    /// Row-by-sample counts.
    pub counts: MatrixView<'a, f64>,
    /// Sample-by-coefficient design matrix.
    pub design: MatrixView<'a, f64>,
    /// One negative-binomial dispersion per row.
    pub dispersions: &'a [f64],
    /// Offsets of the linear predictor.
    pub offsets: PosteriorOffsets<'a>,
    /// Optional row-by-sample observation weights.
    pub weights: Option<MatrixView<'a, f64>>,
    /// MLE estimates used to adapt the prior scale, as `apeglm` does with `mle`.
    pub mle: Option<CoefficientMle<'a>>,
}
#[derive(Clone, Debug)]
/// Options of an `apeglm`-style negative-binomial posterior fit.
pub struct NegativeBinomialPosteriorOptions {
    /// The Cauchy coefficient prior (scale adapted when `mle` is supplied).
    pub prior: CauchyCoefficientPrior,
    /// Coefficient whose sign and threshold probabilities are reported.
    pub coefficient: Option<usize>,
    /// Positive magnitude in natural log units; None omits FSOS calculations.
    pub threshold: Option<f64>,
    /// Flip the sign convention of the reported coefficient, as `apeglm` does.
    pub flip_sign: bool,
    /// Multiplier applied to the reported coefficient.
    pub multiplier: f64,
    /// Return MAP modes only, without Laplace standard deviations.
    pub map_only: bool,
    /// Restart from the opposite sign and flag rows whose modes disagree.
    pub check_opposite_start: bool,
    /// Compiled CR zero-row policy; MAP-only methods ignore this flag.
    pub log_link: bool,
    /// Refinement of the MAP by the compiled (Newton-type) step.
    pub refinement: super::compiled_refinement::NegativeBinomialRefinement,
    /// Finite-difference step for the numerical Hessian.
    pub hessian_step: f64,
    /// Which posterior summaries to compute.
    pub analysis: super::analysis::PosteriorAnalysisOptions,
}
impl Default for NegativeBinomialPosteriorOptions {
    fn default() -> Self {
        Self {
            prior: Default::default(),
            coefficient: None,
            threshold: None,
            flip_sign: true,
            multiplier: 1.0,
            map_only: false,
            check_opposite_start: true,
            log_link: true,
            refinement: Default::default(),
            hessian_step: 1e-3,
            analysis: Default::default(),
        }
    }
}
/// Per-row native convergence, separate from posterior-validity decisions.
#[derive(Clone, Debug)]
pub struct PosteriorMapDiagnostic {
    /// Optimizer convergence code (0 = converged), when the row was fitted.
    pub convergence: Option<i32>,
    /// Objective evaluations.
    pub function_evaluations: usize,
    /// Whether the opposite-sign restart reached a different mode.
    pub opposite_start_disagrees: bool,
    /// Whether the Laplace approximation at the mode is valid.
    pub laplace_valid: bool,
}
#[derive(Clone, Debug)]
/// Output of a negative-binomial posterior batch.
pub struct NegativeBinomialPosterior {
    /// Coefficient owning probability summaries; only the fitter sets this identity.
    summary_coefficient: Option<usize>,
    /// Row-major modes and SDs in natural log coefficient units.
    pub mode: RowMajorMatrix<f64>,
    /// Row-major Laplace posterior standard deviations (absent for MAP-only fits).
    pub standard_deviation: Option<RowMajorMatrix<f64>>,
    /// Per-row false sign rate of the reported coefficient.
    pub false_sign_rate: Option<Vec<Option<f64>>>,
    /// s-values from the false sign rates.
    pub sign_svalues: Option<Vec<Option<f64>>>,
    /// Per-row false-sign-or-small rate against the threshold.
    pub false_sign_or_small_rate: Option<Vec<Option<f64>>>,
    /// s-values from the false-sign-or-small rates.
    pub threshold_svalues: Option<Vec<Option<f64>>>,
    /// The prior actually used.
    pub prior: CauchyCoefficientPrior,
    /// Per-row optimizer diagnostics.
    pub diagnostics: Vec<PosteriorMapDiagnostic>,
    /// Native fitting warnings in batch/phase order, for frontend presentation.
    pub warnings: Vec<String>,
    /// Advanced random-start state, absent on deterministic fitting routes.
    pub random_state: Option<Vec<i32>>,
    /// Additional posterior summaries.
    pub analysis: super::analysis::PosteriorAnalysis,
}
/// Numerical columns for frontend presentation, with no frontend arithmetic.
#[derive(Clone, Debug)]
pub struct PosteriorCoefficientSummary {
    /// Posterior modes.
    pub mode: Vec<f64>,
    /// Posterior standard deviations.
    pub standard_deviation: Vec<Option<f64>>,
    /// Sign s-values.
    pub sign_svalues: Option<Vec<Option<f64>>>,
    /// Threshold s-values.
    pub threshold_svalues: Option<Vec<Option<f64>>>,
}
impl NegativeBinomialPosterior {
    /// Coefficient for which the stored sign/threshold probabilities were fitted.
    pub fn summary_coefficient(&self) -> Option<usize> {
        self.summary_coefficient
    }

    fn validate_summary_shapes(&self) -> Result<(), Error> {
        const SUMMARY_SHAPE: &str = "posterior summary arrays have inconsistent shapes";
        let genes = self.mode.n_rows();
        let columns = self.mode.n_cols();
        for (array, matrix, expected_columns) in [
            (
                PosteriorArray::StandardDeviation,
                &self.standard_deviation,
                Some(columns),
            ),
            (PosteriorArray::Intervals, &self.analysis.intervals, Some(2)),
            (
                PosteriorArray::ContrastMode,
                &self.analysis.contrast_mode,
                None,
            ),
            (
                PosteriorArray::ContrastStandardDeviation,
                &self.analysis.contrast_standard_deviation,
                None,
            ),
        ] {
            if let Some(matrix) = matrix {
                validate_posterior_dimension(
                    array,
                    PosteriorArrayAxis::Rows,
                    genes,
                    matrix.n_rows(),
                    SUMMARY_SHAPE,
                )?;
                if let Some(expected) = expected_columns {
                    validate_posterior_dimension(
                        array,
                        PosteriorArrayAxis::Columns,
                        expected,
                        matrix.n_cols(),
                        SUMMARY_SHAPE,
                    )?;
                }
            }
        }
        validate_posterior_dimension(
            PosteriorArray::Diagnostics,
            PosteriorArrayAxis::Rows,
            genes,
            self.diagnostics.len(),
            SUMMARY_SHAPE,
        )?;
        for (array, values) in [
            (PosteriorArray::FalseSignRate, &self.false_sign_rate),
            (PosteriorArray::SignSvalues, &self.sign_svalues),
            (
                PosteriorArray::FalseSignOrSmallRate,
                &self.false_sign_or_small_rate,
            ),
            (PosteriorArray::ThresholdSvalues, &self.threshold_svalues),
        ] {
            if let Some(values) = values {
                validate_posterior_dimension(
                    array,
                    PosteriorArrayAxis::Rows,
                    genes,
                    values.len(),
                    SUMMARY_SHAPE,
                )?;
            }
        }
        if !self.analysis.grid_diagnostics.is_empty() {
            validate_posterior_dimension(
                PosteriorArray::GridDiagnostics,
                PosteriorArrayAxis::Rows,
                genes,
                self.analysis.grid_diagnostics.len(),
                SUMMARY_SHAPE,
            )?;
        }
        if let (Some(mode), Some(sd)) = (
            &self.analysis.contrast_mode,
            &self.analysis.contrast_standard_deviation,
        ) {
            validate_posterior_dimension(
                PosteriorArray::ContrastStandardDeviation,
                PosteriorArrayAxis::Columns,
                mode.n_cols(),
                sd.n_cols(),
                "posterior contrast arrays have inconsistent shapes",
            )?;
        }
        Ok(())
    }

    // Invariant: `coefficient` is checked against the mode shape above.
    #[allow(clippy::unwrap_used, clippy::expect_used)]
    /// Mode, SD and s-value columns of one coefficient.
    pub fn coefficient_summary(
        &self,
        coefficient: usize,
        scale: CoefficientScale,
    ) -> Result<PosteriorCoefficientSummary, Error> {
        if coefficient >= self.mode.n_cols() {
            return Err(scoped_invalid_options(
                ERROR_SCOPE,
                "coefficient is outside posterior matrix",
            ));
        }
        self.validate_summary_shapes()?;
        if (self.sign_svalues.is_some() || self.threshold_svalues.is_some())
            && self.summary_coefficient != Some(coefficient)
        {
            return Err(Error::PosteriorSummaryCoefficientMismatch {
                requested: coefficient,
                fitted: self.summary_coefficient,
            });
        }
        let factor = match scale {
            CoefficientScale::Natural => 1.0,
            CoefficientScale::Log2 => 1.0 / std::f64::consts::LN_2,
        };
        let mode = (0..self.mode.n_rows())
            .map(|i| self.mode.get(i, coefficient).unwrap() * factor)
            .collect();
        let standard_deviation = (0..self.mode.n_rows())
            .map(|i| {
                self.standard_deviation.as_ref().and_then(|m| {
                    let value = *m.get(i, coefficient).unwrap();
                    (!value.is_nan()).then_some(value * factor)
                })
            })
            .collect();
        Ok(PosteriorCoefficientSummary {
            mode,
            standard_deviation,
            sign_svalues: self.sign_svalues.clone(),
            threshold_svalues: self.threshold_svalues.clone(),
        })
    }
}

fn validate_posterior_dimension(
    array: PosteriorArray,
    axis: PosteriorArrayAxis,
    expected: usize,
    actual: usize,
    reason: &'static str,
) -> Result<(), Error> {
    if actual != expected {
        return Err(Error::PosteriorShapeMismatch {
            array,
            axis,
            expected,
            actual,
            reason,
        });
    }
    Ok(())
}

fn optimize(row: &Row<'_>, start: Vec<f64>) -> Result<rcompat_lbfgsb::OptimResult, Error> {
    let p = start.len();
    let scale = row.objective_gradient(&vec![0.0; p]).0.max(1.0);
    let (value, gradient) = row.objective_gradient(&start);
    if !scale.is_finite() || !value.is_finite() || gradient.iter().any(|v| !v.is_finite()) {
        return Err(scoped_invalid_options(
            ERROR_SCOPE,
            "initial row posterior exceeds the representable range",
        ));
    }
    Ok(super::unbounded_lbfgs::minimize(
        |b| {
            let (value, gradient) = row.objective_gradient(b);
            (
                value / scale + 10.0,
                gradient.into_iter().map(|g| g / scale).collect(),
            )
        },
        start,
    ))
}
fn refine(
    row: &Row<'_>,
    start: Vec<f64>,
    options: &NegativeBinomialPosteriorOptions,
) -> Result<rcompat_lbfgsb::OptimResult, Error> {
    super::compiled_refinement::optimize(row, start, &options.refinement)
}

/// Fit `apeglm`-style negative-binomial posteriors (MAP with a Cauchy
/// prior, Laplace standard deviations, false sign rates) for a batch of rows.
pub fn fit_negative_binomial_posterior(
    input: NegativeBinomialPosteriorInput<'_>,
    options: &NegativeBinomialPosteriorOptions,
) -> Result<NegativeBinomialPosterior, Error> {
    fit_batch(input, options, None, None)
}
/// MAP-only random-start fit. One shared N(0,0.5^2) vector is drawn for
/// the entire batch. On a later fitting error the borrowed RNG remains advanced.
pub fn fit_negative_binomial_posterior_with_rng(
    input: NegativeBinomialPosteriorInput<'_>,
    options: &NegativeBinomialPosteriorOptions,
    rng: &mut super::random::PosteriorRng,
) -> Result<NegativeBinomialPosterior, Error> {
    let mut options = options.clone();
    options.map_only = true;
    options.check_opposite_start = false;
    fit_batch(input, &options, None, Some(rng))
}
type RegressionRequest<'a, 'config> = Option<(
    Option<MatrixView<'a, f64>>,
    &'config super::regression::RegressionPosteriorOptions,
)>;

enum PreparedOffsets<'a> {
    None,
    Log(MatrixView<'a, f64>),
    Factors(MatrixView<'a, f64>),
    LogSizeFactors(Vec<f64>),
}

impl PreparedOffsets<'_> {
    fn log_offset(&self, gene: usize, sample: usize) -> f64 {
        match self {
            Self::None => 0.0,
            Self::Log(matrix) => *matrix.at(gene, sample),
            Self::Factors(matrix) => matrix.at(gene, sample).ln_p(),
            Self::LogSizeFactors(offsets) => offsets[sample],
        }
    }
}

struct PreparedBatch<'a> {
    prior: CauchyCoefficientPrior,
    analysis: super::analysis::PosteriorAnalysis,
    offsets: PreparedOffsets<'a>,
}

impl<'a> PreparedBatch<'a> {
    fn new(
        input: &NegativeBinomialPosteriorInput<'a>,
        options: &NegativeBinomialPosteriorOptions,
        regression: RegressionRequest<'_, '_>,
    ) -> Result<Self, Error> {
        let g = input.counts.n_rows();
        let n = input.counts.n_cols();
        let p = input.design.n_cols();
        if g == 0 || n == 0 || p == 0 || input.design.n_rows() != n || input.dispersions.len() != g
        {
            return Err(scoped_invalid_options(
                ERROR_SCOPE,
                "nonempty counts, sample-aligned design and gene-aligned dispersions are required",
            ));
        }
        if input.design.as_slice().iter().any(|v| !v.is_finite())
            || input
                .counts
                .as_slice()
                .iter()
                .any(|v| !v.is_finite() || (regression.is_none() && *v < 0.0))
        {
            return Err(scoped_invalid_options(
                ERROR_SCOPE,
                "counts must be finite nonnegative numbers and design must be finite",
            ));
        }
        if let Some(m) = input.weights
            && (m.n_rows() != g
                || m.n_cols() != n
                || m.as_slice().iter().any(|v| !v.is_finite() || *v < 0.0))
        {
            return Err(scoped_invalid_options(
                ERROR_SCOPE,
                "finite nonnegative weights must match counts shape",
            ));
        }
        let offsets = match input.offsets {
            PosteriorOffsets::SizeFactors(s) => {
                if s.len() != n || s.iter().any(|v| !v.is_finite() || *v <= 0.0) {
                    return Err(scoped_invalid_options(
                        ERROR_SCOPE,
                        "positive sample-aligned size factors are required",
                    ));
                }
                PreparedOffsets::LogSizeFactors(s.iter().map(|v| v.ln_p()).collect())
            }
            PosteriorOffsets::Log(m) | PosteriorOffsets::Factors(m) => {
                let factors = matches!(input.offsets, PosteriorOffsets::Factors(_));
                if m.n_rows() != g
                    || m.n_cols() != n
                    || m.as_slice()
                        .iter()
                        .any(|v| !v.is_finite() || (factors && *v <= 0.0))
                {
                    return Err(scoped_invalid_options(
                        ERROR_SCOPE,
                        "finite log offsets or positive normalization factors must match counts",
                    ));
                }
                if factors {
                    PreparedOffsets::Factors(m)
                } else {
                    PreparedOffsets::Log(m)
                }
            }
            PosteriorOffsets::None => PreparedOffsets::None,
        };
        if options.coefficient.is_some_and(|c| c >= p)
            || !options.multiplier.is_finite()
            || options.multiplier <= 0.0
            || !options.hessian_step.is_finite()
            || options.hessian_step <= 0.0
            || options
                .threshold
                .is_some_and(|t| !t.is_finite() || (options.flip_sign && t <= 0.0))
        {
            return Err(scoped_invalid_options(
                ERROR_SCOPE,
                "invalid coefficient, multiplier, Hessian step or threshold",
            ));
        }
        options.analysis.validate(
            g,
            p,
            options.coefficient,
            options.threshold,
            options.map_only,
        )?;
        if !options.map_only
            && options.analysis.parameter_standard_errors.is_some()
            && regression.is_some_and(|(parameters, _)| parameters.is_none_or(|m| m.n_cols() != 1))
        {
            return Err(scoped_invalid_options(
                ERROR_SCOPE,
                "parameter uncertainty requires a scalar likelihood parameter",
            ));
        }
        let analysis = super::analysis::PosteriorAnalysis::new(
            &options.analysis,
            g,
            options.coefficient,
            options.map_only,
        )?;
        let mut prior = options.prior.clone();
        if let Some(mle) = input.mle {
            let coefficient = options.coefficient.ok_or_else(|| {
                scoped_invalid_options(ERROR_SCOPE, "MLE prior adaptation requires a coefficient")
            })?;
            if mle.estimates.len() != g || mle.standard_errors.len() != g {
                return Err(scoped_invalid_options(
                    ERROR_SCOPE,
                    "MLE population must match count rows",
                ));
            }
            let factor = match mle.scale {
                CoefficientScale::Natural => 1.0,
                CoefficientScale::Log2 => std::f64::consts::LN_2,
            };
            let x = mle.estimates.iter().map(|v| v * factor).collect::<Vec<_>>();
            let s = mle
                .standard_errors
                .iter()
                .map(|v| v * factor)
                .collect::<Vec<_>>();
            let variance = estimate_cauchy_prior_variance(&x, &s)?;
            prior.scale = (options.multiplier * variance.sqrt()).min(1.0);
            prior.adapted_variance = Some(variance);
            prior.nuisance = (0..p).filter(|&j| j != coefficient).collect();
        }
        prior.validate(p)?;
        Ok(Self {
            prior,
            analysis,
            offsets,
        })
    }
}

enum RowLikelihood<'a> {
    CompiledNegativeBinomial(&'a Row<'a>),
    NegativeBinomialRegression(&'a Row<'a>, &'a RegressionPosteriorOptions),
    General(GeneralRow<'a>),
    BetaBinomial(beta_binomial::Row<'a>, &'a RegressionPosteriorOptions),
}

impl<'a> RowLikelihood<'a> {
    fn new(
        row: &'a Row<'a>,
        regression: RegressionRequest<'_, 'a>,
        gene: usize,
    ) -> Result<Self, Error> {
        let Some((parameters, config)) = regression else {
            return Ok(Self::CompiledNegativeBinomial(row));
        };
        match config.method {
            RegressionMethod::General => Ok(Self::General(GeneralRow {
                y: &row.counts,
                design: row.design,
                parameters: parameters.map(|m| (0..m.n_cols()).map(|j| *m.at(gene, j)).collect()),
                weights: &row.weights,
                offsets: &row.offsets,
                prior: row.prior,
                options: config,
                error: std::cell::RefCell::new(None),
            })),
            RegressionMethod::NegativeBinomialR => {
                Ok(Self::NegativeBinomialRegression(row, config))
            }
            RegressionMethod::BetaBinomialR
            | RegressionMethod::BetaBinomialCR
            | RegressionMethod::BetaBinomialC
            | RegressionMethod::BetaBinomialCStar => {
                let parameters = parameters.ok_or_else(|| {
                    crate::errors::invariant("validated beta-binomial parameters")
                })?;
                Ok(Self::BetaBinomial(
                    beta_binomial::Row {
                        design: row.design,
                        successes: &row.counts,
                        trials: (1..parameters.n_cols())
                            .map(|j| *parameters.at(gene, j))
                            .collect(),
                        weights: &row.weights,
                        concentration: *parameters.at(gene, 0),
                        prior: row.prior,
                    },
                    config,
                ))
            }
        }
    }

    fn optimize(
        &self,
        p: usize,
        options: &NegativeBinomialPosteriorOptions,
        shared_initial: Option<&[f64]>,
        warnings: &mut PosteriorWarnings,
    ) -> Result<MapFit, Error> {
        let alternating_start = |magnitude: f64| {
            (0..p)
                .map(|j| if j % 2 == 0 { magnitude } else { -magnitude })
                .collect::<Vec<_>>()
        };
        let compiled_start =
            || shared_initial.map_or_else(|| alternating_start(0.1), <[f64]>::to_vec);
        let optimization = match self {
            Self::CompiledNegativeBinomial(row) => {
                let initial = compiled_start();
                let mut optimization = optimize(row, initial.clone())?;
                if optimization.convergence < 0 {
                    warnings.primary.push(optimization.message.clone());
                }
                let opposite_start_disagrees = if options.check_opposite_start {
                    let second = optimize(row, initial.into_iter().map(|v| -v).collect())?;
                    if second.convergence < 0 {
                        warnings.opposite.push(second.message.clone());
                    }
                    optimization
                        .par
                        .iter()
                        .zip(second.par)
                        .any(|(a, b)| (a - b).abs() > 0.01)
                } else {
                    false
                };
                let refined = !options.map_only
                    && (optimization.convergence != 0 || opposite_start_disagrees);
                if refined {
                    optimization = refine(row, optimization.par, options)?;
                }
                return Ok(MapFit {
                    optimization,
                    opposite_start_disagrees,
                    refined,
                });
            }
            Self::General(row) => row.optimize(alternating_start(1.0))?,
            Self::BetaBinomial(row, config) => {
                let optimization = if config.method == RegressionMethod::BetaBinomialR {
                    row.optimize(vec![0.0; p], config)?
                } else {
                    row.prefit(compiled_start())?
                };
                if optimization.convergence < 0 {
                    warnings.primary.push(optimization.message.clone());
                }
                optimization
            }
            Self::NegativeBinomialRegression(row, config) => {
                let mut start = alternating_start(1.0);
                if config.log_link {
                    let baseline = (0..row.counts.len())
                        .filter(|&i| {
                            (0..p).filter(|&j| *row.design.at(i, j) == 0.0).count() == p - 1
                        })
                        .collect::<Vec<_>>();
                    let mean = if baseline.is_empty() {
                        row.counts.iter().sum::<f64>() / row.counts.len() as f64
                    } else {
                        baseline.iter().map(|&i| row.counts[i]).sum::<f64>() / baseline.len() as f64
                    };
                    start[0] = if mean == 0.0 { 0.0 } else { mean.ln_p() };
                }
                regression::optimize_nb_r(row, config, start)?
            }
        };
        Ok(MapFit {
            optimization,
            opposite_start_disagrees: false,
            refined: false,
        })
    }
}

#[derive(Default)]
struct PosteriorWarnings {
    primary: Vec<String>,
    opposite: Vec<String>,
}

struct MapFit {
    optimization: OptimResult,
    opposite_start_disagrees: bool,
    refined: bool,
}

struct LaplaceFit {
    map: MapFit,
    covariance: Option<DMatrix<f64>>,
}

fn invert_hessian(mut hessian: DMatrix<f64>) -> Result<Option<DMatrix<f64>>, Error> {
    let p = hessian.nrows();
    match p {
        0..=3 => Ok(hessian.try_inverse_mut().then_some(hessian)),
        4 => {
            let mut small = Matrix4::from_column_slice(hessian.as_slice());
            if !small.try_inverse_mut() {
                return Ok(None);
            }
            hessian.as_mut_slice().copy_from_slice(small.as_slice());
            Ok(Some(hessian))
        }
        _ => {
            let mut inverse = super::try_dense_matrix(p, p, "posterior covariance", |_, _| 0.0)?;
            Ok(try_invert_to(hessian, &mut inverse).then_some(inverse))
        }
    }
}

fn invert_finite_hessian(hessian: DMatrix<f64>) -> Result<Option<DMatrix<f64>>, Error> {
    if hessian.iter().all(|v| v.is_finite()) {
        invert_hessian(hessian)
    } else {
        Ok(None)
    }
}

fn has_valid_covariance_diagonal(matrix: &DMatrix<f64>) -> bool {
    (0..matrix.nrows()).all(|i| {
        let variance = matrix[(i, i)];
        variance.is_finite() && variance > 0.0
    })
}

fn fit_laplace(
    likelihood: &RowLikelihood<'_>,
    options: &NegativeBinomialPosteriorOptions,
    mut map: MapFit,
) -> Result<LaplaceFit, Error> {
    let result = &mut map.optimization;
    let covariance = match likelihood {
        RowLikelihood::General(general) => {
            let hessian = general.hessian(&result.par, options.hessian_step)?;
            if let Some(error) = general.take_error() {
                return Err(error);
            }
            invert_finite_hessian(hessian)?
        }
        RowLikelihood::BetaBinomial(beta_binomial, config) => {
            let hessian = beta_binomial.hessian(&result.par, options.hessian_step)?;
            let compiled_refinement = config.method == RegressionMethod::BetaBinomialCR;
            let finite = hessian.iter().all(|v| v.is_finite());
            // The compiled policy tests the inverse's diagonal before discarding
            // nonfinite Hessians. Keep that retry decision and reuse a valid inverse.
            let inverse = if compiled_refinement || finite {
                invert_hessian(hessian)?
            } else {
                None
            };
            if compiled_refinement
                && (result.convergence != 0
                    || !inverse.as_ref().is_some_and(has_valid_covariance_diagonal))
            {
                *result = beta_binomial.optimize(std::mem::take(&mut result.par), config)?;
                invert_finite_hessian(beta_binomial.hessian(&result.par, options.hessian_step)?)?
            } else {
                inverse.filter(|_| finite)
            }
        }
        RowLikelihood::NegativeBinomialRegression(row, _) => {
            invert_finite_hessian(row.hessian(&result.par, options.hessian_step)?)?
        }
        RowLikelihood::CompiledNegativeBinomial(row) => {
            let covariance =
                invert_finite_hessian(row.hessian(&result.par, options.hessian_step)?)?;
            // Compiled prefits with nonpositive variance use deferred refinement.
            if !map.refined
                && covariance
                    .as_ref()
                    .is_some_and(|m| (0..m.nrows()).any(|i| m[(i, i)] <= 0.0))
            {
                *result = refine(row, std::mem::take(&mut result.par), options)?;
                invert_finite_hessian(row.hessian(&result.par, options.hessian_step)?)?
            } else {
                covariance
            }
        }
    };
    if covariance.is_none() && !matches!(likelihood, RowLikelihood::CompiledNegativeBinomial(_)) {
        return Err(scoped_invalid_options(
            ERROR_SCOPE,
            "general posterior Hessian is singular or nonfinite",
        ));
    }
    Ok(LaplaceFit {
        map,
        covariance: covariance.filter(has_valid_covariance_diagonal),
    })
}

struct PosteriorRates {
    false_sign: Option<Vec<Option<f64>>>,
    false_sign_or_small: Option<Vec<Option<f64>>>,
}

impl PosteriorRates {
    fn new(options: &NegativeBinomialPosteriorOptions, genes: usize) -> Self {
        Self {
            false_sign: options.coefficient.map(|_| vec![None; genes]),
            false_sign_or_small: options
                .threshold
                .filter(|_| options.coefficient.is_some())
                .map(|_| vec![None; genes]),
        }
    }

    fn record<'a>(
        &mut self,
        analysis: &mut PosteriorAnalysis,
        query: AnalysisRow<'a>,
        grid_density: impl FnOnce() -> GridPosteriorRow<'a>,
    ) -> Result<(), Error> {
        let gene = query.gene;
        if let Some(coefficient) = query.coefficient {
            let mean = query.mode[coefficient];
            let error = query.sd[coefficient];
            self.set_false_sign(gene, Some(normal_upper_tail(mean.abs() / error)))?;
            if let Some(rates) = &mut self.false_sign_or_small {
                let threshold = query.threshold.ok_or_else(|| {
                    crate::errors::invariant("threshold set when threshold rates are requested")
                })?;
                let cutoff = if query.flip_sign {
                    (if mean == 0.0 { 0.0 } else { mean.signum() }) * threshold
                } else {
                    threshold
                };
                let z = (cutoff - mean) / error;
                rates[gene] = Some(if query.flip_sign && mean < 0.0 {
                    normal_upper_tail(z)
                } else {
                    normal_upper_tail(-z)
                });
            }
        }
        if query.options.interval.is_some() || query.options.contrasts.is_some() {
            let work = &query.options.work;
            let density = grid_density();
            let mut grid_evaluations = 0usize;
            if let Some(summary) = analysis.row(query, |beta, shift| {
                if grid_evaluations.is_multiple_of(256) {
                    work.check()?;
                }
                grid_evaluations += 1;
                density.value(beta, shift)
            })? {
                self.set_false_sign(gene, summary.false_sign_rate)?;
                if let Some(rates) = &mut self.false_sign_or_small {
                    rates[gene] = summary.threshold_rate;
                }
            }
        }
        Ok(())
    }

    fn set_false_sign(&mut self, gene: usize, rate: Option<f64>) -> Result<(), Error> {
        self.false_sign.as_mut().ok_or_else(|| {
            crate::errors::invariant("false-sign rates allocated for a selected coefficient")
        })?[gene] = rate;
        Ok(())
    }
}

pub(super) fn fit_batch(
    input: NegativeBinomialPosteriorInput<'_>,
    options: &NegativeBinomialPosteriorOptions,
    regression: RegressionRequest<'_, '_>,
    mut random: Option<&mut super::random::PosteriorRng>,
) -> Result<NegativeBinomialPosterior, Error> {
    let mut prepared = PreparedBatch::new(&input, options, regression)?;
    let prior = &prepared.prior;
    let g = input.counts.n_rows();
    let n = input.counts.n_cols();
    let p = input.design.n_cols();
    let na = super::super::location_missing::NA;
    let mut modes = RowMajorMatrix::from_elem(g, p, na)?.into_values();
    let mut sd = if options.map_only {
        None
    } else {
        Some(RowMajorMatrix::from_elem(g, p, na)?.into_values())
    };
    let mut rates = PosteriorRates::new(options, g);
    let mut diagnostics = Vec::with_capacity(g);
    let mut warnings = PosteriorWarnings::default();
    let shared_initial = random.as_deref_mut().map(|rng| rng.start_vector(p));
    let skip_zero = regression.is_none()
        || regression.is_some_and(|(_, config)| {
            config.method == RegressionMethod::NegativeBinomialR && config.log_link
        });
    for gene in 0..g {
        options.analysis.work.check()?;
        let counts = (0..n)
            .map(|i| *input.counts.at(gene, i))
            .collect::<Vec<_>>();
        if skip_zero && counts.iter().all(|&v| v == 0.0) {
            if regression.is_none() && !options.map_only && !options.log_link {
                return Err(Error::PosteriorMissingPrefit { gene });
            }
            diagnostics.push(PosteriorMapDiagnostic {
                convergence: None,
                function_evaluations: 0,
                opposite_start_disagrees: false,
                laplace_valid: false,
            });
            continue;
        }
        let dispersion = input.dispersions[gene];
        if !dispersion.is_finite() || dispersion <= 0.0 || !(1.0 / dispersion).is_finite() {
            return Err(scoped_invalid_options(
                ERROR_SCOPE,
                "positive finite dispersions are required for nonzero rows",
            ));
        }
        let row = Row {
            design: input.design,
            counts,
            weights: (0..n)
                .map(|i| input.weights.map(|m| *m.at(gene, i)).unwrap_or(1.0))
                .collect(),
            offsets: (0..n)
                .map(|sample| prepared.offsets.log_offset(gene, sample))
                .collect(),
            size: 1.0 / dispersion,
            prior,
        };
        let likelihood = RowLikelihood::new(&row, regression, gene)?;
        let mut map = likelihood.optimize(p, options, shared_initial.as_deref(), &mut warnings)?;
        let mut laplace_valid = false;
        if let Some(sd) = &mut sd {
            let laplace = fit_laplace(&likelihood, options, map)?;
            map = laplace.map;
            if let Some(covariance) = laplace.covariance {
                laplace_valid = true;
                let row_sd = &mut sd[gene * p..(gene + 1) * p];
                for (j, value) in row_sd.iter_mut().enumerate() {
                    *value = covariance[(j, j)].sqrt();
                }
                rates.record(
                    &mut prepared.analysis,
                    AnalysisRow {
                        gene,
                        mode: &map.optimization.par,
                        covariance: &covariance,
                        sd: row_sd,
                        coefficient: options.coefficient,
                        options: &options.analysis,
                        threshold: options.threshold,
                        flip_sign: options.flip_sign,
                    },
                    || GridPosteriorRow::new(&row, regression, gene, dispersion, &options.analysis),
                )?;
            }
        }
        modes[gene * p..(gene + 1) * p].copy_from_slice(&map.optimization.par);
        diagnostics.push(PosteriorMapDiagnostic {
            convergence: Some(map.optimization.convergence),
            function_evaluations: map.optimization.counts.function,
            opposite_start_disagrees: map.opposite_start_disagrees,
            laplace_valid,
        });
    }
    let sign_svalues = rates.false_sign.as_ref().map(|v| svalues(v)).transpose()?;
    let threshold_svalues = rates
        .false_sign_or_small
        .as_ref()
        .map(|v| svalues(v))
        .transpose()?;
    warnings.primary.extend(warnings.opposite);
    Ok(NegativeBinomialPosterior {
        summary_coefficient: options.coefficient,
        mode: RowMajorMatrix::from_row_major(g, p, modes)?,
        standard_deviation: sd
            .map(|v| RowMajorMatrix::from_row_major(g, p, v))
            .transpose()?,
        false_sign_rate: rates.false_sign,
        sign_svalues,
        false_sign_or_small_rate: rates.false_sign_or_small,
        threshold_svalues,
        prior: prepared.prior,
        analysis: prepared.analysis,
        diagnostics,
        warnings: warnings.primary,
        random_state: random.as_ref().map(|rng| rng.r_seed()),
    })
}
