//! Batch orchestration; likelihood, prior and optimizers have no R dependency.
use super::{
    negative_binomial::Row,
    prior::{CauchyCoefficientPrior, estimate_cauchy_prior_variance, invalid},
};
use crate::{
    DeseqError,
    matrix::{MatrixView, RowMajorMatrix},
    probability::normal_upper_tail,
    shrinkage::svalues,
};
use rcompat_lbfgsb::{Bounds, OptimControl, optim_lbfgsb_with_gradient};

/// MLE scale supplied for empirical prior adaptation.
#[derive(Clone, Copy, Debug)]
pub enum CoefficientScale {
    Natural,
    Log2,
}
impl CoefficientScale {
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
    pub estimates: &'a [f64],
    pub standard_errors: &'a [f64],
    pub scale: CoefficientScale,
}
/// Log offsets or normalization values; scale conversion stays in the core.
#[derive(Clone, Copy, Debug, Default)]
pub enum PosteriorOffsets<'a> {
    #[default]
    None,
    Log(MatrixView<'a, f64>),
    Factors(MatrixView<'a, f64>),
    SizeFactors(&'a [f64]),
}
/// Borrowed primitive data for a negative-binomial posterior batch.
#[derive(Clone, Copy, Debug)]
pub struct NegativeBinomialPosteriorInput<'a> {
    pub counts: MatrixView<'a, f64>,
    pub design: MatrixView<'a, f64>,
    pub dispersions: &'a [f64],
    pub offsets: PosteriorOffsets<'a>,
    pub weights: Option<MatrixView<'a, f64>>,
    pub mle: Option<CoefficientMle<'a>>,
}
#[derive(Clone, Debug)]
pub struct NegativeBinomialPosteriorOptions {
    pub prior: CauchyCoefficientPrior,
    pub coefficient: Option<usize>,
    /// Positive magnitude in natural log units; None omits FSOS calculations.
    pub threshold: Option<f64>,
    pub flip_sign: bool,
    pub multiplier: f64,
    pub map_only: bool,
    pub check_opposite_start: bool,
    pub hessian_step: f64,
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
            hessian_step: 1e-3,
        }
    }
}
/// Per-row native convergence, separate from posterior-validity decisions.
#[derive(Clone, Debug)]
pub struct PosteriorMapDiagnostic {
    pub convergence: Option<i32>,
    pub function_evaluations: usize,
    pub opposite_start_disagrees: bool,
    pub laplace_valid: bool,
}
#[derive(Clone, Debug)]
pub struct NegativeBinomialPosterior {
    /// Row-major modes and SDs in natural log coefficient units.
    pub mode: RowMajorMatrix<f64>,
    pub standard_deviation: Option<RowMajorMatrix<f64>>,
    pub false_sign_rate: Option<Vec<Option<f64>>>,
    pub sign_svalues: Option<Vec<Option<f64>>>,
    pub false_sign_or_small_rate: Option<Vec<Option<f64>>>,
    pub threshold_svalues: Option<Vec<Option<f64>>>,
    pub prior: CauchyCoefficientPrior,
    pub diagnostics: Vec<PosteriorMapDiagnostic>,
}
/// Numerical columns for frontend presentation, with no frontend arithmetic.
#[derive(Clone, Debug)]
pub struct PosteriorCoefficientSummary {
    pub mode: Vec<f64>,
    pub standard_deviation: Vec<Option<f64>>,
    pub sign_svalues: Option<Vec<Option<f64>>>,
    pub threshold_svalues: Option<Vec<Option<f64>>>,
}
impl NegativeBinomialPosterior {
    pub fn coefficient_summary(
        &self,
        coefficient: usize,
        scale: CoefficientScale,
    ) -> Result<PosteriorCoefficientSummary, DeseqError> {
        if coefficient >= self.mode.n_cols() {
            return Err(invalid("coefficient is outside posterior matrix"));
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
fn optimize(row: &Row<'_>, start: Vec<f64>) -> Result<rcompat_lbfgsb::OptimResult, DeseqError> {
    let p = start.len();
    let scale = row.objective_gradient(&vec![0.0; p]).0.max(1.0);
    let (value, gradient) = row.objective_gradient(&start);
    if !scale.is_finite() || !value.is_finite() || gradient.iter().any(|v| !v.is_finite()) {
        return Err(invalid(
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
fn refine(row: &Row<'_>, start: Vec<f64>) -> Result<rcompat_lbfgsb::OptimResult, DeseqError> {
    let p = start.len();
    let scale = row.objective_gradient(&vec![0.0; p]).0.max(1.0);
    let mut control = OptimControl::default_for_dimension(p);
    control.maxit = 1000;
    control.factr = 1.0;
    control.pgtol = 1e-10;
    let bounds = Bounds::new(vec![f64::NEG_INFINITY; p], vec![f64::INFINITY; p])
        .map_err(|_| invalid("invalid unbounded MAP optimization shape"))?;
    optim_lbfgsb_with_gradient(
        start,
        bounds,
        |b| row.objective_gradient(b).0 / scale,
        |b| {
            row.objective_gradient(b)
                .1
                .into_iter()
                .map(|g| g / scale)
                .collect()
        },
        control,
    )
    .map_err(|e| DeseqError::UnsupportedFeature {
        feature: format!("negative-binomial MAP optimization failed: {e}"),
    })
}
/// Fit Cauchy-shrunk negative-binomial coefficients and Laplace summaries.
pub fn fit_negative_binomial_posterior(
    input: NegativeBinomialPosteriorInput<'_>,
    options: &NegativeBinomialPosteriorOptions,
) -> Result<NegativeBinomialPosterior, DeseqError> {
    fit_batch(input, options, None)
}
pub(super) fn fit_batch(
    input: NegativeBinomialPosteriorInput<'_>,
    options: &NegativeBinomialPosteriorOptions,
    regression: Option<(
        Option<MatrixView<'_, f64>>,
        &super::regression::RegressionPosteriorOptions,
    )>,
) -> Result<NegativeBinomialPosterior, DeseqError> {
    let g = input.counts.n_rows();
    let n = input.counts.n_cols();
    let p = input.design.n_cols();
    if g == 0 || input.design.n_rows() != n || input.dispersions.len() != g {
        return Err(invalid(
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
        return Err(invalid(
            "counts must be finite nonnegative numbers and design must be finite",
        ));
    }
    if let Some(m) = input.weights
        && (m.n_rows() != g
            || m.n_cols() != n
            || m.as_slice().iter().any(|v| !v.is_finite() || *v < 0.0))
    {
        return Err(invalid(
            "finite nonnegative weights must match counts shape",
        ));
    }
    let size_offsets = match input.offsets {
        PosteriorOffsets::SizeFactors(s) => {
            if s.len() != n || s.iter().any(|v| !v.is_finite() || *v <= 0.0) {
                return Err(invalid("positive sample-aligned size factors are required"));
            }
            Some(s.iter().map(|v| v.ln()).collect::<Vec<_>>())
        }
        PosteriorOffsets::Log(m) | PosteriorOffsets::Factors(m) => {
            let factors = matches!(input.offsets, PosteriorOffsets::Factors(_));
            if m.n_rows() != g
                || m.n_cols() != n
                || m.as_slice()
                    .iter()
                    .any(|v| !v.is_finite() || (factors && *v <= 0.0))
            {
                return Err(invalid(
                    "finite log offsets or positive normalization factors must match counts",
                ));
            }
            None
        }
        PosteriorOffsets::None => None,
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
        return Err(invalid(
            "invalid coefficient, multiplier, Hessian step or threshold",
        ));
    }
    let mut prior = options.prior.clone();
    if let Some(mle) = input.mle {
        let coefficient = options
            .coefficient
            .ok_or_else(|| invalid("MLE prior adaptation requires a coefficient"))?;
        if mle.estimates.len() != g || mle.standard_errors.len() != g {
            return Err(invalid("MLE population must match count rows"));
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
    let na = super::super::location_missing::NA;
    let mut modes = RowMajorMatrix::from_elem(g, p, na)?.into_values();
    let mut sd = if options.map_only {
        None
    } else {
        Some(RowMajorMatrix::from_elem(g, p, na)?.into_values())
    };
    let mut fsr = options.coefficient.map(|_| vec![None; g]);
    let mut fsos = options
        .threshold
        .filter(|_| options.coefficient.is_some())
        .map(|_| vec![None; g]);
    let mut diagnostics = Vec::with_capacity(g);
    for gene in 0..g {
        let counts = (0..n)
            .map(|i| *input.counts.at(gene, i))
            .collect::<Vec<_>>();
        let nb_r = regression.is_some_and(|(_, o)| {
            o.method == super::regression::RegressionMethod::NegativeBinomialR
        });
        let skip_zero = regression.is_none()
            || regression.is_some_and(|(_, o)| {
                o.method == super::regression::RegressionMethod::NegativeBinomialR && o.log_link
            });
        if skip_zero && counts.iter().all(|&v| v == 0.0) {
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
            return Err(invalid(
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
                .map(|i| match input.offsets {
                    PosteriorOffsets::None => 0.0,
                    PosteriorOffsets::Log(m) => *m.at(gene, i),
                    PosteriorOffsets::Factors(m) => m.at(gene, i).ln(),
                    PosteriorOffsets::SizeFactors(_) => size_offsets.as_ref().unwrap()[i],
                })
                .collect(),
            size: 1.0 / dispersion,
            prior: &prior,
        };
        let initial = (0..p)
            .map(|j| if j % 2 == 0 { 0.1 } else { -0.1 })
            .collect::<Vec<_>>();
        let general = regression
            .filter(|(_, o)| o.method == super::regression::RegressionMethod::General)
            .map(|(parameters, config)| super::regression::GeneralRow {
                y: &row.counts,
                design: input.design,
                parameters: parameters.map(|m| (0..m.n_cols()).map(|j| *m.at(gene, j)).collect()),
                weights: &row.weights,
                offsets: &row.offsets,
                prior: &prior,
                options: config,
                error: std::cell::RefCell::new(None),
            });
        let mut result = if let Some(general) = &general {
            general.optimize(
                (0..p)
                    .map(|j| if j % 2 == 0 { 1.0 } else { -1.0 })
                    .collect(),
            )?
        } else {
            if nb_r {
                let (_, config) = regression.unwrap();
                let mut start = (0..p)
                    .map(|j| if j % 2 == 0 { 1.0 } else { -1.0 })
                    .collect::<Vec<_>>();
                if config.log_link {
                    let baseline = (0..n)
                        .filter(|&i| {
                            (0..p).filter(|&j| *input.design.at(i, j) == 0.0).count() == p - 1
                        })
                        .collect::<Vec<_>>();
                    let mean = if baseline.is_empty() {
                        row.counts.iter().sum::<f64>() / n as f64
                    } else {
                        baseline.iter().map(|&i| row.counts[i]).sum::<f64>() / baseline.len() as f64
                    };
                    start[0] = if mean == 0.0 { 0.0 } else { mean.ln() };
                }
                super::regression::optimize_nb_r(&row, config, start)?
            } else {
                optimize(&row, initial.clone())?
            }
        };
        let disagreement = if regression.is_none() && options.check_opposite_start {
            let second = optimize(&row, initial.into_iter().map(|v| -v).collect())?;
            result
                .par
                .iter()
                .zip(second.par)
                .any(|(a, b)| (a - b).abs() > 0.01)
        } else {
            false
        };
        if regression.is_none() && !options.map_only && (result.convergence != 0 || disagreement) {
            result = refine(&row, result.par)?;
        }
        for (j, &value) in result.par.iter().enumerate() {
            modes.as_mut_slice()[gene * p + j] = value;
        }
        let mut laplace_valid = false;
        if let Some(sd) = &mut sd {
            let h = if let Some(general) = &general {
                general.hessian(&result.par, options.hessian_step)
            } else {
                row.hessian(&result.par, options.hessian_step)
            };
            if let Some(e) = general.as_ref().and_then(|g| g.take_error()) {
                return Err(e);
            }
            if regression.is_some()
                && (!h.iter().all(|v| v.is_finite()) || h.clone().try_inverse().is_none())
            {
                return Err(invalid(
                    "general posterior Hessian is singular or nonfinite",
                ));
            }
            if h.iter().all(|v| v.is_finite())
                && let Some(covariance) = h.try_inverse()
                && covariance
                    .diagonal()
                    .iter()
                    .all(|v| v.is_finite() && *v > 0.0)
            {
                laplace_valid = true;
                for j in 0..p {
                    sd.as_mut_slice()[gene * p + j] = covariance[(j, j)].sqrt();
                }
                if let Some(c) = options.coefficient {
                    let mean = result.par[c];
                    let error = sd.as_slice()[gene * p + c];
                    fsr.as_mut().unwrap()[gene] = Some(normal_upper_tail(mean.abs() / error));
                    if let Some(rates) = &mut fsos {
                        let cutoff = if options.flip_sign {
                            (if mean == 0.0 { 0.0 } else { mean.signum() })
                                * options.threshold.unwrap()
                        } else {
                            options.threshold.unwrap()
                        };
                        let z = (cutoff - mean) / error;
                        rates[gene] = Some(if options.flip_sign && mean < 0.0 {
                            normal_upper_tail(z)
                        } else {
                            normal_upper_tail(-z)
                        });
                    }
                }
            }
        }
        diagnostics.push(PosteriorMapDiagnostic {
            convergence: Some(result.convergence),
            function_evaluations: result.counts.function,
            opposite_start_disagrees: disagreement,
            laplace_valid,
        });
    }
    let sign_svalues = fsr.as_ref().map(|v| svalues(v)).transpose()?;
    let threshold_svalues = fsos.as_ref().map(|v| svalues(v)).transpose()?;
    Ok(NegativeBinomialPosterior {
        mode: RowMajorMatrix::from_row_major(g, p, modes)?,
        standard_deviation: sd
            .map(|v| RowMajorMatrix::from_row_major(g, p, v))
            .transpose()?,
        false_sign_rate: fsr,
        sign_svalues,
        false_sign_or_small_rate: fsos,
        threshold_svalues,
        prior,
        diagnostics,
    })
}
