//! Optional posterior queries share row covariance; no G-by-p-by-p retention.
use super::{grid, prior::ERROR_SCOPE, regression::RegressionLikelihood};
use crate::errors::{
    ErrorKind, PosteriorArray, PosteriorArrayAxis, PosteriorControl, missing_prerequisite,
    scoped_invalid_options,
};
use crate::{Error, matrix::RowMajorMatrix, probability::normal_upper_quantile};
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
/// How credible intervals are computed.
pub enum PosteriorIntervalKind {
    #[default]
    /// From the Laplace (normal) approximation.
    Laplace,
    /// Highest-density interval on a posterior grid.
    HighestDensity,
    /// Equal-tailed interval on a posterior grid.
    EqualTail,
}
#[derive(Clone, Debug)]
/// Credible-interval controls.
pub struct PosteriorIntervalOptions {
    /// Interval type.
    pub kind: PosteriorIntervalKind,
    /// Credible level, for example 0.95.
    pub level: f64,
    /// Grid cells per side for grid intervals.
    pub grid_intervals: usize,
    /// Grid half-width in posterior standard deviations.
    pub grid_standard_deviations: f64,
    /// Integration points per nuisance coefficient.
    pub nuisance_points: usize,
    /// Nuisance integration half-width in standard deviations.
    pub nuisance_standard_deviations: f64,
}
impl Default for PosteriorIntervalOptions {
    fn default() -> Self {
        Self {
            kind: Default::default(),
            level: 0.95,
            grid_intervals: 50,
            grid_standard_deviations: 5.0,
            nuisance_points: 5,
            nuisance_standard_deviations: 2.0,
        }
    }
}
#[derive(Clone, Debug, Default)]
/// Optional posterior summaries.
pub struct PosteriorAnalysisOptions {
    /// Batch-wide resource controls; numerical grid policy is unchanged.
    pub work: super::work::PosteriorWorkControl,
    /// Credible intervals, if wanted.
    pub interval: Option<PosteriorIntervalOptions>,
    /// Coefficients by contrasts, in natural coefficient units.
    pub contrasts: Option<RowMajorMatrix<f64>>,
    /// One uncertainty SD per observation row, for a scalar likelihood parameter.
    pub parameter_standard_errors: Option<Vec<f64>>,
    /// An explicit grid likelihood may differ from the specialized MAP fitter.
    pub grid_likelihood: Option<RegressionLikelihood>,
}
#[derive(Clone, Debug)]
/// Diagnostics of one row's posterior grid.
pub struct PosteriorGridDiagnostic {
    /// Whether the grid density decreases away from the mode.
    pub monotone: bool,
    /// Probability mass estimated beyond the left end.
    pub outside_left: Option<f64>,
    /// Probability mass estimated beyond the right end.
    pub outside_right: Option<f64>,
    /// Likelihood evaluations.
    pub evaluations: usize,
}
#[derive(Clone, Debug)]
/// Optional posterior summaries.
pub struct PosteriorAnalysis {
    /// Lower and upper bounds, in natural coefficient units.
    pub intervals: Option<RowMajorMatrix<f64>>,
    /// Posterior modes of the requested contrasts.
    pub contrast_mode: Option<RowMajorMatrix<f64>>,
    /// Their Laplace standard deviations.
    pub contrast_standard_deviation: Option<RowMajorMatrix<f64>>,
    /// Per-row grid diagnostics, when a grid was used.
    pub grid_diagnostics: Vec<Option<PosteriorGridDiagnostic>>,
}
impl PosteriorAnalysisOptions {
    /// Preflight work for an active grid query, including optional parameter integration.
    /// Laplace/MAP-only or unselected queries perform no grid evaluations.
    pub fn estimated_grid_evaluations(
        &self,
        genes: usize,
        coefficients: usize,
        coefficient: Option<usize>,
        map_only: bool,
    ) -> Result<usize, Error> {
        let Some(o) = self.interval.as_ref().filter(|o| {
            o.kind != PosteriorIntervalKind::Laplace && !map_only && coefficient.is_some()
        }) else {
            return Ok(0);
        };
        if coefficients == 0 || coefficient.is_some_and(|c| c >= coefficients) {
            return Err(scoped_invalid_options(
                ERROR_SCOPE,
                "posterior coefficient is outside the model",
            ));
        }
        let overflow =
            |operation, reason| self.grid_overflow(o, genes, coefficients, operation, reason);
        let count_overflow =
            |operation| overflow(operation, "posterior grid evaluation count overflows");
        let power = u32::try_from(coefficients - 1)
            .map_err(|_| overflow("nuisance exponent", "posterior dimensions overflow"))?;
        let points = o
            .grid_intervals
            .checked_add(1)
            .ok_or_else(|| count_overflow("interval grid points"))?;
        let nuisance = o
            .nuisance_points
            .checked_pow(power)
            .ok_or_else(|| count_overflow("nuisance grid points"))?;
        self.row_grid_evaluations(o, points, nuisance)
            .ok_or_else(|| count_overflow("row evaluations"))?
            .checked_mul(genes)
            .ok_or_else(|| count_overflow("batch evaluations"))
    }

    pub(super) fn validate(
        &self,
        genes: usize,
        p: usize,
        coefficient: Option<usize>,
        threshold: Option<f64>,
        map_only: bool,
    ) -> Result<(), Error> {
        if let Some(c) = &self.contrasts
            && c.n_rows() != p
        {
            return Err(Error::PosteriorShapeMismatch {
                array: PosteriorArray::Contrasts,
                axis: PosteriorArrayAxis::Rows,
                expected: p,
                actual: c.n_rows(),
                reason: "posterior contrasts must have one row per coefficient",
            });
        }
        if let Some(c) = &self.contrasts {
            c.validate_finite("posterior contrasts")?;
        }
        if let Some(s) = &self.parameter_standard_errors {
            if let Some(threshold) = threshold {
                return Err(Error::InvalidPosteriorControl {
                    control: PosteriorControl::ThresholdWithUncertainty(threshold),
                });
            }
            const UNCERTAINTIES: &str =
                "parameter uncertainties must be finite nonnegative and gene aligned";
            if s.len() != genes {
                return Err(Error::PosteriorShapeMismatch {
                    array: PosteriorArray::ParameterStandardErrors,
                    axis: PosteriorArrayAxis::Rows,
                    expected: genes,
                    actual: s.len(),
                    reason: UNCERTAINTIES,
                });
            }
            if let Some((gene, &value)) = s
                .iter()
                .enumerate()
                .find(|(_, value)| !value.is_finite() || **value < 0.0)
            {
                return Err(Error::InvalidGeneValue {
                    domain: ErrorKind::Options,
                    operation: "posterior parameter standard error",
                    gene,
                    sample: None,
                    value,
                    reason: format!("{ERROR_SCOPE}: {UNCERTAINTIES}"),
                });
            }
            if !map_only
                && coefficient.is_some()
                && self
                    .interval
                    .as_ref()
                    .is_none_or(|o| o.kind == PosteriorIntervalKind::Laplace)
            {
                return Err(missing_prerequisite(
                    ErrorKind::Options,
                    format!("{ERROR_SCOPE}: parameter uncertainty"),
                    "a grid interval",
                ));
            }
        }
        if let Some(o) = &self.interval {
            if !o.level.is_finite() || !(0.0..=1.0).contains(&o.level) {
                return Err(Error::InvalidPosteriorControl {
                    control: PosteriorControl::IntervalLevel(o.level),
                });
            }
            if !map_only && coefficient.is_some() && o.kind != PosteriorIntervalKind::Laplace {
                if let Some(control) = invalid_grid_control(o) {
                    return Err(Error::InvalidPosteriorControl { control });
                }
                let overflow =
                    |operation, reason| self.grid_overflow(o, genes, p, operation, reason);
                let points = o.grid_intervals.checked_add(1).ok_or_else(|| {
                    overflow("interval grid points", "posterior grid dimensions overflow")
                })?;
                let power = p
                    .checked_sub(1)
                    .and_then(|n| u32::try_from(n).ok())
                    .ok_or_else(|| {
                        overflow("nuisance exponent", "posterior dimensions overflow")
                    })?;
                let nuisance = o.nuisance_points.checked_pow(power).ok_or_else(|| {
                    overflow(
                        "nuisance grid points",
                        "posterior nuisance grid dimensions overflow",
                    )
                })?;
                self.row_grid_evaluations(o, points, nuisance)
                    .ok_or_else(|| {
                        overflow(
                            "row evaluations",
                            "posterior grid evaluation count overflows",
                        )
                    })?;
            }
        }
        self.work
            .validate_estimate(self.estimated_grid_evaluations(
                genes,
                p,
                coefficient,
                map_only,
            )?)?;
        Ok(())
    }

    fn row_grid_evaluations(
        &self,
        interval: &PosteriorIntervalOptions,
        points: usize,
        nuisance: usize,
    ) -> Option<usize> {
        points
            .checked_mul(nuisance)?
            .checked_mul(if self.parameter_standard_errors.is_some() {
                interval.nuisance_points
            } else {
                1
            })
    }

    fn grid_overflow(
        &self,
        interval: &PosteriorIntervalOptions,
        genes: usize,
        coefficients: usize,
        operation: &'static str,
        reason: &'static str,
    ) -> Error {
        Error::PosteriorGridOverflow {
            operation,
            grid_intervals: interval.grid_intervals,
            nuisance_points: interval.nuisance_points,
            coefficients,
            genes,
            parameter_uncertainty: self.parameter_standard_errors.is_some(),
            reason,
        }
    }
}

fn invalid_grid_control(o: &PosteriorIntervalOptions) -> Option<PosteriorControl> {
    if o.level == 0.0 || o.level == 1.0 {
        Some(PosteriorControl::GridLevel(o.level))
    } else if o.grid_intervals < 2 || !o.grid_intervals.is_multiple_of(2) {
        Some(PosteriorControl::GridIntervals(o.grid_intervals))
    } else if !o.grid_standard_deviations.is_finite() || o.grid_standard_deviations <= 0.0 {
        Some(PosteriorControl::GridStandardDeviations(
            o.grid_standard_deviations,
        ))
    } else if o.nuisance_points == 0 {
        Some(PosteriorControl::NuisancePoints(o.nuisance_points))
    } else if !o.nuisance_standard_deviations.is_finite() || o.nuisance_standard_deviations < 0.0 {
        Some(PosteriorControl::NuisanceStandardDeviations(
            o.nuisance_standard_deviations,
        ))
    } else {
        None
    }
}
impl PosteriorAnalysis {
    pub(super) fn new(
        options: &PosteriorAnalysisOptions,
        genes: usize,
        coefficient: Option<usize>,
        map_only: bool,
    ) -> Result<Self, Error> {
        let na = f64::from_bits(0x7ff0_0000_0000_07a2);
        let intervals = if options.interval.is_some() && coefficient.is_some() {
            Some(RowMajorMatrix::from_row_major(
                genes,
                2,
                vec![
                    na;
                    genes.checked_mul(2).ok_or_else(|| scoped_invalid_options(
                        ERROR_SCOPE,
                        "posterior interval dimensions overflow"
                    ))?
                ],
            )?)
        } else {
            None
        };
        let contrast = options
            .contrasts
            .as_ref()
            .map(|c| {
                RowMajorMatrix::from_row_major(
                    genes,
                    c.n_cols(),
                    vec![
                        na;
                        genes
                            .checked_mul(c.n_cols())
                            .ok_or_else(|| scoped_invalid_options(
                                ERROR_SCOPE,
                                "posterior contrast dimensions overflow"
                            ))?
                    ],
                )
            })
            .transpose()?;
        Ok(Self {
            intervals,
            contrast_standard_deviation: contrast.clone(),
            contrast_mode: contrast,
            grid_diagnostics: if !map_only
                && options
                    .interval
                    .as_ref()
                    .is_some_and(|o| o.kind != PosteriorIntervalKind::Laplace)
            {
                vec![None; genes]
            } else {
                vec![]
            },
        })
    }
    // Invariant: grid columns come from the posterior coefficient matrix shape.
    #[allow(clippy::unwrap_used, clippy::expect_used)]
    pub(super) fn row(
        &mut self,
        row: AnalysisRow<'_>,
        mut log_posterior: impl FnMut(&[f64], f64) -> Result<f64, Error>,
    ) -> Result<Option<grid::GridSummary>, Error> {
        let AnalysisRow {
            gene,
            mode,
            covariance,
            sd,
            coefficient,
            options,
            threshold,
            flip_sign,
        } = row;
        if let Some(c) = &options.contrasts {
            let k = c.n_cols();
            for column in 0..k {
                let weights = (0..mode.len())
                    .map(|j| *c.get(j, column).unwrap())
                    .collect::<Vec<_>>();
                let mean = mode.iter().zip(&weights).map(|(m, w)| m * w).sum::<f64>();
                let variance = (0..mode.len())
                    .map(|j| {
                        weights[j]
                            * (0..mode.len())
                                .map(|i| weights[i] * covariance[(i, j)])
                                .sum::<f64>()
                    })
                    .sum::<f64>();
                self.contrast_mode
                    .as_mut()
                    .ok_or_else(|| crate::errors::invariant("contrast outputs allocated"))?
                    .as_mut_slice()[gene * k + column] = mean;
                self.contrast_standard_deviation
                    .as_mut()
                    .ok_or_else(|| crate::errors::invariant("contrast outputs allocated"))?
                    .as_mut_slice()[gene * k + column] = variance.sqrt();
            }
        }
        let Some(coefficient) = coefficient else {
            return Ok(None);
        };
        let Some(interval) = &options.interval else {
            return Ok(None);
        };
        if interval.kind == PosteriorIntervalKind::Laplace {
            let z = normal_upper_quantile((1.0 - interval.level) / 2.0);
            self.intervals
                .as_mut()
                .ok_or_else(|| crate::errors::invariant("interval outputs allocated"))?
                .as_mut_slice()[gene * 2..gene * 2 + 2]
                .copy_from_slice(&[
                    mode[coefficient] - z * sd[coefficient],
                    mode[coefficient] + z * sd[coefficient],
                ]);
            return Ok(None);
        }
        let summary = grid::evaluate(
            grid::GridInput {
                mode,
                sd,
                coefficient,
                options: interval,
                parameter_sd: options.parameter_standard_errors.as_ref().map(|s| s[gene]),
                threshold,
                flip_sign,
            },
            &mut log_posterior,
        )?;
        self.intervals
            .as_mut()
            .ok_or_else(|| crate::errors::invariant("interval outputs allocated"))?
            .as_mut_slice()[gene * 2..gene * 2 + 2]
            .copy_from_slice(&summary.interval);
        self.grid_diagnostics[gene] = Some(summary.diagnostic.clone());
        Ok(Some(summary))
    }
}

pub(super) struct AnalysisRow<'a> {
    pub gene: usize,
    pub mode: &'a [f64],
    pub covariance: &'a nalgebra::DMatrix<f64>,
    pub sd: &'a [f64],
    pub coefficient: Option<usize>,
    pub options: &'a PosteriorAnalysisOptions,
    pub threshold: Option<f64>,
    pub flip_sign: bool,
}
/// Only grid queries construct the generic log density. MAP/Laplace arithmetic
/// stays on its existing specialized path, including the original prior policy.
pub(super) struct GridPosteriorRow<'a> {
    row: &'a super::negative_binomial::Row<'a>,
    options: super::regression::RegressionPosteriorOptions,
    parameters: Option<Vec<f64>>,
}
impl<'a> GridPosteriorRow<'a> {
    pub fn new(
        row: &'a super::negative_binomial::Row<'a>,
        regression: Option<(
            Option<crate::matrix::MatrixView<'_, f64>>,
            &super::regression::RegressionPosteriorOptions,
        )>,
        gene: usize,
        dispersion: f64,
        analysis: &PosteriorAnalysisOptions,
    ) -> Self {
        // Copy only density policy. Cloning the full fitting options would
        // copy the gene-aligned parameter-SD vector once per gene (O(G^2)).
        let mut options = super::regression::RegressionPosteriorOptions::default();
        if let Some((_, original)) = regression {
            options.likelihood = original.likelihood.clone();
            options.prior_df = original.prior_df;
            options.prior_noncentrality = original.prior_noncentrality;
            options.nuisance_mean = original.nuisance_mean;
        }
        if let Some(likelihood) = &analysis.grid_likelihood {
            options.likelihood = likelihood.clone();
        }
        let parameters = if let Some((parameters, _)) = regression {
            parameters.map(|m| (0..m.n_cols()).map(|j| *m.at(gene, j)).collect::<Vec<_>>())
        } else {
            Some(vec![dispersion])
        };
        Self {
            row,
            options,
            parameters,
        }
    }
    pub fn value(&self, beta: &[f64], shift: f64) -> Result<f64, Error> {
        let general = super::regression::GeneralRow {
            y: &self.row.counts,
            design: self.row.design,
            parameters: None,
            weights: &self.row.weights,
            offsets: &self.row.offsets,
            prior: self.row.prior,
            options: &self.options,
            error: std::cell::RefCell::new(None),
        };
        let shifted = if shift != 0.0 {
            Some([self
                .parameters
                .as_ref()
                .ok_or_else(|| scoped_invalid_options(ERROR_SCOPE, "missing grid parameter"))?[0]
                + shift])
        } else {
            None
        };
        Ok(-general.value_at_checked(
            beta,
            shifted
                .as_ref()
                .map(|v| &v[..])
                .or(self.parameters.as_deref()),
        )?)
    }
}
