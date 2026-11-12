//! Optional posterior queries share row covariance; no G-by-p-by-p retention.
use super::{grid, prior::invalid, regression::RegressionLikelihood};
use crate::{DeseqError, matrix::RowMajorMatrix, probability::normal_upper_quantile};
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum PosteriorIntervalKind {
    #[default]
    Laplace,
    HighestDensity,
    EqualTail,
}
#[derive(Clone, Debug)]
pub struct PosteriorIntervalOptions {
    pub kind: PosteriorIntervalKind,
    pub level: f64,
    pub grid_intervals: usize,
    pub grid_standard_deviations: f64,
    pub nuisance_points: usize,
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
pub struct PosteriorAnalysisOptions {
    /// Batch-wide resource controls; numerical grid policy is unchanged.
    pub work: super::work::PosteriorWorkControl,
    pub interval: Option<PosteriorIntervalOptions>,
    /// Coefficients by contrasts, in natural coefficient units.
    pub contrasts: Option<RowMajorMatrix<f64>>,
    /// One uncertainty SD per observation row, for a scalar likelihood parameter.
    pub parameter_standard_errors: Option<Vec<f64>>,
    /// An explicit grid likelihood may differ from the specialized MAP fitter.
    pub grid_likelihood: Option<RegressionLikelihood>,
}
#[derive(Clone, Debug)]
pub struct PosteriorGridDiagnostic {
    pub monotone: bool,
    pub outside_left: Option<f64>,
    pub outside_right: Option<f64>,
    pub evaluations: usize,
}
#[derive(Clone, Debug)]
pub struct PosteriorAnalysis {
    /// Lower and upper bounds, in natural coefficient units.
    pub intervals: Option<RowMajorMatrix<f64>>,
    pub contrast_mode: Option<RowMajorMatrix<f64>>,
    pub contrast_standard_deviation: Option<RowMajorMatrix<f64>>,
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
    ) -> Result<usize, DeseqError> {
        let Some(o) = self.interval.as_ref().filter(|o| {
            o.kind != PosteriorIntervalKind::Laplace && !map_only && coefficient.is_some()
        }) else {
            return Ok(0);
        };
        if coefficients == 0 || coefficient.is_some_and(|c| c >= coefficients) {
            return Err(invalid("posterior coefficient is outside the model"));
        }
        let power = u32::try_from(coefficients - 1)
            .map_err(|_| invalid("posterior dimensions overflow"))?;
        o.grid_intervals
            .checked_add(1)
            .and_then(|points| {
                o.nuisance_points
                    .checked_pow(power)
                    .and_then(|n| points.checked_mul(n))
            })
            .and_then(|n| {
                n.checked_mul(if self.parameter_standard_errors.is_some() {
                    o.nuisance_points
                } else {
                    1
                })
            })
            .and_then(|n| n.checked_mul(genes))
            .ok_or_else(|| invalid("posterior grid evaluation count overflows"))
    }

    pub(super) fn validate(
        &self,
        genes: usize,
        p: usize,
        coefficient: Option<usize>,
        threshold: Option<f64>,
        map_only: bool,
    ) -> Result<(), DeseqError> {
        if let Some(c) = &self.contrasts
            && c.n_rows() != p
        {
            return Err(invalid(
                "posterior contrasts must have one row per coefficient",
            ));
        }
        if let Some(c) = &self.contrasts {
            c.validate_finite("posterior contrasts")?;
        }
        if let Some(s) = &self.parameter_standard_errors {
            if threshold.is_some() {
                return Err(invalid(
                    "parameter uncertainty and a threshold cannot be combined",
                ));
            }
            if s.len() != genes || s.iter().any(|v| !v.is_finite() || *v < 0.0) {
                return Err(invalid(
                    "parameter uncertainties must be finite nonnegative and gene aligned",
                ));
            }
            if !map_only
                && coefficient.is_some()
                && self
                    .interval
                    .as_ref()
                    .is_none_or(|o| o.kind == PosteriorIntervalKind::Laplace)
            {
                return Err(invalid("parameter uncertainty requires a grid interval"));
            }
        }
        if let Some(o) = &self.interval {
            if !o.level.is_finite() || !(0.0..=1.0).contains(&o.level) {
                return Err(invalid(
                    "posterior interval level must lie between zero and one",
                ));
            }
            if !map_only && coefficient.is_some() && o.kind != PosteriorIntervalKind::Laplace {
                if o.level == 0.0
                    || o.level == 1.0
                    || o.grid_intervals < 2
                    || o.grid_intervals % 2 != 0
                    || !o.grid_standard_deviations.is_finite()
                    || o.grid_standard_deviations <= 0.0
                    || o.nuisance_points == 0
                    || !o.nuisance_standard_deviations.is_finite()
                    || o.nuisance_standard_deviations < 0.0
                {
                    return Err(invalid("invalid posterior integration grid controls"));
                }
                let points = o
                    .grid_intervals
                    .checked_add(1)
                    .ok_or_else(|| invalid("posterior grid dimensions overflow"))?;
                let nuisance = o
                    .nuisance_points
                    .checked_pow(
                        (p - 1)
                            .try_into()
                            .map_err(|_| invalid("posterior dimensions overflow"))?,
                    )
                    .ok_or_else(|| invalid("posterior nuisance grid dimensions overflow"))?;
                points
                    .checked_mul(nuisance)
                    .and_then(|n| {
                        n.checked_mul(if self.parameter_standard_errors.is_some() {
                            o.nuisance_points
                        } else {
                            1
                        })
                    })
                    .ok_or_else(|| invalid("posterior grid evaluation count overflows"))?;
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
}
impl PosteriorAnalysis {
    pub(super) fn new(
        options: &PosteriorAnalysisOptions,
        genes: usize,
        coefficient: Option<usize>,
        map_only: bool,
    ) -> Result<Self, DeseqError> {
        let na = f64::from_bits(0x7ff0_0000_0000_07a2);
        let intervals = if options.interval.is_some() && coefficient.is_some() {
            Some(RowMajorMatrix::from_row_major(
                genes,
                2,
                vec![
                    na;
                    genes
                        .checked_mul(2)
                        .ok_or_else(|| invalid("posterior interval dimensions overflow"))?
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
                            .ok_or_else(|| invalid("posterior contrast dimensions overflow"))?
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
    pub(super) fn row(
        &mut self,
        row: AnalysisRow<'_>,
        mut log_posterior: impl FnMut(&[f64], f64) -> Result<f64, DeseqError>,
    ) -> Result<Option<grid::GridSummary>, DeseqError> {
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
                self.contrast_mode.as_mut().unwrap().as_mut_slice()[gene * k + column] = mean;
                self.contrast_standard_deviation
                    .as_mut()
                    .unwrap()
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
            self.intervals.as_mut().unwrap().as_mut_slice()[gene * 2..gene * 2 + 2]
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
        self.intervals.as_mut().unwrap().as_mut_slice()[gene * 2..gene * 2 + 2]
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
    pub fn value(&self, beta: &[f64], shift: f64) -> Result<f64, DeseqError> {
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
                .ok_or_else(|| invalid("missing grid parameter"))?[0]
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
