//! Primitive, shared built-in likelihood parameters with explicit row selection.
use super::{
    location_missing::NA,
    log_f,
    truncated_standard::{log_phi, scaled_mean, truncated_second},
};
use crate::{
    DeseqError,
    matrix::{MatrixOrder, MatrixView},
};
use std::sync::Arc;
#[derive(Clone, Debug)]
enum Parameters {
    Count(super::count_likelihood::CountLikelihood),
    NormalMixture {
        weights: Vec<f64>,
        scales: Vec<f64>,
        rows: usize,
    },
    LogF {
        first: Vec<f64>,
        second: Vec<f64>,
    },
}
/// Built-in continuous observation likelihood. Parameter vectors/matrices are
/// shared across selected rows, never copied per posterior or passed back to R.
#[derive(Clone, Debug)]
pub struct BuiltinLocationLikelihood {
    parameters: Arc<Parameters>,
    row: usize,
}
fn invalid(reason: &str) -> DeseqError {
    DeseqError::InvalidOptions {
        reason: format!("built-in likelihood: {reason}"),
    }
}
impl BuiltinLocationLikelihood {
    pub fn poisson(
        y: Vec<f64>,
        scale: Vec<f64>,
        link: super::PoissonLink,
    ) -> Result<Self, DeseqError> {
        Ok(Self {
            parameters: Arc::new(Parameters::Count(
                super::count_likelihood::CountLikelihood::poisson(y, scale, link)?,
            )),
            row: 0,
        })
    }
    pub fn binomial(
        y: Vec<f64>,
        trials: Vec<f64>,
        link: super::BinomialLink,
    ) -> Result<Self, DeseqError> {
        Ok(Self {
            parameters: Arc::new(Parameters::Count(
                super::count_likelihood::CountLikelihood::binomial(y, trials, link)?,
            )),
            row: 0,
        })
    }
    /// Check imported stored log/logit constructor captures without executing R.
    pub fn validate_count_captures(
        &self,
        y: &[f64],
        trials: Option<&[f64]>,
    ) -> Result<(), DeseqError> {
        if let Parameters::Count(c) = self.parameters.as_ref() {
            c.validate_captures(y, trials)
        } else {
            Err(invalid("adjusted captures require a count likelihood"))
        }
    }
    pub(super) fn grid_data(&self, n: usize, se: &[f64]) -> Option<(Vec<f64>, Vec<f64>)> {
        if let Parameters::Count(c) = self.parameters.as_ref() {
            Some(c.grid_data(n, se))
        } else {
            None
        }
    }
    pub(super) fn mode_bounds(
        &self,
        n: usize,
        range: (f64, f64),
        seed: Option<f64>,
    ) -> Result<Option<(f64, f64)>, DeseqError> {
        self.validate_observations(n)?;
        Ok(if let Parameters::Count(c) = self.parameters.as_ref() {
            Some(c.mode_bounds(n, range, seed))
        } else {
            None
        })
    }
    pub(super) fn is_count(&self) -> bool {
        matches!(self.parameters.as_ref(), Parameters::Count(_))
    }
    pub(super) fn support(&self) -> Option<(f64, f64)> {
        if let Parameters::Count(c) = self.parameters.as_ref() {
            c.support()
        } else {
            None
        }
    }

    /// Import canonical primitive capture layouts. Equal total lengths flatten
    /// the captures in column order; unequal SD matrices recycle components.
    pub fn normal_mixture_from_captures(
        weights: MatrixView<'_, f64>,
        scales: MatrixView<'_, f64>,
        scales_are_matrix: bool,
    ) -> Result<Self, DeseqError> {
        super::normal_mixture_captures::import(weights, scales, scales_are_matrix)
    }
    pub fn normal_mixture(weights: Vec<f64>, scales: Vec<f64>) -> Result<Self, DeseqError> {
        let matrix = MatrixView::new(1, scales.len(), &scales, MatrixOrder::RowMajor)?;
        Self::normal_mixture_by_observation(weights, matrix)
    }
    /// One row of SDs per observation, or one shared row; columns match weights.
    pub fn normal_mixture_by_observation(
        weights: Vec<f64>,
        scales: MatrixView<'_, f64>,
    ) -> Result<Self, DeseqError> {
        if weights.is_empty()
            || weights.iter().any(|x| !x.is_finite() || *x < 0.0)
            || !weights.iter().any(|x| *x > 0.0)
            || scales.n_cols() != weights.len()
            || scales.n_rows() == 0
        {
            return Err(invalid(
                "nonnegative nonzero weights and matching noise SD columns are required",
            ));
        }
        let mut values = Vec::with_capacity(
            scales
                .n_rows()
                .checked_mul(scales.n_cols())
                .ok_or_else(|| invalid("noise matrix size overflows"))?,
        );
        for i in 0..scales.n_rows() {
            for j in 0..scales.n_cols() {
                let s = *scales.get(i, j)?;
                if !s.is_finite() || s < 0.0 || !(s * s).is_finite() {
                    return Err(invalid(
                        "noise SDs must be finite nonnegative and squareable",
                    ));
                }
                values.push(s);
            }
        }
        Ok(Self {
            parameters: Arc::new(Parameters::NormalMixture {
                weights,
                scales: values,
                rows: scales.n_rows(),
            }),
            row: 0,
        })
    }
    /// Positive degrees of freedom, shared scalars or per-observation vectors.
    pub fn log_f(first: Vec<f64>, second: Vec<f64>) -> Result<Self, DeseqError> {
        if first.is_empty()
            || second.is_empty()
            || first.iter().chain(&second).any(|x| x.is_nan() || *x <= 0.0)
        {
            return Err(invalid(
                "log-F degrees of freedom must be positive and nonmissing",
            ));
        }
        Ok(Self {
            parameters: Arc::new(Parameters::LogF { first, second }),
            row: 0,
        })
    }
    pub fn name(&self) -> &'static str {
        match self.parameters.as_ref() {
            Parameters::Count(c) => c.name(),
            Parameters::NormalMixture { .. } => "normalmix",
            Parameters::LogF { .. } => "logF",
        }
    }
    pub(super) fn validate_observations(&self, n: usize) -> Result<(), DeseqError> {
        let valid = |length| length == 1 || length == n;
        let ok = match self.parameters.as_ref() {
            Parameters::Count(c) => {
                let (a, b) = c.lengths();
                valid(a) && valid(b)
            }
            Parameters::NormalMixture { rows, .. } => valid(*rows),
            Parameters::LogF { first, second } => valid(first.len()) && valid(second.len()),
        };
        if !ok {
            return Err(invalid(
                "parameter rows must be one or match the complete observation count",
            ));
        }
        Ok(())
    }
    pub fn for_observation(&self, row: usize, n: usize) -> Result<Self, DeseqError> {
        self.validate_observations(n)?;
        if row >= n {
            return Err(invalid("observation index is out of range"));
        }
        Ok(Self {
            parameters: Arc::clone(&self.parameters),
            row,
        })
    }
    fn df(&self, values: &[f64]) -> f64 {
        values[if values.len() == 1 { 0 } else { self.row }]
    }
    fn mixture(&self, f: impl Fn(f64) -> f64) -> f64 {
        let Parameters::NormalMixture {
            weights,
            scales,
            rows,
        } = self.parameters.as_ref()
        else {
            unreachable!()
        };
        let start = if *rows == 1 {
            0
        } else {
            self.row * weights.len()
        };
        weights
            .iter()
            .zip(&scales[start..start + weights.len()])
            .map(|(&w, &s)| w * f(s))
            .sum()
    }
    pub fn log_cdf(&self, x: f64) -> f64 {
        match self.parameters.as_ref() {
            Parameters::Count(c) => c.log_cdf(x, self.row),
            Parameters::NormalMixture { .. } => self
                .mixture(|s| {
                    if s == 0.0 {
                        if x.is_nan() {
                            f64::NAN
                        } else if x < 0.0 {
                            0.0
                        } else {
                            1.0
                        }
                    } else {
                        crate::math::normal_log_cdf(x / s).exp()
                    }
                })
                .ln(),
            Parameters::LogF { first, second } => {
                log_f::log_cdf(x, self.df(first), self.df(second))
            }
        }
    }
    pub fn log_density(&self, x: f64) -> f64 {
        match self.parameters.as_ref() {
            Parameters::Count(c) => c.log_density(x, self.row),
            Parameters::NormalMixture { .. } => self
                .mixture(|s| {
                    if s == 0.0 {
                        if x.is_nan() {
                            f64::NAN
                        } else if x == 0.0 {
                            f64::INFINITY
                        } else {
                            0.0
                        }
                    } else {
                        (log_phi(x / s) - s.ln()).exp()
                    }
                })
                .ln(),
            Parameters::LogF { first, second } => {
                log_f::log_density(x, self.df(first), self.df(second))
            }
        }
    }
    pub fn truncated_moment(&self, a: f64, b: f64, order: u8) -> f64 {
        if a.is_nan() || b.is_nan() || a > b {
            return if super::location_missing::is_na(a) || super::location_missing::is_na(b) {
                NA
            } else {
                f64::NAN
            };
        }
        if order != 1 && order != 2 {
            return f64::NAN;
        }
        match self.parameters.as_ref() {
            Parameters::Count(c) => c.moment(a, b, order, self.row),
            Parameters::NormalMixture { .. } => self.mixture(|s| {
                if order == 1 {
                    scaled_mean(a, b, 0.0, s)
                } else if s == 0.0 {
                    0.0_f64.clamp(a, b).powi(2)
                } else {
                    s * s * truncated_second(a / s, b / s)
                }
            }),
            Parameters::LogF { first, second } => {
                if order == 2 {
                    NA
                } else {
                    log_f::truncated_mean(a, b, self.df(first), self.df(second))
                }
            }
        }
    }
}
