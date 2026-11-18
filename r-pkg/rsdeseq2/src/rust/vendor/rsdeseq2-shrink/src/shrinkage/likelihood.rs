//! Shared location likelihood selection for scalar components and batch methods.
use super::{
    RestrictedLikelihood, student_t,
    truncated_standard::{log_phi, truncated_mean, truncated_second},
};
use crate::{DeseqError, math::normal_log_cdf};
use std::sync::Arc;
/// Location-scale observation likelihood, independent of the calling frontend.
#[derive(Clone, Debug)]
pub enum LocationLikelihood {
    Normal,
    StudentT(f64),
    Builtin(super::BuiltinLocationLikelihood),
    /// Validated elementwise expressions and captured scalar primitives.
    Restricted(Arc<RestrictedLikelihood>),
}
impl LocationLikelihood {
    pub(super) fn validate_observations(&self, n: usize) -> Result<(), DeseqError> {
        if let Self::Builtin(b) = self {
            b.validate_observations(n)
        } else {
            Ok(())
        }
    }
    pub(super) fn at(
        &self,
        row: usize,
        n: usize,
    ) -> Result<std::borrow::Cow<'_, Self>, DeseqError> {
        if let Self::Builtin(b) = self {
            Ok(std::borrow::Cow::Owned(Self::Builtin(
                b.for_observation(row, n)?,
            )))
        } else {
            Ok(std::borrow::Cow::Borrowed(self))
        }
    }
    pub(super) fn validated(self) -> Result<Self, DeseqError> {
        match self {
            Self::StudentT(df) if df.is_nan() || df <= 0.0 => Err(DeseqError::InvalidOptions {
                reason: "Student-t degrees of freedom must be positive".into(),
            }),
            Self::StudentT(df) if df == f64::INFINITY => Ok(Self::Normal),
            other => Ok(other),
        }
    }
    pub(super) fn log_cdf(&self, x: f64) -> f64 {
        match self {
            Self::Normal => normal_log_cdf(x),
            Self::StudentT(df) => student_t::log_cdf(x, *df),
            Self::Builtin(b) => b.log_cdf(x),
            Self::Restricted(f) => f.log_cdf(x),
        }
    }
    pub(super) fn log_density(&self, x: f64) -> f64 {
        match self {
            Self::Normal => log_phi(x),
            Self::StudentT(df) => student_t::log_density(x, *df),
            Self::Builtin(b) => b.log_density(x),
            Self::Restricted(f) => f.log_density(x),
        }
    }
    pub(super) fn moment(&self, a: f64, b: f64, order: u8) -> f64 {
        match self {
            Self::Normal => {
                if order == 1 {
                    truncated_mean(a, b)
                } else {
                    truncated_second(a, b)
                }
            }
            Self::StudentT(df) => student_t::truncated_moment(a, b, *df, i32::from(order)),
            Self::Builtin(l) => l.truncated_moment(a, b, order),
            Self::Restricted(f) => f.truncated_moment(a, b, order).unwrap_or(f64::NAN),
        }
    }
}

impl super::MixtureGrid {
    /// Construct a fixed prior under the requested observation likelihood.
    pub fn prior_with_likelihood(
        &self,
        weights: &[f64],
        likelihood: &LocationLikelihood,
    ) -> Result<super::FixedMixture, DeseqError> {
        match likelihood {
            LocationLikelihood::Normal => self.prior(weights),
            LocationLikelihood::StudentT(df) if *df == f64::INFINITY => self.prior(weights),
            LocationLikelihood::Restricted(f)
                if f.name() == "normal"
                    && !matches!(self.components, super::MixtureComponents::Uniform { .. }) =>
            {
                self.prior(weights)
            }
            _ => match &self.components {
                super::MixtureComponents::Uniform { lower, upper } => {
                    Ok(super::FixedMixture::Uniform(
                        super::UniformMixture::new(weights.to_vec(), lower.clone(), upper.clone())?
                            .with_likelihood(likelihood.clone())?,
                    ))
                }
                _ => Err(DeseqError::InvalidOptions {
                    reason:
                        "normal/truncated-normal mixtures require a normal observation likelihood"
                            .into(),
                }),
            },
        }
    }
}
