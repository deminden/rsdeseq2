//! Common dispersion trend fit representation.

use super::local::LocalDispersionTrendFit;
use super::mean::MeanDispersionTrendFit;
use super::parametric::ParametricDispersionTrendFit;
use crate::errors::DeseqError;

/// Fitted dispersion trend selected by a DESeq2-style `fitType`.
#[derive(Clone, Debug, PartialEq)]
pub enum DispersionTrendFit {
    /// Parametric `asymptDisp + extraPois / mean` trend.
    Parametric(ParametricDispersionTrendFit),
    /// Pure-Rust local log-dispersion trend.
    Local(LocalDispersionTrendFit),
    /// Constant mean dispersion trend.
    Mean(MeanDispersionTrendFit),
}

impl DispersionTrendFit {
    /// Stable DESeq2-style fit type label for this fitted trend.
    pub fn fit_type_label(&self) -> &'static str {
        match self {
            Self::Parametric(_) => "parametric",
            Self::Local(_) => "local",
            Self::Mean(_) => "mean",
        }
    }

    /// Fitted dispersion values for every input row; missing rows are `NaN`.
    pub fn disp_fit(&self) -> &[f64] {
        match self {
            Self::Parametric(fit) => &fit.disp_fit,
            Self::Local(fit) => &fit.disp_fit,
            Self::Mean(fit) => &fit.disp_fit,
        }
    }

    /// Evaluate the fitted trend for all finite positive means, returning
    /// `NaN` for missing rows.
    pub fn evaluate_many_allow_missing(&self, means: &[f64]) -> Result<Vec<f64>, DeseqError> {
        match self {
            Self::Parametric(fit) => fit.trend.evaluate_many_allow_missing(means),
            Self::Local(fit) => fit.trend.evaluate_many_allow_missing(means),
            Self::Mean(fit) => fit.trend.evaluate_many_allow_missing(means),
        }
    }

    /// Rows retained by the selected trend fit rule.
    ///
    /// For `fitType="mean"`, this is the preliminary viability mask using
    /// `dispGeneEst > 100 * minDisp`; use [`Self::use_for_mean`] to inspect the
    /// trimmed-mean mask.
    pub fn use_for_fit(&self) -> &[bool] {
        match self {
            Self::Parametric(fit) => &fit.use_for_fit,
            Self::Local(fit) => &fit.use_for_fit,
            Self::Mean(fit) => &fit.use_for_fit,
        }
    }

    /// Rows used by the `fitType="mean"` trimmed-mean calculation.
    ///
    /// Returns `None` for parametric and local trends.
    pub fn use_for_mean(&self) -> Option<&[bool]> {
        match self {
            Self::Mean(fit) => Some(&fit.use_for_mean),
            Self::Parametric(_) | Self::Local(_) => None,
        }
    }

    /// Number of rows used by the `fitType="mean"` trimmed-mean calculation.
    ///
    /// Returns `None` for parametric and local trends.
    pub fn genes_used_for_mean(&self) -> Option<usize> {
        match self {
            Self::Mean(fit) => Some(fit.genes_used_for_mean),
            Self::Parametric(_) | Self::Local(_) => None,
        }
    }

    /// Whether `fitType="local"` used DESeq2's all-near-minimum dispersion floor.
    ///
    /// Returns `None` for parametric and mean trends.
    pub fn used_min_disp_floor(&self) -> Option<bool> {
        match self {
            Self::Local(fit) => Some(fit.used_min_disp_floor),
            Self::Parametric(_) | Self::Mean(_) => None,
        }
    }

    /// Number of rows retained by the selected trend fit rule.
    pub fn genes_used(&self) -> usize {
        match self {
            Self::Parametric(fit) => fit.genes_used,
            Self::Local(fit) => fit.genes_used,
            Self::Mean(fit) => fit.genes_used,
        }
    }
}
