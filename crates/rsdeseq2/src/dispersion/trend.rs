//! Dispersion trend models, fitting, and row selection.

mod local;
mod mean;
mod model;
mod numerics;
mod parametric;
mod selection;
#[cfg(test)]
mod tests;

pub use parametric::{
    ParametricDispersionTrend, ParametricDispersionTrendFit, ParametricDispersionTrendOptions,
    fit_parametric_dispersion_trend,
};

pub use local::{
    LocalDispersionTrend, LocalDispersionTrendFit, LocalDispersionTrendOptions,
    fit_local_dispersion_trend,
};

pub use mean::{
    MeanDispersionTrend, MeanDispersionTrendFit, MeanDispersionTrendOptions,
    fit_mean_dispersion_trend,
};

pub use model::DispersionTrendFit;

pub use selection::{
    fit_dispersion_trend, local_trend_use_for_fit, mean_trend_use_for_fit, mean_trend_use_for_mean,
    parametric_trend_use_for_fit,
};
