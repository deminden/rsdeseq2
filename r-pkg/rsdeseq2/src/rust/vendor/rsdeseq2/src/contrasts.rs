//! Contrast selection and all-zero rules.

mod aliases;
mod all_zero;
mod counts_view;
mod model_frame;
mod resolution;
mod types;
mod zero_policy;

pub use types::{
    ContrastSpec, FactorLevelContrast, ResolvedResultsContrast, ResultsContrast,
    ResultsContrastAllZero,
};

pub use model_frame::factor_level_contrast_from_model_frame;

pub use resolution::{resolve_coefficient_index, resolve_contrast, resolve_results_contrast};

pub use all_zero::{contrast_all_zero_factor_levels, contrast_all_zero_numeric};

pub use counts_view::ContrastCounts;
pub use zero_policy::ContrastZeroPolicy;
