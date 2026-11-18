//! Small numerical helpers following R's arithmetic.

#[doc(hidden)]
pub mod elementary;
/// Medians as R's `median()` computes them.
pub mod median;
/// Normal-distribution log tail probabilities.
pub mod normal;
/// A small deterministic bounded optimizer.
pub mod optim;
/// Matrix rank by Gaussian elimination and by R's pivoted QR rule.
pub mod rank;
/// The trigamma function.
pub mod special;

pub use median::{median, median_finite};
pub use normal::{normal_log_cdf, normal_log_upper_tail};
pub use special::trigamma;
