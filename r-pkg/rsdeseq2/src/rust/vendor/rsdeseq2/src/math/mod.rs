//! Small numerical helpers used by the core algorithms.

pub mod distributions;
pub(crate) mod elementary;
pub mod median;
pub mod normal;
pub mod optim;
pub mod qr;
pub mod special;
pub(crate) mod sum;

pub use distributions::{
    NegativeBinomialHelpers, negative_binomial_helpers, negative_binomial_log_likelihood,
    negative_binomial_log_likelihood_matrix, negative_binomial_log_likelihood_weighted,
    negative_binomial_log_pmf, negative_binomial_negative_twice_log_likelihood,
};
pub use median::{median, median_finite};
pub use normal::{normal_log_cdf, normal_log_upper_tail};
pub use special::trigamma;
