//! Small numerical helpers used by the core algorithms.

pub mod distributions;
pub mod median;
pub mod optim;
pub mod qr;
pub mod special;

pub use distributions::{
    negative_binomial_helpers, negative_binomial_log_likelihood,
    negative_binomial_log_likelihood_matrix, negative_binomial_log_likelihood_weighted,
    negative_binomial_log_pmf, negative_binomial_negative_twice_log_likelihood,
    NegativeBinomialHelpers,
};
pub use median::{median, median_finite};
pub use special::trigamma;
