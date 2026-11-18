//! Small numerical helpers used by the core algorithms.

/// DESeq2-parameterized negative-binomial likelihood helpers.
pub mod distributions;
#[doc(hidden)]
pub use rcompat_numeric::math::elementary;
// `median` names both the module and the function.
pub use rcompat_numeric::math::{
    median, median_finite, normal, normal_log_cdf, normal_log_upper_tail, optim, rank, trigamma,
};

/// Former name of [`rank`], which computes ranks and has no QR factorization.
#[deprecated(note = "renamed to `rsdeseq2::math::rank`")]
pub mod qr {
    pub use super::rank::*;
}
pub use rcompat_numeric::math::special;

pub use distributions::{
    NegativeBinomialHelpers, negative_binomial_helpers, negative_binomial_log_likelihood,
    negative_binomial_log_likelihood_matrix, negative_binomial_log_likelihood_weighted,
    negative_binomial_log_pmf, negative_binomial_negative_twice_log_likelihood,
};
