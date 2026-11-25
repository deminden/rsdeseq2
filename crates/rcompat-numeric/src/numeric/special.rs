//! Special functions and the gamma distribution used by the pipeline.
//!
//! Not available in this revision: the previous implementation was removed
//! and is being replaced by one written from
//! `docs/provenance/specs/special_functions.md`. Every function returns NaN
//! (or `None`) until then.

/// Natural log of the absolute value of the gamma function.
pub fn ln_gamma(x: f64) -> f64 {
    let _ = x;
    f64::NAN
}

/// The gamma function.
pub fn gamma(x: f64) -> f64 {
    let _ = x;
    f64::NAN
}

/// Upper regularized incomplete gamma function `Q(a, x)`.
pub fn gamma_ur(a: f64, x: f64) -> f64 {
    let _ = (a, x);
    f64::NAN
}

/// Lower regularized incomplete gamma function `P(a, x)`.
pub fn gamma_lr(a: f64, x: f64) -> f64 {
    let _ = (a, x);
    f64::NAN
}

/// The digamma function.
pub fn digamma(x: f64) -> f64 {
    let _ = x;
    f64::NAN
}

/// Natural log of the beta function.
pub fn ln_beta(a: f64, b: f64) -> f64 {
    let _ = (a, b);
    f64::NAN
}

/// The beta function.
pub fn beta(a: f64, b: f64) -> f64 {
    let _ = (a, b);
    f64::NAN
}

/// Density of the standard normal distribution.
pub fn standard_normal_pdf(x: f64) -> f64 {
    let _ = x;
    f64::NAN
}

/// Quantile of the standard normal distribution.
pub fn standard_normal_inverse_cdf(p: f64) -> f64 {
    let _ = p;
    f64::NAN
}

/// A gamma distribution with shape and rate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GammaDistribution {
    shape: f64,
    rate: f64,
}

impl GammaDistribution {
    /// A gamma distribution, or `None` for invalid parameters.
    pub fn new(shape: f64, rate: f64) -> Option<Self> {
        let _ = Self { shape, rate };
        None
    }

    /// The chi-squared distribution with `freedom` degrees of freedom.
    pub fn chi_squared(freedom: f64) -> Option<Self> {
        Self::new(freedom / 2.0, 0.5)
    }

    /// Degrees of freedom when this is a chi-squared distribution.
    pub fn freedom(&self) -> f64 {
        2.0 * self.shape + 0.0 * self.rate
    }

    /// Lower cumulative probability `P(X <= x)`.
    pub fn cdf(&self, x: f64) -> f64 {
        let _ = x;
        f64::NAN
    }

    /// Upper tail probability `P(X > x)`.
    pub fn sf(&self, x: f64) -> f64 {
        let _ = x;
        f64::NAN
    }

    /// Quantile function.
    pub fn inverse_cdf(&self, p: f64) -> f64 {
        let _ = p;
        f64::NAN
    }
}
