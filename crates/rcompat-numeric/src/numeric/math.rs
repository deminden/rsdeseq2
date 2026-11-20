//! Elementary functions that follow the numeric profile.
//!
//! Under the R profiles each method calls the host C
//! library exactly as `f64::ln` and friends do, so results match R built on
//! that host. Under `Portable` they call the pure-Rust `libm` crate, whose
//! results do not depend on the operating system or architecture.
//!
//! `sqrt`, `mul_add` and `powi` are left on `f64`: the first two are
//! correctly rounded everywhere and `powi` is a fixed multiplication chain.

use super::profile::host_libm;

/// Profile-following elementary functions on `f64`.
///
/// Each method mirrors the `f64` method of the same name without the `_p`
/// suffix.
pub trait ProfileMath: Copy {
    /// Natural logarithm.
    fn ln_p(self) -> f64;
    /// Exponential.
    fn exp_p(self) -> f64;
    /// Real power.
    fn powf_p(self, exponent: f64) -> f64;
    /// Base-10 logarithm.
    fn log10_p(self) -> f64;
    /// Base-2 logarithm.
    fn log2_p(self) -> f64;
    /// `exp(x) - 1`.
    fn exp_m1_p(self) -> f64;
    /// `ln(1 + x)`.
    fn ln_1p_p(self) -> f64;
    /// Sine.
    fn sin_p(self) -> f64;
    /// Cosine.
    fn cos_p(self) -> f64;
    /// Tangent.
    fn tan_p(self) -> f64;
    /// Hyperbolic sine.
    fn sinh_p(self) -> f64;
    /// `sqrt(x^2 + y^2)` without undue overflow.
    fn hypot_p(self, other: f64) -> f64;
    /// Inverse hyperbolic sine.
    fn asinh_p(self) -> f64;
}

impl ProfileMath for f64 {
    #[inline]
    fn ln_p(self) -> f64 {
        if host_libm() {
            self.ln()
        } else {
            libm::log(self)
        }
    }

    #[inline]
    fn exp_p(self) -> f64 {
        if host_libm() {
            self.exp()
        } else {
            libm::exp(self)
        }
    }

    #[inline]
    fn powf_p(self, exponent: f64) -> f64 {
        if host_libm() {
            self.powf(exponent)
        } else {
            libm::pow(self, exponent)
        }
    }

    #[inline]
    fn log10_p(self) -> f64 {
        if host_libm() {
            self.log10()
        } else {
            libm::log10(self)
        }
    }

    #[inline]
    fn log2_p(self) -> f64 {
        if host_libm() {
            self.log2()
        } else {
            libm::log2(self)
        }
    }

    #[inline]
    fn exp_m1_p(self) -> f64 {
        if host_libm() {
            self.exp_m1()
        } else {
            libm::expm1(self)
        }
    }

    #[inline]
    fn ln_1p_p(self) -> f64 {
        if host_libm() {
            self.ln_1p()
        } else {
            libm::log1p(self)
        }
    }

    #[inline]
    fn sin_p(self) -> f64 {
        if host_libm() {
            self.sin()
        } else {
            libm::sin(self)
        }
    }

    #[inline]
    fn cos_p(self) -> f64 {
        if host_libm() {
            self.cos()
        } else {
            libm::cos(self)
        }
    }

    #[inline]
    fn tan_p(self) -> f64 {
        if host_libm() {
            self.tan()
        } else {
            libm::tan(self)
        }
    }

    #[inline]
    fn sinh_p(self) -> f64 {
        if host_libm() {
            self.sinh()
        } else {
            libm::sinh(self)
        }
    }

    #[inline]
    fn hypot_p(self, other: f64) -> f64 {
        if host_libm() {
            self.hypot(other)
        } else {
            libm::hypot(self, other)
        }
    }

    #[inline]
    fn asinh_p(self) -> f64 {
        if host_libm() {
            self.asinh()
        } else {
            libm::asinh(self)
        }
    }
}
