//! Special functions and distributions used by the pipeline, following the
//! numeric profile.
//!
//! These are ports of the `statrs` 0.19.1 implementations the crate used
//! before (Lanczos `ln_gamma`, AS 103 `digamma`, the Cephes-style incomplete
//! gamma ratios, Boost's `erf_inv` rational approximations and `statrs`'s
//! safeguarded Newton quantile), with every elementary function routed
//! through [`ProfileMath`]. Under the R profiles the arithmetic is therefore
//! exactly that of `statrs` on the host; under `Portable` it is the same
//! algorithm on pure-Rust maths, identical on every host.
//!
//! statrs is Copyright (c) 2016 Michael Ma and distributed under the MIT
//! license.

// Coefficient tables are copied verbatim from statrs so the port stays bit-identical.
#![allow(clippy::excessive_precision, clippy::unreadable_literal)]

use std::f64::consts as f64_consts;

use super::math::ProfileMath;

const LN_PI: f64 = 1.1447298858494001741434273513530587116472948129153;
const LN_2_SQRT_E_OVER_PI: f64 = 0.6207822376352452223455184457816472122518527279025978;
const TWO_SQRT_E_OVER_PI: f64 = 1.8603827342052657173362492472666631120594218414085755;
const SQRT_2PI: f64 = 2.5066282746310005024157652848110452530069867406099;

/// Auxiliary variable when evaluating the `gamma_ln` function.
const GAMMA_R: f64 = 10.900511;

/// Polynomial coefficients for approximating the `gamma_ln` function.
const GAMMA_DK: &[f64] = &[
    2.48574089138753565546e-5,
    1.05142378581721974210,
    -3.45687097222016235469,
    4.51227709466894823700,
    -2.98285225323576655721,
    1.05639711577126713077,
    -1.95428773191645869583e-1,
    1.70970543404441224307e-2,
    -5.71926117404305781283e-4,
    4.63399473359905636708e-6,
    -2.71994908488607703910e-9,
];

/// Natural log of the gamma function (Pugh 2004 Lanczos approximation).
pub fn ln_gamma(x: f64) -> f64 {
    if x < 0.5 {
        let s = GAMMA_DK
            .iter()
            .enumerate()
            .skip(1)
            .fold(GAMMA_DK[0], |s, t| s + t.1 / (t.0 as f64 - x));

        LN_PI
            - (f64_consts::PI * x).sin_p().ln_p()
            - s.ln_p()
            - LN_2_SQRT_E_OVER_PI
            - (0.5 - x) * ((0.5 - x + GAMMA_R) / f64_consts::E).ln_p()
    } else {
        let s = GAMMA_DK
            .iter()
            .enumerate()
            .skip(1)
            .fold(GAMMA_DK[0], |s, t| s + t.1 / (x + t.0 as f64 - 1.0));

        s.ln_p() + LN_2_SQRT_E_OVER_PI + (x - 0.5) * ((x - 0.5 + GAMMA_R) / f64_consts::E).ln_p()
    }
}

/// The gamma function (Pugh 2004 Lanczos approximation).
pub fn gamma(x: f64) -> f64 {
    if x < 0.5 {
        let s = GAMMA_DK
            .iter()
            .enumerate()
            .skip(1)
            .fold(GAMMA_DK[0], |s, t| s + t.1 / (t.0 as f64 - x));

        f64_consts::PI
            / ((f64_consts::PI * x).sin_p()
                * s
                * TWO_SQRT_E_OVER_PI
                * ((0.5 - x + GAMMA_R) / f64_consts::E).powf_p(0.5 - x))
    } else {
        let s = GAMMA_DK
            .iter()
            .enumerate()
            .skip(1)
            .fold(GAMMA_DK[0], |s, t| s + t.1 / (x + t.0 as f64 - 1.0));

        s * TWO_SQRT_E_OVER_PI * ((x - 0.5 + GAMMA_R) / f64_consts::E).powf_p(x - 0.5)
    }
}

/// `statrs`'s `ulps_eq!(x, 0.0, epsilon = 0.0)` with its default of 5 ULPs.
fn within_five_ulps_of_zero(x: f64) -> bool {
    x == 0.0 || (x > 0.0 && x.to_bits() <= 5) || (x < 0.0 && (-x).to_bits() <= 5)
}

/// Upper regularized incomplete gamma function `Q(a, x)`.
///
/// Panics, as `statrs::function::gamma::gamma_ur` does, unless `a` and `x`
/// are positive and finite (NaN arguments give NaN).
pub fn gamma_ur(a: f64, x: f64) -> f64 {
    if a.is_nan() || x.is_nan() {
        return f64::NAN;
    }
    assert!(a > 0.0 && a != f64::INFINITY, "gamma_ur: invalid a");
    assert!(x > 0.0 && x != f64::INFINITY, "gamma_ur: invalid x");

    let eps = 0.000000000000001;
    let big = 4503599627370496.0;
    let big_inv = 2.22044604925031308085e-16;

    if x < 1.0 || x <= a {
        return 1.0 - gamma_lr(a, x);
    }

    let mut ax = a * x.ln_p() - x - ln_gamma(a);
    if ax < -709.78271289338399 {
        return if a < x { 0.0 } else { 1.0 };
    }

    ax = ax.exp_p();
    let mut y = 1.0 - a;
    let mut z = x + y + 1.0;
    let mut c = 0.0;
    let mut pkm2 = 1.0;
    let mut qkm2 = x;
    let mut pkm1 = x + 1.0;
    let mut qkm1 = z * x;
    let mut ans = pkm1 / qkm1;
    loop {
        y += 1.0;
        z += 2.0;
        c += 1.0;
        let yc = y * c;
        let pk = pkm1 * z - pkm2 * yc;
        let qk = qkm1 * z - qkm2 * yc;

        pkm2 = pkm1;
        pkm1 = pk;
        qkm2 = qkm1;
        qkm1 = qk;

        if pk.abs() > big {
            pkm2 *= big_inv;
            pkm1 *= big_inv;
            qkm2 *= big_inv;
            qkm1 *= big_inv;
        }

        if qk != 0.0 {
            let r = pk / qk;
            let t = ((ans - r) / r).abs();
            ans = r;

            if t <= eps {
                break;
            }
        }
    }
    ans * ax
}

/// Lower regularized incomplete gamma function `P(a, x)`.
///
/// Panics, as `statrs::function::gamma::gamma_lr` does, unless `a` and `x`
/// are positive and finite (NaN arguments give NaN).
pub fn gamma_lr(a: f64, x: f64) -> f64 {
    if a.is_nan() || x.is_nan() {
        return f64::NAN;
    }
    assert!(a > 0.0 && a != f64::INFINITY, "gamma_lr: invalid a");
    assert!(x > 0.0 && x != f64::INFINITY, "gamma_lr: invalid x");

    let eps = 0.000000000000001;
    let big = 4503599627370496.0;
    let big_inv = 2.22044604925031308085e-16;

    if within_five_ulps_of_zero(a) {
        return 1.0;
    }
    if within_five_ulps_of_zero(x) {
        return 0.0;
    }

    let ax = a * x.ln_p() - x - ln_gamma(a);
    if ax < -709.78271289338399 {
        if a < x {
            return 1.0;
        }
        return 0.0;
    }
    if x <= 1.0 || x <= a {
        let mut r2 = a;
        let mut c2 = 1.0;
        let mut ans2 = 1.0;
        loop {
            r2 += 1.0;
            c2 *= x / r2;
            ans2 += c2;

            if c2 / ans2 <= eps {
                break;
            }
        }
        return ax.exp_p() * ans2 / a;
    }

    let mut y = 1.0 - a;
    let mut z = x + y + 1.0;
    let mut c = 0;

    let mut p3 = 1.0;
    let mut q3 = x;
    let mut p2 = x + 1.0;
    let mut q2 = z * x;
    let mut ans = p2 / q2;

    loop {
        y += 1.0;
        z += 2.0;
        c += 1;
        let yc = y * f64::from(c);

        let p = p2 * z - p3 * yc;
        let q = q2 * z - q3 * yc;

        p3 = p2;
        p2 = p;
        q3 = q2;
        q2 = q;

        if p.abs() > big {
            p3 *= big_inv;
            p2 *= big_inv;
            q3 *= big_inv;
            q2 *= big_inv;
        }

        if q != 0.0 {
            let nextans = p / q;
            let error = ((ans - nextans) / nextans).abs();
            ans = nextans;

            if error <= eps {
                break;
            }
        }
    }
    1.0 - ax.exp_p() * ans
}

/// The digamma function (Algorithm AS 103, Bernardo 1976).
pub fn digamma(x: f64) -> f64 {
    let c = 12.0;
    let d1 = -0.57721566490153286;
    let d2 = 1.6449340668482264365;
    let s = 1e-6;
    let s3 = 1.0 / 12.0;
    let s4 = 1.0 / 120.0;
    let s5 = 1.0 / 252.0;
    let s6 = 1.0 / 240.0;
    let s7 = 1.0 / 132.0;

    if x == f64::NEG_INFINITY || x.is_nan() {
        return f64::NAN;
    }
    if x <= 0.0 && x.floor() == x {
        return f64::NEG_INFINITY;
    }
    if x < 0.0 {
        return digamma(1.0 - x) + f64_consts::PI / (-f64_consts::PI * x).tan_p();
    }
    if x <= s {
        return d1 - 1.0 / x + d2 * x;
    }

    let mut result = 0.0;
    let mut z = x;
    while z < c {
        result -= 1.0 / z;
        z += 1.0;
    }

    if z >= c {
        let mut r = 1.0 / z;
        result += z.ln_p() - 0.5 * r;
        r *= r;

        result -= r * (s3 - r * (s4 - r * (s5 - r * (s6 - r * s7))));
    }
    result
}

/// Natural log of the beta function. Panics unless `a` and `b` are positive.
pub fn ln_beta(a: f64, b: f64) -> f64 {
    assert!(a > 0.0, "ln_beta: a must be positive");
    assert!(b > 0.0, "ln_beta: b must be positive");
    ln_gamma(a) + ln_gamma(b) - ln_gamma(a + b)
}

/// The beta function. Panics unless `a` and `b` are positive.
pub fn beta(a: f64, b: f64) -> f64 {
    ln_beta(a, b).exp_p()
}

/// Density of the standard normal distribution.
pub fn standard_normal_pdf(x: f64) -> f64 {
    let d = (x - 0.0) / 1.0;
    (-0.5 * d * d).exp_p() / (SQRT_2PI * 1.0)
}

/// Quantile of the standard normal distribution, as `statrs`'s
/// `Normal::inverse_cdf`. Panics outside `[0, 1]`.
pub fn standard_normal_inverse_cdf(p: f64) -> f64 {
    assert!((0.0..=1.0).contains(&p), "x must be in [0, 1]");
    0.0 - (1.0 * f64_consts::SQRT_2 * erfc_inv(2.0 * p))
}

fn erfc_inv(x: f64) -> f64 {
    if x <= 0.0 {
        f64::INFINITY
    } else if x >= 2.0 {
        f64::NEG_INFINITY
    } else if x > 1.0 {
        erf_inv_impl(-1.0 + x, 2.0 - x, -1.0)
    } else {
        erf_inv_impl(1.0 - x, x, 1.0)
    }
}

fn polynomial(z: f64, coeff: &[f64]) -> f64 {
    coeff.iter().rev().fold(0_f64, |acc, val| acc * z + val)
}

fn erf_inv_impl(p: f64, q: f64, s: f64) -> f64 {
    let result = if p <= 0.5 {
        let y = 0.0891314744949340820313;
        let g = p * (p + 10.0);
        let r = polynomial(p, ERF_INV_IMPL_AN) / polynomial(p, ERF_INV_IMPL_AD);
        g * y + g * r
    } else if q >= 0.25 {
        let y = 2.249481201171875;
        let g = (-2.0 * q.ln_p()).sqrt();
        let xs = q - 0.25;
        let r = polynomial(xs, ERF_INV_IMPL_BN) / polynomial(xs, ERF_INV_IMPL_BD);
        g / (y + r)
    } else {
        let x = (-q.ln_p()).sqrt();
        if x < 3.0 {
            let y = 0.807220458984375;
            let xs = x - 1.125;
            let r = polynomial(xs, ERF_INV_IMPL_CN) / polynomial(xs, ERF_INV_IMPL_CD);
            y * x + r * x
        } else if x < 6.0 {
            let y = 0.93995571136474609375;
            let xs = x - 3.0;
            let r = polynomial(xs, ERF_INV_IMPL_DN) / polynomial(xs, ERF_INV_IMPL_DD);
            y * x + r * x
        } else if x < 18.0 {
            let y = 0.98362827301025390625;
            let xs = x - 6.0;
            let r = polynomial(xs, ERF_INV_IMPL_EN) / polynomial(xs, ERF_INV_IMPL_ED);
            y * x + r * x
        } else if x < 44.0 {
            let y = 0.99714565277099609375;
            let xs = x - 18.0;
            let r = polynomial(xs, ERF_INV_IMPL_FN) / polynomial(xs, ERF_INV_IMPL_FD);
            y * x + r * x
        } else {
            let y = 0.99941349029541015625;
            let xs = x - 44.0;
            let r = polynomial(xs, ERF_INV_IMPL_GN) / polynomial(xs, ERF_INV_IMPL_GD);
            y * x + r * x
        }
    };
    s * result
}

/// A gamma distribution with shape and rate, as `statrs::distribution::Gamma`
/// (only the parts the crate uses).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GammaDistribution {
    shape: f64,
    rate: f64,
}

impl GammaDistribution {
    /// `None` for the parameters `statrs` rejects.
    pub fn new(shape: f64, rate: f64) -> Option<Self> {
        if shape.is_nan() || shape <= 0.0 || rate.is_nan() || rate <= 0.0 {
            return None;
        }
        if shape.is_infinite() && rate.is_infinite() {
            return None;
        }
        Some(Self { shape, rate })
    }

    /// The chi-squared distribution with `freedom` degrees of freedom.
    pub fn chi_squared(freedom: f64) -> Option<Self> {
        Self::new(freedom / 2.0, 0.5)
    }

    /// Degrees of freedom when this is a chi-squared distribution.
    pub fn freedom(&self) -> f64 {
        self.shape * 2.0
    }

    pub fn cdf(&self, x: f64) -> f64 {
        if x <= 0.0 {
            0.0
        } else if x == self.shape && self.rate.is_infinite() {
            1.0
        } else if self.rate.is_infinite() {
            0.0
        } else if x.is_infinite() {
            1.0
        } else {
            gamma_lr(self.shape, x * self.rate)
        }
    }

    pub fn sf(&self, x: f64) -> f64 {
        if x <= 0.0 {
            1.0
        } else if x == self.shape && self.rate.is_infinite() {
            0.0
        } else if self.rate.is_infinite() {
            1.0
        } else if x.is_infinite() {
            0.0
        } else {
            gamma_ur(self.shape, x * self.rate)
        }
    }

    fn ln_pdf(&self, x: f64) -> f64 {
        if x < 0.0 {
            f64::NEG_INFINITY
        } else if self.shape == 1.0 {
            self.rate.ln_p() - self.rate * x
        } else if x.is_infinite() {
            f64::NEG_INFINITY
        } else {
            self.shape * self.rate.ln_p() + (self.shape - 1.0) * x.ln_p()
                - self.rate * x
                - ln_gamma(self.shape)
        }
    }

    fn pdf(&self, x: f64) -> f64 {
        if x < 0.0 {
            0.0
        } else if self.shape == 1.0 {
            self.rate * (-self.rate * x).exp_p()
        } else if self.shape > 160.0 {
            self.ln_pdf(x).exp_p()
        } else if x.is_infinite() {
            0.0
        } else {
            self.rate.powf_p(self.shape) * x.powf_p(self.shape - 1.0) * (-self.rate * x).exp_p()
                / gamma(self.shape)
        }
    }

    /// Quantile by `statrs`'s safeguarded Newton iteration. Panics outside
    /// `[0, 1]`.
    pub fn inverse_cdf(&self, p: f64) -> f64 {
        assert!(
            (0.0..=1.0).contains(&p),
            "default inverse_cdf implementation should be provided probability on [0,1]"
        );
        if p == 0.0 {
            return 0.0;
        }
        if p == 1.0 {
            return f64::INFINITY;
        }
        let mut high = 2.0;
        let mut low = 1.0;
        while self.cdf(low) > p {
            low /= 2.0;
        }
        while self.cdf(high) < p {
            high *= 2.0;
        }
        const MAX_ITERATIONS: usize = 100;
        let mut x = (high + low) / 2.0;
        for _ in 0..MAX_ITERATIONS {
            let value = self.cdf(x);
            if value >= p {
                high = x;
            } else {
                low = x;
            }
            let mut next = x - (value - p) / self.pdf(x);
            if !next.is_finite() || next <= low || next >= high {
                next = (high + low) / 2.0;
            }
            if converged(&mut x, next) {
                break;
            }
        }
        x
    }
}

/// `statrs::prec::convergence`: `approx::relative_eq!` with epsilon 1e-9 and
/// maximum relative difference 1e-14, then advance.
fn converged(x: &mut f64, next: f64) -> bool {
    let a = *x;
    let b = next;
    let equal = if a == b {
        true
    } else if a.is_infinite() || b.is_infinite() {
        false
    } else {
        let difference = (a - b).abs();
        if difference <= 1e-9 {
            true
        } else {
            let largest = a.abs().max(b.abs());
            difference <= largest * 1e-14
        }
    };
    *x = next;
    equal
}

/// Polynomial coefficients for a numerator of `erf_inv_impl`
/// in the interval [0, 0.5].
const ERF_INV_IMPL_AN: &[f64] = &[
    -0.000508781949658280665617,
    -0.00836874819741736770379,
    0.0334806625409744615033,
    -0.0126926147662974029034,
    -0.0365637971411762664006,
    0.0219878681111168899165,
    0.00822687874676915743155,
    -0.00538772965071242932965,
];

/// Polynomial coefficients for a denominator of `erf_inv_impl`
/// in the interval [0, 0.5].
const ERF_INV_IMPL_AD: &[f64] = &[
    1.0,
    -0.970005043303290640362,
    -1.56574558234175846809,
    1.56221558398423026363,
    0.662328840472002992063,
    -0.71228902341542847553,
    -0.0527396382340099713954,
    0.0795283687341571680018,
    -0.00233393759374190016776,
    0.000886216390456424707504,
];

/// Polynomial coefficients for a numerator of `erf_inv_impl`
/// in the interval [0.5, 0.75].
const ERF_INV_IMPL_BN: &[f64] = &[
    -0.202433508355938759655,
    0.105264680699391713268,
    8.37050328343119927838,
    17.6447298408374015486,
    -18.8510648058714251895,
    -44.6382324441786960818,
    17.445385985570866523,
    21.1294655448340526258,
    -3.67192254707729348546,
];

/// Polynomial coefficients for a denominator of `erf_inv_impl`
/// in the interval [0.5, 0.75].
const ERF_INV_IMPL_BD: &[f64] = &[
    1.0,
    6.24264124854247537712,
    3.9713437953343869095,
    -28.6608180499800029974,
    -20.1432634680485188801,
    48.5609213108739935468,
    10.8268667355460159008,
    -22.6436933413139721736,
    1.72114765761200282724,
];

/// Polynomial coefficients for a numerator of `erf_inv_impl`
/// in the interval [0.75, 1] with x less than 3.
const ERF_INV_IMPL_CN: &[f64] = &[
    -0.131102781679951906451,
    -0.163794047193317060787,
    0.117030156341995252019,
    0.387079738972604337464,
    0.337785538912035898924,
    0.142869534408157156766,
    0.0290157910005329060432,
    0.00214558995388805277169,
    -0.679465575181126350155e-6,
    0.285225331782217055858e-7,
    -0.681149956853776992068e-9,
];

/// Polynomial coefficients for a denominator of `erf_inv_impl`
/// in the interval [0.75, 1] with x less than 3.
const ERF_INV_IMPL_CD: &[f64] = &[
    1.0,
    3.46625407242567245975,
    5.38168345707006855425,
    4.77846592945843778382,
    2.59301921623620271374,
    0.848854343457902036425,
    0.152264338295331783612,
    0.01105924229346489121,
];

/// Polynomial coefficients for a numerator of `erf_inv_impl`
/// in the interval [0.75, 1] with x between 3 and 6.
const ERF_INV_IMPL_DN: &[f64] = &[
    -0.0350353787183177984712,
    -0.00222426529213447927281,
    0.0185573306514231072324,
    0.00950804701325919603619,
    0.00187123492819559223345,
    0.000157544617424960554631,
    0.460469890584317994083e-5,
    -0.230404776911882601748e-9,
    0.266339227425782031962e-11,
];

/// Polynomial coefficients for a denominator of `erf_inv_impl`
/// in the interval [0.75, 1] with x between 3 and 6.
const ERF_INV_IMPL_DD: &[f64] = &[
    1.0,
    1.3653349817554063097,
    0.762059164553623404043,
    0.220091105764131249824,
    0.0341589143670947727934,
    0.00263861676657015992959,
    0.764675292302794483503e-4,
];

/// Polynomial coefficients for a numerator of `erf_inv_impl`
/// in the interval [0.75, 1] with x between 6 and 18.
const ERF_INV_IMPL_EN: &[f64] = &[
    -0.0167431005076633737133,
    -0.00112951438745580278863,
    0.00105628862152492910091,
    0.000209386317487588078668,
    0.149624783758342370182e-4,
    0.449696789927706453732e-6,
    0.462596163522878599135e-8,
    -0.281128735628831791805e-13,
    0.99055709973310326855e-16,
];

/// Polynomial coefficients for a denominator of `erf_inv_impl`
/// in the interval [0.75, 1] with x between 6 and 18.
const ERF_INV_IMPL_ED: &[f64] = &[
    1.0,
    0.591429344886417493481,
    0.138151865749083321638,
    0.0160746087093676504695,
    0.000964011807005165528527,
    0.275335474764726041141e-4,
    0.282243172016108031869e-6,
];

/// Polynomial coefficients for a numerator of `erf_inv_impl`
/// in the interval [0.75, 1] with x between 18 and 44.
const ERF_INV_IMPL_FN: &[f64] = &[
    -0.0024978212791898131227,
    -0.779190719229053954292e-5,
    0.254723037413027451751e-4,
    0.162397777342510920873e-5,
    0.396341011304801168516e-7,
    0.411632831190944208473e-9,
    0.145596286718675035587e-11,
    -0.116765012397184275695e-17,
];

/// Polynomial coefficients for a denominator of `erf_inv_impl`
/// in the interval [0.75, 1] with x between 18 and 44.
const ERF_INV_IMPL_FD: &[f64] = &[
    1.0,
    0.207123112214422517181,
    0.0169410838120975906478,
    0.000690538265622684595676,
    0.145007359818232637924e-4,
    0.144437756628144157666e-6,
    0.509761276599778486139e-9,
];

/// Polynomial coefficients for a numerator of `erf_inv_impl`
/// in the interval [0.75, 1] with x greater than 44.
const ERF_INV_IMPL_GN: &[f64] = &[
    -0.000539042911019078575891,
    -0.28398759004727721098e-6,
    0.899465114892291446442e-6,
    0.229345859265920864296e-7,
    0.225561444863500149219e-9,
    0.947846627503022684216e-12,
    0.135880130108924861008e-14,
    -0.348890393399948882918e-21,
];

/// Polynomial coefficients for a denominator of `erf_inv_impl`
/// in the interval [0.75, 1] with x greater than 44.
const ERF_INV_IMPL_GD: &[f64] = &[
    1.0,
    0.0845746234001899436914,
    0.00282092984726264681981,
    0.468292921940894236786e-4,
    0.399968812193862100054e-6,
    0.161809290887904476097e-8,
    0.231558608310259605225e-11,
];

#[cfg(test)]
mod tests {
    use super::*;
    use crate::numeric::profile::host_libm;

    fn grid() -> Vec<f64> {
        let mut values = vec![
            1e-300, 1e-12, 1e-6, 0.1, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0,
        ];
        let mut x = 0.013;
        while x < 2_000.0 {
            values.push(x);
            x *= 1.37;
        }
        values
    }

    // Under the R profiles every port is the statrs function, bit for bit.
    #[test]
    fn ports_match_statrs_under_host_libm() {
        if !host_libm() {
            return;
        }
        use statrs::distribution::{ChiSquared, Continuous, ContinuousCDF, Normal};
        for &x in &grid() {
            assert_eq!(
                ln_gamma(x).to_bits(),
                statrs::function::gamma::ln_gamma(x).to_bits(),
                "ln_gamma {x}"
            );
            assert_eq!(
                ln_gamma(-x + 0.25).to_bits(),
                statrs::function::gamma::ln_gamma(-x + 0.25).to_bits(),
                "ln_gamma {x}"
            );
            assert_eq!(
                gamma(x.min(150.0)).to_bits(),
                statrs::function::gamma::gamma(x.min(150.0)).to_bits(),
                "gamma {x}"
            );
            assert_eq!(
                digamma(x).to_bits(),
                statrs::function::gamma::digamma(x).to_bits(),
                "digamma {x}"
            );
            assert_eq!(
                digamma(-x - 0.3).to_bits(),
                statrs::function::gamma::digamma(-x - 0.3).to_bits(),
                "digamma {x}"
            );
            assert_eq!(
                standard_normal_pdf(x - 3.0).to_bits(),
                Normal::standard().pdf(x - 3.0).to_bits()
            );
            for &a in &[0.3, 1.0, 2.5, 7.0, 40.0, 400.0] {
                assert_eq!(
                    gamma_lr(a, x).to_bits(),
                    statrs::function::gamma::gamma_lr(a, x).to_bits(),
                    "gamma_lr {a} {x}"
                );
                assert_eq!(
                    gamma_ur(a, x).to_bits(),
                    statrs::function::gamma::gamma_ur(a, x).to_bits(),
                    "gamma_ur {a} {x}"
                );
                assert_eq!(
                    ln_beta(a, x).to_bits(),
                    statrs::function::beta::ln_beta(a, x).to_bits()
                );
                assert_eq!(
                    beta(a, x.min(50.0)).to_bits(),
                    statrs::function::beta::beta(a, x.min(50.0)).to_bits()
                );
                let ours = GammaDistribution::chi_squared(a).unwrap();
                let theirs = ChiSquared::new(a).unwrap();
                assert_eq!(
                    ours.sf(x).to_bits(),
                    theirs.sf(x).to_bits(),
                    "chisq sf {a} {x}"
                );
                assert_eq!(
                    ours.cdf(x).to_bits(),
                    theirs.cdf(x).to_bits(),
                    "chisq cdf {a} {x}"
                );
            }
        }
        for i in 1..200 {
            let p = f64::from(i) / 200.0;
            assert_eq!(
                standard_normal_inverse_cdf(p).to_bits(),
                Normal::standard().inverse_cdf(p).to_bits(),
                "qnorm {p}"
            );
            for &df in &[1.0, 2.0, 5.0, 30.0] {
                assert_eq!(
                    GammaDistribution::chi_squared(df)
                        .unwrap()
                        .inverse_cdf(p)
                        .to_bits(),
                    ChiSquared::new(df).unwrap().inverse_cdf(p).to_bits(),
                    "qchisq {df} {p}"
                );
            }
        }
        for p in [1e-300, 1e-100, 1e-20, 1e-8, 0.999_999_9] {
            assert_eq!(
                standard_normal_inverse_cdf(p).to_bits(),
                Normal::standard().inverse_cdf(p).to_bits()
            );
        }
    }
}
