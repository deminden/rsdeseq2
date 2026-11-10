//! Heavy-tailed coefficient prior and empirical variance adaptation.
use crate::DeseqError;
pub(super) fn invalid(reason: &str) -> DeseqError {
    DeseqError::InvalidOptions {
        reason: format!("heavy-tailed shrinkage: {reason}"),
    }
}
/// Independent Cauchy priors, with selected Gaussian nuisance coefficients.
#[derive(Clone, Debug)]
pub struct CauchyCoefficientPrior {
    /// Zero-based coefficient indices with Gaussian rather than Cauchy priors.
    pub nuisance: Vec<usize>,
    pub scale: f64,
    pub nuisance_scale: f64,
    /// Empirical variance, when adapted from coefficient MLEs and SEs.
    pub adapted_variance: Option<f64>,
}
impl Default for CauchyCoefficientPrior {
    fn default() -> Self {
        Self {
            nuisance: vec![0],
            scale: 1.0,
            nuisance_scale: 15.0,
            adapted_variance: None,
        }
    }
}
impl CauchyCoefficientPrior {
    pub(super) fn validate(&self, coefficients: usize) -> Result<(), DeseqError> {
        if !self.scale.is_finite()
            || self.scale <= 0.0
            || !self.nuisance_scale.is_finite()
            || self.nuisance_scale <= 0.0
            || !(self.scale * self.scale).is_finite()
            || self.scale * self.scale == 0.0
            || !(self.nuisance_scale * self.nuisance_scale).is_finite()
            || self.nuisance_scale * self.nuisance_scale == 0.0
            || self.nuisance.iter().any(|&i| i >= coefficients)
        {
            return Err(invalid(
                "finite positive prior scales and valid nuisance indices are required",
            ));
        }
        Ok(())
    }
    pub(super) fn contribution(&self, beta: &[f64], gradient: &mut [f64]) -> f64 {
        let scale2 = self.scale * self.scale;
        let nuisance2 = self.nuisance_scale * self.nuisance_scale;
        beta.iter()
            .zip(gradient)
            .enumerate()
            .map(|(j, (&b, g))| {
                if self.nuisance.contains(&j) {
                    *g += b / nuisance2;
                    b * b / (2.0 * nuisance2)
                } else {
                    *g += (2.0 * b) / (scale2 + b * b);
                    (b * b / scale2).ln_1p()
                }
            })
            .sum()
    }
}

// Bracket-preserving Brent interpolation, with the reference default absolute
// tolerance EPS^(1/4). This is a standalone implementation of the root algorithm.
fn root(f: impl Fn(f64) -> f64, mut a: f64, mut b: f64, tolerance: f64) -> Result<f64, DeseqError> {
    let mut fa = f(a);
    let mut fb = f(b);
    if fa == 0.0 {
        return Ok(a);
    }
    if fb == 0.0 {
        return Ok(b);
    }
    if !fa.is_finite() || !fb.is_finite() || fa.signum() == fb.signum() {
        return Err(invalid(
            "prior variance root must have finite opposite endpoint signs",
        ));
    }
    let mut c = a;
    let mut fc = fa;
    let mut d = b - a;
    let mut e = d;
    for _ in 0..1000 {
        if fb.signum() == fc.signum() {
            c = a;
            fc = fa;
            d = b - a;
            e = d;
        }
        if fc.abs() < fb.abs() {
            a = b;
            b = c;
            c = a;
            fa = fb;
            fb = fc;
            fc = fa;
        }
        let tol = 2.0 * f64::EPSILON * b.abs() + 0.5 * tolerance;
        let half = 0.5 * (c - b);
        if half.abs() <= tol || fb == 0.0 {
            return Ok(b);
        }
        if e.abs() >= tol && fa.abs() > fb.abs() {
            let ratio = fb / fa;
            let (mut p, mut q) = if a == c {
                (2.0 * half * ratio, 1.0 - ratio)
            } else {
                let q = fa / fc;
                let r = fb / fc;
                (
                    ratio * (2.0 * half * q * (q - r) - (b - a) * (r - 1.0)),
                    (q - 1.0) * (r - 1.0) * (ratio - 1.0),
                )
            };
            if p > 0.0 {
                q = -q;
            } else {
                p = -p;
            }
            if 2.0 * p < (3.0 * half * q - (tol * q).abs()).min((e * q).abs()) {
                e = d;
                d = p / q;
            } else {
                d = half;
                e = d;
            }
        } else {
            d = half;
            e = d;
        }
        a = b;
        fa = fb;
        b += if d.abs() > tol { d } else { tol.copysign(half) };
        fb = f(b);
        if !fb.is_finite() {
            return Err(invalid("nonfinite prior variance objective"));
        }
    }
    Err(invalid("prior variance root did not converge"))
}
/// Empirical Gaussian variance underlying adaptive Cauchy scale selection.
/// Missing MLEs are omitted. Missing/infinite SEs of retained MLEs are errors.
/// Returns the reference lower boundary or a bracketed root on [1e-6,400].
pub fn estimate_cauchy_prior_variance(
    estimates: &[f64],
    standard_errors: &[f64],
) -> Result<f64, DeseqError> {
    if estimates.len() != standard_errors.len() {
        return Err(invalid("MLEs and SEs must have matching lengths"));
    }
    let mut squared = vec![];
    for (&x, &s) in estimates.iter().zip(standard_errors) {
        if x.is_nan() {
            continue;
        }
        if !x.is_finite() || !s.is_finite() || !(x * x).is_finite() || !(s * s).is_finite() {
            return Err(invalid(
                "retained MLEs and SEs must be finite and squareable",
            ));
        }
        squared.push((x * x, s * s));
    }
    if squared.is_empty() {
        return Err(invalid("at least one nonmissing MLE is required"));
    }
    let objective = |a: f64| {
        let mut numerator = 0.0;
        let mut denominator = 0.0;
        for &(s, d) in &squared {
            let weight = 1.0 / (2.0 * (a + d).powi(2));
            numerator += (s - d) * weight;
            denominator += weight;
        }
        numerator / denominator - a
    };
    let lower = 1e-6;
    if objective(lower) < 0.0 {
        Ok(lower)
    } else {
        root(objective, lower, 400.0, f64::EPSILON.sqrt().sqrt())
    }
}
