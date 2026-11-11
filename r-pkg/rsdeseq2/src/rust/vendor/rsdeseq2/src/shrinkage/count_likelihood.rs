//! Count-likelihood kernels expressed as gamma/beta integrals. No R callbacks.
use super::{beta_normalization::log_regularized_beta, gamma_probability};
use crate::DeseqError;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PoissonLink {
    Identity,
    Log,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BinomialLink {
    Identity,
    Logit,
}
#[derive(Clone, Debug)]
pub(super) enum CountLikelihood {
    Poisson {
        y: Vec<f64>,
        scale: Vec<f64>,
        link: PoissonLink,
    },
    Binomial {
        y: Vec<f64>,
        trials: Vec<f64>,
        link: BinomialLink,
    },
}
fn invalid(s: &str) -> DeseqError {
    DeseqError::InvalidOptions {
        reason: format!("count likelihood: {s}"),
    }
}
fn at(v: &[f64], row: usize) -> f64 {
    v[if v.len() == 1 { 0 } else { row }]
}
fn vector(v: &[f64], positive: bool) -> bool {
    !v.is_empty()
        && v.iter()
            .all(|x| x.is_finite() && if positive { *x > 0.0 } else { *x >= 0.0 })
}
fn beta_cdf(x: f64, a: f64, b: f64) -> f64 {
    if x > (a + 1.0) / (a + b + 2.0) {
        (0.5 - log_regularized_beta(b, a, 1.0 - x).exp()) + 0.5
    } else {
        log_regularized_beta(a, b, x).exp()
    }
}
fn gamma_density(x: f64, a: f64, rate: f64) -> f64 {
    if x.is_nan() {
        return f64::NAN;
    }
    if x < 0.0 || x.is_infinite() {
        return f64::NEG_INFINITY;
    }
    if x == 0.0 {
        return if a < 1.0 {
            f64::INFINITY
        } else if a == 1.0 {
            rate.ln()
        } else {
            f64::NEG_INFINITY
        };
    }
    if a >= 1e5 {
        return super::large_shape::gamma_log_density(x, a, rate);
    }
    a * rate.ln() + (a - 1.0) * x.ln() - rate * x - statrs::function::gamma::ln_gamma(a)
}
fn beta_density(x: f64, a: f64, b: f64) -> f64 {
    if x.is_nan() {
        return f64::NAN;
    }
    if !(0.0..=1.0).contains(&x) {
        return f64::NEG_INFINITY;
    }
    if x == 0.0 {
        return if a < 1.0 {
            f64::INFINITY
        } else if a == 1.0 {
            b.ln()
        } else {
            f64::NEG_INFINITY
        };
    }
    if x == 1.0 {
        return if b < 1.0 {
            f64::INFINITY
        } else if b == 1.0 {
            a.ln()
        } else {
            f64::NEG_INFINITY
        };
    }
    if super::large_shape::beta_large(a, b) {
        return super::large_shape::beta_log_density(x, a, b);
    }
    let normalization = if a.max(b) >= 1e5 {
        super::large_shape::log_beta(a, b)
    } else {
        statrs::function::beta::ln_beta(a, b)
    };
    (a - 1.0) * x.ln() + (b - 1.0) * (-x).ln_1p() - normalization
}
// Source zero-mass intervals select the denser endpoint (right on a tie).
fn moment(
    a: f64,
    b: f64,
    order: u8,
    factor: f64,
    cdf: impl Fn(f64, u8) -> f64,
    density: impl Fn(f64) -> f64,
) -> f64 {
    let mass = cdf(b, 0) - cdf(a, 0);
    if mass == 0.0 {
        let x = if density(a) > density(b) { a } else { b };
        return x.powi(i32::from(order));
    }
    if a == b {
        return a.powi(i32::from(order));
    }
    factor * (cdf(b, order) - cdf(a, order)) / mass
}
impl CountLikelihood {
    pub(super) fn poisson(
        y: Vec<f64>,
        scale: Vec<f64>,
        link: PoissonLink,
    ) -> Result<Self, DeseqError> {
        if !vector(&y, false) || !vector(&scale, true) {
            return Err(invalid(
                "finite nonnegative counts and positive scales are required",
            ));
        }
        Ok(Self::Poisson { y, scale, link })
    }
    pub(super) fn binomial(
        y: Vec<f64>,
        trials: Vec<f64>,
        link: BinomialLink,
    ) -> Result<Self, DeseqError> {
        if !vector(&y, false) || !vector(&trials, false) {
            return Err(invalid(
                "finite nonnegative counts and trial numbers are required",
            ));
        }
        let n = y.len().max(trials.len());
        if (y.len() != 1 && y.len() != n)
            || (trials.len() != 1 && trials.len() != n)
            || (0..n).any(|i| at(&y, i) > at(&trials, i))
        {
            return Err(invalid("aligned counts must not exceed trials"));
        }
        Ok(Self::Binomial { y, trials, link })
    }
    pub(super) fn grid_data(&self, n: usize, se: &[f64]) -> (Vec<f64>, Vec<f64>) {
        match self {
            Self::Poisson { y, scale, .. } => {
                let x = (0..n).map(|i| at(y, i) / at(scale, i)).collect::<Vec<_>>();
                let s = x
                    .iter()
                    .enumerate()
                    .map(|(i, x)| x.sqrt() / at(scale, i))
                    .collect();
                (x, s)
            }
            Self::Binomial { y, .. } => ((0..n).map(|i| at(y, i)).collect(), se.to_vec()),
        }
    }
    pub(super) fn mode_bounds(&self, n: usize, range: (f64, f64), seed: Option<f64>) -> (f64, f64) {
        match self {
            Self::Poisson { y, scale, link } => {
                let mut lo = f64::INFINITY;
                let mut hi = f64::NEG_INFINITY;
                for i in 0..n {
                    let v = at(y, i) / at(scale, i);
                    let v = if *link == PoissonLink::Log {
                        (v + 0.01).ln()
                    } else {
                        v
                    };
                    lo = lo.min(v);
                    hi = hi.max(v);
                }
                if *link == PoissonLink::Identity
                    && let Some(seed) = seed
                {
                    lo = lo.min(seed);
                    hi = hi.max(seed);
                }
                (lo, hi)
            }
            Self::Binomial { y, trials, link } => {
                if *link == BinomialLink::Identity {
                    return (range.0.max(0.0), range.1.min(1.0));
                }
                let mut lo = f64::INFINITY;
                let mut hi = f64::NEG_INFINITY;
                for i in 0..n {
                    let p = (at(y, i) + 0.01) / (at(trials, i) + 0.02);
                    let v = (p / (1.0 - p)).ln();
                    lo = lo.min(v);
                    hi = hi.max(v);
                }
                (lo, hi)
            }
        }
    }
    pub(super) fn validate_captures(&self, y: &[f64], n: Option<&[f64]>) -> Result<(), DeseqError> {
        let (observed, trials) = match self {
            Self::Poisson {
                y,
                link: PoissonLink::Log,
                ..
            } => (y, None),
            Self::Binomial {
                y,
                trials,
                link: BinomialLink::Logit,
            } => (y, Some(trials)),
            _ => return Err(invalid("adjusted captures require a log/logit likelihood")),
        };
        if y.len() != observed.len() || !y.iter().zip(observed).all(|(a, b)| *a == *b + 1e-5) {
            return Err(invalid("modified adjusted count capture"));
        }
        match (n, trials) {
            (None, None) => Ok(()),
            (Some(a), Some(b))
                if a.len() == b.len() && a.iter().zip(b).all(|(a, b)| *a == *b + 2e-5) =>
            {
                Ok(())
            }
            _ => Err(invalid("modified adjusted trials capture")),
        }
    }
    pub(super) fn lengths(&self) -> (usize, usize) {
        match self {
            Self::Poisson { y, scale, .. } => (y.len(), scale.len()),
            Self::Binomial { y, trials, .. } => (y.len(), trials.len()),
        }
    }
    pub(super) fn name(&self) -> &'static str {
        match self {
            Self::Poisson { .. } => "pois",
            Self::Binomial { .. } => "binom",
        }
    }
    pub(super) fn support(&self) -> Option<(f64, f64)> {
        match self {
            Self::Poisson {
                link: PoissonLink::Identity,
                ..
            } => Some((0.0, f64::INFINITY)),
            Self::Binomial {
                link: BinomialLink::Identity,
                ..
            } => Some((0.0, 1.0)),
            _ => None,
        }
    }
    pub(super) fn log_cdf(&self, x: f64, row: usize) -> f64 {
        match self {
            Self::Poisson { y, scale, link } => {
                let y = at(y, row);
                let scale = at(scale, row);
                let (x, a, normalization) = match link {
                    PoissonLink::Identity => (x.abs(), y + 1.0, scale.ln()),
                    PoissonLink::Log => ((-x).exp(), y + 1e-5, (y + 1e-5).ln()),
                };
                gamma_probability::log_gamma_probability(a, x * scale, true) - normalization
            }
            Self::Binomial { y, trials, link } => {
                let y = at(y, row);
                let n = at(trials, row);
                let (x, a, b, normalization) = match link {
                    BinomialLink::Identity => (x.abs(), y + 1.0, n - y + 1.0, -(n + 1.0).ln()),
                    BinomialLink::Logit => {
                        let a = y + 1e-5;
                        let n = n + 2e-5;
                        (1.0 / (1.0 + (-x).exp()), a, n - a, (n / (a * (n - a))).ln())
                    }
                };
                log_regularized_beta(a, b, x) + normalization
            }
        }
    }
    pub(super) fn log_density(&self, x: f64, row: usize) -> f64 {
        match self {
            Self::Poisson { y, scale, link } => {
                let y = at(y, row);
                let r = at(scale, row);
                match link {
                    PoissonLink::Identity => gamma_density(x.abs(), y + 1.0, r) - r.ln(),
                    PoissonLink::Log => gamma_density((-x).exp(), y + 1e-5, r) - (y + 1e-5).ln(),
                }
            }
            Self::Binomial { y, trials, link } => {
                let y = at(y, row);
                let n = at(trials, row);
                match link {
                    BinomialLink::Identity => {
                        beta_density(x.abs(), y + 1.0, n - y + 1.0) - (n + 1.0).ln()
                    }
                    BinomialLink::Logit => {
                        let a = y + 1e-5;
                        let n = n + 2e-5;
                        beta_density(1.0 / (1.0 + (-x).exp()), a, n - a) + (n / (a * (n - a))).ln()
                    }
                }
            }
        }
    }
    pub(super) fn moment(&self, a: f64, b: f64, order: u8, row: usize) -> f64 {
        let value = match self {
            Self::Poisson { y, scale, link } => {
                let y = at(y, row);
                let r = at(scale, row);
                let (a, b, shape) = match link {
                    PoissonLink::Identity => (-b, -a, y + 1.0),
                    PoissonLink::Log => ((-b).exp(), (-a).exp(), y + 1e-5),
                };
                let factor = if order == 1 {
                    shape / r
                } else {
                    shape * (shape + 1.0) / (r * r)
                };
                moment(
                    a,
                    b,
                    order,
                    factor,
                    |x, k| gamma_probability::gamma_probability(shape + f64::from(k), x * r, true),
                    |x| gamma_density(x, shape, r),
                )
            }
            Self::Binomial { y, trials, link } => {
                let y = at(y, row);
                let n = at(trials, row);
                let (a, b, alpha, beta) = match link {
                    BinomialLink::Identity => (-b, -a, y + 1.0, n - y + 1.0),
                    BinomialLink::Logit => {
                        let alpha = y + 1e-5;
                        (
                            1.0 / (1.0 + (-b).exp()),
                            1.0 / (1.0 + (-a).exp()),
                            alpha,
                            (n + 2e-5) - alpha,
                        )
                    }
                };
                let total = alpha + beta;
                let factor = if order == 1 {
                    alpha / total
                } else {
                    alpha * (alpha + 1.0) / (total * (total + 1.0))
                };
                moment(
                    a,
                    b,
                    order,
                    factor,
                    |x, k| beta_cdf(x, alpha + f64::from(k), beta),
                    |x| beta_density(x, alpha, beta),
                )
            }
        };
        if order == 1 { -value } else { value }
    }
}
