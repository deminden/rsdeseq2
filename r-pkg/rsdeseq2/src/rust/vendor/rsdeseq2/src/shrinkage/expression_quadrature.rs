//! Adaptive Gauss integration for imported elementwise density expressions.
use crate::probability::quadrature::GAUSS_LEGENDRE_16;
fn panel(f: &impl Fn(f64) -> f64, a: f64, b: f64) -> f64 {
    let middle = (a + b) / 2.0;
    let half = (b - a) / 2.0;
    GAUSS_LEGENDRE_16
        .iter()
        .map(|&(x, w)| half * w * (f(middle - half * x) + f(middle + half * x)))
        .sum()
}
struct Interval {
    a: f64,
    b: f64,
    value: f64,
    error: f64,
}
fn interval(f: &impl Fn(f64) -> f64, a: f64, b: f64) -> Interval {
    let middle = (a + b) / 2.0;
    let whole = panel(f, a, b);
    let split = panel(f, a, middle) + panel(f, middle, b);
    Interval {
        a,
        b,
        value: split,
        error: (split - whole).abs(),
    }
}
pub(super) fn integrate(f: impl Fn(f64) -> f64, a: f64, b: f64) -> Option<f64> {
    if a.is_nan() || b.is_nan() {
        return None;
    }
    if a > b {
        return integrate(f, b, a).map(|v| -v);
    }
    if a == b {
        return Some(0.0);
    }
    // Nodes lie strictly inside the interval, including under infinite-bound maps.
    let transformed = |t: f64| {
        if a.is_finite() && b.is_finite() {
            return f(t);
        }
        if a.is_infinite() && b.is_infinite() {
            let q = (std::f64::consts::PI * (t - 0.5)).tan();
            return f(q) * std::f64::consts::PI * (1.0 + q * q);
        }
        let q = t / (1.0 - t);
        let x = if a.is_finite() { a + q } else { b - q };
        f(x) / (1.0 - t).powi(2)
    };
    let (a, b) = if a.is_finite() && b.is_finite() {
        (a, b)
    } else {
        (0.0, 1.0)
    };
    let mut intervals = vec![interval(&transformed, a, b)];
    for _ in 0..100 {
        let value: f64 = intervals.iter().map(|p| p.value).sum();
        let error: f64 = intervals.iter().map(|p| p.error).sum();
        if !value.is_finite() || !error.is_finite() {
            return None;
        }
        if error <= 1e-11_f64.max(1e-11 * value.abs()) {
            return Some(value);
        }
        let index = intervals
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.error.total_cmp(&b.1.error))?
            .0;
        let old = intervals.swap_remove(index);
        let middle = (old.a + old.b) / 2.0;
        if middle == old.a || middle == old.b {
            return None;
        }
        intervals.push(interval(&transformed, old.a, middle));
        intervals.push(interval(&transformed, middle, old.b));
    }
    None
}
