//! Fifteen-point Gauss-Kronrod integration of infinite log-F moment intervals.
//! Nodes/weights are mathematical quadrature constants; the QUADPACK reference
//! is https://www.netlib.org/quadpack/dqk15i.f . No R routine is called or linked.
const NODES: [f64; 8] = [
    0.9914553711208126,
    0.9491079123427585,
    0.8648644233597691,
    0.7415311855993945,
    0.5860872354676911,
    0.4058451513773972,
    0.2077849550078985,
    0.0,
];
const KRONROD: [f64; 8] = [
    0.02293532201052922,
    0.06309209262997855,
    0.1047900103222502,
    0.1406532597155259,
    0.1690047266392679,
    0.1903505780647854,
    0.2044329400752989,
    0.2094821410847278,
];
const GAUSS: [f64; 4] = [
    0.1294849661688697,
    0.2797053914892767,
    0.3818300505051189,
    0.4179591836734694,
];
struct Panel {
    a: f64,
    b: f64,
    value: f64,
    error: f64,
}
fn panel(f: &impl Fn(f64) -> f64, a: f64, b: f64) -> Option<Panel> {
    let middle = (a + b) / 2.0;
    let half = (b - a) / 2.0;
    let center = f(middle);
    if !center.is_finite() {
        return None;
    }
    let mut values = [(0.0, 0.0); 7];
    let mut k = KRONROD[7] * center;
    let mut g = GAUSS[3] * center;
    let mut absolute = KRONROD[7] * center.abs();
    for i in 0..7 {
        let offset = half * NODES[i];
        let first = f(middle - offset);
        let second = f(middle + offset);
        if !first.is_finite() || !second.is_finite() {
            return None;
        }
        values[i] = (first, second);
        k += KRONROD[i] * (first + second);
        absolute += KRONROD[i] * (first.abs() + second.abs());
        if i % 2 == 1 {
            g += GAUSS[i / 2] * (first + second);
        }
    }
    let mean = k / 2.0;
    let mut ascending = KRONROD[7] * (center - mean).abs();
    for i in 0..7 {
        ascending += KRONROD[i] * ((values[i].0 - mean).abs() + (values[i].1 - mean).abs());
    }
    let value = k * half;
    ascending *= half;
    absolute *= half;
    let mut error = ((k - g) * half).abs();
    if ascending != 0.0 && error != 0.0 {
        error = ascending * (200.0 * error / ascending).powf(1.5).min(1.0);
    }
    error = error.max(50.0 * f64::EPSILON * absolute);
    Some(Panel { a, b, value, error })
}
pub(super) fn integrate_infinite(f: impl Fn(f64) -> f64, a: f64, b: f64) -> Option<f64> {
    let transformed = |t: f64| {
        let x = (1.0 - t) / t;
        let value = if a.is_infinite() && b.is_infinite() {
            f(x) + f(-x)
        } else if a.is_finite() {
            f(a + x)
        } else {
            f(b - x)
        };
        value / (t * t)
    };
    let mut panels = vec![panel(&transformed, 0.0, 1.0)?];
    let tolerance = f64::EPSILON.sqrt().sqrt();
    for _ in 0..100 {
        let value: f64 = panels.iter().map(|p| p.value).sum();
        let error: f64 = panels.iter().map(|p| p.error).sum();
        if error <= tolerance.max(tolerance * value.abs()) {
            return Some(value);
        }
        let i = panels
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.error.total_cmp(&b.1.error))?
            .0;
        let p = panels.swap_remove(i);
        let mid = (p.a + p.b) / 2.0;
        panels.push(panel(&transformed, p.a, mid)?);
        panels.push(panel(&transformed, mid, p.b)?);
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn kronrod_constants_integrate_polynomial_and_infinite_exponential_truths() {
        for degree in 0..=22 {
            let p = panel(&|x| x.powi(degree), -1.0, 1.0).unwrap();
            let truth = if degree % 2 == 0 {
                2.0 / f64::from(degree + 1)
            } else {
                0.0
            };
            assert!((p.value - truth).abs() < 5e-15);
        }
        let value = integrate_infinite(|x| (-x).exp(), 0.0, f64::INFINITY).unwrap();
        assert!((value - 1.0).abs() < 1e-4);
        assert!(integrate_infinite(|_| f64::NAN, 0.0, f64::INFINITY).is_none());
    }
}
