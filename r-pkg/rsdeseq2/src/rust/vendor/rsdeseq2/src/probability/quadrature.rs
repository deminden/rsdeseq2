//! Fixed Gauss-Legendre quadrature for scaled distribution tails.

// Positive nodes/weights of the 16-point Gauss-Legendre rule, generated using
// 90-digit Legendre recurrence/Newton iteration, then rounded to f64.
pub(crate) const GAUSS_LEGENDRE_16: [(f64, f64); 8] = [
    (0.9894009349916499, 0.027152459411754096),
    (0.9445750230732326, 0.062253523938647894),
    (0.8656312023878318, 0.09515851168249279),
    (0.755404408355003, 0.12462897125553388),
    (0.6178762444026438, 0.14959598881657674),
    (0.45801677765722737, 0.16915651939500254),
    (0.2816035507792589, 0.18260341504492358),
    (0.09501250983763744, 0.1894506104550685),
];

// Integrate a nonnegative, rapidly decreasing tail in a caller-scaled
// coordinate, using successively doubled panels and compensated summation.
pub(crate) fn integrate_positive_tail(density: impl Fn(f64) -> f64) -> f64 {
    let mut integral = 0.0;
    let mut compensation = 0.0;
    let mut left = 0.0;
    for right in [1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0] {
        let center = 0.5 * (left + right);
        let half_width = 0.5 * (right - left);
        for (node, weight) in GAUSS_LEGENDRE_16 {
            let offset = half_width * node;
            let contribution =
                half_width * weight * (density(center - offset) + density(center + offset));
            let delta = contribution - compensation;
            let next = integral + delta;
            compensation = (next - integral) - delta;
            integral = next;
        }
        left = right;
    }
    integral
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gauss_rule_reproduces_analytic_polynomial_moments_through_degree_31() {
        for degree in 0..=31 {
            let actual = GAUSS_LEGENDRE_16
                .iter()
                .map(|&(x, w)| w * (x.powi(degree) + (-x).powi(degree)))
                .sum::<f64>();
            let expected = if degree % 2 == 0 {
                2.0 / f64::from(degree + 1)
            } else {
                0.0
            };
            assert!(
                (actual - expected).abs() <= 4.0 * f64::EPSILON,
                "degree {degree}: {actual} != {expected}"
            );
        }
    }

    #[test]
    fn scaled_tail_integrator_matches_exponential_and_gaussian_integrals() {
        let exponential = integrate_positive_tail(|x| (-x).exp());
        let gaussian = integrate_positive_tail(|x| 0.3989422804014327 * (-0.5 * x * x).exp());
        assert!((exponential - 1.0).abs() <= 2.0 * f64::EPSILON);
        assert!((gaussian - 0.5).abs() <= f64::EPSILON);
    }
}
