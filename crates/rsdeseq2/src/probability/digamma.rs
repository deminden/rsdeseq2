//! Positive digamma via DLMF 5.5 recurrence and 5.11 Bernoulli expansion.
pub(crate) fn digamma_positive(mut x: f64) -> f64 {
    if x == f64::INFINITY {
        return x;
    }
    if !x.is_finite() || x <= 0.0 {
        return f64::NAN;
    }
    let mut correction = 0.0;
    while x < 32.0 {
        correction -= 1.0 / x;
        x += 1.0;
    }
    let inv = 1.0 / x;
    let square = inv * inv;
    correction + x.ln()
        - 0.5 * inv
        - square
            * (1.0 / 12.0
                + square
                    * (-1.0 / 120.0
                        + square
                            * (1.0 / 252.0
                                + square
                                    * (-1.0 / 240.0
                                        + square * (1.0 / 132.0 - square * (691.0 / 32760.0))))))
}
