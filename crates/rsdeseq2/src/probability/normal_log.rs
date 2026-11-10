//! Rational normal approximations from Cody (Math. Comp., 1969), evaluated
//! directly in log space. Splitting the square preserves low-order information
//! needed when posterior moments subtract nearby log probabilities.
#![allow(clippy::excessive_precision)]

fn polynomial(x: f64, ascending: &[f64]) -> f64 {
    let mut value = *ascending.last().unwrap();
    for &coefficient in ascending[..ascending.len() - 1].iter().rev() {
        value = value * x + coefficient;
    }
    value
}
fn tail_factor(q: f64) -> f64 {
    if q <= 32.0_f64.sqrt() {
        let numerator = [
            9842.7148383839780218,
            11602.651437647350124,
            6848.1904505362823326,
            2494.5375852903726711,
            597.27027639480026226,
            93.506656132177855979,
            8.8831497943883759412,
            0.39894151208813466764,
            1.0765576773720192317e-8,
        ];
        let denominator = [
            19685.429676859990727,
            38912.003286093271411,
            34900.952721145977266,
            18615.571640885098091,
            6485.558298266760755,
            1519.377599407554805,
            235.38790178262499861,
            22.266688044328115691,
            1.0,
        ];
        polynomial(q, &numerator) / polynomial(q, &denominator)
    } else {
        let reciprocal_square = 1.0 / (q * q);
        let numerator = [
            2.9112874951168792e-5,
            0.001421619193227893466,
            0.022235277870649807,
            0.1274011611602473639,
            0.21589853405795699,
            0.02307344176494017303,
        ];
        let denominator = [
            7.29751555083966205e-5,
            0.00378239633202758244,
            0.0659881378689285515,
            0.468238212480865118,
            1.28426009614491121,
            1.0,
        ];
        (0.39894228040143267794
            - reciprocal_square * polynomial(reciprocal_square, &numerator)
                / polynomial(reciprocal_square, &denominator))
            / q
    }
}
pub(crate) fn normal_log_upper_tail(q: f64) -> f64 {
    if q.is_nan() {
        return q;
    }
    let magnitude = q.abs();
    if magnitude <= 0.67448975 {
        let square = q * q;
        let numerator = [
            18154.981253343561249,
            1067.6894854603709582,
            161.02823106855587881,
            2.2352520354606839287,
            0.065682337918207449113,
        ];
        let denominator = [
            45507.789335026729956,
            10260.932208618978205,
            976.09855173777669322,
            47.20258190468824187,
            1.0,
        ];
        return (0.5 - q * polynomial(square, &numerator) / polynomial(square, &denominator)).ln();
    }
    if magnitude > 1e155 {
        return if q > 0.0 { f64::NEG_INFINITY } else { -0.0 };
    }
    let coarse = (magnitude * 16.0).trunc() / 16.0;
    let leading = -coarse * (coarse / 2.0);
    let remainder = -((magnitude - coarse) * (magnitude + coarse)) / 2.0;
    let factor = tail_factor(magnitude);
    if q > 0.0 {
        leading + remainder + factor.ln()
    } else {
        (-leading.exp() * remainder.exp() * factor).ln_1p()
    }
}
