//! R scalar reductions using IEEE binary128 on Linux arm64.

mod wide;
use wide::Wide;

fn sum(values: impl Iterator<Item = f64>) -> Wide {
    values.map(Wide::from_double).fold(Wide::ZERO, Wide::add)
}

pub(super) fn sum_as_double(values: impl Iterator<Item = f64>) -> f64 {
    sum(values).to_double()
}

pub(in crate::design) fn sum_argument(values: impl Iterator<Item = f64>) -> f64 {
    let total = sum(values);
    let narrowed = total.to_double();
    if total.exceeds_double() {
        f64::INFINITY.copysign(narrowed)
    } else {
        narrowed
    }
}

pub(super) fn column_mean(values: impl Iterator<Item = f64>, count: usize) -> f64 {
    sum(values).divide(count).to_double()
}

pub(super) fn mean(values: &[f64]) -> f64 {
    if values.is_empty() {
        return f64::NAN;
    }
    let n = values.len();
    let total = sum(values.iter().copied());
    let finite = total.to_double().is_finite();
    let initial = if finite {
        total.divide(n)
    } else {
        sum(values.iter().map(|x| x / n as f64))
    };
    let narrowed = initial.to_double();
    if !narrowed.is_finite() {
        return narrowed;
    }
    let corrections = values
        .iter()
        .map(|&x| Wide::from_double(x).subtract(initial));
    let correction = if finite {
        corrections.fold(Wide::ZERO, Wide::add).divide(n)
    } else {
        corrections.map(|x| x.divide(n)).fold(Wide::ZERO, Wide::add)
    };
    initial.add(correction).to_double()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rounds_each_operation_at_binary128_precision() {
        let wide = 2.0_f64.powi(112);
        assert_eq!(sum_as_double([wide, 1.0, -wide].into_iter()), 1.0);
        let wide = 2.0_f64.powi(113);
        assert_eq!(sum_as_double([wide, 1.0, -wide].into_iter()), 0.0);
    }

    #[test]
    fn narrows_subnormals_and_preserves_r_sum_special_values() {
        let smallest = f64::from_bits(1);
        assert_eq!(column_mean([smallest].into_iter(), 2).to_bits(), 0);
        assert_eq!(column_mean([3.0 * smallest].into_iter(), 2).to_bits(), 2);
        let near_limit = [f64::MAX, 2.0_f64.powi(960)];
        assert_eq!(sum_as_double(near_limit.into_iter()), f64::MAX);
        assert_eq!(sum_argument(near_limit.into_iter()), f64::INFINITY);
        let r_na = f64::from_bits(0x7ff0_0000_0000_07a2);
        assert_eq!(
            sum_as_double([r_na, f64::NAN].into_iter()).to_bits(),
            0x7ff8_0000_0000_07a2
        );
    }
}
