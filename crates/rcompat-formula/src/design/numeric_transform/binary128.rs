//! R scalar reductions using IEEE binary128 on Linux arm64.

use rustc_apfloat::ieee::{Double, Quad};
use rustc_apfloat::{Float, FloatConvert};

fn wide(x: f64) -> Quad {
    let bits = if x.is_nan() {
        x.to_bits() | (1 << 51)
    } else {
        x.to_bits()
    };
    Double::from_bits(u128::from(bits))
        .convert(&mut false)
        .value
}

fn narrow(x: Quad) -> f64 {
    let value: Double = x.convert(&mut false).value;
    f64::from_bits(value.to_bits() as u64)
}

fn sum(values: impl Iterator<Item = f64>) -> Quad {
    values.map(wide).fold(Quad::ZERO, add)
}

pub(super) fn sum_as_double(values: impl Iterator<Item = f64>) -> f64 {
    narrow(sum(values))
}

pub(in crate::design) fn sum_argument(values: impl Iterator<Item = f64>) -> f64 {
    let total = sum(values);
    let limit = wide(f64::MAX);
    if total > limit {
        f64::INFINITY
    } else if total < -limit {
        f64::NEG_INFINITY
    } else {
        narrow(total)
    }
}

pub(super) fn column_mean(values: impl Iterator<Item = f64>, count: usize) -> f64 {
    narrow((sum(values) / Quad::from_u128(count as u128).value).value)
}

pub(super) fn mean(values: &[f64]) -> f64 {
    if values.is_empty() {
        return f64::NAN;
    }
    let n = Quad::from_u128(values.len() as u128).value;
    let total = sum(values.iter().copied());
    let finite = narrow(total).is_finite();
    let initial = if finite {
        (total / n).value
    } else {
        sum(values.iter().map(|x| x / values.len() as f64))
    };
    if !narrow(initial).is_finite() {
        return narrow(initial);
    }
    let corrections = values.iter().map(|&x| (wide(x) - initial).value);
    let correction = if finite {
        (corrections.fold(Quad::ZERO, add) / n).value
    } else {
        corrections.map(|x| (x / n).value).fold(Quad::ZERO, add)
    };
    narrow((initial + correction).value)
}

fn add(a: Quad, b: Quad) -> Quad {
    // GCC's binary128 addition preserves the first NaN operand. Opposite
    // infinities create a positive quiet NaN with the payload bits set.
    if a.is_nan() {
        return a;
    }
    if b.is_nan() {
        return b;
    }
    let sum = (a + b).value;
    if sum.is_nan() {
        wide(f64::from_bits(0x7fff_ffff_ffff_ffff))
    } else {
        sum
    }
}
