//! Vector reductions in a fixed evaluation order.
//!
//! R's `qr()` reduces columns with the dot product, Euclidean norm and scaled
//! vector update of the BLAS it is linked against. These functions perform
//! the same operations in the order and with the rounding recorded for the
//! reference host (`docs/provenance/specs/reductions.md`), so that results
//! built on them, such as the `poly()` and `contr.poly()` bases, are
//! identical on every platform.

/// Number of interleaved accumulators in the blocked part of [`dot`].
const LANES: usize = 16;
/// Length of one block of the blocked part of [`dot`].
const BLOCK: usize = 2 * LANES;

/// Dot product of two equal-length slices.
///
/// The longest prefix whose length is a multiple of 32 is accumulated into
/// sixteen interleaved fused multiply-add accumulators (position `p` goes to
/// accumulator `p mod 16`). The accumulators are combined by a fixed pairwise
/// tree, odd-numbered and even-numbered lanes separately and then added. The
/// remaining tail is added in order with fused multiply-adds. For fewer than
/// 32 elements this is a left-to-right fused accumulation from zero.
///
/// Only the common prefix is used if the lengths differ; callers pass slices
/// of equal length.
pub fn dot(x: &[f64], y: &[f64]) -> f64 {
    let n = x.len().min(y.len());
    let blocked = BLOCK * (n / BLOCK);
    let mut lanes = [0.0_f64; LANES];
    for (p, (&xp, &yp)) in x[..blocked].iter().zip(&y[..blocked]).enumerate() {
        let lane = &mut lanes[p % LANES];
        *lane = xp.mul_add(yp, *lane);
    }
    let a = lanes;
    let odd = ((a[0] + a[2]) + (a[4] + a[6])) + ((a[8] + a[10]) + (a[12] + a[14]));
    let even = ((a[1] + a[3]) + (a[5] + a[7])) + ((a[9] + a[11]) + (a[13] + a[15]));
    let mut sum = odd + even;
    for (&xp, &yp) in x[blocked..n].iter().zip(&y[blocked..n]) {
        sum = xp.mul_add(yp, sum);
    }
    sum
}

/// Euclidean norm by a one-pass scaled sum of squares.
///
/// The scan keeps `scale`, the largest magnitude seen so far, and `ssq`, the
/// sum of squares in units of `scale` (Blue, 1978; Higham, 2002), so entries
/// near the underflow or overflow thresholds do not lose the result. Zero
/// entries are skipped and the result depends on the order of the entries.
///
/// Rounding of the updates (chosen because it reproduces every row of the
/// recorded observations, `dnrm2.tsv` and `dnrm2_random.tsv`, bit for bit,
/// while the other candidate groupings and fusings do not):
///
/// - a new largest magnitude `a`: `r = scale / a`, then
///   `ssq = 1 + ssq * (r * r)` with every operation rounded separately;
/// - any other entry `a`: `r = a / scale`, then `ssq = fma(r, r, ssq)`;
/// - the result is `scale * sqrt(ssq)`.
pub fn norm(x: &[f64]) -> f64 {
    let mut scale = 0.0_f64;
    let mut ssq = 1.0_f64;
    for &value in x {
        if value == 0.0 {
            continue;
        }
        let magnitude = value.abs();
        if magnitude > scale {
            let ratio = scale / magnitude;
            ssq = 1.0 + ssq * (ratio * ratio);
            scale = magnitude;
        } else {
            let ratio = magnitude / scale;
            ssq = ratio.mul_add(ratio, ssq);
        }
    }
    scale * ssq.sqrt()
}

/// Scaled vector update `y <- alpha * x + y` with one rounding per element.
///
/// Every element becomes `fma(alpha, x[i], y[i])`. Only the common prefix is
/// updated if the lengths differ; callers pass slices of equal length.
pub fn axpy(alpha: f64, x: &[f64], y: &mut [f64]) {
    for (yi, &xi) in y.iter_mut().zip(x) {
        *yi = alpha.mul_add(xi, *yi);
    }
}

#[cfg(test)]
mod tests {
    use super::{axpy, dot, norm};

    const DDOT_RANDOM: &str = include_str!(
        "../../../../docs/provenance/observations/reductions_r462_macos/ddot_random.tsv"
    );
    const DDOT_FMA: &str =
        include_str!("../../../../docs/provenance/observations/reductions_r462_macos/ddot_fma.tsv");
    const DDOT_FMA_PAIRS: &str = include_str!(
        "../../../../docs/provenance/observations/reductions_r462_macos/ddot_fma_pairs.tsv"
    );
    const DAXPY_FMA: &str = include_str!(
        "../../../../docs/provenance/observations/reductions_r462_macos/daxpy_fma.tsv"
    );
    const DNRM2: &str =
        include_str!("../../../../docs/provenance/observations/reductions_r462_macos/dnrm2.tsv");
    const DNRM2_RANDOM: &str = include_str!(
        "../../../../docs/provenance/observations/reductions_r462_macos/dnrm2_random.tsv"
    );

    // Product a * b = 1 - 2^-60 is inexact in binary64 (it rounds to 1).
    const A: f64 = 1.0 + 1.0 / 1_073_741_824.0;
    const B: f64 = 1.0 - 1.0 / 1_073_741_824.0;

    fn bits(text: &str) -> f64 {
        f64::from_bits(u64::from_str_radix(text, 16).unwrap())
    }

    fn vector(text: &str) -> Vec<f64> {
        text.split(',').map(bits).collect()
    }

    fn rows(table: &str) -> impl Iterator<Item = Vec<&str>> {
        table
            .lines()
            .skip(1)
            .filter(|line| !line.is_empty())
            .map(|line| line.split('\t').collect())
    }

    #[test]
    fn dot_reproduces_recorded_random_products() {
        let mut count = 0;
        for row in rows(DDOT_RANDOM) {
            let (x, y) = (vector(row[0]), vector(row[1]));
            let expected = bits(row[2]);
            assert_eq!(x.len(), y.len());
            let got = dot(&x, &y);
            assert_eq!(got.to_bits(), expected.to_bits(), "n = {}", x.len());
            count += 1;
        }
        assert_eq!(count, 300);
    }

    #[test]
    fn dot_fuses_every_product_into_its_accumulator() {
        let mut count = 0;
        for row in rows(DDOT_FMA) {
            let n: usize = row[0].parse().unwrap();
            let position: usize = row[1].parse().unwrap();
            let (mut x, mut y) = (vec![0.0; n], vec![0.0; n]);
            x[position - 1] = A;
            y[position - 1] = B;
            let other = if position == 1 { n } else { 1 };
            if n > 1 {
                x[other - 1] = -1.0;
                y[other - 1] = 1.0;
            }
            assert_eq!(dot(&x, &y).to_bits(), bits(row[2]).to_bits(), "{row:?}");
            count += 1;
        }
        assert!(count > 0);
        for row in rows(DDOT_FMA_PAIRS) {
            let n: usize = row[0].parse().unwrap();
            let previous: usize = row[1].parse().unwrap();
            let position: usize = row[2].parse().unwrap();
            let (mut x, mut y) = (vec![0.0; n], vec![0.0; n]);
            x[previous - 1] = -1.0;
            y[previous - 1] = 1.0;
            x[position - 1] = A;
            y[position - 1] = B;
            assert_eq!(dot(&x, &y).to_bits(), bits(row[3]).to_bits(), "{row:?}");
        }
    }

    #[test]
    fn axpy_rounds_each_element_once() {
        let mut count = 0;
        for row in rows(DAXPY_FMA) {
            let n: usize = row[0].parse().unwrap();
            let expected = vector(row[1]);
            let x = vec![A; n];
            let mut y = vec![-1.0; n];
            axpy(B, &x, &mut y);
            assert_eq!(expected.len(), n);
            for (got, want) in y.iter().zip(&expected) {
                assert_eq!(got.to_bits(), want.to_bits(), "n = {n}");
            }
            count += 1;
        }
        assert!(count > 0);
    }

    #[test]
    fn norm_reproduces_recorded_cases() {
        let mut count = 0;
        for row in rows(DNRM2) {
            let x = vector(row[2]);
            assert_eq!(x.len(), row[1].parse::<usize>().unwrap());
            assert_eq!(norm(&x).to_bits(), bits(row[3]).to_bits(), "{}", row[0]);
            count += 1;
        }
        assert!(count > 0);
        let mut count = 0;
        for row in rows(DNRM2_RANDOM) {
            let x = vector(row[0]);
            assert_eq!(
                norm(&x).to_bits(),
                bits(row[1]).to_bits(),
                "n = {}",
                x.len()
            );
            count += 1;
        }
        assert_eq!(count, 300);
    }

    #[test]
    fn empty_and_zero_inputs() {
        assert_eq!(dot(&[], &[]), 0.0);
        assert_eq!(norm(&[]), 0.0);
        assert_eq!(norm(&[0.0, 0.0]), 0.0);
        assert_eq!(norm(&[0.0, -3.0, 0.0, 4.0]), 5.0);
    }
}
