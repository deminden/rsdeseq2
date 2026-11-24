//! Fixed binary64 arithmetic for polynomial QR.
//!
//! These reductions express the rounding observed in the R reference
//! profiles. The block size determines where the sequential tail starts.

use crate::numeric::profile::{r_linux_arm64, r_linux_x86};

pub(super) fn dot(x: &[f64], y: &[f64]) -> f64 {
    if r_linux_x86() {
        lane_dot::<16>(x, y, false)
    } else if r_linux_arm64() {
        two_lane_dot(x, y)
    } else {
        lane_dot::<32>(x, y, true)
    }
}

fn lane_dot<const BLOCK: usize>(x: &[f64], y: &[f64], fused_tail: bool) -> f64 {
    let (x_blocks, x_tail) = x.as_chunks::<BLOCK>();
    let (y_blocks, y_tail) = y.as_chunks::<BLOCK>();
    let mut lane_totals = [0.0; 16];
    for (x_block, y_block) in x_blocks.iter().zip(y_blocks) {
        for (left, right) in x_block
            .as_chunks::<16>()
            .0
            .iter()
            .zip(y_block.as_chunks::<16>().0)
        {
            for (total, (&a, &b)) in lane_totals.iter_mut().zip(left.iter().zip(right)) {
                *total = a.mul_add(b, *total);
            }
        }
    }
    let groups = lane_totals.as_chunks::<4>().0;
    let parity = |i| {
        let partial = |four: &[f64; 4]| four[i] + four[i + 2];
        (partial(&groups[0]) + partial(&groups[1])) + (partial(&groups[2]) + partial(&groups[3]))
    };
    let sum = parity(0) + parity(1);
    x_tail.iter().zip(y_tail).fold(sum, |total, (&a, &b)| {
        if fused_tail {
            a.mul_add(b, total)
        } else {
            total + a * b
        }
    })
}

fn two_lane_dot(x: &[f64], y: &[f64]) -> f64 {
    let end = x.len() - x.len() % 4;
    let pairs = x[..end]
        .as_chunks::<2>()
        .0
        .iter()
        .zip(y[..end].as_chunks::<2>().0);
    let totals = pairs.fold([0.0; 2], |mut totals, (left, right)| {
        for (sum, (&a, &b)) in totals.iter_mut().zip(left.iter().zip(right)) {
            *sum += a * b;
        }
        totals
    });
    x[end..]
        .iter()
        .zip(&y[end..])
        .fold(totals[0] + totals[1], |total, (&a, &b)| a.mul_add(b, total))
}

pub(super) fn norm(x: &[f64]) -> f64 {
    if r_linux_x86() {
        return super::extended::norm(x);
    }
    // Keep the accumulated sum of squares in units of the largest magnitude
    // seen so far. Rescale it when a larger magnitude arrives.
    let (scale, squares) = x
        .iter()
        .copied()
        .fold((0.0_f64, 1.0), |(scale, squares), value| {
            let magnitude = value.abs();
            if magnitude == 0.0 {
                return (scale, squares);
            }
            if magnitude > scale {
                let relative = scale / magnitude;
                (magnitude, 1.0 + squares * (relative * relative))
            } else {
                let relative = magnitude / scale;
                (scale, relative.mul_add(relative, squares))
            }
        });
    scale * squares.sqrt()
}

pub(super) fn axpy(weight: f64, x: &[f64], y: &mut [f64]) {
    y.iter_mut().zip(x).for_each(|(target, source)| {
        *target = weight.mul_add(*source, *target);
    });
}
