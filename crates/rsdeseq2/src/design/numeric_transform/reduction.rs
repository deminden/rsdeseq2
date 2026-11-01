// Arithmetic schedules adapted from OpenBLAS 0.3.26.
// See OPENBLAS-LICENSE in this directory for its BSD license.
//! Explicit, portable reduction order for the small polynomial QR problem.
//!
//! These operations specify rounding rather than delegating it to whichever
//! BLAS the calling frontend happens to load. Rust and R use the same schedule.
pub(super) fn dot(x: &[f64], y: &[f64]) -> f64 {
    #[cfg(all(target_os = "linux", target_arch = "x86_64"))]
    {
        dot_reference_haswell(x, y)
    }
    #[cfg(not(all(target_os = "linux", target_arch = "x86_64")))]
    dot_binary64(x, y)
}
#[cfg(not(all(target_os = "linux", target_arch = "x86_64")))]
fn dot_binary64(x: &[f64], y: &[f64]) -> f64 {
    let blocks = x.len() / 32;
    let mut lanes = [0.0; 16];
    for block in 0..blocks {
        let start = block * 32;
        for half in [0, 16] {
            for lane in 0..16 {
                lanes[lane] = x[start + half + lane].mul_add(y[start + half + lane], lanes[lane]);
            }
        }
    }
    let mut sum = if blocks == 0 {
        0.0
    } else {
        let left_even = (lanes[0] + lanes[2]) + (lanes[4] + lanes[6]);
        let left_odd = (lanes[1] + lanes[3]) + (lanes[5] + lanes[7]);
        let right_even = (lanes[8] + lanes[10]) + (lanes[12] + lanes[14]);
        let right_odd = (lanes[9] + lanes[11]) + (lanes[13] + lanes[15]);
        (left_even + right_even) + (left_odd + right_odd)
    };
    for index in blocks * 32..x.len() {
        sum = x[index].mul_add(y[index], sum);
    }
    sum
}
pub(super) fn norm(x: &[f64]) -> f64 {
    #[cfg(all(target_os = "linux", target_arch = "x86_64"))]
    {
        super::extended::norm(x)
    }
    #[cfg(not(all(target_os = "linux", target_arch = "x86_64")))]
    norm_binary64(x)
}
#[cfg(not(all(target_os = "linux", target_arch = "x86_64")))]
fn norm_binary64(x: &[f64]) -> f64 {
    let mut scale = 0.0;
    let mut squares = 1.0;
    for value in x.iter().map(|v| v.abs()).filter(|&v| v != 0.0) {
        if value > scale {
            let ratio = scale / value;
            squares = 1.0 + squares * (ratio * ratio);
            scale = value;
        } else {
            let ratio = value / scale;
            squares = ratio.mul_add(ratio, squares);
        }
    }
    scale * squares.sqrt()
}

#[cfg(all(target_os = "linux", target_arch = "x86_64"))]
fn dot_reference_haswell(x: &[f64], y: &[f64]) -> f64 {
    let mut lanes = [0.0; 16];
    let end = x.len() / 16 * 16;
    for start in (0..end).step_by(16) {
        for i in 0..16 {
            lanes[i] = x[start + i].mul_add(y[start + i], lanes[i]);
        }
    }
    let mut pair = [0.0; 8];
    for group in 0..4 {
        for i in 0..2 {
            pair[group * 2 + i] = lanes[group * 4 + i] + lanes[group * 4 + i + 2];
        }
    }
    let even = (pair[0] + pair[2]) + (pair[4] + pair[6]);
    let odd = (pair[1] + pair[3]) + (pair[5] + pair[7]);
    let mut total = even + odd;
    for i in end..x.len() {
        total += x[i] * y[i];
    }
    total
}
#[cfg(all(target_os = "linux", target_arch = "x86_64"))]
pub(super) fn axpy(weight: f64, x: &[f64], y: &mut [f64]) {
    let end = x.len() / 16 * 16;
    for (i, (v, &x)) in y.iter_mut().zip(x).enumerate() {
        *v = if cfg!(all(target_os = "linux", target_arch = "x86_64")) && i >= end {
            *v + weight * x
        } else {
            weight.mul_add(x, *v)
        };
    }
}
