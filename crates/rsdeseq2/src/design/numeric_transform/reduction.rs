//! Explicit, portable reduction order for the small polynomial QR problem.
//!
//! These operations specify rounding rather than delegating it to whichever
//! BLAS the calling frontend happens to load. Rust and R use the same schedule.
pub(super) fn dot(x: &[f64], y: &[f64]) -> f64 {
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
