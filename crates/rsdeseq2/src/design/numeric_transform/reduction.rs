// Arithmetic schedules adapted from OpenBLAS 0.3.26.
// See OPENBLAS-LICENSE in this directory for its BSD license.
//! Explicit, portable reduction order for the small polynomial QR problem.
//!
//! These operations specify rounding rather than delegating it to whichever
//! BLAS the calling frontend happens to load. Rust and R use the same schedule.
pub(super) fn dot(x: &[f64], y: &[f64]) -> f64 {
    #[cfg(all(target_os = "linux", target_arch = "x86_64"))]
    {
        dot_skylake(x, y)
    }
    #[cfg(not(all(target_os = "linux", target_arch = "x86_64")))]
    dot_binary64(x, y)
}
#[cfg(not(all(target_os = "linux", target_arch = "x86_64")))]
fn dot_binary64(x: &[f64], y: &[f64]) -> f64 {
    let block_width = if cfg!(all(target_os = "linux", target_arch = "x86_64")) {
        16
    } else {
        32
    };
    let blocks = x.len() / block_width;
    let mut lanes = [0.0; 16];
    for block in 0..blocks {
        let start = block * block_width;
        for half in (0..block_width).step_by(16) {
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
    for index in blocks * block_width..x.len() {
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
fn dot_skylake(x: &[f64], y: &[f64]) -> f64 {
    let mut wide = [0.0; 32];
    let n32 = x.len() / 32 * 32;
    for start in (0..n32).step_by(32) {
        for i in 0..32 {
            wide[i] = x[start + i].mul_add(y[start + i], wide[i]);
        }
    }
    let mut lanes = [0.0; 16];
    for block in 0..4 {
        for i in 0..4 {
            lanes[block * 4 + i] = wide[block * 8 + i] + wide[block * 8 + i + 4];
        }
    }
    let n16 = x.len() / 16 * 16;
    for start in (n32..n16).step_by(16) {
        for i in 0..16 {
            lanes[i] = x[start + i].mul_add(y[start + i], lanes[i]);
        }
    }
    let mut reduced = [0.0; 4];
    for i in 0..4 {
        reduced[i] = ((lanes[i] + lanes[i + 4]) + lanes[i + 8]) + lanes[i + 12];
    }
    let mut total = (reduced[0] + reduced[2]) + (reduced[1] + reduced[3]);
    for i in n16..x.len() {
        total = x[i].mul_add(y[i], total);
    }
    total
}
