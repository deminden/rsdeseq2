//! Centered polynomial QR with relative-norm column deferral.
//! This operation does not change the rank policy used by model fitting.
#[cfg(all(target_os = "linux", target_arch = "x86_64"))]
use super::reduction::axpy;
use super::reduction::{dot, norm};

// Caller supplies a finite, nonempty column-major problem with width <= rows.
pub(in crate::design) fn basis(mut columns: Vec<Vec<f64>>) -> (Vec<Vec<f64>>, usize) {
    let mut tails = columns.iter().map(|c| norm(c)).collect::<Vec<_>>();
    let mut original = tails.clone();
    let n = columns[0].len();
    let width = columns.len();
    let mut reflectors = vec![0.0; width];
    let mut active = width;
    for step in 0..width {
        while step < active && tails[step] < original[step] * 1e-7 {
            columns[step..].rotate_left(1);
            tails[step..].rotate_left(1);
            original[step..].rotate_left(1);
            active -= 1;
        }
        if step == n - 1 {
            break;
        }
        let magnitude = norm(&columns[step][step..]);
        if magnitude == 0.0 {
            continue;
        }
        let diagonal = columns[step][step];
        let signed = if diagonal == 0.0 {
            magnitude
        } else {
            magnitude.copysign(diagonal)
        };
        let reciprocal = 1.0 / signed;
        for value in &mut columns[step][step..] {
            *value *= reciprocal;
        }
        columns[step][step] += 1.0;
        for other in step + 1..width {
            let inner = dot(&columns[step][step..], &columns[other][step..]);
            let weight = -inner / columns[step][step];
            if weight != 0.0 {
                let (previous, remaining) = columns.split_at_mut(other);
                #[cfg(all(target_os = "linux", target_arch = "x86_64"))]
                axpy(weight, &previous[step][step..], &mut remaining[0][step..]);
                #[cfg(not(all(target_os = "linux", target_arch = "x86_64")))]
                for (value, &reflector) in
                    remaining[0][step..].iter_mut().zip(&previous[step][step..])
                {
                    *value = weight.mul_add(reflector, *value);
                }
            }
            if tails[other] != 0.0 {
                let fraction = columns[other][step].abs() / tails[other];
                let remaining = if cfg!(all(target_os = "linux", target_arch = "x86_64")) {
                    (1.0 - fraction * fraction).max(0.0)
                } else {
                    (-fraction).mul_add(fraction, 1.0).max(0.0)
                };
                tails[other] = if remaining >= 1e-6 {
                    tails[other] * remaining.sqrt()
                } else {
                    norm(&columns[other][step + 1..])
                };
            }
        }
        reflectors[step] = columns[step][step];
        columns[step][step] = -signed;
    }
    let reflections = active.min(n - 1);
    let mut output = Vec::with_capacity(width);
    for column in 0..width {
        let mut q = vec![0.0; n];
        q[column] = columns[column][column];
        for step in (0..reflections).rev() {
            if reflectors[step] == 0.0 {
                continue;
            }
            let diagonal = columns[step][step];
            columns[step][step] = reflectors[step];
            let inner = dot(&columns[step][step..], &q[step..]);
            let weight = -inner / reflectors[step];
            if weight != 0.0 {
                #[cfg(all(target_os = "linux", target_arch = "x86_64"))]
                axpy(weight, &columns[step][step..], &mut q[step..]);
                #[cfg(not(all(target_os = "linux", target_arch = "x86_64")))]
                for r in step..n {
                    q[r] = weight.mul_add(columns[step][r], q[r]);
                }
            }
            columns[step][step] = diagonal;
        }
        output.push(q);
    }
    (output, active)
}
