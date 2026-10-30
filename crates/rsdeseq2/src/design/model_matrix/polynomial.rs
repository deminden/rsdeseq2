//! Polynomial contrast construction with relative-norm column deferral.
//!
//! This is separate from fitting rank checks: changing the order/number of
//! applied reflections changes the basis even when its span is unchanged.
use super::polynomial_reduction::{dot, norm};
pub(super) fn contrasts(n: usize) -> Vec<f64> {
    let center = (n as f64 + 1.0) / 2.0;
    let mut columns = (0..n)
        .map(|degree| {
            (0..n)
                .map(|r| {
                    let x = r as f64 + 1.0 - center;
                    match degree {
                        0 => 1.0,
                        1 => x,
                        2 => x * x,
                        _ => x.powf(degree as f64),
                    }
                })
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let mut tails = columns.iter().map(|c| norm(c)).collect::<Vec<_>>();
    let mut original = tails.clone();
    let mut reflectors = vec![0.0; n];
    let mut active = n;
    for step in 0..n {
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
        for other in step + 1..n {
            let inner = dot(&columns[step][step..], &columns[other][step..]);
            let weight = -inner / columns[step][step];
            if weight != 0.0 {
                let (previous, remaining) = columns.split_at_mut(other);
                for (value, &reflector) in
                    remaining[0][step..].iter_mut().zip(&previous[step][step..])
                {
                    *value = weight.mul_add(reflector, *value);
                }
            }
            if tails[other] != 0.0 {
                let fraction = columns[other][step].abs() / tails[other];
                let remaining = (-fraction).mul_add(fraction, 1.0).max(0.0);
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
    let mut output = vec![0.0; n * (n - 1)];
    for column in 1..n {
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
                for r in step..n {
                    q[r] = weight.mul_add(columns[step][r], q[r]);
                }
            }
            columns[step][step] = diagonal;
        }
        let length = q.iter().map(|v| v * v).sum::<f64>().sqrt();
        for row in 0..n {
            output[row * (n - 1) + column - 1] = q[row] / length;
        }
    }
    output
}
