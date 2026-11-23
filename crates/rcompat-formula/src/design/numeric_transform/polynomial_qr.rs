//! Centered polynomial QR with relative-norm column deferral.
//! This operation does not change the rank policy used by model fitting.
use super::reduction::axpy;
use super::reduction::{dot, norm};
use crate::errors::Error;
use crate::matrix::try_alloc_matrix;

pub(in crate::design) fn try_columns(
    rows: usize,
    width: usize,
    value: impl Fn(usize, usize) -> f64,
) -> Result<Vec<Vec<f64>>, Error> {
    let mut columns = try_alloc_matrix(width, 1, "polynomial QR columns")?;
    for column in 0..width {
        let mut values = try_alloc_matrix(rows, 1, "polynomial QR column values")?;
        values.extend((0..rows).map(|row| value(row, column)));
        columns.push(values);
    }
    Ok(columns)
}

// Caller supplies a finite, nonempty column-major problem with width <= rows.
pub(in crate::design) fn visit_basis(
    mut columns: Vec<Vec<f64>>,
    mut visit: impl FnMut(usize, &[f64]),
) -> Result<usize, Error> {
    let n = columns[0].len();
    let width = columns.len();
    let mut tails = try_alloc_matrix(width, 1, "polynomial QR column norms")?;
    tails.extend(columns.iter().map(|c| norm(c)));
    let mut original = try_alloc_matrix(width, 1, "polynomial QR original norms")?;
    original.extend_from_slice(&tails);
    let mut reflectors = try_alloc_matrix(width, 1, "polynomial QR reflectors")?;
    reflectors.resize(width, 0.0);
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
                axpy(weight, &previous[step][step..], &mut remaining[0][step..]);
            }
            if tails[other] != 0.0 {
                let fraction = columns[other][step].abs() / tails[other];
                let remaining = if crate::numeric::profile::r_linux_x86() {
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
    let mut q = try_alloc_matrix(n, 1, "polynomial QR basis column")?;
    q.resize(n, 0.0);
    for column in 0..width {
        q.fill(0.0);
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
                axpy(weight, &columns[step][step..], &mut q[step..]);
            }
            columns[step][step] = diagonal;
        }
        visit(column, &q);
    }
    Ok(active)
}
