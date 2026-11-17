use super::*;
use crate::numeric::checked::{self, CheckedExt};

pub(super) fn lowess_sorted(
    x: &[f64],
    y: &[f64],
    span_fraction: f64,
    robust_iterations: usize,
    delta: f64,
) -> Result<Vec<f64>, DeseqError> {
    let n = x.len();
    let span = ((span_fraction * n as f64 + 1e-7) as usize).clamp(2, n);
    let mut robustness_weights = vec![1.0; n];
    let mut fitted = vec![0.0; n];

    for iteration in 0..=robust_iterations {
        fitted = lowess_delta_pass(x, y, span, &robustness_weights, delta)?;
        if iteration == robust_iterations {
            break;
        }
        let residuals = y
            .iter()
            .copied()
            .zip(fitted.iter().copied())
            .enumerate()
            .map(|(idx, (observed, fitted))| {
                checked::try_sub(observed, fitted)
                    .nonfinite_at("lowess residual", Some(idx))
                    .map(f64::abs)
            })
            .collect::<Result<Vec<_>, DeseqError>>()?;
        let mean_absolute_residual =
            checked::mean(&residuals).nonfinite_at("lowess mean absolute residual", None)?;
        let scale = six_mad(&residuals)?;
        if scale
            < checked::try_mul(1e-7, mean_absolute_residual)
                .nonfinite_at("lowess residual tolerance", Some(0))?
        {
            break;
        }
        for (weight, residual) in robustness_weights.iter_mut().zip(residuals.iter().copied()) {
            let absolute = residual.abs();
            if absolute
                <= checked::try_mul(0.001, scale)
                    .nonfinite_at("lowess robustness low cutoff", Some(0))?
            {
                *weight = 1.0;
            } else if absolute
                >= checked::try_mul(0.999, scale)
                    .nonfinite_at("lowess robustness high cutoff", Some(0))?
            {
                *weight = 0.0;
            } else {
                let ratio = checked::try_div_finite_result(absolute, scale)
                    .nonfinite_at("lowess robustness ratio", Some(0))?;
                let ratio_sq = checked::try_mul(ratio, ratio)
                    .nonfinite_at("lowess robustness ratio square", Some(0))?;
                let one_minus_ratio_sq = checked::try_sub(1.0, ratio_sq)
                    .nonfinite_at("lowess robustness weight base", Some(0))?;
                *weight = checked::try_mul(one_minus_ratio_sq, one_minus_ratio_sq)
                    .nonfinite_at("lowess robustness weight", Some(0))?;
            }
        }
    }
    Ok(fitted)
}

fn lowess_delta_pass(
    x: &[f64],
    y: &[f64],
    span: usize,
    robustness_weights: &[f64],
    delta: f64,
) -> Result<Vec<f64>, DeseqError> {
    let n = x.len();
    let mut fitted = vec![0.0; n];
    let mut n_left = 0_usize;
    let mut n_right = span - 1;
    let mut last_estimated: Option<usize> = None;
    let mut idx = 0_usize;
    loop {
        while n_right < n - 1 {
            let left_radius = checked::try_sub(x[idx], x[n_left])
                .nonfinite_at("lowess left radius", Some(idx))?;
            let next_right_radius = checked::try_sub(x[n_right + 1], x[idx])
                .nonfinite_at("lowess right radius", Some(idx))?;
            if left_radius > next_right_radius {
                n_left += 1;
                n_right += 1;
            } else {
                break;
            }
        }

        let fitted_value =
            lowest_at(x, y, idx, n_left, n_right, robustness_weights)?.unwrap_or(y[idx]);
        fitted[idx] = fitted_value;

        if let Some(last) = last_estimated
            && last < idx - 1
        {
            let denominator = checked::try_sub(x[idx], x[last])
                .nonfinite_at("lowess interpolation denominator", Some(idx))?;
            for point in last + 1..idx {
                fitted[point] = if denominator.abs() <= f64::EPSILON {
                    fitted[last]
                } else {
                    let fraction = checked::try_div_finite_result(
                        checked::try_sub(x[point], x[last])
                            .nonfinite_at("lowess interpolation x delta", Some(point))?,
                        denominator,
                    )
                    .nonfinite_at("lowess interpolation fraction", Some(point))?;
                    checked::try_add(
                        checked::try_mul(fraction, fitted[idx])
                            .nonfinite_at("lowess interpolation upper term", Some(point))?,
                        checked::try_mul(
                            checked::try_sub(1.0, fraction)
                                .nonfinite_at("lowess interpolation complement", Some(point))?,
                            fitted[last],
                        )
                        .nonfinite_at("lowess interpolation lower term", Some(point))?,
                    )
                    .nonfinite_at("lowess interpolation fitted value", Some(point))?
                };
            }
        }
        let mut last = idx;
        let cutoff = checked::try_add(x[last], delta.max(0.0))
            .nonfinite_at("lowess delta cutoff", Some(last))?;
        let mut next = last + 1;
        while next < n {
            if x[next] > cutoff {
                break;
            }
            if x[next] == x[last] {
                fitted[next] = fitted[last];
                last = next;
            }
            next += 1;
        }
        last_estimated = Some(last);
        if last >= n - 1 {
            break;
        }
        idx = (next.saturating_sub(1)).max(last + 1);
    }
    Ok(fitted)
}

fn lowest_at(
    x: &[f64],
    y: &[f64],
    idx: usize,
    n_left: usize,
    n_right: usize,
    robustness_weights: &[f64],
) -> Result<Option<f64>, DeseqError> {
    let target_x = x[idx];
    let range = checked::try_sub(x[x.len() - 1], x[0]).nonfinite_at("lowess x range", Some(idx))?;
    let left_bandwidth =
        checked::try_sub(target_x, x[n_left]).nonfinite_at("lowess left bandwidth", Some(idx))?;
    let right_bandwidth =
        checked::try_sub(x[n_right], target_x).nonfinite_at("lowess right bandwidth", Some(idx))?;
    let bandwidth = left_bandwidth.max(right_bandwidth);
    let high =
        checked::try_mul(0.999, bandwidth).nonfinite_at("lowess high bandwidth", Some(idx))?;
    let low = checked::try_mul(0.001, bandwidth).nonfinite_at("lowess low bandwidth", Some(idx))?;
    let mut weights = Vec::new();
    let mut sum_weights = 0.0;
    let mut point = n_left;
    while point < x.len() {
        let distance = checked::try_sub(x[point], target_x)
            .nonfinite_at("lowess distance", Some(point))?
            .abs();
        if distance <= high {
            let proximity = if distance <= low {
                1.0
            } else {
                let ratio = checked::try_div_finite_result(distance, bandwidth)
                    .nonfinite_at("lowess distance ratio", Some(point))?;
                let ratio_sq = checked::try_mul(ratio, ratio)
                    .nonfinite_at("lowess ratio square", Some(point))?;
                let ratio_cubed = checked::try_mul(ratio_sq, ratio)
                    .nonfinite_at("lowess ratio cube", Some(point))?;
                let one_minus_ratio_cubed = checked::try_sub(1.0, ratio_cubed)
                    .nonfinite_at("lowess tricube base", Some(point))?;
                checked::try_mul(
                    checked::try_mul(one_minus_ratio_cubed, one_minus_ratio_cubed)
                        .nonfinite_at("lowess tricube square", Some(point))?,
                    one_minus_ratio_cubed,
                )
                .nonfinite_at("lowess tricube proximity", Some(point))?
            };
            let weight = checked::try_mul(proximity, robustness_weights[point])
                .nonfinite_at("lowess tricube robustness weight", Some(point))?;
            weights.push((point, weight));
            sum_weights = checked::try_add(sum_weights, weight)
                .nonfinite_at("lowess weight sum", Some(point))?;
        } else if x[point] > target_x {
            break;
        }
        point += 1;
    }
    if sum_weights <= 0.0 {
        return Ok(None);
    }

    for (_, weight) in &mut weights {
        *weight = checked::try_div_finite_result(*weight, sum_weights)
            .nonfinite_at("lowess normalized weight", Some(0))?;
    }

    if bandwidth > 0.0 {
        let center = weights
            .iter()
            .map(|(point, weight)| {
                checked::try_mul(*weight, x[*point]).nonfinite_at("lowess weighted x", Some(*point))
            })
            .collect::<Result<Vec<_>, DeseqError>>()?;
        let center = checked::try_sum(center).nonfinite("lowess center")?;
        let mut variance = 0.0;
        for (point, weight) in &weights {
            let residual = checked::try_sub(x[*point], center)
                .nonfinite_at("lowess centered x", Some(*point))?;
            let residual_sq = checked::try_mul(residual, residual)
                .nonfinite_at("lowess centered x square", Some(*point))?;
            variance = checked::try_add(
                variance,
                checked::try_mul(*weight, residual_sq)
                    .nonfinite_at("lowess weighted variance term", Some(*point))?,
            )
            .nonfinite_at("lowess variance", Some(*point))?;
        }
        if variance.sqrt()
            > checked::try_mul(0.001, range)
                .nonfinite_at("lowess slope variance threshold", Some(idx))?
        {
            let slope_factor = checked::try_div_finite_result(
                checked::try_sub(target_x, center)
                    .nonfinite_at("lowess slope target delta", Some(idx))?,
                variance,
            )
            .nonfinite_at("lowess slope factor", Some(idx))?;
            for (point, weight) in &mut weights {
                let centered = checked::try_sub(x[*point], center)
                    .nonfinite_at("lowess slope centered x", Some(*point))?;
                let adjustment = checked::try_add(
                    checked::try_mul(slope_factor, centered)
                        .nonfinite_at("lowess slope weight adjustment", Some(*point))?,
                    1.0,
                )
                .nonfinite_at("lowess slope weight factor", Some(*point))?;
                *weight = checked::try_mul(*weight, adjustment)
                    .nonfinite_at("lowess slope weight", Some(*point))?;
            }
        }
    }

    let mut fitted = 0.0;
    for (point, weight) in weights {
        fitted = checked::try_add(
            fitted,
            checked::try_mul(weight, y[point])
                .nonfinite_at("lowess fitted weighted y", Some(point))?,
        )
        .nonfinite_at("lowess fitted value", Some(point))?;
    }
    Ok(Some(fitted))
}
