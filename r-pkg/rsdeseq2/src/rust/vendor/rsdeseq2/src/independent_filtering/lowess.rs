use super::*;

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
                checked_sub(observed, fitted, idx, "lowess residual").map(f64::abs)
            })
            .collect::<Result<Vec<_>, DeseqError>>()?;
        let mean_absolute_residual =
            checked_scaled_mean(&residuals, "lowess mean absolute residual")?;
        let scale = six_mad(&residuals)?;
        if scale < checked_mul(1e-7, mean_absolute_residual, 0, "lowess residual tolerance")? {
            break;
        }
        for (weight, residual) in robustness_weights.iter_mut().zip(residuals.iter().copied()) {
            let absolute = residual.abs();
            if absolute <= checked_mul(0.001, scale, 0, "lowess robustness low cutoff")? {
                *weight = 1.0;
            } else if absolute >= checked_mul(0.999, scale, 0, "lowess robustness high cutoff")? {
                *weight = 0.0;
            } else {
                let ratio = checked_div(absolute, scale, 0, "lowess robustness ratio")?;
                let ratio_sq = checked_mul(ratio, ratio, 0, "lowess robustness ratio square")?;
                let one_minus_ratio_sq =
                    checked_sub(1.0, ratio_sq, 0, "lowess robustness weight base")?;
                *weight = checked_mul(
                    one_minus_ratio_sq,
                    one_minus_ratio_sq,
                    0,
                    "lowess robustness weight",
                )?;
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
            let left_radius = checked_sub(x[idx], x[n_left], idx, "lowess left radius")?;
            let next_right_radius =
                checked_sub(x[n_right + 1], x[idx], idx, "lowess right radius")?;
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
            let denominator =
                checked_sub(x[idx], x[last], idx, "lowess interpolation denominator")?;
            for point in last + 1..idx {
                fitted[point] = if denominator.abs() <= f64::EPSILON {
                    fitted[last]
                } else {
                    let fraction = checked_div(
                        checked_sub(x[point], x[last], point, "lowess interpolation x delta")?,
                        denominator,
                        point,
                        "lowess interpolation fraction",
                    )?;
                    checked_add(
                        checked_mul(
                            fraction,
                            fitted[idx],
                            point,
                            "lowess interpolation upper term",
                        )?,
                        checked_mul(
                            checked_sub(1.0, fraction, point, "lowess interpolation complement")?,
                            fitted[last],
                            point,
                            "lowess interpolation lower term",
                        )?,
                        point,
                        "lowess interpolation fitted value",
                    )?
                };
            }
        }
        let mut last = idx;
        let cutoff = checked_add(x[last], delta.max(0.0), last, "lowess delta cutoff")?;
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
    let range = checked_sub(x[x.len() - 1], x[0], idx, "lowess x range")?;
    let left_bandwidth = checked_sub(target_x, x[n_left], idx, "lowess left bandwidth")?;
    let right_bandwidth = checked_sub(x[n_right], target_x, idx, "lowess right bandwidth")?;
    let bandwidth = left_bandwidth.max(right_bandwidth);
    let high = checked_mul(0.999, bandwidth, idx, "lowess high bandwidth")?;
    let low = checked_mul(0.001, bandwidth, idx, "lowess low bandwidth")?;
    let mut weights = Vec::new();
    let mut sum_weights = 0.0;
    let mut point = n_left;
    while point < x.len() {
        let distance = checked_sub(x[point], target_x, point, "lowess distance")?.abs();
        if distance <= high {
            let proximity = if distance <= low {
                1.0
            } else {
                let ratio = checked_div(distance, bandwidth, point, "lowess distance ratio")?;
                let ratio_sq = checked_mul(ratio, ratio, point, "lowess ratio square")?;
                let ratio_cubed = checked_mul(ratio_sq, ratio, point, "lowess ratio cube")?;
                let one_minus_ratio_cubed =
                    checked_sub(1.0, ratio_cubed, point, "lowess tricube base")?;
                checked_mul(
                    checked_mul(
                        one_minus_ratio_cubed,
                        one_minus_ratio_cubed,
                        point,
                        "lowess tricube square",
                    )?,
                    one_minus_ratio_cubed,
                    point,
                    "lowess tricube proximity",
                )?
            };
            let weight = checked_mul(
                proximity,
                robustness_weights[point],
                point,
                "lowess tricube robustness weight",
            )?;
            weights.push((point, weight));
            sum_weights = checked_add(sum_weights, weight, point, "lowess weight sum")?;
        } else if x[point] > target_x {
            break;
        }
        point += 1;
    }
    if sum_weights <= 0.0 {
        return Ok(None);
    }

    for (_, weight) in &mut weights {
        *weight = checked_div(*weight, sum_weights, 0, "lowess normalized weight")?;
    }

    if bandwidth > 0.0 {
        let center = weights
            .iter()
            .map(|(point, weight)| checked_mul(*weight, x[*point], *point, "lowess weighted x"))
            .collect::<Result<Vec<_>, DeseqError>>()?;
        let center = checked_sum(center, "lowess center")?;
        let mut variance = 0.0;
        for (point, weight) in &weights {
            let residual = checked_sub(x[*point], center, *point, "lowess centered x")?;
            let residual_sq = checked_mul(residual, residual, *point, "lowess centered x square")?;
            variance = checked_add(
                variance,
                checked_mul(
                    *weight,
                    residual_sq,
                    *point,
                    "lowess weighted variance term",
                )?,
                *point,
                "lowess variance",
            )?;
        }
        if variance.sqrt() > checked_mul(0.001, range, idx, "lowess slope variance threshold")? {
            let slope_factor = checked_div(
                checked_sub(target_x, center, idx, "lowess slope target delta")?,
                variance,
                idx,
                "lowess slope factor",
            )?;
            for (point, weight) in &mut weights {
                let centered = checked_sub(x[*point], center, *point, "lowess slope centered x")?;
                let adjustment = checked_add(
                    checked_mul(
                        slope_factor,
                        centered,
                        *point,
                        "lowess slope weight adjustment",
                    )?,
                    1.0,
                    *point,
                    "lowess slope weight factor",
                )?;
                *weight = checked_mul(*weight, adjustment, *point, "lowess slope weight")?;
            }
        }
    }

    let mut fitted = 0.0;
    for (point, weight) in weights {
        fitted = checked_add(
            fitted,
            checked_mul(weight, y[point], point, "lowess fitted weighted y")?,
            point,
            "lowess fitted value",
        )?;
    }
    Ok(Some(fitted))
}
