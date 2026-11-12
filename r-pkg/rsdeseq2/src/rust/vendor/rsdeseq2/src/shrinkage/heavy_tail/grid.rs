//! Streaming Cartesian quadrature and posterior interval queries.
//! Simpson normalization, rectangular probability queries and exponential edge
//! interpolation retain the established source-grid approximation policies.
use super::{
    analysis::{PosteriorGridDiagnostic, PosteriorIntervalKind, PosteriorIntervalOptions},
    prior::invalid,
};
use crate::DeseqError;
#[derive(Clone, Debug)]
pub(super) struct GridSummary {
    pub interval: [f64; 2],
    pub false_sign_rate: Option<f64>,
    pub threshold_rate: Option<f64>,
    pub diagnostic: PosteriorGridDiagnostic,
}
fn sequence(from: f64, to: f64, index: usize, intervals: usize) -> f64 {
    if index == intervals {
        to
    } else {
        from + index as f64 * ((to - from) / intervals as f64)
    }
}
fn integral(values: &[f64], width: f64) -> f64 {
    let n = values.len() - 1;
    width / 3.0
        * (values[0]
            + values[n]
            + 2.0 * (2..n).step_by(2).map(|i| values[i]).sum::<f64>()
            + 4.0 * (1..n).step_by(2).map(|i| values[i]).sum::<f64>())
}
fn solve_area(x0: f64, x1: f64, y0: f64, y1: f64, area: f64) -> f64 {
    let slope = (y1.ln() - y0.ln()) / (x1 - x0);
    x0 + (slope / y0 * area).ln_1p() / slope
}
fn allocate(n: usize) -> Result<Vec<f64>, DeseqError> {
    let mut out = vec![];
    out.try_reserve_exact(n)
        .map_err(|_| invalid("posterior grid allocation exceeds available memory"))?;
    out.resize(n, 0.0);
    Ok(out)
}
fn marginal(
    input: &GridInput<'_>,
    betas: &[f64],
    shift: f64,
    log_posterior: &mut impl FnMut(&[f64], f64) -> Result<f64, DeseqError>,
    evaluations: &mut usize,
) -> Result<Vec<f64>, DeseqError> {
    let GridInput {
        mode,
        sd,
        coefficient,
        options: o,
        ..
    } = *input;
    let points = if o.nuisance_points == 1 {
        1
    } else {
        o.nuisance_points
            .checked_pow(
                (mode.len() - 1)
                    .try_into()
                    .map_err(|_| invalid("posterior dimension overflows"))?,
            )
            .ok_or_else(|| invalid("posterior grid overflows"))?
    };
    let mut sums = allocate(betas.len())?;
    let mut peak = f64::NEG_INFINITY;
    let mut candidate = mode.to_vec();
    for row in 0..points {
        let mut digits = row;
        for j in 0..mode.len() {
            if j != coefficient {
                candidate[j] = if o.nuisance_points == 1 {
                    mode[j]
                } else {
                    let index = digits % o.nuisance_points;
                    digits /= o.nuisance_points;
                    mode[j]
                        + sd[j]
                            * sequence(
                                -o.nuisance_standard_deviations,
                                o.nuisance_standard_deviations,
                                index,
                                o.nuisance_points - 1,
                            )
                };
            }
        }
        for (j, &beta) in betas.iter().enumerate() {
            candidate[coefficient] = beta;
            let value = log_posterior(&candidate, shift)?;
            *evaluations += 1;
            if !value.is_finite() {
                return Err(invalid("posterior grid log density must be finite"));
            }
            if value > peak {
                let ratio = (peak - value).exp();
                for sum in &mut sums {
                    *sum *= ratio;
                }
                peak = value;
            }
            sums[j] += (value - peak).exp();
        }
    }
    Ok(sums)
}
#[derive(Clone, Copy)]
struct Edge {
    area: f64,
    zero: f64,
    slope: f64,
}
fn edge(first: f64, second: f64, x: f64, width: f64, left: bool) -> Result<Edge, DeseqError> {
    if first <= 0.0 || second <= 0.0 {
        return Ok(Edge {
            area: 0.0,
            zero: 0.0,
            slope: 0.0,
        });
    }
    let slope = if left {
        (second.ln() - first.ln()) / width
    } else {
        (first.ln() - second.ln()) / width
    };
    if (left && slope <= 0.0) || (!left && slope >= 0.0) || !slope.is_finite() {
        return Err(invalid(
            "posterior grid boundary does not have a decaying tail",
        ));
    }
    let amplitude = (first.ln() - slope * x).exp();
    let endpoint = (slope * x).exp();
    if !amplitude.is_finite() || !endpoint.is_finite() {
        return Ok(Edge {
            area: 0.0,
            zero: 0.0,
            slope,
        });
    }
    let zero = if left {
        amplitude / slope
    } else {
        -amplitude / slope
    };
    Ok(Edge {
        area: zero * endpoint,
        zero,
        slope,
    })
}
pub(super) struct GridInput<'a> {
    pub mode: &'a [f64],
    pub sd: &'a [f64],
    pub coefficient: usize,
    pub options: &'a PosteriorIntervalOptions,
    pub parameter_sd: Option<f64>,
    pub threshold: Option<f64>,
    pub flip_sign: bool,
}
pub(super) fn evaluate(
    input: GridInput<'_>,
    log_posterior: &mut impl FnMut(&[f64], f64) -> Result<f64, DeseqError>,
) -> Result<GridSummary, DeseqError> {
    let GridInput {
        mode,
        sd,
        coefficient,
        options: o,
        parameter_sd,
        threshold,
        flip_sign,
    } = input;
    let n = o.grid_intervals;
    let lower = mode[coefficient] - o.grid_standard_deviations * sd[coefficient];
    let upper = mode[coefficient] + o.grid_standard_deviations * sd[coefficient];
    let width = (upper - lower) / n as f64;
    if !lower.is_finite() || !upper.is_finite() || width <= 0.0 || !width.is_finite() {
        return Err(invalid("posterior grid bounds must be finite and distinct"));
    }
    let mut betas = allocate(n + 1)?;
    for (i, x) in betas.iter_mut().enumerate() {
        *x = sequence(lower, upper, i, n);
    }
    let mut evaluations = 0;
    let mut values = if let Some(error) = parameter_sd.filter(|_| o.nuisance_points != 1) {
        let mut mixed = allocate(n + 1)?;
        for j in 0..o.nuisance_points {
            let z = sequence(
                -o.nuisance_standard_deviations / 2.0,
                o.nuisance_standard_deviations / 2.0,
                j,
                o.nuisance_points - 1,
            );
            let part = marginal(&input, &betas, error * z, log_posterior, &mut evaluations)?;
            let area = integral(&part, width);
            let weight = (-0.5 * z * z).exp() * 0.3989422804014327;
            for (sum, v) in mixed.iter_mut().zip(part) {
                *sum += v / area * weight;
            }
        }
        mixed
    } else {
        marginal(&input, &betas, 0.0, log_posterior, &mut evaluations)?
    };
    let monotone_side = |left: bool| {
        let side = betas
            .iter()
            .zip(&values)
            .filter(|(b, _)| {
                if left {
                    **b < mode[coefficient]
                } else {
                    **b > mode[coefficient]
                }
            })
            .map(|(_, v)| *v)
            .collect::<Vec<_>>();
        side.windows(2)
            .all(|p| if left { p[1] >= p[0] } else { p[1] <= p[0] })
    };
    let na = f64::from_bits(0x7ff0_0000_0000_07a2);
    if !monotone_side(true) || !monotone_side(false) {
        return Ok(GridSummary {
            interval: [na; 2],
            false_sign_rate: None,
            threshold_rate: None,
            diagnostic: PosteriorGridDiagnostic {
                monotone: false,
                outside_left: None,
                outside_right: None,
                evaluations,
            },
        });
    }
    let left = edge(values[0], values[1], lower, width, true)?;
    let right = edge(values[n], values[n - 1], upper, width, false)?;
    let total = integral(&values, width) + left.area + right.area;
    if !total.is_finite() || total <= 0.0 {
        return Err(invalid("posterior grid normalization is invalid"));
    }
    for value in &mut values {
        *value /= total;
    }
    let outside_left = left.area / total;
    let outside_right = right.area / total;
    let mut cumulative = allocate(n + 1)?;
    let mut sum = outside_left;
    for (c, &v) in cumulative.iter_mut().zip(&values) {
        sum += width * v;
        *c = sum;
    }
    let quantile = |prob: f64| -> Result<f64, DeseqError> {
        let first = cumulative
            .iter()
            .position(|&v| v > prob)
            .ok_or_else(|| invalid("posterior grid quantile exceeds grid mass"))?;
        if first == 0 {
            // The source's zero R index drops the lower endpoint and recycles
            // the upper bound. Invert the already-normalized exponential tail
            // instead; this also handles mass in the first rectangular cell.
            if outside_left > 0.0 && left.slope > 0.0 {
                return Ok(lower + (prob / outside_left).ln() / left.slope);
            }
            return Err(invalid("posterior grid quantile has no left-tail mass"));
        }
        let j = first - 1;
        Ok(solve_area(
            betas[j],
            betas[j + 1],
            values[j],
            values[j + 1],
            prob - cumulative[j],
        ))
    };
    let alpha = (1.0 - o.level) / 2.0;
    let interval = if outside_left + outside_right >= 1.0 - o.level {
        let lo = if outside_left < alpha {
            quantile(alpha)?
        } else {
            solve_area(betas[0], betas[1], values[1], values[0], alpha)
        };
        let hi = if outside_right < alpha {
            quantile(1.0 - alpha)?
        } else {
            solve_area(betas[n - 1], betas[n], values[n - 1], values[n], alpha)
        };
        [lo, hi]
    } else if o.kind == PosteriorIntervalKind::EqualTail {
        [quantile(alpha)?, quantile(1.0 - alpha)?]
    } else {
        let mut sorted = values.clone();
        sorted.sort_by(|a, b| b.total_cmp(a));
        let mut mass = 0.0;
        let mut height = None;
        for v in sorted {
            mass += width * v;
            if mass >= o.level {
                height = Some(v);
                break;
            }
        }
        let height = height.ok_or_else(|| invalid("posterior HPD mass exceeds grid"))?;
        let indices = values
            .iter()
            .enumerate()
            .filter(|(_, v)| **v >= height)
            .map(|(i, _)| i)
            .collect::<Vec<_>>();
        if indices.len() < 2 || indices.windows(2).any(|v| v[1] != v[0] + 1) {
            return Err(invalid(
                "posterior HPD region is not one contiguous interval",
            ));
        }
        let a = indices[0];
        let b = *indices.last().unwrap();
        let rem = (mass - o.level) / 2.0;
        [
            solve_area(betas[a], betas[a + 1], values[a], values[a + 1], rem),
            solve_area(betas[b - 1], betas[b], values[b - 1], values[b], rem),
        ]
    };
    let lower_probability = |x: f64| {
        let last = betas.iter().rposition(|&v| v < x);
        if let Some(last) = last {
            if last == n {
                return None;
            }
            Some(
                width * values[..=last].iter().sum::<f64>() - (betas[last + 1] - x) * values[last]
                    + outside_left,
            )
        } else {
            Some(left.zero * (left.slope * x).exp() / total)
        }
    };
    let fsr = if parameter_sd.is_some() {
        None
    } else if mode[coefficient] > 0.0 {
        lower_probability(0.0)
    } else {
        if let Some(first) = betas.iter().position(|&v| v > 0.0) {
            if first == 0 {
                None
            } else {
                Some(
                    width * values[first..].iter().sum::<f64>()
                        + betas[first - 1] * values[first]
                        + outside_right,
                )
            }
        } else {
            Some(right.zero / total)
        }
    };
    let threshold_rate = threshold
        .and_then(|t| {
            lower_probability(if flip_sign && mode[coefficient] < 0.0 {
                -t
            } else {
                t
            })
        })
        .map(|p| {
            if flip_sign && mode[coefficient] < 0.0 {
                1.0 - p
            } else {
                p
            }
        });
    Ok(GridSummary {
        interval,
        false_sign_rate: fsr,
        threshold_rate,
        diagnostic: PosteriorGridDiagnostic {
            monotone: true,
            outside_left: Some(outside_left),
            outside_right: Some(outside_right),
            evaluations,
        },
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn fixed_grid_queries_match_original_reference_without_optimizer_drift() {
        let reference: serde_json::Value = serde_json::from_str(include_str!(
            "../../../tests/data/posterior_analysis_r462.json"
        ))
        .unwrap();
        let mut failures = vec![];
        for c in reference["grid_leaves"].as_array().unwrap() {
            let mean = c["mean"].as_f64().unwrap();
            let options = PosteriorIntervalOptions {
                kind: if c["kind"] == "HPD" {
                    PosteriorIntervalKind::HighestDensity
                } else {
                    PosteriorIntervalKind::EqualTail
                },
                level: 0.9,
                grid_standard_deviations: c["nsd"].as_f64().unwrap(),
                ..Default::default()
            };
            let mut density = |b: &[f64], _: f64| {
                Ok(-0.5 * ((b[1] - mean) / 0.7).powi(2)
                    - 0.5 * (b[0] / 1e8).powi(2)
                    - 0.5 * (b[1] / 1e8).powi(2))
            };
            let result = evaluate(
                GridInput {
                    mode: &[0.0, mean],
                    sd: &[1.0, 0.7],
                    coefficient: 1,
                    options: &options,
                    parameter_sd: None,
                    threshold: Some(0.2),
                    flip_sign: true,
                },
                &mut density,
            )
            .unwrap();
            let expected =
                |field: &str, i: usize| c[field][i].as_str().unwrap().parse::<f64>().unwrap();
            for (field, actual) in [
                ("fsr", result.false_sign_rate),
                ("threshold", result.threshold_rate),
            ] {
                let Some(actual) = actual else {
                    assert_eq!(c[field][0], "NA", "{} {field}", c["name"]);
                    continue;
                };
                if (actual - expected(field, 0)).abs() > 2e-10 {
                    failures.push(format!(
                        "{} {field}: {actual} != {}",
                        c["name"],
                        expected(field, 0)
                    ));
                }
            }
            for (i, &interval) in result.interval.iter().enumerate() {
                assert!(
                    (interval - expected("interval", i)).abs()
                        <= 2.0 * options.grid_standard_deviations * 0.7 / 50.0 + 1e-10,
                    "{} interval",
                    c["name"]
                );
            }
            assert!(result.diagnostic.monotone);
            assert_eq!(result.diagnostic.evaluations, 51 * 5);
            assert!(
                (result.diagnostic.outside_left.unwrap() - expected("outside", 0)).abs() < 2e-10
            );
            assert!(
                (result.diagnostic.outside_right.unwrap() - expected("outside", 1)).abs() < 2e-10
            );
        }
        assert!(failures.is_empty(), "{}", failures.join("\n"));
    }
}
