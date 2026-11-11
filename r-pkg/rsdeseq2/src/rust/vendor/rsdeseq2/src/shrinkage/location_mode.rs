//! Standalone likelihood profiling for an estimated unimodal-prior location.
use super::location_missing::{NA, excluded, is_na, log_marginal};
use super::{LocationLikelihood, LocationShrinkageOptions, location_fit::fit_location};
use crate::DeseqError;

/// Search bounds, specified on the original estimate scale, before alpha scaling.
#[derive(Clone, Copy, Debug)]
pub enum PriorModeSearch {
    ObservationRange,
    /// Expand the data range to include one additional coordinate.
    ExpandedObservationRange(f64),
    /// Explicit native bounds, independent of the data range.
    Interval(f64, f64),
}
/// Profile maximum and numerical work, with warnings in evaluation order.
#[derive(Clone, Debug)]
pub struct LocationModeEstimate {
    pub mode: f64,
    pub log_likelihood: f64,
    pub evaluations: usize,
    pub bounds: (f64, f64),
    pub warnings: Vec<String>,
}
fn invalid(reason: &str) -> DeseqError {
    DeseqError::InvalidOptions {
        reason: format!("prior mode search: {reason}"),
    }
}

// Golden-section fallback and inverse parabolic interpolation of three finite
// function values. The bracket always contains the current best point. Arithmetic
// and tolerances follow the published Brent method, without an R implementation.
fn minimize(
    mut objective: impl FnMut(f64) -> Result<f64, DeseqError>,
    mut lower: f64,
    mut upper: f64,
    tolerance: f64,
) -> Result<(f64, f64, usize), DeseqError> {
    let golden = (3.0 - 5.0_f64.sqrt()) / 2.0;
    let mut best = lower + golden * (upper - lower);
    let mut previous = best;
    let mut older = best;
    let mut value = objective(best)?;
    if !value.is_finite() {
        return Err(invalid("profile likelihood must be finite"));
    }
    let mut previous_value = value;
    let mut older_value = value;
    let mut step = 0.0_f64;
    let mut history = 0.0_f64;
    for evaluations in 1..=1000 {
        let midpoint = (lower + upper) * 0.5;
        let spacing = f64::EPSILON.sqrt() * best.abs() + tolerance / 3.0;
        if (best - midpoint).abs() <= 2.0 * spacing - (upper - lower) * 0.5 {
            return Ok((best, value, evaluations));
        }
        if evaluations == 1000 {
            break;
        }
        let mut interpolate = false;
        if history.abs() > spacing {
            let first = (best - previous) * (value - older_value);
            let second = (best - older) * (value - previous_value);
            let mut numerator = (best - older) * second - (best - previous) * first;
            let mut denominator = 2.0 * (second - first);
            if denominator > 0.0 {
                numerator = -numerator;
            }
            denominator = denominator.abs();
            let saved = history;
            history = step;
            if numerator.abs() < (0.5 * denominator * saved).abs()
                && numerator > denominator * (lower - best)
                && numerator < denominator * (upper - best)
            {
                step = numerator / denominator;
                let trial = best + step;
                if trial - lower < 2.0 * spacing || upper - trial < 2.0 * spacing {
                    step = spacing.copysign(midpoint - best);
                }
                interpolate = true;
            }
        }
        if !interpolate {
            history = if best >= midpoint {
                lower - best
            } else {
                upper - best
            };
            step = golden * history;
        }
        let trial = best
            + if step.abs() >= spacing {
                step
            } else {
                spacing.copysign(step)
            };
        let trial_value = objective(trial)?;
        if !trial_value.is_finite() {
            return Err(invalid("profile likelihood must be finite"));
        }
        if trial_value <= value {
            if trial >= best {
                lower = best;
            } else {
                upper = best;
            }
            older = previous;
            older_value = previous_value;
            previous = best;
            previous_value = value;
            best = trial;
            value = trial_value;
        } else {
            if trial < best {
                lower = trial;
            } else {
                upper = trial;
            }
            if trial_value <= previous_value || previous == best {
                older = previous;
                older_value = previous_value;
                previous = trial;
                previous_value = trial_value;
            } else if trial_value <= older_value || older == best || older == previous {
                older = trial;
                older_value = trial_value;
            }
        }
    }
    Err(invalid("bounded profile search exhausted 1000 evaluations"))
}

/// Profile the unpenalized marginal likelihood, refitting mixture weights at
/// every trial location. Observation weights and Dirichlet priors affect those
/// fits; the outer likelihood itself remains unweighted, as in original ashr.
/// No posterior moment, sign-rate, pruning or CDF query is evaluated here.
pub fn estimate_location_mode(
    estimates: &[f64],
    standard_errors: &[f64],
    options: &LocationShrinkageOptions,
    search: PriorModeSearch,
) -> Result<LocationModeEstimate, DeseqError> {
    if options.supplied_grid.is_some() || options.fixed_prior {
        return Err(invalid(
            "an estimated mode requires an unfixed generated prior",
        ));
    }
    let (lower, upper) = match search {
        PriorModeSearch::Interval(a, b) => (a, b),
        PriorModeSearch::ObservationRange | PriorModeSearch::ExpandedObservationRange(_) => {
            if estimates.is_empty() || estimates.iter().any(|x| !x.is_finite()) {
                return Err(invalid(
                    "finite nonmissing estimates are required for data-derived bounds",
                ));
            }
            let mut lower = estimates.iter().copied().fold(f64::INFINITY, f64::min);
            let mut upper = estimates.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            if let PriorModeSearch::ExpandedObservationRange(seed) = search {
                if !seed.is_finite() {
                    return Err(invalid("the range expansion must be finite"));
                }
                lower = lower.min(seed);
                upper = upper.max(seed);
            }
            (lower, upper)
        }
    };
    if !lower.is_finite()
        || !upper.is_finite()
        || lower >= upper
        || !(upper - lower).is_finite()
        || !(lower + upper).is_finite()
    {
        return Err(invalid("finite increasing search bounds are required"));
    }
    let tolerance = (upper - lower) * f64::EPSILON.sqrt().sqrt();
    if tolerance == 0.0 {
        return Err(invalid("the search interval is too small"));
    }
    let mut trial_options = options.clone();
    trial_options.thresholds.clear();
    let mut warnings = vec![];
    let mut warning_cache = vec![];
    let (mode, objective, evaluations) = minimize(
        |mode| {
            trial_options.grid.mode = mode;
            let fitted = fit_location(estimates, standard_errors, &trial_options)?;
            let mut trial_warnings = fitted.warnings;
            let prior = fitted
                .grid
                .prior_with_likelihood(&fitted.fitted_weights, &trial_options.likelihood)?;
            if let LocationLikelihood::Restricted(f) = &trial_options.likelihood {
                trial_warnings.extend(super::likelihood_diagnostics::likelihood(
                    f,
                    &fitted.grid.components,
                    &fitted.x,
                    &fitted.se,
                ));
            }
            warnings.extend_from_slice(&trial_warnings);
            warning_cache.push((mode.to_bits(), trial_warnings));
            let mut ll = 0.0;
            for i in 0..estimates.len() {
                let x = fitted.x[i];
                let se = fitted.se[i];
                let value = if se > 0.0 && se.is_finite() && x.is_finite() {
                    log_marginal(
                        &prior.log_likelihoods_at(x, se, i, estimates.len())?,
                        &fitted.fitted_weights,
                    )
                } else {
                    excluded(
                        &fitted.grid.components,
                        &fitted.fitted_weights,
                        x,
                        se,
                        &[],
                        0,
                        0,
                    )
                    .1
                };
                ll += value
                    - options.alpha
                        * standard_errors[if standard_errors.len() == 1 { 0 } else { i }].ln();
            }
            if fitted.x.iter().chain(&fitted.se).any(|x| is_na(*x)) {
                ll = NA;
            }
            if !ll.is_finite() {
                return Err(invalid("profile likelihood must be finite"));
            }
            Ok(-ll)
        },
        lower,
        upper,
        tolerance,
    )?;
    // stats::optimize re-evaluates its selected point to populate $objective.
    // That deterministic fit is already cached: reproduce its warnings without
    // another fit. The evaluation count reports actual native fits only.
    if let Some((_, messages)) = warning_cache
        .iter()
        .find(|(bits, _)| *bits == mode.to_bits())
    {
        warnings.extend_from_slice(messages);
    }
    Ok(LocationModeEstimate {
        mode,
        log_likelihood: -objective,
        evaluations,
        bounds: (lower, upper),
        warnings,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn bounded_search_has_independent_interior_boundary_and_flat_solutions() {
        let (x, f, _) = minimize(|x| Ok((x - 0.37).powi(2)), -2.0, 3.0, 1e-10).unwrap();
        assert!((x - 0.37).abs() < 1e-8 && f < 1e-16);
        let (x, _, _) = minimize(|x| Ok(-x), -2.0, 3.0, 1e-10).unwrap();
        assert!((x - 3.0).abs() < 1e-7);
        let (x, f, _) = minimize(|_| Ok(4.0), -2.0, 3.0, 1e-10).unwrap();
        assert!((-2.0..3.0).contains(&x) && f == 4.0);
        assert!(
            minimize(
                |_| Err(invalid("injected objective failure")),
                0.0,
                1.0,
                1e-6
            )
            .is_err()
        );
    }
}
