//! Grid preparation and component fitting shared by normal-likelihood shrinkage.
use super::{
    FixedMixture, MixtureFit, MixtureFitOptions, NormalMixture, TruncatedNormalMixture,
    UniformMixture, fit_mixture_proportions,
};
use crate::{
    DeseqError,
    matrix::{MatrixOrder, MatrixView},
};

/// Component grid families accepted by the normal-likelihood shrinkage engine.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MixtureFamily {
    Normal,
    Uniform,
    HalfUniform,
    PositiveUniform,
    NegativeUniform,
    HalfNormal,
}
/// Dirichlet specification for a generated grid. Parameters below one require
/// a variational method and are outside this convex weight-fitting operation.
#[derive(Clone, Debug)]
pub enum MixturePriorPenalty {
    Uniform,
    /// Dirichlet 1/K parameters for variational inference.
    Unit,
    NullBiased(f64),
    Parameters(Vec<f64>),
}
/// Original ashr method shorthand, distinct from individual grid controls.
#[derive(Clone, Copy, Debug)]
pub enum AdaptiveMethod {
    Fdr,
    Shrink,
}
/// Preparation controls. Observation-scale transformations belong to the caller.
#[derive(Clone, Debug)]
pub struct MixtureGridOptions {
    pub family: MixtureFamily,
    pub scales: Option<Vec<f64>>,
    pub point_mass: bool,
    pub grid_multiplier: f64,
    pub mode: f64,
    pub range: (f64, f64),
    pub penalty: MixturePriorPenalty,
    /// Source EM initialization, before asymmetric component expansion.
    pub concentrated_initialization: bool,
    pub variational_prior: bool,
}
impl Default for MixtureGridOptions {
    fn default() -> Self {
        Self {
            family: MixtureFamily::Normal,
            scales: None,
            point_mass: false,
            grid_multiplier: std::f64::consts::SQRT_2,
            mode: 0.0,
            range: (f64::NEG_INFINITY, f64::INFINITY),
            penalty: MixturePriorPenalty::Uniform,
            concentrated_initialization: false,
            variational_prior: false,
        }
    }
}
/// Component parameters retain their generated order for result metadata.
#[derive(Clone, Debug)]
pub enum MixtureComponents {
    Normal {
        means: Vec<f64>,
        sds: Vec<f64>,
    },
    Uniform {
        lower: Vec<f64>,
        upper: Vec<f64>,
    },
    TruncatedNormal {
        means: Vec<f64>,
        sds: Vec<f64>,
        lower: Vec<f64>,
        upper: Vec<f64>,
    },
}
/// Prepared grid, initial probabilities, penalty and explicit warning messages.
#[derive(Clone, Debug)]
pub struct MixtureGrid {
    pub components: MixtureComponents,
    pub initial_weights: Vec<f64>,
    pub dirichlet: Vec<f64>,
    pub warnings: Vec<String>,
}
fn invalid(reason: &str) -> DeseqError {
    DeseqError::InvalidOptions {
        reason: format!("mixture grid: {reason}"),
    }
}

/// Construct a supplied-scale or geometric automatic grid. Missing observations
/// and zero/infinite SEs do not choose the automatic scale range. This operation
/// neither runs an R function nor fits the mixture weights.
pub fn prepare_mixture_grid(
    estimates: &[f64],
    standard_errors: &[f64],
    options: &MixtureGridOptions,
) -> Result<MixtureGrid, DeseqError> {
    if standard_errors.len() != estimates.len() && standard_errors.len() != 1 {
        return Err(invalid("SEs must have length one or match estimates"));
    }
    if !options.mode.is_finite()
        || !options.grid_multiplier.is_finite()
        || options.grid_multiplier <= 1.0
        || options.range.0.is_nan()
        || options.range.1.is_nan()
        || options.range.0 > options.range.1
    {
        return Err(invalid(
            "finite mode, multiplier > 1 and ordered nonmissing range are required",
        ));
    }
    let mut scales = match &options.scales {
        Some(scales) => scales.clone(),
        None => {
            let mut minimum = f64::INFINITY;
            let mut signal = f64::NEG_INFINITY;
            for (i, &x) in estimates.iter().enumerate() {
                let s = standard_errors[if standard_errors.len() == 1 { 0 } else { i }];
                if x.is_nan() || s.is_nan() || s == 0.0 || s == f64::INFINITY {
                    continue;
                }
                if !x.is_finite() || !s.is_finite() || s < 0.0 {
                    return Err(invalid(
                        "automatic grid requires finite observations and nonnegative SEs",
                    ));
                }
                minimum = minimum.min(s);
                signal = signal.max((x - options.mode).powi(2) - s * s);
            }
            let minimum = minimum / 10.0;
            let maximum = if signal <= 0.0 {
                8.0 * minimum
            } else {
                2.0 * signal.sqrt()
            };
            if minimum <= 0.0 || !minimum.is_finite() || !maximum.is_finite() {
                return Err(invalid("no representable automatic scale range"));
            }
            let steps = ((maximum / minimum).log2() / options.grid_multiplier.log2()).ceil();
            if !steps.is_finite() || !(0.0..=1_000_000.0).contains(&steps) {
                return Err(invalid(
                    "automatic grid size is not representable or exceeds one million intervals",
                ));
            }
            (-(steps as i32)..=0)
                .map(|exponent| options.grid_multiplier.powf(f64::from(exponent)) * maximum)
                .collect()
        }
    };
    if scales.is_empty() || scales.iter().any(|v| !v.is_finite() || *v < 0.0) {
        return Err(invalid("scales must be nonempty, finite and nonnegative"));
    }
    if options.point_mass {
        scales.insert(0, 0.0);
    }
    let null = scales
        .iter()
        .enumerate()
        .min_by(|a, b| a.1.total_cmp(b.1))
        .unwrap()
        .0;
    let mut weights = vec![1.0; scales.len()];
    if options.concentrated_initialization {
        let k = scales.len();
        let n = estimates.len();
        if k < n {
            weights.fill(1.0 / n as f64);
            weights[null] = (n - k + 1) as f64 / n as f64;
        } else {
            weights.fill(1.0 / k as f64);
        }
        let sum = weights.iter().sum::<f64>();
        for w in &mut weights {
            *w /= sum;
        }
    }
    let mut dirichlet = match &options.penalty {
        MixturePriorPenalty::Uniform => vec![1.0; scales.len()],
        MixturePriorPenalty::Unit if options.variational_prior => {
            vec![1.0 / scales.len() as f64; scales.len()]
        }
        MixturePriorPenalty::NullBiased(weight) => {
            let mut p = vec![1.0; scales.len()];
            p[null] = *weight;
            p
        }
        MixturePriorPenalty::Parameters(values) if values.len() == scales.len() => values.clone(),
        _ => return Err(invalid("penalty length must match the initial scale grid")),
    };
    if dirichlet
        .iter()
        .any(|v| !v.is_finite() || *v <= 0.0 || (!options.variational_prior && *v < 1.0))
    {
        return Err(invalid(
            "convex fitting requires finite Dirichlet parameters >= 1",
        ));
    }
    let mut warnings = vec![];
    let mut family = options.family;
    if matches!(
        family,
        MixtureFamily::Uniform
            | MixtureFamily::HalfUniform
            | MixtureFamily::PositiveUniform
            | MixtureFamily::NegativeUniform
    ) {
        if options.range.0 == options.mode {
            family = MixtureFamily::PositiveUniform;
        } else if options.range.1 == options.mode {
            family = MixtureFamily::NegativeUniform;
        }
    }
    if options.family == MixtureFamily::HalfUniform
        && !matches!(options.penalty, MixturePriorPenalty::NullBiased(_))
    {
        warnings.push("Use of halfuniform without nullbiased prior can lead to misleading local false sign rates, and so is not recommended".into());
    }
    let mode = options.mode;
    let components = match family {
        MixtureFamily::Normal => MixtureComponents::Normal {
            means: vec![mode; scales.len()],
            sds: scales,
        },
        MixtureFamily::HalfNormal | MixtureFamily::HalfUniform => {
            let duplicate = (0..scales.len())
                .filter(|&i| scales[null] > 0.0 || i != null)
                .collect::<Vec<_>>();
            let divisor = if scales[null] > 0.0 {
                2.0
            } else {
                2.0 - weights[null]
            };
            let mut expanded = scales.clone();
            for &i in &duplicate {
                expanded.push(scales[i]);
                weights.push(weights[i]);
                dirichlet.push(dirichlet[i]);
            }
            for w in &mut weights {
                *w /= divisor;
            }
            let mut lower = vec![
                if family == MixtureFamily::HalfNormal {
                    f64::NEG_INFINITY
                } else {
                    mode
                };
                scales.len()
            ];
            let mut upper = vec![
                if family == MixtureFamily::HalfNormal {
                    0.0
                } else {
                    mode
                };
                scales.len()
            ];
            if family == MixtureFamily::HalfNormal {
                lower.extend(vec![0.0; duplicate.len()]);
                upper.extend(vec![f64::INFINITY; duplicate.len()]);
                MixtureComponents::TruncatedNormal {
                    means: vec![mode; expanded.len()],
                    sds: expanded,
                    lower,
                    upper,
                }
            } else {
                for (a, &sd) in lower.iter_mut().zip(&scales) {
                    *a = mode - sd;
                }
                lower.extend(vec![mode; duplicate.len()]);
                upper.extend(duplicate.iter().map(|&i| mode + scales[i]));
                MixtureComponents::Uniform { lower, upper }
            }
        }
        _ => {
            let lower = scales
                .iter()
                .map(|s| {
                    if family == MixtureFamily::PositiveUniform {
                        mode
                    } else {
                        mode - s
                    }
                })
                .collect();
            let upper = scales
                .iter()
                .map(|s| {
                    if family == MixtureFamily::NegativeUniform {
                        mode
                    } else {
                        mode + s
                    }
                })
                .collect();
            MixtureComponents::Uniform { lower, upper }
        }
    };
    let components = if let MixtureComponents::Uniform { lower, upper } = components {
        let mut a = vec![];
        let mut b = vec![];
        let mut w = vec![];
        let mut p = vec![];
        for i in 0..lower.len() {
            let left = lower[i].max(options.range.0);
            let right = upper[i].min(options.range.1);
            if !a.iter().zip(&b).any(|(&aa, &bb)| aa == left && bb == right) {
                a.push(left);
                b.push(right);
                w.push(weights[i]);
                p.push(dirichlet[i]);
            }
        }
        let total: f64 = w.iter().sum();
        for x in &mut w {
            *x /= total;
        }
        weights = w;
        dirichlet = p;
        MixtureComponents::Uniform { lower: a, upper: b }
    } else {
        if options.range != (f64::NEG_INFINITY, f64::INFINITY) {
            warnings
                .push("Can't constrain grange for normal/halfnormal mixture prior case.".into());
        }
        components
    };
    let grid = MixtureGrid {
        components,
        initial_weights: weights,
        dirichlet,
        warnings,
    };
    grid.prior(&grid.initial_weights)?;
    Ok(grid)
}
impl MixtureGrid {
    /// Validate a supplied component grid with the requested prior penalty.
    /// Supplied mixtures use their first component as the null component.
    pub fn from_supplied(
        components: MixtureComponents,
        weights: Vec<f64>,
        penalty: &MixturePriorPenalty,
    ) -> Result<Self, DeseqError> {
        Self::from_supplied_with_policy(components, weights, penalty, false)
    }
    pub fn from_supplied_variational(
        components: MixtureComponents,
        weights: Vec<f64>,
        penalty: &MixturePriorPenalty,
    ) -> Result<Self, DeseqError> {
        Self::from_supplied_with_policy(components, weights, penalty, true)
    }
    fn from_supplied_with_policy(
        components: MixtureComponents,
        weights: Vec<f64>,
        penalty: &MixturePriorPenalty,
        variational: bool,
    ) -> Result<Self, DeseqError> {
        let mut dirichlet = vec![1.0; weights.len()];
        match penalty {
            MixturePriorPenalty::Uniform => {}
            MixturePriorPenalty::Unit if variational => dirichlet.fill(1.0 / weights.len() as f64),
            MixturePriorPenalty::Unit => {
                return Err(invalid("unit prior requires variational fitting"));
            }
            MixturePriorPenalty::NullBiased(w) => {
                if let Some(first) = dirichlet.first_mut() {
                    *first = *w;
                }
            }
            MixturePriorPenalty::Parameters(p) => dirichlet = p.clone(),
        }
        if dirichlet.len() != weights.len()
            || dirichlet
                .iter()
                .any(|v| !v.is_finite() || *v <= 0.0 || (!variational && *v < 1.0))
        {
            return Err(invalid(
                "convex fitting requires matching finite Dirichlet parameters >=1",
            ));
        }
        let grid = Self {
            components,
            initial_weights: weights,
            dirichlet,
            warnings: vec![],
        };
        grid.prior(&grid.initial_weights)?;
        Ok(grid)
    }

    /// Construct the same component prior with explicitly supplied weights.
    pub fn prior(&self, weights: &[f64]) -> Result<FixedMixture, DeseqError> {
        Ok(match &self.components {
            MixtureComponents::Normal { means, sds } => FixedMixture::Normal(NormalMixture::new(
                weights.to_vec(),
                means.clone(),
                sds.clone(),
            )?),
            MixtureComponents::Uniform { lower, upper } => FixedMixture::Uniform(
                UniformMixture::new(weights.to_vec(), lower.clone(), upper.clone())?,
            ),
            MixtureComponents::TruncatedNormal {
                means,
                sds,
                lower,
                upper,
            } => FixedMixture::TruncatedNormal(TruncatedNormalMixture::new(
                weights.to_vec(),
                means.clone(),
                sds.clone(),
                lower.clone(),
                upper.clone(),
            )?),
        })
    }
    /// Fit component probabilities to finite observations using the shared
    /// native optimizer. The caller selects/excludes rows before this operation.
    pub fn fit(
        &self,
        estimates: &[f64],
        standard_errors: &[f64],
        observation_weights: Option<&[f64]>,
        options: MixtureFitOptions,
    ) -> Result<MixtureFit, DeseqError> {
        if estimates.len() != standard_errors.len() {
            return Err(invalid("fit estimates and SEs must have matching lengths"));
        }
        estimates
            .len()
            .checked_mul(self.initial_weights.len())
            .ok_or_else(|| invalid("likelihood matrix shape overflows"))?;
        let prior = self.prior(&self.initial_weights)?;
        let mut likelihoods = vec![];
        for (&x, &s) in estimates.iter().zip(standard_errors) {
            likelihoods.extend(prior.log_likelihoods(x, s)?);
        }
        fit_mixture_proportions(
            MatrixView::new(
                estimates.len(),
                self.initial_weights.len(),
                &likelihoods,
                MatrixOrder::RowMajor,
            )?,
            observation_weights,
            &self.dirichlet,
            Some(&self.initial_weights),
            options,
        )
    }
}

impl std::str::FromStr for MixtureFamily {
    type Err = DeseqError;
    fn from_str(name: &str) -> Result<Self, Self::Err> {
        match name {
            "normal" => Ok(Self::Normal),
            "uniform" => Ok(Self::Uniform),
            "halfuniform" => Ok(Self::HalfUniform),
            "+uniform" => Ok(Self::PositiveUniform),
            "-uniform" => Ok(Self::NegativeUniform),
            "halfnormal" => Ok(Self::HalfNormal),
            _ => Err(invalid("unknown normal-likelihood component family")),
        }
    }
}

impl MixtureGridOptions {
    /// Apply method defaults while preserving explicit-argument conflicts.
    pub fn with_ashr_method(
        mut self,
        method: Option<AdaptiveMethod>,
        point_mass: Option<bool>,
        penalty: Option<MixturePriorPenalty>,
        nullweight: f64,
    ) -> Result<Self, DeseqError> {
        if method.is_some() && (point_mass.is_some() || penalty.is_some()) {
            return Err(invalid(
                "specify either method or pointmass/prior, not both",
            ));
        }
        if matches!(method, Some(AdaptiveMethod::Shrink)) {
            self.point_mass = false;
            self.penalty = MixturePriorPenalty::Uniform;
        } else {
            self.point_mass = point_mass.unwrap_or(true);
            self.penalty = penalty.unwrap_or(MixturePriorPenalty::NullBiased(nullweight));
        }
        Ok(self)
    }
}
