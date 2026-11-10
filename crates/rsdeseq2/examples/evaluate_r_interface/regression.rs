//! Data-only transport for the shared general-regression posterior engine.
use super::{
    bits, heavy_tail::Matrix, parse_numeric_input as number, parse_numeric_inputs as numbers,
};
use rsdeseq2::{
    matrix::{MatrixOrder, MatrixView},
    shrinkage::{
        CauchyCoefficientPrior, CoefficientMle, CoefficientScale, NegativeBinomialPosterior,
        NegativeBinomialPosteriorInput, NegativeBinomialPosteriorOptions, PosteriorAnalysisOptions,
        PosteriorIntervalKind, PosteriorIntervalOptions, PosteriorOffsets, PosteriorRng,
        RegressionLikelihood, RegressionMethod, RegressionOptimizer, RegressionPosteriorInput,
        RegressionPosteriorOptions, RestrictedRegressionLikelihood,
        fit_negative_binomial_posterior, fit_negative_binomial_posterior_with_rng,
        fit_regression_posterior, fit_regression_posterior_with_rng,
    },
};
use serde::Deserialize;
use serde_json::{Value, json};
use std::{error::Error, sync::Arc};
#[derive(Deserialize)]
pub struct Input {
    pub regression_cases: Vec<Case>,
    #[serde(default)]
    pub density_cases: Vec<Leaf>,
}
#[derive(Deserialize)]
pub struct Case {
    pub name: String,
    pub observations: Matrix,
    pub design: Matrix,
    pub parameters: Option<Matrix>,
    pub offsets: Option<Matrix>,
    pub weights: Option<Matrix>,
    pub mle: Option<Matrix>,
    pub coefficient: usize,
    pub nuisance: Vec<usize>,
    pub scale: String,
    pub nuisance_scale: String,
    pub prior_df: String,
    pub prior_mean: String,
    pub nuisance_mean: String,
    pub method: String,
    pub optimizer: String,
    pub bounds: [String; 2],
    pub log_link: bool,
    pub threshold: Option<String>,
    pub analysis: Option<Analysis>,
    pub likelihood: Option<Value>,
    pub rng_state: Option<Vec<String>>,
    pub flip_sign: Option<bool>,
}
#[derive(Deserialize)]
pub struct Analysis {
    pub kind: String,
    pub level: String,
    pub grid_intervals: usize,
    pub grid_standard_deviations: String,
    pub nuisance_points: usize,
    pub nuisance_standard_deviations: String,
    pub contrasts: Option<Matrix>,
    pub parameter_standard_errors: Option<Vec<String>>,
}
impl Analysis {
    fn options(&self) -> Result<PosteriorAnalysisOptions, Box<dyn Error>> {
        Ok(PosteriorAnalysisOptions {
            interval: Some(PosteriorIntervalOptions {
                kind: match self.kind.as_str() {
                    "laplace" => PosteriorIntervalKind::Laplace,
                    "HPD" => PosteriorIntervalKind::HighestDensity,
                    "credible" => PosteriorIntervalKind::EqualTail,
                    _ => return Err("invalid interval kind".into()),
                },
                level: number(&self.level)?,
                grid_intervals: self.grid_intervals,
                grid_standard_deviations: number(&self.grid_standard_deviations)?,
                nuisance_points: self.nuisance_points,
                nuisance_standard_deviations: number(&self.nuisance_standard_deviations)?,
            }),
            contrasts: self
                .contrasts
                .as_ref()
                .map(|m| {
                    let v = numbers(&m.values)?;
                    let view = m.view(&v)?;
                    Ok::<_, Box<dyn Error>>(rsdeseq2::matrix::RowMajorMatrix::from_row_major(
                        m.rows,
                        m.columns,
                        (0..m.rows)
                            .flat_map(|i| (0..m.columns).map(move |j| *view.get(i, j).unwrap()))
                            .collect(),
                    )?)
                })
                .transpose()?,
            parameter_standard_errors: self
                .parameter_standard_errors
                .as_ref()
                .map(|s| numbers(s))
                .transpose()?,
            ..Default::default()
        })
    }
}
#[derive(Deserialize)]
pub struct Leaf {
    pub name: String,
    pub likelihood: Value,
    pub observations: Vec<String>,
    pub parameters: Vec<String>,
}
pub fn likelihood(v: &Value) -> Result<RestrictedRegressionLikelihood, Box<dyn Error>> {
    let (arguments, expression, constants) = super::likelihood::parts(v)?;
    Ok(RestrictedRegressionLikelihood::new(
        arguments, expression, constants,
    )?)
}
pub fn leaf(c: &Leaf) -> Result<Vec<f64>, Box<dyn Error>> {
    let y = numbers(&c.observations)?;
    let parameters = numbers(&c.parameters)?;
    let n = y.len();
    let design = vec![1.0; n];
    Ok(likelihood(&c.likelihood)?.evaluate(
        &y,
        MatrixView::new(n, 1, &design, MatrixOrder::RowMajor)?,
        &[0.0],
        Some(&parameters),
        &vec![0.0; n],
    )?)
}
pub fn run(c: &Case) -> Result<NegativeBinomialPosterior, Box<dyn Error>> {
    let y = numbers(&c.observations.values)?;
    let x = numbers(&c.design.values)?;
    let parameters = c
        .parameters
        .as_ref()
        .map(|m| numbers(&m.values))
        .transpose()?;
    let offsets = c.offsets.as_ref().map(|m| numbers(&m.values)).transpose()?;
    let weights = c.weights.as_ref().map(|m| numbers(&m.values)).transpose()?;
    let mle = c.mle.as_ref().map(|m| numbers(&m.values)).transpose()?;
    let input = RegressionPosteriorInput {
        observations: c.observations.view(&y)?,
        design: c.design.view(&x)?,
        parameters: c
            .parameters
            .as_ref()
            .zip(parameters.as_ref())
            .map(|(m, v)| m.view(v))
            .transpose()?,
        offsets: c
            .offsets
            .as_ref()
            .zip(offsets.as_ref())
            .map(|(m, v)| m.view(v))
            .transpose()?
            .map(PosteriorOffsets::Log)
            .unwrap_or_default(),
        weights: c
            .weights
            .as_ref()
            .zip(weights.as_ref())
            .map(|(m, v)| m.view(v))
            .transpose()?,
        mle: c
            .mle
            .as_ref()
            .zip(mle.as_ref())
            .map(|(m, v)| {
                if m.columns != 2 || m.rows != c.observations.rows {
                    return Err("invalid MLE fixture shape");
                }
                Ok(CoefficientMle {
                    estimates: &v[..m.rows],
                    standard_errors: &v[m.rows..],
                    scale: CoefficientScale::Natural,
                })
            })
            .transpose()?,
    };
    let mut options = RegressionPosteriorOptions {
        posterior: NegativeBinomialPosteriorOptions {
            prior: CauchyCoefficientPrior {
                nuisance: c.nuisance.clone(),
                scale: number(&c.scale)?,
                nuisance_scale: number(&c.nuisance_scale)?,
                adapted_variance: None,
            },
            coefficient: Some(c.coefficient),
            threshold: c.threshold.as_ref().map(|v| number(v)).transpose()?,
            analysis: c
                .analysis
                .as_ref()
                .map(Analysis::options)
                .transpose()?
                .unwrap_or_default(),
            flip_sign: c.flip_sign.unwrap_or(true),
            ..Default::default()
        },
        likelihood: if c.method.starts_with("betabin") {
            RegressionLikelihood::BetaBinomial
        } else {
            match &c.likelihood {
                None => RegressionLikelihood::NegativeBinomial,
                Some(v) => RegressionLikelihood::Restricted(Arc::new(likelihood(v)?)),
            }
        },
        method: match c.method.as_str() {
            "general" | "nbinomC*" | "nbinomCR" | "nbinomC" => RegressionMethod::General,
            "nbinomR" => RegressionMethod::NegativeBinomialR,
            "betabinR" => RegressionMethod::BetaBinomialR,
            "betabinCR" => RegressionMethod::BetaBinomialCR,
            "betabinC" => RegressionMethod::BetaBinomialC,
            "betabinC*" => RegressionMethod::BetaBinomialCStar,
            _ => return Err("invalid fixture method".into()),
        },
        optimizer: match c.optimizer.as_str() {
            "BFGS" => RegressionOptimizer::Bfgs,
            "L-BFGS-B" => RegressionOptimizer::Lbfgsb,
            _ => return Err("invalid fixture optimizer".into()),
        },
        prior_df: number(&c.prior_df)?,
        prior_noncentrality: number(&c.prior_mean)?,
        nuisance_mean: number(&c.nuisance_mean)?,
        log_link: c.log_link,
        bounds: (number(&c.bounds[0])?, number(&c.bounds[1])?),
    };
    if c.method != "general" && c.analysis.as_ref().is_some_and(|a| a.kind != "laplace") {
        options.posterior.analysis.grid_likelihood = c
            .likelihood
            .as_ref()
            .map(|v| likelihood(v).map(|f| RegressionLikelihood::Restricted(Arc::new(f))))
            .transpose()?;
    }
    if c.method.starts_with("nbinomC") {
        let parameter = input.parameters.ok_or("missing NB dispersions")?;
        if parameter.n_cols() != 1 {
            return Err("NB dispersions require one column".into());
        }
        let dispersions = (0..parameter.n_rows())
            .map(|i| parameter.get(i, 0).copied())
            .collect::<Result<Vec<_>, _>>()?;
        let input = NegativeBinomialPosteriorInput {
            counts: input.observations,
            design: input.design,
            dispersions: &dispersions,
            offsets: input.offsets,
            weights: input.weights,
            mle: input.mle,
        };
        let mut posterior = options.posterior;
        posterior.check_opposite_start = c.method == "nbinomCR";
        posterior.map_only = c.method != "nbinomCR";
        return Ok(if let Some(state) = &c.rng_state {
            let mut rng = PosteriorRng::from_r_seed(&super::random::state(state)?)?;
            fit_negative_binomial_posterior_with_rng(input, &posterior, &mut rng)?
        } else {
            fit_negative_binomial_posterior(input, &posterior)?
        });
    }
    Ok(if let Some(state) = &c.rng_state {
        let mut rng = PosteriorRng::from_r_seed(&super::random::state(state)?)?;
        fit_regression_posterior_with_rng(input, &options, &mut rng)?
    } else {
        fit_regression_posterior(input, &options)?
    })
}
pub fn evaluate(input: Input) -> Result<Value, Box<dyn Error>> {
    let cases = input
        .regression_cases
        .iter()
        .map(|c| match run(c) {
            Ok(r) => json!({"name": c.name, "value": super::heavy_tail::encode(r, c.coefficient)}),
            Err(e) => json!({"name": c.name, "error": e.to_string()}),
        })
        .collect::<Vec<_>>();
    let leaves = input.density_cases.iter().map(|c| {
        Ok(json!({"name": c.name, "values": leaf(c)?.into_iter().map(bits).collect::<Vec<_>>()}))
    }).collect::<Result<Vec<Value>, Box<dyn Error>>>()?;
    Ok(json!({"cases": cases, "leaves": leaves}))
}
