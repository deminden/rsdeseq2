//! Negative-binomial apeglm argument transport; all arithmetic is in core.
use super::ashr::{choice, field, flag, number};
use super::posterior::{numbers, scalar};
use super::{Data, Value, fields, get, named, text};
use rsdeseq2::{
    matrix::{MatrixOrder, MatrixView},
    shrinkage::{
        CauchyCoefficientPrior, CoefficientMle, CoefficientScale, NegativeBinomialPosteriorInput,
        NegativeBinomialPosteriorOptions, PosteriorOffsets, RegressionLikelihood, RegressionMethod,
        RegressionOptimizer, RegressionPosteriorInput, RegressionPosteriorOptions,
        RestrictedRegressionLikelihood, fit_negative_binomial_posterior,
        fit_negative_binomial_posterior_with_rng, fit_regression_posterior,
        fit_regression_posterior_with_rng,
    },
};
fn matrix(v: &Value) -> Result<MatrixView<'_, f64>, String> {
    let [r, c] = v.dimensions.ok_or("apeglm requires a numeric matrix")?;
    MatrixView::new(r, c, numbers(v)?, MatrixOrder::RowMajor).map_err(|e| e.to_string())
}
fn present(v: Option<&Value>) -> Option<&Value> {
    v.filter(|v| !matches!(v.data, Data::Null))
}
fn coefficient(v: &Value, p: usize) -> Result<usize, String> {
    let v = scalar(v)?;
    if !v.is_finite() || v.fract() != 0.0 || v < 1.0 || v > p as f64 {
        return Err("coef must be a one-based coefficient index".into());
    }
    Ok(v as usize - 1)
}
fn indices(v: &Value, p: usize) -> Result<Vec<usize>, String> {
    numbers(v)?
        .iter()
        .map(|&x| coefficient(&Value::numbers(vec![x]), p))
        .collect()
}
fn optional(v: Option<Vec<Option<f64>>>) -> Value {
    v.map(|v| {
        Value::numbers(
            v.into_iter()
                .map(|v| v.unwrap_or(f64::from_bits(0x7ff00000000007a2)))
                .collect(),
        )
    })
    .unwrap_or_else(|| Value::new(Data::Null))
}
pub(super) fn evaluate(request: &Value) -> Result<Value, String> {
    let args = get(request, "arguments")?;
    let (names, _) = fields(args)?;
    let allowed = [
        "Y",
        "x",
        "log.lik",
        "param",
        "coef",
        "mle",
        "no.shrink",
        "interval.type",
        "interval.level",
        "threshold",
        "contrasts",
        "weights",
        "offset",
        "flip.sign",
        "prior.control",
        "multiplier",
        "ngrid",
        "nsd",
        "ngrid.nuis",
        "nsd.nuis",
        "log.link",
        "param.sd",
        "method",
        "optim.method",
        "bounds",
    ];
    for name in names {
        if !allowed.contains(&name.as_str()) {
            return Err(format!("unused apeglm argument '{name}'"));
        }
    }
    let method = choice(
        field(args, "method"),
        &[
            "general",
            "nbinomR",
            "nbinomCR",
            "nbinomC",
            "nbinomC*",
            "betabinR",
            "betabinCR",
            "betabinC",
            "betabinC*",
        ],
    )?;
    if choice(
        field(args, "interval.type"),
        &["laplace", "HPD", "credible"],
    )? != "laplace"
    {
        return Err("native apeglm grid intervals are not implemented yet".into());
    }
    let general = matches!(
        method.as_str(),
        "general" | "nbinomR" | "betabinR" | "betabinCR" | "betabinC" | "betabinC*"
    );
    let optimizer = choice(field(args, "optim.method"), &["BFGS", "L-BFGS-B"])?;
    let log_link = field(args, "log.link")
        .map(flag)
        .transpose()?
        .unwrap_or(true);
    if !general && (optimizer != "BFGS" || !log_link) {
        return Err(
            "native apeglm alternate optimizers/log-link policy are not implemented yet".into(),
        );
    }
    if present(field(args, "contrasts")).is_some() || present(field(args, "param.sd")).is_some() {
        return Err(
            "native apeglm contrasts and parameter uncertainty are not implemented yet".into(),
        );
    }
    let counts = matrix(get(args, "Y")?)?;
    let design = matrix(get(args, "x")?)?;
    let p = design.n_cols();
    let target = present(field(args, "coef"))
        .map(|v| coefficient(v, p))
        .transpose()?;
    if method.starts_with("nbinom") && target == Some(0) {
        return Err("nbinom apeglm requires a coefficient after the intercept".into());
    }
    let mut prior = CauchyCoefficientPrior::default();
    let supplied_prior = present(field(args, "prior.control"));
    let mut prior_df = 1.0;
    let mut prior_ncp = 0.0;
    let mut nuisance_mean = 0.0;
    if let Some(pc) = supplied_prior {
        prior_ncp = scalar(get(pc, "prior.mean")?)?;
        prior_df = scalar(get(pc, "prior.df")?)?;
        nuisance_mean = scalar(get(pc, "prior.no.shrink.mean")?)?;
        if method != "general" && (prior_ncp != 0.0 || prior_df != 1.0 || nuisance_mean != 0.0) {
            return Err("nbinom apeglm requires zero-centered Cauchy priors".into());
        }
        prior.nuisance = indices(get(pc, "no.shrink")?, p)?;
        prior.scale = scalar(get(pc, "prior.scale")?)?;
        prior.nuisance_scale = scalar(get(pc, "prior.no.shrink.scale")?)?;
        prior.adapted_variance = present(field(pc, "prior.var")).map(scalar).transpose()?;
    }
    let all_nuisance = field(args, "no.shrink")
        .map(flag)
        .transpose()?
        .unwrap_or(false);
    if all_nuisance {
        prior.nuisance = (0..p).collect();
    }
    let weights = present(field(args, "weights")).map(matrix).transpose()?;
    let offsets = if let Some(offset) = present(field(args, "offset")) {
        if let Some(v) = field(offset, "size_factors") {
            PosteriorOffsets::SizeFactors(numbers(v)?)
        } else if let Some(v) = field(offset, "factors") {
            PosteriorOffsets::Factors(matrix(v)?)
        } else {
            PosteriorOffsets::Log(matrix(offset)?)
        }
    } else {
        PosteriorOffsets::None
    };
    let mut mle_x = vec![];
    let mut mle_s = vec![];
    let mle = if let Some(v) = present(field(args, "mle")) {
        let scale = if let Some(x) = field(v, "estimates") {
            mle_x = numbers(x)?.to_vec();
            mle_s = numbers(get(v, "standard_errors")?)?.to_vec();
            match text(get(v, "scale")?)? {
                "log2" => CoefficientScale::Log2,
                "natural" => CoefficientScale::Natural,
                _ => return Err("unknown MLE scale".into()),
            }
        } else {
            let m = matrix(v)?;
            if m.n_cols() != 2 {
                return Err("mle requires two columns".into());
            }
            for i in 0..m.n_rows() {
                mle_x.push(*m.get(i, 0).map_err(|e| e.to_string())?);
                mle_s.push(*m.get(i, 1).map_err(|e| e.to_string())?);
            }
            CoefficientScale::Natural
        };
        Some(CoefficientMle {
            estimates: &mle_x,
            standard_errors: &mle_s,
            scale,
        })
    } else {
        None
    };
    let threshold = present(field(args, "threshold"))
        .map(|v| {
            if let Some(value) = field(v, "value") {
                let scale = match text(get(v, "scale")?)? {
                    "log2" => CoefficientScale::Log2,
                    "natural" => CoefficientScale::Natural,
                    _ => return Err("unknown threshold scale".into()),
                };
                scalar(value).map(|v| scale.to_natural(v))
            } else {
                scalar(v)
            }
        })
        .transpose()?;
    let posterior = NegativeBinomialPosteriorOptions {
        prior,
        coefficient: target,
        threshold,
        flip_sign: field(args, "flip.sign")
            .map(flag)
            .transpose()?
            .unwrap_or(true),
        multiplier: number(args, "multiplier", 1.0)?,
        map_only: matches!(
            method.as_str(),
            "nbinomC" | "betabinC" | "nbinomC*" | "betabinC*"
        ),
        check_opposite_start: method == "nbinomCR",
        ..Default::default()
    };
    let mut rng = if method.ends_with('*') {
        Some(super::posterior_rng::restore(field(request, "rng_state"))?)
    } else {
        None
    };
    let result = if general {
        let parameter = get(args, "param")?;
        let parameters = if parameter.dimensions.is_some() {
            matrix(parameter)?
        } else {
            MatrixView::new(
                counts.n_rows(),
                1,
                numbers(parameter)?,
                MatrixOrder::RowMajor,
            )
            .map_err(|e| e.to_string())?
        };
        let likelihood = if method.starts_with("betabin") {
            RegressionLikelihood::BetaBinomial
        } else if method == "nbinomR" {
            RegressionLikelihood::NegativeBinomial
        } else {
            let model = get(args, "log.lik")?;
            match text(get(model, "type")?)? {
                "negative_binomial" => RegressionLikelihood::NegativeBinomial,
                "restricted" => {
                    let (a, e, c) =
                        super::location_shrinkage::function_parts(get(model, "model")?)?;
                    RegressionLikelihood::Restricted(std::sync::Arc::new(
                        RestrictedRegressionLikelihood::new(a, e, c).map_err(|e| e.to_string())?,
                    ))
                }
                _ => return Err("unsupported general likelihood descriptor".into()),
            }
        };
        let bounds = if optimizer == "L-BFGS-B" {
            if let Some(v) = present(field(args, "bounds")) {
                let v = numbers(v)?;
                if v.len() != 2 {
                    return Err("bounds requires two values".into());
                }
                (v[0], v[1])
            } else {
                (f64::NEG_INFINITY, f64::INFINITY)
            }
        } else {
            (f64::NEG_INFINITY, f64::INFINITY)
        };
        let input = RegressionPosteriorInput {
            observations: counts,
            design,
            parameters: Some(parameters),
            offsets,
            weights,
            mle,
        };
        let options = RegressionPosteriorOptions {
            posterior,
            likelihood,
            prior_df,
            prior_noncentrality: prior_ncp,
            nuisance_mean,
            optimizer: if optimizer == "BFGS" {
                RegressionOptimizer::Bfgs
            } else {
                RegressionOptimizer::Lbfgsb
            },
            method: match method.as_str() {
                "nbinomR" => RegressionMethod::NegativeBinomialR,
                "betabinR" => RegressionMethod::BetaBinomialR,
                "betabinCR" => RegressionMethod::BetaBinomialCR,
                "betabinC" => RegressionMethod::BetaBinomialC,
                "betabinC*" => RegressionMethod::BetaBinomialCStar,
                _ => RegressionMethod::General,
            },
            log_link,
            bounds,
        };
        if let Some(rng) = rng.as_mut() {
            fit_regression_posterior_with_rng(input, &options, rng)
        } else {
            fit_regression_posterior(input, &options)
        }
    } else {
        let input = NegativeBinomialPosteriorInput {
            counts,
            design,
            dispersions: numbers(get(args, "param")?)?,
            offsets,
            weights,
            mle,
        };
        if let Some(rng) = rng.as_mut() {
            fit_negative_binomial_posterior_with_rng(input, &posterior, rng)
        } else {
            fit_negative_binomial_posterior(input, &posterior)
        }
    };
    let result = match result {
        Ok(result) => result,
        Err(error) => {
            if let Some(rng) = rng.as_ref() {
                let state = rng.has_drawn().then(|| rng.r_seed());
                return Ok(named(
                    &["error", "random_state"],
                    vec![
                        Value::new(Data::Text(vec![Some(error.to_string())])),
                        super::posterior_rng::encode(state.as_deref()),
                    ],
                ));
            }
            return Err(error.to_string());
        }
    };
    let target = target.ok_or("wrapper apeglm request requires a coefficient")?;
    let summary = result
        .coefficient_summary(target, CoefficientScale::Log2)
        .map_err(|e| e.to_string())?;
    // Preserve extra primitive prior metadata; replace numerical fields with
    // the core's final values, including empirical adaptation.
    let mut pc = if let Some(v) = supplied_prior {
        v.clone()
    } else {
        named(
            &[
                "no.shrink",
                "prior.mean",
                "prior.scale",
                "prior.df",
                "prior.no.shrink.mean",
                "prior.no.shrink.scale",
            ],
            vec![
                Value::numbers(vec![1.0]),
                Value::numbers(vec![0.0]),
                Value::numbers(vec![1.0]),
                Value::numbers(vec![1.0]),
                Value::numbers(vec![0.0]),
                Value::numbers(vec![15.0]),
            ],
        )
    };
    let mut replace = |name: &str, value: Value| -> Result<(), String> {
        let Data::List(values) = &mut pc.data else {
            return Err("prior.control requires a primitive list".into());
        };
        let names = pc.names.as_mut().ok_or("prior.control requires names")?;
        if let Some(i) = names.iter().position(|n| n == name) {
            values[i] = value;
        } else {
            names.push(name.into());
            values.push(value);
        }
        Ok(())
    };
    if mle.is_some() || all_nuisance {
        replace("no.shrink", super::indices(result.prior.nuisance.clone()))?;
    }
    if mle.is_some() {
        replace("prior.scale", Value::numbers(vec![result.prior.scale]))?;
        if let Some(v) = result.prior.adapted_variance {
            replace("prior.var", Value::numbers(vec![v]))?;
        }
    }
    let convergence = result
        .diagnostics
        .iter()
        .map(|d| d.convergence.map(|v| v as f64))
        .collect::<Vec<_>>();
    Ok(named(
        &[
            "posterior_mode",
            "posterior_sd",
            "svalue",
            "threshold_svalue",
            "prior_control",
            "convergence",
            "warnings",
            "random_state",
        ],
        vec![
            Value::numbers(summary.mode),
            optional(Some(summary.standard_deviation)),
            optional(summary.sign_svalues),
            optional(summary.threshold_svalues),
            pc,
            optional(Some(convergence)),
            Value::new(Data::Text(result.warnings.into_iter().map(Some).collect())),
            super::posterior_rng::encode(result.random_state.as_deref()),
        ],
    ))
}
