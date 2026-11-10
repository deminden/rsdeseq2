use rsdeseq2::shrinkage::{FixedMixture, TruncatedNormalMixture, UniformMixture};
use serde::Deserialize;
use serde_json::{Value, json};
use std::error::Error;
#[derive(Deserialize)]
pub struct Input {
    #[serde(alias = "truncated_cases")]
    uniform_cases: Vec<Case>,
}
#[derive(Deserialize)]
struct Case {
    name: String,
    weights: Vec<String>,
    means: Option<Vec<String>>,
    sds: Option<Vec<String>>,
    lower: Vec<String>,
    upper: Vec<String>,
    estimates: Vec<String>,
    standard_errors: Vec<String>,
    thresholds: Vec<String>,
}
pub fn evaluate(input: Input) -> Result<Value, Box<dyn Error>> {
    let mut cases = vec![];
    for c in input.uniform_cases {
        let weights = super::parse_numeric_inputs(&c.weights)?;
        let lower = super::parse_numeric_inputs(&c.lower)?;
        let upper = super::parse_numeric_inputs(&c.upper)?;
        let prior = match (c.means, c.sds) {
            (Some(means), Some(sds)) => FixedMixture::TruncatedNormal(TruncatedNormalMixture::new(
                weights,
                super::parse_numeric_inputs(&means)?,
                super::parse_numeric_inputs(&sds)?,
                lower,
                upper,
            )?),
            (None, None) => FixedMixture::Uniform(UniformMixture::new(weights, lower, upper)?),
            _ => return Err("truncated mixture requires both means and SDs".into()),
        };
        let estimates = super::parse_numeric_inputs(&c.estimates)?;
        let errors = super::parse_numeric_inputs(&c.standard_errors)?;
        if estimates.len() != errors.len() {
            return Err("uniform input lengths differ".into());
        }
        let thresholds = super::parse_numeric_inputs(&c.thresholds)?;
        let mut likelihoods = vec![];
        let mut means = vec![];
        let mut sds = vec![];
        let mut weights = vec![];
        let mut cdf = vec![];
        let posteriors = prior.posteriors(&estimates, &errors)?;
        for ((&estimate, &se), p) in estimates.iter().zip(&errors).zip(posteriors) {
            likelihoods.extend(
                prior
                    .log_likelihoods(estimate, se)?
                    .into_iter()
                    .map(super::bits),
            );
            means.push(super::bits(p.mean()));
            sds.push(super::bits(p.standard_deviation()));
            weights.extend(p.weights().iter().copied().map(super::bits));
            cdf.extend(thresholds.iter().map(|c| super::bits(p.cdf(*c))));
        }
        cases.push(json!({"name":c.name,"log_likelihoods":likelihoods,"posterior_mean":means,"posterior_sd":sds,"posterior_weights":weights,"cdf":cdf}));
    }
    Ok(json!({"cases":cases}))
}
