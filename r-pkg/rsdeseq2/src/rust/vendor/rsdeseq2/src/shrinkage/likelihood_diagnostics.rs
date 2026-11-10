//! Preserve optional-method diagnostics without repeating statistical work.
//! Cached warning sequences follow the original full-summary operation graph.
use super::{MixtureComponents, RestrictedLikelihood};
#[derive(Default)]
pub(super) struct Diagnostics {
    pub likelihood: Vec<String>,
    pub mean: Vec<String>,
    pub second: Vec<String>,
    pub cdf: Vec<String>,
    pub null: Vec<String>,
}
fn coordinates(x: &[f64], s: &[f64], ends: &[f64], reverse: bool) -> Vec<f64> {
    ends.iter()
        .flat_map(|&v| {
            x.iter()
                .zip(s)
                .map(move |(&x, &s)| if reverse { (v - x) / s } else { (x - v) / s })
        })
        .collect()
}
pub(super) fn likelihood(
    f: &RestrictedLikelihood,
    components: &MixtureComponents,
    x: &[f64],
    s: &[f64],
) -> Vec<String> {
    let MixtureComponents::Uniform { lower, upper } = components else {
        return vec![];
    };
    let mut warnings = f.cdf_warnings(&coordinates(x, s, lower, false));
    warnings.extend(f.cdf_warnings(&coordinates(x, s, upper, false)));
    warnings.extend(f.density_warnings(&coordinates(x, s, upper, false)));
    warnings
}
pub(super) fn cdf(
    f: &RestrictedLikelihood,
    components: &MixtureComponents,
    x: &[f64],
    s: &[f64],
    cutoff: f64,
) -> Vec<String> {
    let MixtureComponents::Uniform { lower, upper } = components else {
        return vec![];
    };
    let keep = lower
        .iter()
        .zip(upper)
        .enumerate()
        .filter_map(|(i, (&a, &b))| (a <= cutoff && cutoff < b).then_some(i))
        .collect::<Vec<_>>();
    if keep.is_empty() {
        return vec![];
    }
    let a = keep.iter().map(|&i| lower[i]).collect::<Vec<_>>();
    let b = keep.iter().map(|&i| upper[i]).collect::<Vec<_>>();
    let mut warnings = f.cdf_warnings(&coordinates(x, s, &a, false));
    warnings.extend(f.cdf_warnings(&coordinates(x, s, &vec![cutoff; keep.len()], false)));
    warnings.extend(f.cdf_warnings(&coordinates(x, s, &b, false)));
    warnings
}
impl Diagnostics {
    pub(super) fn new(
        f: &RestrictedLikelihood,
        components: &MixtureComponents,
        x: &[f64],
        s: &[f64],
    ) -> Self {
        let mut d = Self {
            likelihood: likelihood(f, components, x, s),
            cdf: cdf(f, components, x, s, 0.0),
            null: f.density_warnings(&x.iter().zip(s).map(|(x, s)| x / s).collect::<Vec<_>>()),
            ..Default::default()
        };
        if let MixtureComponents::Uniform { lower, upper } = components {
            d.mean = f.moment_warnings(
                &coordinates(x, s, upper, false),
                &coordinates(x, s, lower, false),
                1,
            );
            let a = coordinates(x, s, lower, true);
            let b = coordinates(x, s, upper, true);
            d.second = f.moment_warnings(&a, &b, 1);
            d.second.extend(f.moment_warnings(&a, &b, 2));
        }
        d
    }
    pub(super) fn append_full_summary(
        &self,
        out: &mut Vec<String>,
        svalue: &[String],
        qvalue: &[String],
    ) {
        let mut probability = self.likelihood.clone();
        probability.extend_from_slice(&self.cdf);
        let mut negative = probability;
        negative.extend_from_slice(&self.likelihood);
        let mut sign = negative.clone();
        sign.extend_from_slice(&self.likelihood);
        // NegativeProb, PositiveProb, lfsr, svalue, lfdr, qvalue, PM, PSD.
        out.extend_from_slice(&negative);
        out.extend_from_slice(&sign);
        out.extend_from_slice(&sign);
        out.extend_from_slice(&sign);
        out.extend_from_slice(svalue);
        out.extend_from_slice(&self.likelihood);
        out.extend_from_slice(&self.likelihood);
        out.extend_from_slice(qvalue);
        out.extend_from_slice(&self.likelihood);
        out.extend_from_slice(&self.mean);
        out.extend_from_slice(&self.likelihood);
        out.extend_from_slice(&self.second);
        out.extend_from_slice(&self.likelihood);
        out.extend_from_slice(&self.mean);
    }
}
