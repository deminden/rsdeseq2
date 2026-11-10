use rsdeseq2::shrinkage::{
    MixtureComponents, MixtureFamily, MixtureGridOptions, MixturePriorPenalty, prepare_mixture_grid,
};
use serde::Deserialize;
#[derive(Deserialize)]
struct Reference {
    cases: Vec<Case>,
}
#[derive(Deserialize)]
struct Case {
    name: String,
    family: String,
    estimates: Vec<String>,
    standard_errors: Vec<String>,
    scales: Option<Vec<String>>,
    point_mass: bool,
    mode: String,
    range: Vec<String>,
    nullweight: Option<String>,
    component_class: String,
    initial_weights: Vec<String>,
    dirichlet: Vec<String>,
    means: Option<Vec<String>>,
    sds: Option<Vec<String>>,
    lower: Option<Vec<String>>,
    upper: Option<Vec<String>>,
    warnings: Vec<String>,
}
fn number(v: &str) -> f64 {
    if let Some(bits) = v.strip_prefix("bits:") {
        return f64::from_bits(u64::from_str_radix(bits, 16).unwrap());
    }
    match v {
        "NA" | "NaN" => f64::NAN,
        "Inf" => f64::INFINITY,
        "-Inf" => f64::NEG_INFINITY,
        _ => v.parse().unwrap(),
    }
}
fn numbers(v: &[String]) -> Vec<f64> {
    v.iter().map(|v| number(v)).collect()
}
fn close(actual: &[f64], expected: &[String], context: &str) {
    assert_eq!(actual.len(), expected.len(), "{context}: lengths");
    for (i, (&a, b)) in actual.iter().zip(expected).enumerate() {
        let b = number(b);
        assert!(
            a == b || (a.is_finite() && b.is_finite() && (a - b).abs() <= 2e-14 + 2e-14 * b.abs()),
            "{context}[{i}]: {a:.17e} != {b:.17e}"
        );
    }
}
#[test]
fn all_six_normal_likelihood_grids_match_original_preparation() {
    let r: Reference = serde_json::from_str(include_str!("data/mixture_grid_r462.json")).unwrap();
    assert_eq!(r.cases.len(), 24);
    for c in r.cases {
        let o = MixtureGridOptions {
            family: c.family.parse::<MixtureFamily>().unwrap(),
            scales: c.scales.as_ref().map(|v| numbers(v)),
            point_mass: c.point_mass,
            mode: number(&c.mode),
            range: (number(&c.range[0]), number(&c.range[1])),
            penalty: c
                .nullweight
                .as_ref()
                .map(|v| MixturePriorPenalty::NullBiased(number(v)))
                .unwrap_or(MixturePriorPenalty::Uniform),
            ..Default::default()
        };
        let g =
            prepare_mixture_grid(&numbers(&c.estimates), &numbers(&c.standard_errors), &o).unwrap();
        close(
            &g.initial_weights,
            &c.initial_weights,
            &format!("{} weights", c.name),
        );
        close(&g.dirichlet, &c.dirichlet, &format!("{} penalty", c.name));
        assert_eq!(g.warnings, c.warnings, "{} warnings", c.name);
        match g.components {
            MixtureComponents::Normal { means, sds } => {
                assert_eq!(c.component_class, "normalmix");
                close(&means, c.means.as_ref().unwrap(), &c.name);
                close(&sds, c.sds.as_ref().unwrap(), &c.name);
            }
            MixtureComponents::Uniform { lower, upper } => {
                assert_eq!(c.component_class, "unimix");
                close(&lower, c.lower.as_ref().unwrap(), &c.name);
                close(&upper, c.upper.as_ref().unwrap(), &c.name);
            }
            MixtureComponents::TruncatedNormal {
                means,
                sds,
                lower,
                upper,
            } => {
                assert_eq!(c.component_class, "tnormalmix");
                close(&means, c.means.as_ref().unwrap(), &c.name);
                close(&sds, c.sds.as_ref().unwrap(), &c.name);
                close(&lower, c.lower.as_ref().unwrap(), &c.name);
                close(&upper, c.upper.as_ref().unwrap(), &c.name);
            }
        }
    }
}
