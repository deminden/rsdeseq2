//! Frozen original-DESeq2 row decisions with raw p-values supplied exactly.
use rsdeseq2::prelude::*;
use serde::Deserialize;

#[derive(Deserialize)]
struct Reference {
    genes: usize,
    samples: usize,
    counts: Vec<u32>,
    cooks: Vec<String>,
    maximum: Vec<String>,
    cutoff: String,
    cases: Vec<Case>,
}
#[derive(Deserialize)]
struct Factor {
    name: String,
    levels: usize,
    codes: Vec<i32>,
}
#[derive(Deserialize)]
struct Case {
    name: String,
    kind: String,
    variables: Vec<String>,
    coefficients: usize,
    design: Vec<String>,
    factors: Vec<Factor>,
    raw_pvalue: Vec<Option<String>>,
    pvalue: Vec<Option<String>>,
    padj: Vec<Option<String>>,
    observed_eligible: bool,
    masked: Vec<bool>,
}
fn numbers(values: Vec<String>) -> Vec<f64> {
    values.into_iter().map(|v| v.parse().unwrap()).collect()
}
fn probabilities(values: Vec<Option<String>>) -> Vec<Option<f64>> {
    values
        .into_iter()
        .map(|v| v.map(|v| f64::from_bits(u64::from_str_radix(&v, 16).unwrap())))
        .collect()
}
fn bits(values: &[Option<f64>]) -> Vec<Option<u64>> {
    values.iter().map(|v| v.map(f64::to_bits)).collect()
}

#[test]
fn original_deseq2_results_masks_and_adjustments_match_all_formula_and_matrix_origins() {
    let r: Reference = serde_json::from_str(include_str!("data/cooks_policy_r462.json")).unwrap();
    let counts = CountView::from_unsigned(
        MatrixView::new(r.genes, r.samples, &r.counts, MatrixOrder::RowMajor).unwrap(),
    );
    let cooks = numbers(r.cooks);
    let heuristic = LowCountHeuristic::new(
        counts,
        MatrixView::new(r.genes, r.samples, &cooks, MatrixOrder::RowMajor).unwrap(),
    )
    .unwrap();
    let maximum = numbers(r.maximum).into_iter().map(Some).collect::<Vec<_>>();
    assert_eq!(r.cases.len(), 12);
    for case in r.cases {
        let origin = if case.kind == "matrix" {
            DesignOrigin::matrix()
        } else {
            DesignOrigin::formula(case.variables).unwrap()
        };
        let values = numbers(case.design);
        let design =
            MatrixView::new(r.samples, case.coefficients, &values, MatrixOrder::RowMajor).unwrap();
        let eligible = case
            .factors
            .iter()
            .any(|f| low_count_eligible(&origin, design, &f.name, &f.codes, f.levels).unwrap());
        assert_eq!(eligible, case.observed_eligible, "{}", case.name);
        let raw = probabilities(case.raw_pvalue);
        let actual = CooksPolicy::new(
            Some(r.cutoff.parse().unwrap()),
            eligible.then_some(heuristic),
        )
        .unwrap()
        .evaluate(raw.iter().copied(), maximum.iter().copied())
        .unwrap();
        let expected = probabilities(case.pvalue);
        assert_eq!(bits(&actual.pvalue), bits(&expected), "{}", case.name);
        assert_eq!(
            actual
                .pvalue
                .iter()
                .zip(&raw)
                .map(|(p, raw)| p.is_none() && raw.is_some())
                .collect::<Vec<_>>(),
            case.masked,
            "{}",
            case.name
        );
        assert_eq!(
            bits(&bh_adjust(&actual.pvalue)),
            bits(&probabilities(case.padj)),
            "{} adjustment",
            case.name
        );
    }
}
