//! `DESeq(betaPrior = TRUE)` through `run` (hardening 4.3d), against R 4.6.2 /
//! DESeq2 1.52.0 references from
//! scripts/references/generate_beta_prior_references.R.

use std::path::PathBuf;

use rsdeseq2::design::{FormulaFactorColumn, FormulaModelFrame};
use rsdeseq2::io::{read_count_matrix_tsv, read_labeled_design_matrix_tsv};
use rsdeseq2::prelude::*;

mod common;

fn data(case: &str, name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/breadth")
        .join(case)
        .join(name)
}

fn optional(text: &str) -> Option<f64> {
    (text != "NA").then(|| text.parse().unwrap())
}

/// The sample table the R script rebuilds from the design matrix.
fn model_frame(case: &str) -> FormulaModelFrame {
    let labelled = read_labeled_design_matrix_tsv(data(case, "design_matrix.tsv")).unwrap();
    let design = labelled.design;
    let names = design.coefficient_names().unwrap().to_vec();
    let column = |name: &str| names.iter().position(|n| n == name);
    let condition = column("conditionB").unwrap();
    let batches: Vec<(usize, String)> = names
        .iter()
        .enumerate()
        .filter_map(|(i, n)| n.strip_prefix("batch").map(|level| (i, level.to_string())))
        .collect();
    let mut factors = Vec::new();
    if !batches.is_empty() {
        let sample_levels = (0..design.n_samples())
            .map(|s| {
                batches
                    .iter()
                    .find(|(i, _)| design.matrix().get(s, *i).copied() == Some(1.0))
                    .map_or("b1".to_string(), |(_, level)| level.clone())
            })
            .collect();
        let mut levels = vec!["b1".to_string()];
        levels.extend(batches.iter().map(|(_, level)| level.clone()));
        factors.push(FormulaFactorColumn {
            name: "batch".into(),
            sample_levels,
            levels: Some(levels),
            reference: None,
        });
    }
    factors.push(FormulaFactorColumn {
        name: "condition".into(),
        sample_levels: (0..design.n_samples())
            .map(|s| {
                if design.matrix().get(s, condition).copied() == Some(1.0) {
                    "B"
                } else {
                    "A"
                }
                .to_string()
            })
            .collect(),
        levels: Some(vec!["A".into(), "B".into()]),
        reference: None,
    });
    FormulaModelFrame {
        factors,
        numeric_covariates: Vec::new(),
    }
}

/// Largest relative gap per field, rows whose missingness differs, and
/// significance flips at padj < 0.1.
fn compare(case: &str, results: &DeseqResults) -> (f64, usize, usize) {
    let reference = common::read_tsv(&data(case, "reference_results_beta_prior.tsv"));
    assert_eq!(reference.len(), results.rows.len(), "{case}: genes");
    let mut gap = 0.0f64;
    let (mut missing, mut flips) = (0, 0);
    for (row, expected) in results.rows.iter().zip(&reference) {
        for (actual, name) in [
            (row.log2_fold_change, "log2FoldChange"),
            (row.lfc_se, "lfcSE"),
            (row.stat, "stat"),
            (row.pvalue, "pvalue"),
            (row.padj, "padj"),
            (row.dispersion, "dispersion"),
        ] {
            match (actual.filter(|v| v.is_finite()), optional(&expected[name])) {
                (Some(a), Some(e)) => gap = gap.max((a - e).abs() / e.abs().max(1e-12)),
                (None, None) => {}
                _ => missing += 1,
            }
        }
        let significant = |p: Option<f64>| p.is_some_and(|p| p < 0.1);
        if significant(row.padj) != significant(optional(&expected["padj"])) {
            flips += 1;
        }
    }
    eprintln!("{case}: max relative gap {gap:.3e}, missing {missing}, flips {flips}");
    (gap, missing, flips)
}

fn check(case: &str, formula: &str) {
    let counts = read_count_matrix_tsv(data(case, "counts.tsv")).unwrap();
    let builder = DeseqBuilder::new().model_frame(model_frame(case));
    let output = builder
        .run(&Analysis::formula(&counts, formula).beta_prior(true))
        .unwrap();
    let RunOutput::BetaPrior(output) = output else {
        panic!("{case}: expected a beta-prior output");
    };
    // betaPriorVar, in expanded-column order.
    let reference = common::read_tsv(&data(case, "reference_beta_prior_variance.tsv"));
    let names = output.expanded_design.coefficient_names().unwrap();
    assert_eq!(names.len(), reference.len(), "{case}: expanded columns");
    for ((name, actual), expected) in names.iter().zip(&output.beta_prior_variance).zip(&reference) {
        assert_eq!(*name, expected["coefficient"], "{case}: column order");
        let expected: f64 = expected["variance"].parse().unwrap();
        assert!(
            (actual - expected).abs() <= 1e-12 * expected,
            "{case}: prior variance {name}: {actual} vs {expected}"
        );
    }
    let (gap, missing, flips) = compare(case, &output.results);
    assert_eq!(missing, 0, "{case}: missingness differs from R");
    assert_eq!(flips, 0, "{case}: significance flips");
    assert!(gap < 1e-9, "{case}: max relative gap {gap:e}");
}

#[test]
fn three_vs_three_matches_deseq_beta_prior() {
    check("three_vs_three", "~ condition");
}

#[test]
fn batch_and_condition_match_deseq_beta_prior() {
    check("many_samples", "~ batch + condition");
}

#[test]
fn beta_prior_needs_a_formula_wald_run() {
    let counts = read_count_matrix_tsv(data("three_vs_three", "counts.tsv")).unwrap();
    let design = read_labeled_design_matrix_tsv(data("three_vs_three", "design_matrix.tsv"))
        .unwrap()
        .design;
    let builder = DeseqBuilder::new().model_frame(model_frame("three_vs_three"));
    let error = builder
        .run(&Analysis::new(&counts, &design).beta_prior(true))
        .unwrap_err();
    assert!(error.to_string().contains("needs a formula design"), "{error}");
    let error = builder
        .run(
            &Analysis::formula(&counts, "~ condition")
                .lrt_formula("~ 1")
                .beta_prior(true),
        )
        .unwrap_err();
    assert!(error.to_string().contains("Wald test"), "{error}");
}
