//! Reference breadth (hardening 6.6): `DESeq()` + `results()` defaults on
//! three simulated shapes the GTEx panels do not cover, from R 4.6.2 /
//! DESeq2 1.52.0 (scripts/references/generate_breadth_references.R):
//!
//! - `three_vs_three`: 1,999 genes x 6 samples, the smallest common design.
//! - `sparse`: 983 genes x 48 cells, 40 % extra zeros, `sfType = "poscounts"`.
//! - `many_samples`: 400 genes x 500 samples, `~ batch + condition`.
//!
//! The gates are the end-to-end gates of `gtex_end_to_end`: identical
//! missingness, no significance flips, tight bounds for every gene whose
//! dispersion matches R, and at most one MAP knife-edge gene.

use std::path::PathBuf;

use rsdeseq2::io::{read_count_matrix_tsv, read_labeled_design_matrix_tsv};
use rsdeseq2::prelude::*;

mod common;

const DISPERSION_MATCH_LOG_GAP: f64 = 1e-10;
const KNIFE_EDGE_LOG_GAP: f64 = 1e-4;

fn data(case: &str, name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/breadth")
        .join(case)
        .join(name)
}

fn optional(text: &str) -> Option<f64> {
    (text != "NA").then(|| text.parse().unwrap())
}

fn max_gap(
    actual: &[Option<f64>],
    expected: &[Option<f64>],
    keep: &dyn Fn(usize) -> bool,
    log_scale: bool,
) -> (f64, usize) {
    let mut gap = 0.0f64;
    let mut missing = 0;
    for (i, (a, e)) in actual.iter().zip(expected).enumerate() {
        if !keep(i) {
            continue;
        }
        match (a, e) {
            (Some(a), Some(e)) if log_scale && *a > 0.0 && *e > 0.0 => {
                gap = gap.max((a.log10() - e.log10()).abs());
            }
            (Some(a), Some(e)) => gap = gap.max((a - e).abs() / e.abs().max(1.0)),
            (None, None) => {}
            _ => missing += 1,
        }
    }
    (gap, missing)
}

fn check_case(case: &str, builder: DeseqBuilder, matched_bound: f64) {
    let counts = read_count_matrix_tsv(data(case, "counts.tsv")).unwrap();
    let design = read_labeled_design_matrix_tsv(data(case, "design_matrix.tsv"))
        .unwrap()
        .design;
    let cutoff = default_cooks_cutoff(design.n_samples(), design.n_coefficients())
        .unwrap()
        .unwrap();
    let options = CooksReplacementOptions::new(cutoff);
    let (fit, results) = builder
        .run(&Analysis::new(&counts, &design).cooks_replacement(&options))
        .unwrap()
        .into_fit_and_results()
        .unwrap();
    let results = results.unwrap();

    let reference = common::read_tsv(&data(case, "reference_results.tsv"));
    assert_eq!(reference.len(), counts.n_genes(), "{case}: genes");
    let column = |name: &str| -> Vec<Option<f64>> {
        reference.iter().map(|row| optional(&row[name])).collect()
    };
    let reference_size_factors = common::read_tsv(&data(case, "reference_size_factors.tsv"));
    for (actual, row) in fit.size_factors.iter().zip(&reference_size_factors) {
        let expected: f64 = row["size_factor"].parse().unwrap();
        assert!(
            (actual - expected).abs() <= 1e-12 * expected,
            "{case}: size factor {actual} vs {expected}"
        );
    }

    let dispersion = fit
        .dispersion
        .as_ref()
        .unwrap()
        .iter()
        .map(|v| Some(*v).filter(|v| v.is_finite()))
        .collect::<Vec<_>>();
    let reference_dispersion = column("dispersion");
    let log_gap = |i: usize| match (dispersion[i], reference_dispersion[i]) {
        (Some(a), Some(e)) => (a.ln() - e.ln()).abs(),
        (None, None) => 0.0,
        _ => f64::INFINITY,
    };
    let knife_edge = (0..dispersion.len())
        .filter(|&i| log_gap(i) > DISPERSION_MATCH_LOG_GAP)
        .collect::<Vec<_>>();
    for &i in &knife_edge {
        eprintln!(
            "{case}: knife-edge dispersion gene {} log gap {:.3e}",
            reference[i]["gene"],
            log_gap(i)
        );
    }
    assert!(
        knife_edge.len() <= 1,
        "{case}: {} genes left the dispersion match",
        knife_edge.len()
    );
    for &i in &knife_edge {
        assert!(
            log_gap(i) <= KNIFE_EDGE_LOG_GAP,
            "{case}: knife-edge gap {}",
            log_gap(i)
        );
    }
    let strict = |i: usize| !knife_edge.contains(&i);
    let edge = |i: usize| knife_edge.contains(&i);

    let fields: [(&str, Vec<Option<f64>>, bool); 6] = [
        (
            "baseMean",
            results.rows.iter().map(|r| Some(r.base_mean)).collect(),
            false,
        ),
        (
            "log2FoldChange",
            results.rows.iter().map(|r| r.log2_fold_change).collect(),
            false,
        ),
        (
            "lfcSE",
            results.rows.iter().map(|r| r.lfc_se).collect(),
            false,
        ),
        ("stat", results.rows.iter().map(|r| r.stat).collect(), false),
        (
            "pvalue",
            results.rows.iter().map(|r| r.pvalue).collect(),
            true,
        ),
        ("padj", results.rows.iter().map(|r| r.padj).collect(), true),
    ];
    for (name, actual, log_scale) in &fields {
        let expected = column(name);
        let (tight, tight_missing) = max_gap(actual, &expected, &strict, *log_scale);
        let (loose, loose_missing) = max_gap(actual, &expected, &edge, *log_scale);
        eprintln!("{case} {name}: matched {tight:.3e}, knife-edge {loose:.3e}");
        assert_eq!(
            tight_missing + loose_missing,
            0,
            "{case} {name}: missingness differs"
        );
        assert!(
            tight <= matched_bound,
            "{case} {name}: {tight:e} > {matched_bound:e}"
        );
        assert!(loose <= 1e-4, "{case} {name}: knife-edge {loose:e} > 1e-4");
    }
    let padj = &fields[5].1;
    for alpha in [0.05, 0.1] {
        assert_eq!(
            common::significance_flips(padj, &column("padj"), alpha),
            0,
            "{case}: calls flipped at {alpha}"
        );
    }
}

#[test]
fn three_vs_three_matches_deseq2_defaults() {
    check_case("three_vs_three", DeseqBuilder::new(), 1e-8);
}

#[test]
fn sparse_poscounts_matches_deseq2_defaults() {
    check_case(
        "sparse",
        DeseqBuilder::new().size_factor_method(SizeFactorMethod::PosCounts),
        1e-8,
    );
}

#[test]
fn five_hundred_samples_match_deseq2_defaults() {
    check_case("many_samples", DeseqBuilder::new(), 1e-8);
}
