//! `DESeq(betaPrior = TRUE)` through `run` (hardening 4.3d), against R 4.6.2 /
//! DESeq2 1.52.0 references from
//! scripts/references/generate_beta_prior_references.R.

use std::path::PathBuf;

use rsdeseq2::design::{FormulaFactorColumn, FormulaModelFrame, FormulaNumericColumn};
use rsdeseq2::errors::RunConflict;
use rsdeseq2::io::{read_count_matrix_tsv, read_labeled_design_matrix_tsv};
use rsdeseq2::prelude::*;

use crate::support as common;

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
    for ((name, actual), expected) in names
        .iter()
        .zip(&output.beta_prior_variance)
        .zip(&reference)
    {
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
fn formula_sample_offsets_match_cell_factors_at_early_stage_and_beta_prior() {
    use rsdeseq2::normalization::normalization_factors_with_offsets;

    let case = "three_vs_three";
    let counts = read_count_matrix_tsv(data(case, "counts.tsv")).unwrap();
    let mut frame = model_frame(case);
    let offsets = vec![0.0, 0.1, -0.1, 0.2, -0.2, 0.0];
    frame.numeric_covariates.push(FormulaNumericColumn {
        name: "exposure".into(),
        values: offsets.clone(),
    });
    let base = vec![0.85, 1.0, 1.1, 0.9, 1.0, 1.2];
    let effective = normalization_factors_with_offsets(
        counts.n_genes(),
        counts.n_samples(),
        Some(&base),
        None,
        &offsets,
    )
    .unwrap();
    let by_sample = DeseqBuilder::new()
        .model_frame(frame.clone())
        .size_factors(base);
    let by_cell = DeseqBuilder::new()
        .model_frame(frame)
        .normalization_factors(effective.clone());
    let formula_offsets = Analysis::formula(&counts, "~ condition + offset(exposure)");
    let formula_cells = Analysis::formula(&counts, "~ condition");

    for (sample_request, cell_request) in [
        (
            formula_offsets.stop_after(Stage::SizeFactors),
            formula_cells.stop_after(Stage::SizeFactors),
        ),
        (
            formula_offsets
                .stop_after(Stage::GeneWise)
                .means(MeanModel::LinearMu),
            formula_cells
                .stop_after(Stage::GeneWise)
                .means(MeanModel::LinearMu),
        ),
        (
            formula_offsets.beta_prior(true),
            formula_cells.beta_prior(true),
        ),
    ] {
        let sample = by_sample.run(&sample_request).unwrap();
        let cell = by_cell.run(&cell_request).unwrap();
        let fit = sample.fit();
        assert_eq!(fit.size_factors, vec![1.0; counts.n_samples()]);
        assert_eq!(fit.normalization_factors.as_ref(), Some(&effective));
        assert_eq!(fit.base_mean, cell.fit().base_mean);
        assert_eq!(fit.base_var, cell.fit().base_var);
        let as_bits = |values: Option<&[f64]>| {
            values.map(|row| row.iter().map(|value| value.to_bits()).collect::<Vec<_>>())
        };
        assert_eq!(
            as_bits(fit.disp_gene_est.as_deref()),
            as_bits(cell.fit().disp_gene_est.as_deref()),
        );
        if let (Some(sample), Some(cell)) = (sample.results(), cell.results()) {
            assert_eq!(sample.rows.len(), cell.rows.len());
            for (sample, cell) in sample.rows.iter().zip(&cell.rows) {
                assert_eq!(sample, cell);
            }
            assert_eq!(sample.metadata, cell.metadata);
            assert_eq!(sample.independent_filtering, cell.independent_filtering);
        } else {
            assert!(sample.results().is_none() && cell.results().is_none());
        }
    }
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
    std::assert_matches!(
        error,
        DeseqError::InvalidRun {
            conflict: RunConflict::BetaPriorNeedsFormula
        }
    );
    let error = builder
        .run(
            &Analysis::formula(&counts, "~ condition")
                .lrt_formula("~ 1")
                .beta_prior(true),
        )
        .unwrap_err();
    std::assert_matches!(
        error,
        DeseqError::InvalidRun {
            conflict: RunConflict::BetaPriorNeedsWald
        }
    );
    let error = builder
        .run(
            &Analysis::formula(&counts, "~ condition")
                .beta_prior(true)
                .result(ResultColumn::Coefficient(usize::MAX)),
        )
        .unwrap_err();
    std::assert_matches!(
        error,
        DeseqError::InvalidDimensions { context, expected: 1, actual: usize::MAX }
            if context == "beta-prior coefficient index"
    );
}
