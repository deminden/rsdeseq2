//! Property tests for invariances that every DESeq2-compatible stage must keep.
//!
//! Inputs are small random count matrices, so each case runs in milliseconds;
//! the case count is kept low because the fits are deterministic per input.

use crate::support::two_group_design;
use proptest::prelude::*;
use proptest::test_runner::FileFailurePersistence;
use rsdeseq2::prelude::*;

const SAMPLES: usize = 6;

fn config() -> ProptestConfig {
    ProptestConfig {
        cases: 48,
        // The integration binary's root is not lib.rs or main.rs, so name the
        // regression file directly (paths are relative to the crate).
        failure_persistence: Some(Box::new(FileFailurePersistence::Direct(
            "tests/properties.proptest-regressions",
        ))),
        ..ProptestConfig::default()
    }
}

/// Genes x SAMPLES counts with every gene positive in every sample, so ratio
/// size factors are always defined.
fn positive_counts(max_genes: usize) -> impl Strategy<Value = (usize, Vec<u32>)> {
    (4..=max_genes).prop_flat_map(|genes| {
        prop::collection::vec(1u32..5_000, genes * SAMPLES).prop_map(move |values| (genes, values))
    })
}

fn matrix(genes: usize, values: Vec<u32>) -> CountMatrix {
    CountMatrix::from_row_major_u32(genes, SAMPLES, values).unwrap()
}

fn close(a: f64, b: f64, relative: f64) -> bool {
    (a - b).abs() <= relative * a.abs().max(b.abs()).max(f64::MIN_POSITIVE)
}

proptest! {
    #![proptest_config(config())]

    #[test]
    fn size_factors_ignore_gene_order((genes, values) in positive_counts(30), rotation in 0usize..30) {
        let counts = matrix(genes, values.clone());
        let shift = rotation % genes;
        let mut rotated = values[shift * SAMPLES..].to_vec();
        rotated.extend_from_slice(&values[..shift * SAMPLES]);
        let original = estimate_size_factors(&counts, SizeFactorMethod::Ratio).unwrap();
        let permuted = estimate_size_factors(&matrix(genes, rotated), SizeFactorMethod::Ratio).unwrap();
        for (a, b) in original.iter().zip(&permuted) {
            prop_assert!(close(*a, *b, 1e-12), "{a} vs {b}");
        }
    }

    #[test]
    fn size_factors_follow_sample_order((genes, values) in positive_counts(30)) {
        let counts = matrix(genes, values.clone());
        let reversed: Vec<u32> = values
            .chunks(SAMPLES)
            .flat_map(|row| row.iter().rev().copied())
            .collect();
        let original = estimate_size_factors(&counts, SizeFactorMethod::Ratio).unwrap();
        let mirrored = estimate_size_factors(&matrix(genes, reversed), SizeFactorMethod::Ratio).unwrap();
        for (sample, factor) in original.iter().enumerate() {
            prop_assert!(close(*factor, mirrored[SAMPLES - 1 - sample], 1e-12));
        }
    }

    #[test]
    fn scaling_a_sample_scales_its_relative_size_factor(
        (genes, values) in positive_counts(30),
        k in 2u32..5,
    ) {
        let counts = matrix(genes, values.clone());
        let scaled: Vec<u32> = values
            .iter()
            .enumerate()
            .map(|(index, &value)| if index % SAMPLES == 0 { value * k } else { value })
            .collect();
        let before = estimate_size_factors(&counts, SizeFactorMethod::Ratio).unwrap();
        let after = estimate_size_factors(&matrix(genes, scaled), SizeFactorMethod::Ratio).unwrap();
        // Geometric means change too, so compare sample 0 relative to sample 1.
        let ratio_before = before[0] / before[1];
        let ratio_after = after[0] / after[1];
        prop_assert!(close(ratio_after, ratio_before * f64::from(k), 1e-9));
    }

    #[test]
    fn bh_adjustment_is_monotone_bounded_and_at_least_p(
        pvalues in prop::collection::vec(prop::option::of(0.0f64..=1.0), 1..60),
    ) {
        let adjusted = bh_adjust(&pvalues);
        prop_assert_eq!(adjusted.len(), pvalues.len());
        let mut pairs = Vec::new();
        for (p, q) in pvalues.iter().zip(&adjusted) {
            match (p, q) {
                (None, None) => {}
                (Some(p), Some(q)) => {
                    prop_assert!(*q <= 1.0 && *q >= *p, "p={p} q={q}");
                    pairs.push((*p, *q));
                }
                _ => prop_assert!(false, "missingness must be preserved"),
            }
        }
        pairs.sort_by(|a, b| a.0.total_cmp(&b.0));
        for window in pairs.windows(2) {
            prop_assert!(window[0].1 <= window[1].1, "not monotone: {window:?}");
        }
    }

    #[test]
    fn negating_a_contrast_negates_effect_and_keeps_pvalue(
        (genes, values) in positive_counts(12),
        dispersion in 0.01f64..1.0,
    ) {
        let counts = matrix(genes, values);
        let design = two_group_design(SAMPLES, None);
        let dispersions = vec![dispersion; genes];
        let builder = DeseqBuilder::new();
        let (_, forward) = builder
            .fit_fixed_dispersion_wald_contrast(&counts, &design, &dispersions, &[0.0, 1.0])
            .unwrap();
        let (_, backward) = builder
            .fit_fixed_dispersion_wald_contrast(&counts, &design, &dispersions, &[0.0, -1.0])
            .unwrap();
        for (a, b) in forward.rows.iter().zip(&backward.rows) {
            prop_assert_eq!(a.log2_fold_change.map(|v| -v), b.log2_fold_change);
            prop_assert_eq!(a.stat.map(|v| -v), b.stat);
            prop_assert_eq!(a.pvalue, b.pvalue);
        }
    }

    #[test]
    fn results_do_not_depend_on_thread_count(
        (genes, values) in positive_counts(12),
        dispersion in 0.01f64..1.0,
    ) {
        let counts = matrix(genes, values);
        let design = two_group_design(SAMPLES, None);
        let dispersions = vec![dispersion; genes];
        let one = DeseqBuilder::new()
            .threads(1)
            .fit_fixed_dispersion_wald(&counts, &design, &dispersions, 1)
            .unwrap();
        let four = DeseqBuilder::new()
            .threads(4)
            .fit_fixed_dispersion_wald(&counts, &design, &dispersions, 1)
            .unwrap();
        // Debug text is exact for every non-NaN f64 and treats NaN equal to
        // itself, unlike PartialEq on the fit's optional NaN fields.
        prop_assert!(format!("{one:?}") == format!("{four:?}"), "1 vs 4 threads differ");
    }

    #[test]
    fn all_zero_genes_are_missing_not_errors(
        (genes, mut values) in positive_counts(12),
        zero_gene in 0usize..12,
        dispersion in 0.01f64..1.0,
    ) {
        let zero_gene = zero_gene % genes;
        values[zero_gene * SAMPLES..(zero_gene + 1) * SAMPLES].fill(0);
        let counts = matrix(genes, values);
        let design = two_group_design(SAMPLES, None);
        let dispersions = vec![dispersion; genes];
        let (_, results) = DeseqBuilder::new()
            .size_factors(vec![1.0; SAMPLES])
            .fit_fixed_dispersion_wald(&counts, &design, &dispersions, 1)
            .unwrap();
        let row = &results.rows[zero_gene];
        prop_assert_eq!(row.base_mean, 0.0);
        prop_assert!(row.log2_fold_change.is_none() && row.pvalue.is_none() && row.padj.is_none());
    }
}
