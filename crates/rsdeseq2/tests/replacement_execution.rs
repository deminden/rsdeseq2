//! Execution policy from original DESeq2, shared across fitting/selection routes.
use rsdeseq2::prelude::*;
use serde_json::Value;

fn numbers(value: &Value) -> Vec<f64> {
    value
        .as_array()
        .unwrap()
        .iter()
        .map(|v| {
            f64::from_bits(
                u64::from_str_radix(v.as_str().unwrap().strip_prefix("bits:").unwrap(), 16)
                    .unwrap(),
            )
        })
        .collect()
}
fn matrix_bits(m: &Option<RowMajorMatrix<f64>>) -> Option<Vec<u64>> {
    m.as_ref()
        .map(|m| m.as_slice().iter().map(|v| v.to_bits()).collect())
}

#[test]
fn ineligible_samples_do_not_refit_or_change_original_fits_and_results() {
    let r: Value =
        serde_json::from_str(include_str!("data/replacement_execution_r462.json")).unwrap();
    let genes = r["genes"].as_u64().unwrap() as usize;
    let samples = r["samples"].as_u64().unwrap() as usize;
    let counts = CountMatrix::from_row_major_u32(
        genes,
        samples,
        r["counts"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_u64().unwrap() as u32)
            .collect(),
    )
    .unwrap();
    let design = DesignMatrix::from_row_major(samples, 2, numbers(&r["design"]), None).unwrap();
    let reduced = DesignMatrix::from_row_major(samples, 1, vec![1.0; samples], None).unwrap();
    let dispersions = vec![0.1; genes];
    let contrast = [0.0, 1.0];
    let mut options = CooksReplacementOptions::new(
        resolve_cooks_cutoff(CooksCutoff::Default, samples, 2)
            .unwrap()
            .unwrap(),
    );
    options.min_replicates = 7;
    let mut unexpected_refits = Vec::new();
    let cases = r["cases"].as_array().unwrap();
    assert_eq!(cases.len(), 8);
    for c in cases {
        assert_eq!(c["refit_calls"], 0);
        assert_eq!(c["gene_fit_calls"], 1);
        assert_eq!(c["complete_results_identical_to_disabled"], true);
        let lrt = c["test"] == "LRT";
        let mut builder = DeseqBuilder::new()
            .automatic_dispersion_means()
            .fit_type(FitType::Mean)
            .size_factors(numbers(&r["size_factors"]))
            .disable_cooks_cutoff()
            .disable_independent_filtering();
        if lrt {
            builder = builder.test(TestType::Lrt).reduced_design(reduced.clone());
        }
        if c["normalization_factors"].as_bool().unwrap() {
            builder = builder.normalization_factors(
                RowMajorMatrix::from_row_major(
                    genes,
                    samples,
                    numbers(&r["normalization_factors"]),
                )
                .unwrap(),
            );
        }
        if c["weights"].as_bool().unwrap() {
            builder = builder.observation_weights(
                RowMajorMatrix::from_row_major(genes, samples, numbers(&r["weights"])).unwrap(),
            );
        }
        for fixed in [false, true] {
            for use_contrast in [false, true] {
                let output = match (fixed, lrt, use_contrast) {
                    (false, _, false) => builder
                        .fit_with_test_results_with_cooks_replacement(&counts, &design, &options),
                    (false, _, true) => builder
                        .fit_with_test_results_contrast_with_cooks_replacement(
                            &counts, &design, &contrast, &options,
                        ),
                    (true, false, false) => builder
                        .fit_fixed_dispersion_wald_with_cooks_replacement(
                            &counts,
                            &design,
                            &dispersions,
                            1,
                            &options,
                        )
                        .map(CooksReplacementTestOutput::Wald),
                    (true, false, true) => builder
                        .fit_fixed_dispersion_wald_contrast_with_cooks_replacement(
                            &counts,
                            &design,
                            &dispersions,
                            &contrast,
                            &options,
                        )
                        .map(CooksReplacementTestOutput::Wald),
                    (true, true, false) => builder
                        .fit_fixed_dispersion_lrt_with_cooks_replacement(
                            &counts,
                            &design,
                            &reduced,
                            &dispersions,
                            1,
                            &options,
                        )
                        .map(CooksReplacementTestOutput::Lrt),
                    (true, true, true) => builder
                        .fit_fixed_dispersion_lrt_contrast_with_cooks_replacement(
                            &counts,
                            &design,
                            &reduced,
                            &dispersions,
                            &contrast,
                            &options,
                        )
                        .map(CooksReplacementTestOutput::Lrt),
                }
                .unwrap();
                let (original, refit, raw_results, final_results, plan) = match &output {
                    CooksReplacementTestOutput::Wald(o) => (
                        &o.original_fit,
                        &o.refit_fit,
                        &o.original_results,
                        &o.results,
                        &o.refit_plan,
                    ),
                    CooksReplacementTestOutput::Lrt(o) => (
                        &o.original_fit,
                        &o.refit_fit,
                        &o.original_results,
                        &o.results,
                        &o.refit_plan,
                    ),
                };
                assert!(!plan.has_replaceable_samples());
                // A real candidate is essential: a naturally empty plan could pass
                // this test even if every unnecessary fitting branch still ran.
                assert!(
                    plan.n_refit > 0 && plan.should_refit,
                    "{} fixed={fixed}",
                    c["name"]
                );
                assert_eq!(plan.replacement.replaced_counts, counts);
                if refit.is_some() {
                    unexpected_refits.push(format!(
                        "{} fixed={fixed} contrast={use_contrast}",
                        c["name"]
                    ));
                    continue;
                }
                assert_eq!(final_results.rows, raw_results.rows);
                let snapshots = [
                    &original.beta,
                    &original.beta_se,
                    &original.beta_covariance,
                    &original.mu,
                    &original.cooks,
                ]
                .map(matrix_bits);
                let original_zero = original.all_zero.clone();
                let (merged, results, retained) = output.into_merged_fit().unwrap();
                assert_eq!(
                    [
                        &merged.beta,
                        &merged.beta_se,
                        &merged.beta_covariance,
                        &merged.mu,
                        &merged.cooks
                    ]
                    .map(matrix_bits),
                    snapshots
                );
                assert_eq!(merged.all_zero, original_zero);
                let request = ResultRequest::new(
                    ResultsContrast::Numeric(contrast.to_vec()),
                    WaldTestOptions::default(),
                    CooksCutoff::Disabled,
                    IndependentFilteringOptions {
                        enabled: false,
                        ..Default::default()
                    },
                )
                .unwrap();
                let later = merged
                    .results_for_request_with_replacement(&counts, &request, None, &retained)
                    .unwrap();
                for (actual, expected) in later.rows.iter().zip(&results.rows) {
                    assert_eq!(actual.base_mean, expected.base_mean);
                    assert_eq!(actual.stat, expected.stat);
                    assert_eq!(actual.pvalue, expected.pvalue);
                    assert_eq!(actual.padj, expected.padj);
                }
            }
        }
    }
    assert!(
        unexpected_refits.is_empty(),
        "unnecessary refits: {unexpected_refits:?}"
    );
}
