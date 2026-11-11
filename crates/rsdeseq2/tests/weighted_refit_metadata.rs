//! Original replacement counts and raw-weight base metadata, independent of optimizers.
use rsdeseq2::{cooks::prepare_cooks_replacement_refit_with_weights, prelude::*};
use serde_json::Value;
fn number(v: &Value) -> f64 {
    let s = v.as_str().unwrap();
    if let Some(h) = s.strip_prefix("bits:") {
        f64::from_bits(u64::from_str_radix(h, 16).unwrap())
    } else {
        s.parse().unwrap()
    }
}
fn matrix(v: &Value) -> RowMajorMatrix<f64> {
    let r = v["rows"].as_u64().unwrap() as usize;
    let c = v["columns"].as_u64().unwrap() as usize;
    let values = v["values"]
        .as_array()
        .unwrap()
        .iter()
        .map(number)
        .collect::<Vec<_>>();
    RowMajorMatrix::from_row_major(
        r,
        c,
        (0..r)
            .flat_map(|i| {
                (0..c).map({
                    let values = &values;
                    move |j| values[i + j * r]
                })
            })
            .collect(),
    )
    .unwrap()
}
#[test]
fn replacement_metadata_matches_original_deseq2_with_raw_weights_and_both_normalizations() {
    let reference: Value =
        serde_json::from_str(include_str!("data/weighted_refit_metadata_r462.json")).unwrap();
    for c in reference["cases"].as_array().unwrap() {
        let counts = matrix(&c["counts"]);
        let counts = CountMatrix::from_row_major_u32(
            counts.n_rows(),
            counts.n_cols(),
            counts.as_slice().iter().map(|v| *v as u32).collect(),
        )
        .unwrap();
        let m = matrix(&c["design"]);
        let design =
            DesignMatrix::from_row_major(m.n_rows(), m.n_cols(), m.into_values(), None).unwrap();
        let sf = c["size_factors"]
            .as_array()
            .unwrap()
            .iter()
            .map(number)
            .collect::<Vec<_>>();
        let nf =
            (!c["normalization_factors"].is_null()).then(|| matrix(&c["normalization_factors"]));
        let weights = (!c["weights"].is_null()).then(|| matrix(&c["weights"]));
        let normalized = if let Some(nf) = &nf {
            normalized_counts_with_factors(&counts, nf).unwrap()
        } else {
            normalized_counts(&counts, &sf).unwrap()
        };
        let mut options = CooksReplacementOptions::new(number(&c["cutoff"]));
        options.min_replicates = 3;
        options.which_samples = Some(
            c["which_samples"]
                .as_array()
                .unwrap()
                .iter()
                .map(|v| v.as_bool().unwrap())
                .collect(),
        );
        let plan = prepare_cooks_replacement_refit_with_weights(
            &counts,
            &normalized,
            &sf,
            nf.as_ref(),
            &matrix(&c["cooks"]),
            &design,
            &options,
            weights.as_ref(),
        )
        .unwrap();
        assert_eq!(
            plan.replacement
                .replaced_counts
                .as_slice()
                .iter()
                .map(|v| f64::from(*v))
                .collect::<Vec<_>>(),
            matrix(&c["replaced_counts"]).as_slice(),
            "{} counts",
            c["name"]
        );
        assert_eq!(plan.n_refit, c["n_refit"].as_u64().unwrap() as usize);
        assert_eq!(
            plan.replaced_all_zero,
            c["all_zero"]
                .as_array()
                .unwrap()
                .iter()
                .map(|v| v.as_bool().unwrap())
                .collect::<Vec<_>>()
        );
        for (field, actual) in [
            ("base_mean", &plan.replaced_base_mean),
            ("base_var", &plan.replaced_base_var),
        ] {
            for (i, (a, b)) in actual
                .iter()
                .zip(c[field].as_array().unwrap().iter().map(number))
                .enumerate()
            {
                assert!(
                    a.is_finite() && (a - b).abs() <= 1e-10 + 1e-12 * b.abs(),
                    "{} {field} row {i}: {a} != {b}",
                    c["name"]
                );
            }
        }
    }
}

#[test]
fn weighted_failure_state_after_replacement_matches_original_fits_and_later_diagnostics() {
    let reference: Value =
        serde_json::from_str(include_str!("data/weighted_replacement_state_r462.json")).unwrap();
    let numbers = |v: &Value| v.as_array().unwrap().iter().map(number).collect::<Vec<_>>();
    let flags = |v: &Value| {
        v.as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_bool().unwrap())
            .collect::<Vec<_>>()
    };
    let genes = reference["genes"].as_u64().unwrap() as usize;
    let samples = reference["samples"].as_u64().unwrap() as usize;
    let failed = reference["failed_row"].as_u64().unwrap() as usize;
    let design =
        DesignMatrix::from_row_major(samples, 2, numbers(&reference["design"]), None).unwrap();
    let request = ResultRequest::new(
        ResultsContrast::Numeric(vec![0.0, 1.0]),
        WaldTestOptions::default(),
        CooksCutoff::Disabled,
        IndependentFilteringOptions {
            enabled: false,
            ..Default::default()
        },
    )
    .unwrap();
    let mut mismatches = Vec::new();
    for c in reference["cases"].as_array().unwrap() {
        let counts = CountMatrix::from_row_major_u32(
            genes,
            samples,
            c["counts"]
                .as_array()
                .unwrap()
                .iter()
                .map(|v| v.as_u64().unwrap() as u32)
                .collect(),
        )
        .unwrap();
        let mut builder = DeseqBuilder::new()
            .automatic_dispersion_means()
            .fit_type(FitType::Mean)
            .size_factors(numbers(&reference["size_factors"]))
            .observation_weights(
                RowMajorMatrix::from_row_major(genes, samples, numbers(&reference["weights"]))
                    .unwrap(),
            )
            .disable_cooks_cutoff()
            .disable_independent_filtering();
        if c["normalization_factors"].as_bool().unwrap() {
            builder = builder.normalization_factors(
                RowMajorMatrix::from_row_major(genes, samples, numbers(&reference["factors"]))
                    .unwrap(),
            );
        }
        if c["test"] == "LRT" {
            builder = builder.test(TestType::Lrt).reduced_design(
                DesignMatrix::from_row_major(samples, 1, vec![1.0; samples], None).unwrap(),
            );
        }
        let (fit, results, plan) = if c["kind"] == "disabled" {
            let (f, r) = builder.fit_with_results(&counts, &design).unwrap();
            (f, r, None)
        } else {
            let mut options = CooksReplacementOptions::new(
                resolve_cooks_cutoff(CooksCutoff::Default, samples, 2)
                    .unwrap()
                    .unwrap(),
            );
            options.min_replicates = 3;
            let (f, r, p) = builder
                .fit_with_test_results_with_cooks_replacement(&counts, &design, &options)
                .unwrap()
                .into_merged_fit()
                .unwrap();
            (f, r, Some(p))
        };
        assert_eq!(
            plan.as_ref().map_or(0, |p| p.n_refit),
            c["n_refit"].as_u64().unwrap() as usize,
            "{} replacement count",
            c["name"]
        );
        if fit.all_zero != flags(&c["all_zero"]) {
            mismatches.push(c["name"].as_str().unwrap());
        }
        assert_eq!(
            fit.weights_fail.as_ref().unwrap(),
            &flags(&c["weights_fail"])
        );
        for (actual, expected) in [
            (fit.base_mean[failed], number(&c["failed_mean"])),
            (fit.base_var[failed], number(&c["failed_variance"])),
        ] {
            assert!(
                (actual - expected).abs() <= 1e-10 + 1e-12 * expected.abs(),
                "{} failed-row base summary",
                c["name"]
            );
        }
        let later = if let Some(p) = &plan {
            fit.results_for_request_with_replacement(&counts, &request, None, p)
                .unwrap()
        } else {
            fit.results_for_request(&counts, &request, None).unwrap()
        };
        for rows in [&results.rows, &later.rows] {
            assert_eq!(
                rows[failed].stat.is_none(),
                c["failed_stat_missing"].as_bool().unwrap(),
                "{} statistic",
                c["name"]
            );
            assert_eq!(
                rows[failed].pvalue.is_none(),
                c["failed_pvalue_missing"].as_bool().unwrap(),
                "{} p-value",
                c["name"]
            );
            assert_eq!(
                rows[failed].converged.is_none(),
                c["failed_convergence_missing"].as_bool().unwrap(),
                "{} convergence",
                c["name"]
            );
        }
    }
    assert!(
        mismatches.is_empty(),
        "allZero differs from original DESeq2: {mismatches:?}"
    );
}
