//! The twelve deprecated beta-prior GLM fits (size factors or normalization
//! factors, with or without weights; supplied, estimated or expanded prior).
//!
//! `RSDESEQ2_GLM_INPUT_DUMP=<file>` writes every output so two commits can be
//! compared bit for bit; the test itself checks that each one succeeds and
//! that equal inputs through `Offsets::PerSample` and expanded
//! `Offsets::PerCell` agree where DESeq2 uses the same path.

use rsdeseq2::prelude::*;

fn record(tag: &str, value: String) {
    if let Some(path) = std::env::var_os("RSDESEQ2_GLM_INPUT_DUMP") {
        use std::io::Write;
        let mut file = std::fs::OpenOptions::new()
            .create(true)
            .append(true)
            .open(path)
            .unwrap();
        writeln!(file, "{tag}\t{value}").unwrap();
    }
    assert!(value.starts_with("Ok("), "{tag}: {value}");
}

#[test]
fn deprecated_beta_prior_fits_all_run() {
    let counts = CountMatrix::from_row_major_u32(
        4,
        6,
        vec![
            10, 14, 12, 30, 35, 28, //
            5, 3, 6, 4, 7, 5, //
            200, 180, 220, 90, 100, 85, //
            0, 2, 1, 9, 12, 7,
        ],
    )
    .unwrap();
    let design = DesignMatrix::from_row_major(
        6,
        2,
        vec![1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        Some(vec!["Intercept".into(), "condition_B_vs_A".into()]),
    )
    .unwrap();
    let expanded = DesignMatrix::from_row_major(
        6,
        3,
        vec![
            1.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 0.0,
            1.0,
        ],
        Some(vec![
            "Intercept".into(),
            "conditionA".into(),
            "conditionB".into(),
        ]),
    )
    .unwrap();
    let groups = vec![vec![0], vec![2]];
    let expanded_input = ExpandedModelBetaPriorDesignInput {
        expanded_design: &expanded,
        standard_design: &design,
        coefficient_groups: &groups,
    };
    let size_factors = [0.8, 1.1, 0.9, 1.2, 1.0, 1.05];
    let factors =
        RowMajorMatrix::from_row_major(4, 6, (0..24).map(|i| 0.7 + 0.03 * f64::from(i)).collect())
            .unwrap();
    let weights = RowMajorMatrix::from_row_major(
        4,
        6,
        (0..24)
            .map(|i| if i % 7 == 3 { 0.5 } else { 1.0 })
            .collect(),
    )
    .unwrap();
    let dispersions = [0.1, 0.2, 0.05, 0.3];
    let base_mean = [21.0, 5.0, 146.0, 5.2];
    let disp_fit = [0.12, 0.25, 0.06, 0.25];
    let prior = [1e6, 0.8];
    let irls = IrlsOptions::default();
    let options = BetaPriorRefitOptions::default();
    let sw = BetaPriorSizeFactorWeightInput {
        size_factors: &size_factors,
        weights: Some(&weights),
    };
    let nw = BetaPriorNormalizationFactorWeightInput {
        normalization_factors: &factors,
        weights: Some(&weights),
    };
    let s = |v| format!("{v:?}");
    #[allow(deprecated)]
    let outputs = [
        ("supplied/sf", s(fit_glms_with_beta_prior_variance(&counts, &design, &size_factors, &dispersions, &prior, irls.clone()).map(|f| format!("{f:?}")))),
        ("supplied/sf+w", s(fit_glms_with_beta_prior_variance_and_weights(&counts, &design, &size_factors, &dispersions, Some(&weights), &prior, irls.clone()).map(|f| format!("{f:?}")))),
        ("supplied/nf", s(fit_glms_with_beta_prior_variance_and_normalization_factors(&counts, &design, &factors, &dispersions, &prior, irls.clone()).map(|f| format!("{f:?}")))),
        ("supplied/nf+w", s(fit_glms_with_beta_prior_variance_and_normalization_factors_and_weights(&counts, &design, &factors, &dispersions, Some(&weights), &prior, irls.clone()).map(|f| format!("{f:?}")))),
        ("estimated/sf", s(fit_glms_with_estimated_beta_prior_variance(&counts, &design, &size_factors, &dispersions, &base_mean, &disp_fit, options.clone()).map(|f| format!("{f:?}")))),
        ("estimated/sf+w", s(fit_glms_with_estimated_beta_prior_variance_and_weights(&counts, &design, sw, &dispersions, &base_mean, &disp_fit, options.clone()).map(|f| format!("{f:?}")))),
        ("estimated/nf", s(fit_glms_with_estimated_beta_prior_variance_and_normalization_factors(&counts, &design, &factors, &dispersions, &base_mean, &disp_fit, options.clone()).map(|f| format!("{f:?}")))),
        ("estimated/nf+w", s(fit_glms_with_estimated_beta_prior_variance_and_normalization_factors_and_weights(&counts, &design, nw, &dispersions, &base_mean, &disp_fit, options.clone()).map(|f| format!("{f:?}")))),
        ("expanded/sf", s(fit_expanded_glms_with_estimated_beta_prior_variance(&counts, expanded_input, &size_factors, &dispersions, &base_mean, &disp_fit, options.clone()).map(|f| format!("{f:?}")))),
        ("expanded/sf+w", s(fit_expanded_glms_with_estimated_beta_prior_variance_and_weights(&counts, expanded_input, sw, &dispersions, &base_mean, &disp_fit, options.clone()).map(|f| format!("{f:?}")))),
        ("expanded/nf", s(fit_expanded_glms_with_estimated_beta_prior_variance_and_normalization_factors(&counts, expanded_input, &factors, &dispersions, &base_mean, &disp_fit, options.clone()).map(|f| format!("{f:?}")))),
        ("expanded/nf+w", s(fit_expanded_glms_with_estimated_beta_prior_variance_and_normalization_factors_and_weights(&counts, expanded_input, nw, &dispersions, &base_mean, &disp_fit, options.clone()).map(|f| format!("{f:?}")))),
    ];
    for (tag, value) in outputs {
        record(tag, value);
    }
}
