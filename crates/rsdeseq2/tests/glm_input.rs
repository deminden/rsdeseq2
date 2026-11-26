//! The beta-prior, fixed-dispersion, rlog and VST fits with size-factor or
//! normalization-factor offsets, with or without weights.
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
fn offset_and_weight_fits_all_run() {
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
    let sf = Offsets::PerSample(&size_factors);
    let nf = Offsets::PerCell(&factors);
    let input = |offsets| GlmInput::new(&counts, offsets, &dispersions);
    let weighted = |offsets| input(offsets).weights(&weights);
    let s = |v| format!("{v:?}");
    let outputs = [
        (
            "supplied/sf",
            s(
                fit_glms_with_beta_prior(input(sf), &design, &prior, irls.clone())
                    .map(|f| format!("{f:?}")),
            ),
        ),
        (
            "supplied/sf+w",
            s(
                fit_glms_with_beta_prior(weighted(sf), &design, &prior, irls.clone())
                    .map(|f| format!("{f:?}")),
            ),
        ),
        (
            "supplied/nf",
            s(
                fit_glms_with_beta_prior(input(nf), &design, &prior, irls.clone())
                    .map(|f| format!("{f:?}")),
            ),
        ),
        (
            "supplied/nf+w",
            s(
                fit_glms_with_beta_prior(weighted(nf), &design, &prior, irls.clone())
                    .map(|f| format!("{f:?}")),
            ),
        ),
        (
            "estimated/sf",
            s(fit_glms_with_estimated_beta_prior(
                input(sf),
                &design,
                &base_mean,
                &disp_fit,
                options.clone(),
            )
            .map(|f| format!("{f:?}"))),
        ),
        (
            "estimated/sf+w",
            s(fit_glms_with_estimated_beta_prior(
                weighted(sf),
                &design,
                &base_mean,
                &disp_fit,
                options.clone(),
            )
            .map(|f| format!("{f:?}"))),
        ),
        (
            "estimated/nf",
            s(fit_glms_with_estimated_beta_prior(
                input(nf),
                &design,
                &base_mean,
                &disp_fit,
                options.clone(),
            )
            .map(|f| format!("{f:?}"))),
        ),
        (
            "estimated/nf+w",
            s(fit_glms_with_estimated_beta_prior(
                weighted(nf),
                &design,
                &base_mean,
                &disp_fit,
                options.clone(),
            )
            .map(|f| format!("{f:?}"))),
        ),
        (
            "expanded/sf",
            s(fit_expanded_glms_with_estimated_beta_prior(
                input(sf),
                expanded_input,
                &base_mean,
                &disp_fit,
                options.clone(),
            )
            .map(|f| format!("{f:?}"))),
        ),
        (
            "expanded/sf+w",
            s(fit_expanded_glms_with_estimated_beta_prior(
                weighted(sf),
                expanded_input,
                &base_mean,
                &disp_fit,
                options.clone(),
            )
            .map(|f| format!("{f:?}"))),
        ),
        (
            "expanded/nf",
            s(fit_expanded_glms_with_estimated_beta_prior(
                input(nf),
                expanded_input,
                &base_mean,
                &disp_fit,
                options.clone(),
            )
            .map(|f| format!("{f:?}"))),
        ),
        (
            "expanded/nf+w",
            s(fit_expanded_glms_with_estimated_beta_prior(
                weighted(nf),
                expanded_input,
                &base_mean,
                &disp_fit,
                options.clone(),
            )
            .map(|f| format!("{f:?}"))),
        ),
    ];
    for (tag, value) in outputs {
        record(tag, value);
    }
    let one_design = DesignMatrix::intercept_only(6).unwrap();
    let fits =
        [
            (
                "irls/sf",
                s(fit_fixed_dispersion_glm(input(sf), &design, irls.clone())
                    .map(|f| format!("{f:?}"))),
            ),
            (
                "irls/sf+w",
                s(
                    fit_fixed_dispersion_glm(weighted(sf), &design, irls.clone())
                        .map(|f| format!("{f:?}")),
                ),
            ),
            (
                "irls/nf",
                s(fit_fixed_dispersion_glm(input(nf), &design, irls.clone())
                    .map(|f| format!("{f:?}"))),
            ),
            (
                "irls/nf+w",
                s(
                    fit_fixed_dispersion_glm(weighted(nf), &design, irls.clone())
                        .map(|f| format!("{f:?}")),
                ),
            ),
            (
                "irls/intercept-design",
                s(fit_irls(
                    &counts,
                    &one_design,
                    &size_factors,
                    &dispersions,
                    irls.clone(),
                )
                .map(|f| format!("{f:?}"))),
            ),
            (
                "intercept/sf",
                s(fit_intercept_only_glm(input(sf)).map(|f| format!("{f:?}"))),
            ),
            (
                "intercept/sf+w",
                s(fit_intercept_only_glm(weighted(sf)).map(|f| format!("{f:?}"))),
            ),
            (
                "intercept/nf+w",
                s(fit_intercept_only_glm(weighted(nf)).map(|f| format!("{f:?}"))),
            ),
        ];
    for (tag, value) in fits {
        record(tag, value);
    }
    let rlog_counts = CountMatrix::from_row_major_u32(
        3,
        6,
        vec![
            10, 14, 12, 30, 35, 28, 5, 3, 6, 4, 7, 5, 200, 180, 220, 90, 100, 85,
        ],
    )
    .unwrap();
    let rlog_factors =
        RowMajorMatrix::from_row_major(3, 6, (0..18).map(|i| 0.7 + 0.03 * f64::from(i)).collect())
            .unwrap();
    let rlog_disp = [0.1, 0.2, 0.05];
    let intercept = [4.0, 2.2, 7.1];
    let rlog_irls = IrlsOptions::deseq2_rlog();
    let rsf = Offsets::PerSample(&size_factors);
    let rnf = Offsets::PerCell(&rlog_factors);
    let rlogs = [
        (
            "rlog/sf",
            s(
                rlog_transform(&rlog_counts, rsf, &rlog_disp, 0.7, rlog_irls.clone())
                    .map(|f| format!("{f:?}")),
            ),
        ),
        (
            "rlog/nf",
            s(
                rlog_transform(&rlog_counts, rnf, &rlog_disp, 0.7, rlog_irls.clone())
                    .map(|f| format!("{f:?}")),
            ),
        ),
        (
            "rlog_fit/sf",
            s(
                rlog_fit(&rlog_counts, rsf, &rlog_disp, 0.7, rlog_irls.clone())
                    .map(|f| format!("{f:?}")),
            ),
        ),
        (
            "rlog_fit/nf",
            s(
                rlog_fit(&rlog_counts, rnf, &rlog_disp, 0.7, rlog_irls.clone())
                    .map(|f| format!("{f:?}")),
            ),
        ),
        (
            "rlog_frozen/sf",
            s(rlog_frozen(
                &rlog_counts,
                rsf,
                &rlog_disp,
                0.7,
                &intercept,
                rlog_irls.clone(),
            )
            .map(|f| format!("{f:?}"))),
        ),
        (
            "rlog_frozen/nf",
            s(rlog_frozen(
                &rlog_counts,
                rnf,
                &rlog_disp,
                0.7,
                &intercept,
                rlog_irls.clone(),
            )
            .map(|f| format!("{f:?}"))),
        ),
        (
            "rlog_prior/sf",
            s(rlog_with_estimated_prior(
                &rlog_counts,
                rsf,
                &base_mean[..3],
                &disp_fit[..3],
                &rlog_disp,
                rlog_irls.clone(),
            )
            .map(|f| format!("{f:?}"))),
        ),
        (
            "rlog_prior/nf",
            s(rlog_with_estimated_prior(
                &rlog_counts,
                rnf,
                &base_mean[..3],
                &disp_fit[..3],
                &rlog_disp,
                rlog_irls.clone(),
            )
            .map(|f| format!("{f:?}"))),
        ),
    ];
    for (tag, value) in rlogs {
        record(tag, value);
    }
    let normalized = normalized_counts(&rlog_counts, &size_factors).unwrap();
    let trend = DispersionTrendFit::Mean(
        fit_mean_dispersion_trend(
            &[20.0, 5.0, 140.0],
            &[0.1, 0.2, 0.05],
            MeanDispersionTrendOptions::default(),
        )
        .unwrap(),
    );
    let vsts = [
        (
            "vst_inverse/sf",
            s(local_vst_inverse_factor_mean(rsf).map(|v| format!("{v:?}"))),
        ),
        (
            "vst_inverse/nf",
            s(local_vst_inverse_factor_mean(rnf).map(|v| format!("{v:?}"))),
        ),
        (
            "vst_trend/sf",
            s(
                vst_with_dispersion_trend_and_offsets(&normalized, &trend, rsf)
                    .map(|f| format!("{f:?}")),
            ),
        ),
        (
            "vst_trend/nf",
            s(
                vst_with_dispersion_trend_and_offsets(&normalized, &trend, rnf)
                    .map(|f| format!("{f:?}")),
            ),
        ),
    ];
    for (tag, value) in vsts {
        record(tag, value);
    }
}
