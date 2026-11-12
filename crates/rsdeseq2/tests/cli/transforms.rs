use super::*;

#[test]
fn cli_rlog_rejects_incomplete_option_combinations() {
    let dir = temp_dir("rlog-invalid-options");
    let frozen = dir.join("frozen.tsv");
    write_frozen_intercept_fixture(&frozen);
    let cases = [
        (
            vec!["--frozen-intercept", frozen.to_str().unwrap()],
            "--rlog-prior-variance is required with --frozen-intercept",
        ),
        (
            vec!["--rlog-prior-variance", "2.5"],
            "--rlog-prior-variance requires --frozen-intercept",
        ),
        (
            vec!["--blind=false"],
            "invalid dimensions for rlog design path",
        ),
    ];
    let counts = reference_data_path("counts.tsv");
    let output = dir.join("rlog.tsv");
    for (options, error) in cases {
        let mut args = vec![
            "rlog",
            "--counts",
            counts.to_str().unwrap(),
            "--fit-type",
            "mean",
            "--output",
            output.to_str().unwrap(),
        ];
        args.extend(options);
        run_cli_failure_contains(&args, error);
    }
}

#[test]
fn cli_weighted_wald_and_rlog_match_library_numerically_and_differ_from_unweighted() {
    let dir = temp_dir("weighted-routing-numerical");
    let weights_path = dir.join("weights.tsv");
    write_observation_weight_fixture(&weights_path);
    let weights = rsdeseq2::io::read_observation_weights_tsv(&weights_path).unwrap();
    let counts = common::reference_counts();
    let design = common::reference_full_design();
    let base = DeseqBuilder::new()
        .fit_type(FitType::Mean)
        .disable_cooks_cutoff()
        .disable_independent_filtering();
    let weighted = base.clone().observation_weights(weights);
    let (_, expected) = weighted.fit_with_results(&counts, &design).unwrap();
    let (_, unweighted) = base.fit_with_results(&counts, &design).unwrap();
    assert_ne!(
        expected.rows, unweighted.rows,
        "weights must change Wald values"
    );
    let output = dir.join("wald.tsv");
    let oracle = dir.join("library.tsv");
    run_cli(&[
        "wald",
        "--counts",
        reference_data_path("counts.tsv").to_str().unwrap(),
        "--design",
        reference_data_path("design_full.tsv").to_str().unwrap(),
        "--observation-weights",
        weights_path.to_str().unwrap(),
        "--fit-type",
        "mean",
        "--disable-cooks-cutoff",
        "--disable-independent-filtering",
        "--output",
        output.to_str().unwrap(),
    ]);
    rsdeseq2::io::write_deseq_results_tsv(&oracle, &expected).unwrap();
    assert_eq!(fs::read(&output).unwrap(), fs::read(&oracle).unwrap());
    let blind_design =
        DesignMatrix::from_row_major(4, 1, vec![1.0; 4], Some(vec!["Intercept".into()])).unwrap();
    let transformed = weighted
        .clone()
        .automatic_dispersion_means()
        .rlog_glm_mu(&counts, &blind_design)
        .unwrap();
    let unweighted = base
        .automatic_dispersion_means()
        .rlog_glm_mu(&counts, &blind_design)
        .unwrap();
    assert_ne!(
        transformed.transformed, unweighted.transformed,
        "weights must change rlog values"
    );
    let output = dir.join("rlog.tsv");
    run_cli(&[
        "rlog",
        "--counts",
        reference_data_path("counts.tsv").to_str().unwrap(),
        "--observation-weights",
        weights_path.to_str().unwrap(),
        "--fit-type",
        "mean",
        "--output",
        output.to_str().unwrap(),
    ]);
    assert_numeric_tsv(&output, transformed.transformed.as_slice());
}
