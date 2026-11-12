use super::*;

#[test]
fn cli_wald_rejects_additive_beta_prior_when_reported_design_disagrees() {
    let dir = temp_dir("wald-beta-prior-additive-design-mismatch");
    let standard = dir.join("standard.tsv");
    let condition_levels = dir.join("condition.tsv");
    let batch_levels = dir.join("batch.tsv");
    let dispersions = dir.join("dispersions.tsv");
    let base_mean = dir.join("base_mean.tsv");
    let disp_fit = dir.join("disp_fit.tsv");
    let output = dir.join("wald.tsv");
    write_standard_contrast_design_fixture(&standard);
    write_sample_level_fixture(&condition_levels);
    write_batch_sample_level_fixture(&batch_levels);
    write_gene_numeric_fixture(&dispersions, "dispersion", &[0.08, 0.12, 0.05, 0.2]);
    write_gene_numeric_fixture(&base_mean, "baseMean", &[16.5, 3.0, 97.5, 7.5]);
    write_gene_numeric_fixture(&disp_fit, "dispFit", &[0.08, 0.12, 0.05, 0.2]);
    let sample_levels = format!("{},{}", condition_levels.display(), batch_levels.display());

    run_cli_failure_contains(
        &[
            "wald",
            "--counts",
            reference_data_path("counts.tsv").to_str().unwrap(),
            "--design",
            standard.to_str().unwrap(),
            "--fit-type",
            "mean",
            "--disable-cooks-cutoff",
            "--beta-prior-additive-factors",
            "condition,batch",
            "--beta-prior-additive-references",
            "A,X",
            "--beta-prior-additive-sample-levels",
            &sample_levels,
            "--beta-prior-dispersions",
            dispersions.to_str().unwrap(),
            "--beta-prior-base-mean",
            base_mean.to_str().unwrap(),
            "--beta-prior-disp-fit",
            disp_fit.to_str().unwrap(),
            "--coefficient-name",
            "batch_Y_vs_X",
            "--output",
            output.to_str().unwrap(),
        ],
        "invalid options: reported --design does not match the beta-prior additive design generated from sample levels",
    );
}

#[test]
fn cli_wald_rejects_additive_beta_prior_mismatched_factor_lists() {
    let dir = temp_dir("wald-beta-prior-additive-list-mismatch");
    let standard = dir.join("standard.tsv");
    let condition_levels = dir.join("condition.tsv");
    let batch_levels = dir.join("batch.tsv");
    let dispersions = dir.join("dispersions.tsv");
    let base_mean = dir.join("base_mean.tsv");
    let disp_fit = dir.join("disp_fit.tsv");
    let output = dir.join("wald.tsv");
    write_additive_contrast_design_fixture(&standard);
    write_sample_level_fixture(&condition_levels);
    write_batch_sample_level_fixture(&batch_levels);
    write_gene_numeric_fixture(&dispersions, "dispersion", &[0.08, 0.12, 0.05, 0.2]);
    write_gene_numeric_fixture(&base_mean, "baseMean", &[16.5, 3.0, 97.5, 7.5]);
    write_gene_numeric_fixture(&disp_fit, "dispFit", &[0.08, 0.12, 0.05, 0.2]);
    let sample_levels = format!("{},{}", condition_levels.display(), batch_levels.display());

    run_cli_failure_contains(
        &[
            "wald",
            "--counts",
            reference_data_path("counts.tsv").to_str().unwrap(),
            "--design",
            standard.to_str().unwrap(),
            "--fit-type",
            "mean",
            "--disable-cooks-cutoff",
            "--beta-prior-additive-factors",
            "condition,batch",
            "--beta-prior-additive-references",
            "A",
            "--beta-prior-additive-sample-levels",
            &sample_levels,
            "--beta-prior-dispersions",
            dispersions.to_str().unwrap(),
            "--beta-prior-base-mean",
            base_mean.to_str().unwrap(),
            "--beta-prior-disp-fit",
            disp_fit.to_str().unwrap(),
            "--coefficient-name",
            "batch_Y_vs_X",
            "--output",
            output.to_str().unwrap(),
        ],
        "invalid dimensions for beta-prior additive factor inputs",
    );
}

#[test]
fn cli_wald_rejects_factor_beta_prior_when_reported_design_disagrees() {
    let dir = temp_dir("wald-beta-prior-factor-design-mismatch");
    let standard = dir.join("standard.tsv");
    let levels = dir.join("levels.tsv");
    let dispersions = dir.join("dispersions.tsv");
    let base_mean = dir.join("base_mean.tsv");
    let disp_fit = dir.join("disp_fit.tsv");
    let output = dir.join("wald.tsv");
    write_mismatched_standard_contrast_design_fixture(&standard);
    write_sample_level_fixture(&levels);
    write_gene_numeric_fixture(&dispersions, "dispersion", &[0.08, 0.12, 0.05, 0.2]);
    write_gene_numeric_fixture(&base_mean, "baseMean", &[16.5, 3.0, 97.5, 7.5]);
    write_gene_numeric_fixture(&disp_fit, "dispFit", &[0.08, 0.12, 0.05, 0.2]);

    run_cli_failure_contains(
        &[
            "wald",
            "--counts",
            reference_data_path("counts.tsv").to_str().unwrap(),
            "--design",
            standard.to_str().unwrap(),
            "--fit-type",
            "mean",
            "--disable-cooks-cutoff",
            "--beta-prior-factor",
            "condition",
            "--beta-prior-reference",
            "A",
            "--beta-prior-sample-levels",
            levels.to_str().unwrap(),
            "--beta-prior-dispersions",
            dispersions.to_str().unwrap(),
            "--beta-prior-base-mean",
            base_mean.to_str().unwrap(),
            "--beta-prior-disp-fit",
            disp_fit.to_str().unwrap(),
            "--coefficient-name",
            "condition_B_vs_A",
            "--output",
            output.to_str().unwrap(),
        ],
        "invalid options: reported --design does not match the beta-prior factor design generated from sample levels",
    );
}

#[test]
fn cli_wald_rejects_incomplete_expanded_beta_prior_inputs() {
    let dir = temp_dir("wald-beta-prior-expanded-incomplete");
    let standard = dir.join("standard.tsv");
    let expanded = dir.join("expanded.tsv");
    let output = dir.join("wald.tsv");
    write_standard_contrast_design_fixture(&standard);
    write_expanded_factor_design_fixture(&expanded);

    run_cli_failure_contains(
        &[
            "wald",
            "--counts",
            reference_data_path("counts.tsv").to_str().unwrap(),
            "--design",
            standard.to_str().unwrap(),
            "--fit-type",
            "mean",
            "--disable-cooks-cutoff",
            "--beta-prior-expanded-design",
            expanded.to_str().unwrap(),
            "--beta-prior-coefficient-groups",
            "0|2",
            "--output",
            output.to_str().unwrap(),
        ],
        "invalid dimensions for beta-prior expanded Wald inputs",
    );
}

#[test]
fn cli_wald_rejects_incomplete_factor_beta_prior_inputs() {
    let dir = temp_dir("wald-beta-prior-factor-incomplete");
    let standard = dir.join("standard.tsv");
    let levels = dir.join("levels.tsv");
    let dispersions = dir.join("dispersions.tsv");
    let base_mean = dir.join("base_mean.tsv");
    let output = dir.join("wald.tsv");
    write_standard_contrast_design_fixture(&standard);
    write_sample_level_fixture(&levels);
    write_gene_numeric_fixture(&dispersions, "dispersion", &[0.08, 0.12, 0.05, 0.2]);
    write_gene_numeric_fixture(&base_mean, "baseMean", &[16.5, 3.0, 97.5, 7.5]);

    run_cli_failure_contains(
        &[
            "wald",
            "--counts",
            reference_data_path("counts.tsv").to_str().unwrap(),
            "--design",
            standard.to_str().unwrap(),
            "--fit-type",
            "mean",
            "--disable-cooks-cutoff",
            "--beta-prior-factor",
            "condition",
            "--beta-prior-reference",
            "A",
            "--beta-prior-sample-levels",
            levels.to_str().unwrap(),
            "--beta-prior-dispersions",
            dispersions.to_str().unwrap(),
            "--beta-prior-base-mean",
            base_mean.to_str().unwrap(),
            "--coefficient-name",
            "condition_B_vs_A",
            "--output",
            output.to_str().unwrap(),
        ],
        "invalid dimensions for beta-prior factor Wald inputs",
    );
}

#[test]
fn cli_wald_rejects_incomplete_additive_beta_prior_inputs() {
    let dir = temp_dir("wald-beta-prior-additive-incomplete");
    let standard = dir.join("standard.tsv");
    let condition_levels = dir.join("condition.tsv");
    let batch_levels = dir.join("batch.tsv");
    let output = dir.join("wald.tsv");
    write_additive_contrast_design_fixture(&standard);
    write_sample_level_fixture(&condition_levels);
    write_batch_sample_level_fixture(&batch_levels);
    let sample_levels = format!("{},{}", condition_levels.display(), batch_levels.display());

    run_cli_failure_contains(
        &[
            "wald",
            "--counts",
            reference_data_path("counts.tsv").to_str().unwrap(),
            "--design",
            standard.to_str().unwrap(),
            "--fit-type",
            "mean",
            "--disable-cooks-cutoff",
            "--beta-prior-additive-factors",
            "condition,batch",
            "--beta-prior-additive-references",
            "A,X",
            "--beta-prior-additive-sample-levels",
            &sample_levels,
            "--coefficient-name",
            "batch_Y_vs_X",
            "--output",
            output.to_str().unwrap(),
        ],
        "invalid dimensions for beta-prior additive-factor Wald inputs",
    );
}

#[test]
fn cli_wald_rejects_factor_level_contrast_for_beta_prior_cli() {
    let dir = temp_dir("wald-beta-prior-factor-level-contrast");
    let standard = dir.join("standard.tsv");
    let expanded = dir.join("expanded.tsv");
    let dispersions = dir.join("dispersions.tsv");
    let base_mean = dir.join("base_mean.tsv");
    let disp_fit = dir.join("disp_fit.tsv");
    let output = dir.join("wald.tsv");
    write_standard_contrast_design_fixture(&standard);
    write_expanded_factor_design_fixture(&expanded);
    write_gene_numeric_fixture(&dispersions, "dispersion", &[0.08, 0.12, 0.05, 0.2]);
    write_gene_numeric_fixture(&base_mean, "baseMean", &[16.5, 3.0, 97.5, 7.5]);
    write_gene_numeric_fixture(&disp_fit, "dispFit", &[0.08, 0.12, 0.05, 0.2]);

    run_cli_failure_contains(
        &[
            "wald",
            "--counts",
            reference_data_path("counts.tsv").to_str().unwrap(),
            "--design",
            standard.to_str().unwrap(),
            "--fit-type",
            "mean",
            "--disable-cooks-cutoff",
            "--beta-prior-expanded-design",
            expanded.to_str().unwrap(),
            "--beta-prior-coefficient-groups",
            "0|2",
            "--beta-prior-dispersions",
            dispersions.to_str().unwrap(),
            "--beta-prior-base-mean",
            base_mean.to_str().unwrap(),
            "--beta-prior-disp-fit",
            disp_fit.to_str().unwrap(),
            "--contrast-factor",
            "condition",
            "--contrast-numerator",
            "B",
            "--contrast-denominator",
            "A",
            "--output",
            output.to_str().unwrap(),
        ],
        "invalid options: beta-prior expanded Wald CLI currently accepts coefficient, named, list, or numeric contrasts",
    );
}

#[test]
fn cli_wald_rejects_expanded_beta_prior_group_index_outside_expanded_design() {
    let dir = temp_dir("wald-beta-prior-expanded-bad-group");
    let standard = dir.join("standard.tsv");
    let expanded = dir.join("expanded.tsv");
    let dispersions = dir.join("dispersions.tsv");
    let base_mean = dir.join("base_mean.tsv");
    let disp_fit = dir.join("disp_fit.tsv");
    let output = dir.join("wald.tsv");
    write_standard_contrast_design_fixture(&standard);
    write_expanded_factor_design_fixture(&expanded);
    write_gene_numeric_fixture(&dispersions, "dispersion", &[0.08, 0.12, 0.05, 0.2]);
    write_gene_numeric_fixture(&base_mean, "baseMean", &[16.5, 3.0, 97.5, 7.5]);
    write_gene_numeric_fixture(&disp_fit, "dispFit", &[0.08, 0.12, 0.05, 0.2]);

    run_cli_failure_contains(
        &[
            "wald",
            "--counts",
            reference_data_path("counts.tsv").to_str().unwrap(),
            "--design",
            standard.to_str().unwrap(),
            "--fit-type",
            "mean",
            "--disable-cooks-cutoff",
            "--beta-prior-expanded-design",
            expanded.to_str().unwrap(),
            "--beta-prior-coefficient-groups",
            "0|5",
            "--beta-prior-dispersions",
            dispersions.to_str().unwrap(),
            "--beta-prior-base-mean",
            base_mean.to_str().unwrap(),
            "--beta-prior-disp-fit",
            disp_fit.to_str().unwrap(),
            "--coefficient-name",
            "condition_B_vs_A",
            "--output",
            output.to_str().unwrap(),
        ],
        "invalid options: coefficient group index 5 is outside the expanded design columns",
    );
}

#[test]
fn cli_wald_rejects_duplicate_expanded_beta_prior_group_indices() {
    let dir = temp_dir("wald-beta-prior-expanded-duplicate-group");
    let standard = dir.join("standard.tsv");
    let expanded = dir.join("expanded.tsv");
    let dispersions = dir.join("dispersions.tsv");
    let base_mean = dir.join("base_mean.tsv");
    let disp_fit = dir.join("disp_fit.tsv");
    let output = dir.join("wald.tsv");
    write_standard_contrast_design_fixture(&standard);
    write_expanded_factor_design_fixture(&expanded);
    write_gene_numeric_fixture(&dispersions, "dispersion", &[0.08, 0.12, 0.05, 0.2]);
    write_gene_numeric_fixture(&base_mean, "baseMean", &[16.5, 3.0, 97.5, 7.5]);
    write_gene_numeric_fixture(&disp_fit, "dispFit", &[0.08, 0.12, 0.05, 0.2]);

    run_cli_failure_contains(
        &[
            "wald",
            "--counts",
            reference_data_path("counts.tsv").to_str().unwrap(),
            "--design",
            standard.to_str().unwrap(),
            "--fit-type",
            "mean",
            "--disable-cooks-cutoff",
            "--beta-prior-expanded-design",
            expanded.to_str().unwrap(),
            "--beta-prior-coefficient-groups",
            "0|0",
            "--beta-prior-dispersions",
            dispersions.to_str().unwrap(),
            "--beta-prior-base-mean",
            base_mean.to_str().unwrap(),
            "--beta-prior-disp-fit",
            disp_fit.to_str().unwrap(),
            "--coefficient-name",
            "condition_B_vs_A",
            "--output",
            output.to_str().unwrap(),
        ],
        "invalid options: coefficient group index 0 appears more than once",
    );
}

#[test]
fn cli_wald_rejects_expanded_beta_prior_group_count_mismatch() {
    let dir = temp_dir("wald-beta-prior-expanded-group-count");
    let standard = dir.join("standard.tsv");
    let expanded = dir.join("expanded.tsv");
    let dispersions = dir.join("dispersions.tsv");
    let base_mean = dir.join("base_mean.tsv");
    let disp_fit = dir.join("disp_fit.tsv");
    let output = dir.join("wald.tsv");
    write_standard_contrast_design_fixture(&standard);
    write_expanded_factor_design_fixture(&expanded);
    write_gene_numeric_fixture(&dispersions, "dispersion", &[0.08, 0.12, 0.05, 0.2]);
    write_gene_numeric_fixture(&base_mean, "baseMean", &[16.5, 3.0, 97.5, 7.5]);
    write_gene_numeric_fixture(&disp_fit, "dispFit", &[0.08, 0.12, 0.05, 0.2]);

    run_cli_failure_contains(
        &[
            "wald",
            "--counts",
            reference_data_path("counts.tsv").to_str().unwrap(),
            "--design",
            standard.to_str().unwrap(),
            "--fit-type",
            "mean",
            "--disable-cooks-cutoff",
            "--beta-prior-expanded-design",
            expanded.to_str().unwrap(),
            "--beta-prior-coefficient-groups",
            "0|1|2",
            "--beta-prior-dispersions",
            dispersions.to_str().unwrap(),
            "--beta-prior-base-mean",
            base_mean.to_str().unwrap(),
            "--beta-prior-disp-fit",
            disp_fit.to_str().unwrap(),
            "--coefficient-name",
            "condition_B_vs_A",
            "--output",
            output.to_str().unwrap(),
        ],
        "invalid dimensions for beta-prior coefficient groups",
    );
}

#[test]
fn cli_wald_rejects_each_pair_of_beta_prior_input_modes() {
    let dir = temp_dir("beta-prior-exclusive");
    let design = dir.join("design.tsv");
    let expanded = dir.join("expanded.tsv");
    let levels = dir.join("levels.tsv");
    let batch = dir.join("batch.tsv");
    write_standard_contrast_design_fixture(&design);
    write_expanded_factor_design_fixture(&expanded);
    write_sample_level_fixture(&levels);
    write_batch_sample_level_fixture(&batch);
    let additive_levels = format!("{},{}", levels.display(), batch.display());
    let modes = [
        vec![
            "--beta-prior-expanded-design",
            expanded.to_str().unwrap(),
            "--beta-prior-coefficient-groups",
            "0|2",
        ],
        vec![
            "--beta-prior-factor",
            "condition",
            "--beta-prior-reference",
            "A",
            "--beta-prior-sample-levels",
            levels.to_str().unwrap(),
        ],
        vec![
            "--beta-prior-additive-factors",
            "condition,batch",
            "--beta-prior-additive-references",
            "A,X",
            "--beta-prior-additive-sample-levels",
            &additive_levels,
        ],
    ];
    for (a, b) in [(0, 1), (0, 2), (1, 2)] {
        let counts = reference_data_path("counts.tsv");
        let output = dir.join("results.tsv");
        let mut args = vec![
            "wald",
            "--counts",
            counts.to_str().unwrap(),
            "--design",
            design.to_str().unwrap(),
            "--output",
            output.to_str().unwrap(),
        ];
        args.extend(&modes[a]);
        args.extend(&modes[b]);
        run_cli_failure_contains(
            &args,
            "beta-prior expanded matrix, one-factor, and additive-factor inputs are mutually exclusive",
        );
    }
}
