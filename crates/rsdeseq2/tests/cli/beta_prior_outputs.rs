use super::*;

#[test]
fn cli_wald_expanded_beta_prior_writes_cooks_replacement_sidecars() {
    let dir = temp_dir("wald-beta-prior-expanded-replacement");
    let standard = dir.join("standard.tsv");
    let expanded = dir.join("expanded.tsv");
    let dispersions = dir.join("dispersions.tsv");
    let base_mean = dir.join("base_mean.tsv");
    let disp_fit = dir.join("disp_fit.tsv");
    let output = dir.join("wald.tsv");
    let replacement_metadata = dir.join("replacement_metadata.tsv");
    let replacement_row_metadata = dir.join("replacement_rows.tsv");
    let replaced_counts = dir.join("replaced_counts.tsv");
    let candidate_counts = dir.join("candidate_counts.tsv");
    let outlier_cells = dir.join("outlier_cells.tsv");
    write_standard_contrast_design_fixture(&standard);
    write_expanded_factor_design_fixture(&expanded);
    write_gene_numeric_fixture(&dispersions, "dispersion", &[0.08, 0.12, 0.05, 0.2]);
    write_gene_numeric_fixture(&base_mean, "baseMean", &[16.5, 3.0, 97.5, 7.5]);
    write_gene_numeric_fixture(&disp_fit, "dispFit", &[0.08, 0.12, 0.05, 0.2]);

    run_cli(&[
        "wald",
        "--counts",
        reference_data_path("counts.tsv").to_str().unwrap(),
        "--design",
        standard.to_str().unwrap(),
        "--fit-type",
        "mean",
        "--cooks-cutoff",
        "0.1",
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
        "--coefficient-name",
        "condition_B_vs_A",
        "--cooks-replacement-metadata-output",
        replacement_metadata.to_str().unwrap(),
        "--cooks-replacement-row-metadata-output",
        replacement_row_metadata.to_str().unwrap(),
        "--cooks-replaced-counts-output",
        replaced_counts.to_str().unwrap(),
        "--cooks-candidate-replacement-counts-output",
        candidate_counts.to_str().unwrap(),
        "--cooks-outlier-cells-output",
        outlier_cells.to_str().unwrap(),
        "--output",
        output.to_str().unwrap(),
    ]);

    assert_deseq_results_table(&output);
    assert_tsv_starts_with(&replacement_metadata, "name\tvalue\n");
    assert_tsv_contains(&replacement_metadata, "nRefit\t");
    assert_tsv_starts_with(
        &replacement_row_metadata,
        "gene\treplace\trefitReplace\tnewAllZero\treplacedAllZero\t",
    );
    assert_matrix_table(&replaced_counts);
    assert_matrix_table(&candidate_counts);
    assert_tsv_starts_with(&outlier_cells, "gene\tsample1\tsample2\tsample3\tsample4\n");
}

#[test]
fn cli_wald_factor_beta_prior_writes_cooks_replacement_sidecars() {
    let dir = temp_dir("wald-beta-prior-factor-replacement");
    let standard = dir.join("standard.tsv");
    let levels = dir.join("levels.tsv");
    let dispersions = dir.join("dispersions.tsv");
    let base_mean = dir.join("base_mean.tsv");
    let disp_fit = dir.join("disp_fit.tsv");
    let output = dir.join("wald.tsv");
    let replacement_metadata = dir.join("replacement_metadata.tsv");
    let replacement_row_metadata = dir.join("replacement_rows.tsv");
    let replaced_counts = dir.join("replaced_counts.tsv");
    let candidate_counts = dir.join("candidate_counts.tsv");
    let outlier_cells = dir.join("outlier_cells.tsv");
    write_standard_contrast_design_fixture(&standard);
    write_sample_level_fixture(&levels);
    write_gene_numeric_fixture(&dispersions, "dispersion", &[0.08, 0.12, 0.05, 0.2]);
    write_gene_numeric_fixture(&base_mean, "baseMean", &[16.5, 3.0, 97.5, 7.5]);
    write_gene_numeric_fixture(&disp_fit, "dispFit", &[0.08, 0.12, 0.05, 0.2]);

    run_cli(&[
        "wald",
        "--counts",
        reference_data_path("counts.tsv").to_str().unwrap(),
        "--design",
        standard.to_str().unwrap(),
        "--fit-type",
        "mean",
        "--cooks-cutoff",
        "0.1",
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
        "--cooks-replacement-metadata-output",
        replacement_metadata.to_str().unwrap(),
        "--cooks-replacement-row-metadata-output",
        replacement_row_metadata.to_str().unwrap(),
        "--cooks-replaced-counts-output",
        replaced_counts.to_str().unwrap(),
        "--cooks-candidate-replacement-counts-output",
        candidate_counts.to_str().unwrap(),
        "--cooks-outlier-cells-output",
        outlier_cells.to_str().unwrap(),
        "--output",
        output.to_str().unwrap(),
    ]);

    assert_deseq_results_table(&output);
    assert_tsv_starts_with(&replacement_metadata, "name\tvalue\n");
    assert_tsv_contains(&replacement_metadata, "nRefit\t");
    assert_tsv_starts_with(
        &replacement_row_metadata,
        "gene\treplace\trefitReplace\tnewAllZero\treplacedAllZero\t",
    );
    assert_matrix_table(&replaced_counts);
    assert_matrix_table(&candidate_counts);
    assert_tsv_starts_with(&outlier_cells, "gene\tsample1\tsample2\tsample3\tsample4\n");
}

#[test]
fn cli_wald_additive_beta_prior_writes_cooks_replacement_sidecars() {
    let dir = temp_dir("wald-beta-prior-additive-replacement");
    let standard = dir.join("standard.tsv");
    let condition_levels = dir.join("condition.tsv");
    let batch_levels = dir.join("batch.tsv");
    let dispersions = dir.join("dispersions.tsv");
    let base_mean = dir.join("base_mean.tsv");
    let disp_fit = dir.join("disp_fit.tsv");
    let output = dir.join("wald.tsv");
    let replacement_metadata = dir.join("replacement_metadata.tsv");
    let replacement_row_metadata = dir.join("replacement_rows.tsv");
    let replaced_counts = dir.join("replaced_counts.tsv");
    let candidate_counts = dir.join("candidate_counts.tsv");
    let outlier_cells = dir.join("outlier_cells.tsv");
    write_additive_contrast_design_fixture(&standard);
    write_sample_level_fixture(&condition_levels);
    write_batch_sample_level_fixture(&batch_levels);
    write_gene_numeric_fixture(&dispersions, "dispersion", &[0.08, 0.12, 0.05, 0.2]);
    write_gene_numeric_fixture(&base_mean, "baseMean", &[16.5, 3.0, 97.5, 7.5]);
    write_gene_numeric_fixture(&disp_fit, "dispFit", &[0.08, 0.12, 0.05, 0.2]);
    let sample_levels = format!("{},{}", condition_levels.display(), batch_levels.display());

    run_cli(&[
        "wald",
        "--counts",
        reference_data_path("counts.tsv").to_str().unwrap(),
        "--design",
        standard.to_str().unwrap(),
        "--fit-type",
        "mean",
        "--cooks-cutoff",
        "0.1",
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
        "--cooks-replacement-metadata-output",
        replacement_metadata.to_str().unwrap(),
        "--cooks-replacement-row-metadata-output",
        replacement_row_metadata.to_str().unwrap(),
        "--cooks-replaced-counts-output",
        replaced_counts.to_str().unwrap(),
        "--cooks-candidate-replacement-counts-output",
        candidate_counts.to_str().unwrap(),
        "--cooks-outlier-cells-output",
        outlier_cells.to_str().unwrap(),
        "--output",
        output.to_str().unwrap(),
    ]);

    assert_deseq_results_table(&output);
    assert_tsv_starts_with(&replacement_metadata, "name\tvalue\n");
    assert_tsv_contains(&replacement_metadata, "nRefit\t");
    assert_tsv_starts_with(
        &replacement_row_metadata,
        "gene\treplace\trefitReplace\tnewAllZero\treplacedAllZero\t",
    );
    assert_matrix_table(&replaced_counts);
    assert_matrix_table(&candidate_counts);
    assert_tsv_starts_with(&outlier_cells, "gene\tsample1\tsample2\tsample3\tsample4\n");
}
