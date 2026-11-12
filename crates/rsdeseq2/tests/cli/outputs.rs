use super::*;

#[test]
fn cli_wald_writes_result_and_cooks_sidecars() {
    let dir = temp_dir("wald-sidecars");
    let output = dir.join("wald.tsv");
    let cooks = dir.join("cooks.tsv");
    let column_metadata = dir.join("result_columns.tsv");
    let table_metadata = dir.join("result_table_metadata.tsv");
    let filter_metadata = dir.join("filter_metadata.tsv");
    let filter_num_rej = dir.join("filter_num_rej.tsv");
    let filter_lowess = dir.join("filter_lowess.tsv");

    run_cli(&[
        "wald",
        "--counts",
        reference_data_path("counts.tsv").to_str().unwrap(),
        "--design",
        reference_data_path("design_full.tsv").to_str().unwrap(),
        "--fit-type",
        "mean",
        "--coefficient",
        "1",
        "--cooks-distance-output",
        cooks.to_str().unwrap(),
        "--result-column-metadata-output",
        column_metadata.to_str().unwrap(),
        "--result-table-metadata-output",
        table_metadata.to_str().unwrap(),
        "--independent-filter-metadata-output",
        filter_metadata.to_str().unwrap(),
        "--independent-filter-num-rej-output",
        filter_num_rej.to_str().unwrap(),
        "--independent-filter-lowess-output",
        filter_lowess.to_str().unwrap(),
        "--output",
        output.to_str().unwrap(),
    ]);

    assert_deseq_results_table(&output);
    assert_tsv_starts_with(&cooks, "gene\tsample1\tsample2\tsample3\tsample4\n");
    assert_tsv_contains(&cooks, "gene1\t");
    assert_tsv_starts_with(&column_metadata, "name\ttype\tdescription\n");
    assert_tsv_contains(&column_metadata, "baseMean\tresults\t");
    assert_tsv_starts_with(&table_metadata, "name\tvalue\n");
    assert_tsv_contains(&table_metadata, "testType\tWald\n");
    assert_tsv_starts_with(&filter_metadata, "name\tvalue\n");
    assert_tsv_contains(&filter_metadata, "alpha\t");
    assert_tsv_starts_with(&filter_num_rej, "theta\tnumRej\n");
    assert_tsv_starts_with(&filter_lowess, "x\ty\n");
    let marker: serde_json::Value =
        serde_json::from_slice(&fs::read(dir.join("wald.tsv.complete.json")).unwrap()).unwrap();
    let files = marker["files"].as_array().unwrap();
    assert_eq!(files.len(), 7);
    for file in files {
        use sha2::{Digest, Sha256};
        let bytes = fs::read(file["path"].as_str().unwrap()).unwrap();
        assert_eq!(file["bytes"].as_u64().unwrap(), bytes.len() as u64);
        assert_eq!(
            file["sha256"].as_str().unwrap(),
            Sha256::digest(&bytes)
                .iter()
                .map(|byte| format!("{byte:02x}"))
                .collect::<String>()
        );
    }
}

#[test]
fn cli_wald_writes_cooks_replacement_sidecars() {
    let dir = temp_dir("wald-replacement-sidecars");
    let output = dir.join("wald.tsv");
    let replacement_metadata = dir.join("replacement_metadata.tsv");
    let replacement_row_metadata = dir.join("replacement_rows.tsv");
    let replaced_counts = dir.join("replaced_counts.tsv");
    let candidate_counts = dir.join("candidate_counts.tsv");
    let outlier_cells = dir.join("outlier_cells.tsv");

    run_cli(&[
        "wald",
        "--counts",
        reference_data_path("counts.tsv").to_str().unwrap(),
        "--design",
        reference_data_path("design_full.tsv").to_str().unwrap(),
        "--fit-type",
        "mean",
        "--coefficient",
        "1",
        "--cooks-cutoff",
        "0",
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
    assert_tsv_starts_with(
        &replaced_counts,
        "gene\tsample1\tsample2\tsample3\tsample4\n",
    );
    assert_tsv_starts_with(
        &candidate_counts,
        "gene\tsample1\tsample2\tsample3\tsample4\n",
    );
    assert_tsv_starts_with(&outlier_cells, "gene\tsample1\tsample2\tsample3\tsample4\n");
}

#[test]
fn cli_wald_rejects_replacement_sidecars_without_replacement_refit() {
    let dir = temp_dir("wald-replacement-sidecars-disabled");
    let output = dir.join("wald.tsv");
    let replacement_metadata = dir.join("replacement_metadata.tsv");
    fs::write(&output, "previous result\n").unwrap();

    run_cli_failure_contains(
        &[
            "wald",
            "--counts",
            reference_data_path("counts.tsv").to_str().unwrap(),
            "--design",
            reference_data_path("design_full.tsv").to_str().unwrap(),
            "--fit-type",
            "mean",
            "--coefficient",
            "1",
            "--disable-cooks-cutoff",
            "--cooks-replacement-metadata-output",
            replacement_metadata.to_str().unwrap(),
            "--output",
            output.to_str().unwrap(),
        ],
        "invalid options: Cook's replacement sidecar output requires Cook's replacement/refit to run",
    );
    assert_eq!(fs::read_to_string(output).unwrap(), "previous result\n");
    assert!(!replacement_metadata.exists());
    assert!(!dir.join("wald.tsv.complete.json").exists());
    assert_eq!(fs::read_dir(dir).unwrap().count(), 1);
}

#[test]
fn cli_rejects_input_and_sidecar_collisions_before_writing() {
    let dir = temp_dir("path-collisions");
    let counts = dir.join("counts.tsv");
    fs::copy(reference_data_path("counts.tsv"), &counts).unwrap();
    let original = fs::read(&counts).unwrap();
    let output = dir.join("results.tsv");
    for (main, sidecar) in [(&counts, &output), (&output, &output)] {
        run_cli_failure_contains(
            &[
                "wald",
                "--counts",
                counts.to_str().unwrap(),
                "--design",
                reference_data_path("design_full.tsv").to_str().unwrap(),
                "--output",
                main.to_str().unwrap(),
                "--fit-diagnostics-output",
                sidecar.to_str().unwrap(),
            ],
            "path collision",
        );
        assert_eq!(fs::read(&counts).unwrap(), original);
        assert!(!output.exists());
        assert!(!dir.join("results.tsv.complete.json").exists());
    }
    let error = Command::new(env!("CARGO_BIN_EXE_rsdeseq2"))
        .args([
            "size-factors",
            "--counts",
            counts.to_str().unwrap(),
            "--output",
            dir.join("missing/results.tsv").to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(!error.status.success());
    assert_eq!(fs::read_dir(dir).unwrap().count(), 1);
}

#[test]
fn cli_lrt_writes_result_and_cooks_sidecars() {
    let dir = temp_dir("lrt-sidecars");
    let output = dir.join("lrt.tsv");
    let cooks = dir.join("cooks.tsv");
    let column_metadata = dir.join("result_columns.tsv");
    let table_metadata = dir.join("result_table_metadata.tsv");

    run_cli(&[
        "lrt",
        "--counts",
        reference_data_path("counts.tsv").to_str().unwrap(),
        "--design",
        reference_data_path("design_full.tsv").to_str().unwrap(),
        "--reduced-design",
        reference_data_path("design_reduced.tsv").to_str().unwrap(),
        "--fit-type",
        "mean",
        "--coefficient",
        "1",
        "--cooks-distance-output",
        cooks.to_str().unwrap(),
        "--result-column-metadata-output",
        column_metadata.to_str().unwrap(),
        "--result-table-metadata-output",
        table_metadata.to_str().unwrap(),
        "--output",
        output.to_str().unwrap(),
    ]);

    assert_deseq_results_table(&output);
    assert_tsv_starts_with(&cooks, "gene\tsample1\tsample2\tsample3\tsample4\n");
    assert_tsv_contains(&cooks, "gene1\t");
    assert_tsv_starts_with(&column_metadata, "name\ttype\tdescription\n");
    assert_tsv_contains(&column_metadata, "stat\tresults\t");
    assert_tsv_starts_with(&table_metadata, "name\tvalue\n");
    assert_tsv_contains(&table_metadata, "testType\tLRT\n");
}
