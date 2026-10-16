mod common;
use rsdeseq2::prelude::*;

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicUsize, Ordering};

static TEMP_COUNTER: AtomicUsize = AtomicUsize::new(0);

fn temp_dir(label: &str) -> PathBuf {
    let id = TEMP_COUNTER.fetch_add(1, Ordering::Relaxed);
    let path =
        std::env::temp_dir().join(format!("rsdeseq2-cli-{label}-{}-{id}", std::process::id()));
    fs::create_dir_all(&path).unwrap();
    path
}

fn reference_data_path(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/deseq2_reference")
        .join(name)
}

fn write_size_factor_fixture(path: &Path) {
    fs::write(
        path,
        "\
sample\tsize_factor
sample1\t1
sample3\t1.5
sample4\t2
sample2\t0.5
",
    )
    .unwrap();
}

fn write_wald_t_degrees_of_freedom_fixture(path: &Path) {
    fs::write(
        path,
        "\
gene\tdf
gene1\t4
gene2\t4
gene3\t4
gene4\t4
",
    )
    .unwrap();
}

fn write_frozen_intercept_fixture(path: &Path) {
    fs::write(
        path,
        "\
gene\tintercept
gene3\t2.75
gene1\t3.5
gene4\t4.25
gene2\t-0.5
",
    )
    .unwrap();
}

fn write_observation_weight_fixture(path: &Path) {
    let text = fs::read_to_string(reference_data_path("observation_weights.tsv")).unwrap();
    fs::write(
        path,
        text.replace("gene2\t1\t1\t0\t0", "gene2\t1\t1\t0.25\t0.5"),
    )
    .unwrap();
}

fn write_standard_contrast_design_fixture(path: &Path) {
    fs::write(
        path,
        "\
sample\tIntercept\tcondition_B_vs_A
sample3\t1\t1
sample1\t1\t0
sample4\t1\t1
sample2\t1\t0
",
    )
    .unwrap();
}

fn write_additive_contrast_design_fixture(path: &Path) {
    fs::write(
        path,
        "\
sample\tIntercept\tcondition_B_vs_A\tbatch_Y_vs_X
sample3\t1\t1\t1
sample1\t1\t0\t0
sample4\t1\t1\t0
sample2\t1\t0\t1
",
    )
    .unwrap();
}

fn write_mismatched_standard_contrast_design_fixture(path: &Path) {
    fs::write(
        path,
        "\
sample\tIntercept\tcondition_B_vs_A
sample3\t1\t0
sample1\t1\t0
sample4\t1\t1
sample2\t1\t1
",
    )
    .unwrap();
}

fn write_expanded_factor_design_fixture(path: &Path) {
    fs::write(
        path,
        "\
sample\tIntercept\tconditionA\tconditionB
sample3\t1\t0\t1
sample1\t1\t1\t0
sample4\t1\t0\t1
sample2\t1\t1\t0
",
    )
    .unwrap();
}

fn write_gene_numeric_fixture(path: &Path, column: &str, values: &[f64]) {
    let genes = ["gene1", "gene2", "gene3", "gene4"];
    let mut text = format!("gene\t{column}\n");
    for (gene, value) in genes.iter().zip(values.iter()) {
        text.push_str(&format!("{gene}\t{value}\n"));
    }
    fs::write(path, text).unwrap();
}

fn write_ambiguous_coefficient_design_fixture(path: &Path) {
    fs::write(
        path,
        "\
sample\t.Intercept.\tIntercept
sample3\t1\t1
sample1\t1\t0
sample4\t1\t1
sample2\t1\t0
",
    )
    .unwrap();
}

fn write_intercept_design_fixture(path: &Path) {
    fs::write(
        path,
        "\
sample\tIntercept
sample3\t1
sample1\t1
sample4\t1
sample2\t1
",
    )
    .unwrap();
}

fn write_sample_level_fixture(path: &Path) {
    fs::write(
        path,
        "\
sample\tcondition
sample3\tB
sample1\tA
sample4\tB
sample2\tA
",
    )
    .unwrap();
}

fn write_batch_sample_level_fixture(path: &Path) {
    fs::write(
        path,
        "\
sample\tbatch
sample3\tY
sample1\tX
sample4\tX
sample2\tY
",
    )
    .unwrap();
}

fn write_cell_type_contrast_design_fixture(path: &Path) {
    fs::write(
        path,
        "\
sample\tIntercept\tcell.type_B.cell_vs_A.cell
sample3\t1\t1
sample1\t1\t0
sample4\t1\t1
sample2\t1\t0
",
    )
    .unwrap();
}

fn write_ambiguous_cell_type_sample_level_fixture(path: &Path) {
    fs::write(
        path,
        "\
sample\tcell.type
sample3\tB cell
sample1\tA cell
sample4\tB cell
sample2\tA-cell
",
    )
    .unwrap();
}

fn run_cli(args: &[&str]) {
    let status = Command::new(env!("CARGO_BIN_EXE_rsdeseq2"))
        .args(args)
        .status()
        .unwrap();
    assert!(status.success(), "CLI exited with status {status}");
}

fn run_cli_failure_contains(args: &[&str], needle: &str) {
    let output = Command::new(env!("CARGO_BIN_EXE_rsdeseq2"))
        .args(args)
        .output()
        .unwrap();
    assert!(
        !output.status.success(),
        "CLI unexpectedly succeeded with stdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains(needle),
        "expected stderr to contain {needle:?}, got:\n{stderr}"
    );
}

fn assert_deseq_results_table(path: &Path) {
    let table = fs::read_to_string(path).unwrap();
    let mut lines = table.lines();
    let header = lines.next().unwrap();
    assert!(
        header
            == "gene\tbaseMean\tlog2FoldChange\tlfcSE\tstat\tpvalue\tpadj\tdispersion\tconverged"
            || header
                == "gene\tbaseMean\tlog2FoldChange\tlfcSE\tstat\tpvalue\tpadj\tdispersion\tconverged\tfiltered",
        "unexpected result header: {header}"
    );

    let rows = lines.collect::<Vec<_>>();
    assert_eq!(rows.len(), 4);
    for (i, row) in rows.iter().enumerate() {
        let fields = row.split('\t').collect::<Vec<_>>();
        assert_eq!(fields.len(), header.split('\t').count());
        assert_eq!(fields[0], format!("gene{}", i + 1));
        for field in &fields[1..8] {
            if *field != "NA" {
                assert!(
                    field.parse::<f64>().unwrap().is_finite(),
                    "non-finite result {field}"
                );
            }
        }
    }
}

fn assert_tsv_starts_with(path: &Path, header: &str) {
    let table = fs::read_to_string(path).unwrap();
    assert!(
        table.starts_with(header),
        "unexpected TSV header for {}:\n{table}",
        path.display()
    );
}

fn assert_tsv_contains(path: &Path, needle: &str) {
    let table = fs::read_to_string(path).unwrap();
    assert!(
        table.contains(needle),
        "expected {} to contain {needle:?}, got:\n{table}",
        path.display()
    );
}

fn assert_matrix_table(path: &Path) {
    let table = fs::read_to_string(path).unwrap();
    let mut lines = table.lines();
    assert_eq!(
        lines.next().unwrap(),
        "gene\tsample1\tsample2\tsample3\tsample4"
    );
    let rows = lines.collect::<Vec<_>>();
    assert_eq!(rows.len(), 4);
    for row in rows {
        let fields = row.split('\t').collect::<Vec<_>>();
        assert_eq!(fields.len(), 5);
        for field in fields.iter().skip(1) {
            let value = field.parse::<f64>().unwrap();
            assert!(value.is_finite(), "non-finite matrix value {field}");
        }
    }
}

#[test]
fn cli_rlog_accepts_frozen_intercepts_in_blind_mode() {
    let dir = temp_dir("rlog-frozen-blind");
    let frozen_intercept = dir.join("frozen_intercept.tsv");
    let output = dir.join("rlog.tsv");
    write_frozen_intercept_fixture(&frozen_intercept);

    run_cli(&[
        "rlog",
        "--counts",
        reference_data_path("counts.tsv").to_str().unwrap(),
        "--frozen-intercept",
        frozen_intercept.to_str().unwrap(),
        "--rlog-prior-variance",
        "2.5",
        "--fit-type",
        "mean",
        "--output",
        output.to_str().unwrap(),
    ]);

    assert_matrix_table(&output);
}

#[test]
fn cli_rlog_accepts_frozen_intercepts_in_design_aware_mode() {
    let dir = temp_dir("rlog-frozen-design");
    let frozen_intercept = dir.join("frozen_intercept.tsv");
    let output = dir.join("rlog.tsv");
    write_frozen_intercept_fixture(&frozen_intercept);

    run_cli(&[
        "rlog",
        "--counts",
        reference_data_path("counts.tsv").to_str().unwrap(),
        "--design",
        reference_data_path("design_full.tsv").to_str().unwrap(),
        "--blind=false",
        "--frozen-intercept",
        frozen_intercept.to_str().unwrap(),
        "--rlog-prior-variance",
        "2.5",
        "--fit-type",
        "mean",
        "--output",
        output.to_str().unwrap(),
    ]);

    assert_matrix_table(&output);
}

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
}

#[test]
fn cli_wald_accepts_lfc_threshold_alternative() {
    let dir = temp_dir("wald-threshold");
    let output = dir.join("wald.tsv");

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
        "--lfc-threshold",
        "0.5",
        "--alternative",
        "greater",
        "--output",
        output.to_str().unwrap(),
    ]);

    assert_deseq_results_table(&output);
}

#[test]
fn cli_wald_accepts_expanded_beta_prior_coefficient_workflow() {
    let dir = temp_dir("wald-beta-prior-expanded");
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

    run_cli(&[
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
        "--coefficient-name",
        "condition_B_vs_A",
        "--output",
        output.to_str().unwrap(),
    ]);

    assert_deseq_results_table(&output);
}

#[test]
fn cli_wald_accepts_expanded_beta_prior_contrast_workflow() {
    let dir = temp_dir("wald-beta-prior-expanded-contrast");
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

    run_cli(&[
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
        "--contrast-name",
        "condition_B_vs_A",
        "--output",
        output.to_str().unwrap(),
    ]);

    assert_deseq_results_table(&output);
}

#[test]
fn cli_wald_accepts_expanded_beta_prior_normalization_factors() {
    let dir = temp_dir("wald-beta-prior-expanded-nf");
    let standard = dir.join("standard.tsv");
    let expanded = dir.join("expanded.tsv");
    let dispersions = dir.join("dispersions.tsv");
    let base_mean = dir.join("base_mean.tsv");
    let disp_fit = dir.join("disp_fit.tsv");
    let output = dir.join("wald.tsv");
    write_standard_contrast_design_fixture(&standard);
    write_expanded_factor_design_fixture(&expanded);
    write_gene_numeric_fixture(&dispersions, "dispersion", &[0.08, 0.12, 0.05, 0.2]);
    write_gene_numeric_fixture(&base_mean, "baseMean", &[12.375, 2.5, 60.9375, 7.03125]);
    write_gene_numeric_fixture(&disp_fit, "dispFit", &[0.08, 0.12, 0.05, 0.2]);

    run_cli(&[
        "wald",
        "--counts",
        reference_data_path("counts.tsv").to_str().unwrap(),
        "--design",
        standard.to_str().unwrap(),
        "--fit-type",
        "mean",
        "--normalization-factors",
        reference_data_path("normalization_factors.tsv")
            .to_str()
            .unwrap(),
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
        "--coefficient-name",
        "condition_B_vs_A",
        "--output",
        output.to_str().unwrap(),
    ]);

    assert_deseq_results_table(&output);
}

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
fn cli_wald_accepts_factor_expanded_beta_prior_coefficient_workflow() {
    let dir = temp_dir("wald-beta-prior-factor");
    let standard = dir.join("standard.tsv");
    let levels = dir.join("levels.tsv");
    let dispersions = dir.join("dispersions.tsv");
    let base_mean = dir.join("base_mean.tsv");
    let disp_fit = dir.join("disp_fit.tsv");
    let output = dir.join("wald.tsv");
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
    ]);

    assert_deseq_results_table(&output);
}

#[test]
fn cli_wald_accepts_factor_expanded_beta_prior_contrast_workflow() {
    let dir = temp_dir("wald-beta-prior-factor-contrast");
    let standard = dir.join("standard.tsv");
    let levels = dir.join("levels.tsv");
    let dispersions = dir.join("dispersions.tsv");
    let base_mean = dir.join("base_mean.tsv");
    let disp_fit = dir.join("disp_fit.tsv");
    let output = dir.join("wald.tsv");
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
        "--contrast-name",
        "condition_B_vs_A",
        "--output",
        output.to_str().unwrap(),
    ]);

    assert_deseq_results_table(&output);
}

#[test]
fn cli_wald_accepts_factor_expanded_beta_prior_normalization_factors() {
    let dir = temp_dir("wald-beta-prior-factor-nf");
    let standard = dir.join("standard.tsv");
    let levels = dir.join("levels.tsv");
    let dispersions = dir.join("dispersions.tsv");
    let base_mean = dir.join("base_mean.tsv");
    let disp_fit = dir.join("disp_fit.tsv");
    let output = dir.join("wald.tsv");
    write_standard_contrast_design_fixture(&standard);
    write_sample_level_fixture(&levels);
    write_gene_numeric_fixture(&dispersions, "dispersion", &[0.08, 0.12, 0.05, 0.2]);
    write_gene_numeric_fixture(&base_mean, "baseMean", &[12.375, 2.5, 60.9375, 7.03125]);
    write_gene_numeric_fixture(&disp_fit, "dispFit", &[0.08, 0.12, 0.05, 0.2]);

    run_cli(&[
        "wald",
        "--counts",
        reference_data_path("counts.tsv").to_str().unwrap(),
        "--design",
        standard.to_str().unwrap(),
        "--fit-type",
        "mean",
        "--normalization-factors",
        reference_data_path("normalization_factors.tsv")
            .to_str()
            .unwrap(),
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
    ]);

    assert_deseq_results_table(&output);
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
fn cli_wald_accepts_additive_expanded_beta_prior_coefficient_workflow() {
    let dir = temp_dir("wald-beta-prior-additive");
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

    run_cli(&[
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
    ]);

    assert_deseq_results_table(&output);
}

#[test]
fn cli_wald_accepts_additive_expanded_beta_prior_normalization_factors() {
    let dir = temp_dir("wald-beta-prior-additive-nf");
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
    write_gene_numeric_fixture(&base_mean, "baseMean", &[12.375, 2.5, 60.9375, 7.03125]);
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
        "--normalization-factors",
        reference_data_path("normalization_factors.tsv")
            .to_str()
            .unwrap(),
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
        "--contrast-name",
        "batch_Y_vs_X",
        "--output",
        output.to_str().unwrap(),
    ]);

    assert_deseq_results_table(&output);
}

#[test]
fn cli_wald_accepts_additive_expanded_beta_prior_observation_weights() {
    let dir = temp_dir("wald-beta-prior-additive-weights");
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

    run_cli(&[
        "wald",
        "--counts",
        reference_data_path("counts.tsv").to_str().unwrap(),
        "--design",
        standard.to_str().unwrap(),
        "--fit-type",
        "mean",
        "--observation-weights",
        reference_data_path("observation_weights.tsv")
            .to_str()
            .unwrap(),
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
    ]);

    assert_deseq_results_table(&output);
}

#[test]
fn cli_wald_accepts_additive_expanded_beta_prior_list_contrast() {
    let dir = temp_dir("wald-beta-prior-additive-list");
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

    run_cli(&[
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
        "--contrast-positive",
        "condition_B_vs_A,batch_Y_vs_X",
        "--contrast-negative-weight",
        "-0.5",
        "--output",
        output.to_str().unwrap(),
    ]);

    assert_deseq_results_table(&output);
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
fn cli_wald_rejects_ambiguous_r_cleaned_coefficient_name_alias() {
    let dir = temp_dir("wald-coefficient-name-ambiguous");
    let design = dir.join("design.tsv");
    let output = dir.join("wald.tsv");
    write_ambiguous_coefficient_design_fixture(&design);

    run_cli_failure_contains(
        &[
            "wald",
            "--counts",
            reference_data_path("counts.tsv").to_str().unwrap(),
            "--design",
            design.to_str().unwrap(),
            "--fit-type",
            "mean",
            "--coefficient-name",
            "(Intercept)",
            "--output",
            output.to_str().unwrap(),
        ],
        "invalid options: coefficient '(Intercept)' resolves ambiguously after R-style cleanup",
    );
}

#[test]
fn cli_wald_rejects_factor_level_contrast_without_sample_levels() {
    let dir = temp_dir("wald-contrast-factor-missing-levels");
    let design = dir.join("design.tsv");
    let output = dir.join("wald.tsv");
    write_standard_contrast_design_fixture(&design);

    run_cli_failure_contains(
        &[
            "wald",
            "--counts",
            reference_data_path("counts.tsv").to_str().unwrap(),
            "--design",
            design.to_str().unwrap(),
            "--fit-type",
            "mean",
            "--contrast-factor",
            "condition",
            "--contrast-numerator",
            "B",
            "--contrast-denominator",
            "A",
            "--output",
            output.to_str().unwrap(),
        ],
        "invalid options: factor-level results contrast requires --contrast-sample-levels",
    );
}

#[test]
fn cli_wald_rejects_ambiguous_cleaned_factor_level_alias() {
    let dir = temp_dir("wald-contrast-factor-ambiguous-level");
    let design = dir.join("design.tsv");
    let levels = dir.join("levels.tsv");
    let output = dir.join("wald.tsv");
    write_cell_type_contrast_design_fixture(&design);
    write_ambiguous_cell_type_sample_level_fixture(&levels);

    run_cli_failure_contains(
        &[
            "wald",
            "--counts",
            reference_data_path("counts.tsv").to_str().unwrap(),
            "--design",
            design.to_str().unwrap(),
            "--fit-type",
            "mean",
            "--contrast-factor",
            "cell.type",
            "--contrast-numerator",
            "B.cell",
            "--contrast-denominator",
            "A.cell",
            "--contrast-sample-levels",
            levels.to_str().unwrap(),
            "--output",
            output.to_str().unwrap(),
        ],
        "denominator level 'A.cell' resolves ambiguously",
    );
}

#[test]
fn cli_wald_accepts_t_pvalue_options() {
    let dir = temp_dir("wald-t");
    let residual_output = dir.join("wald_t_residual.tsv");
    let scalar_output = dir.join("wald_t_scalar.tsv");
    let per_gene_df = dir.join("wald_t_df.tsv");
    let per_gene_output = dir.join("wald_t_per_gene.tsv");
    write_wald_t_degrees_of_freedom_fixture(&per_gene_df);

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
        "--use-t",
        "--output",
        residual_output.to_str().unwrap(),
    ]);
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
        "--t-degrees-of-freedom",
        "4",
        "--output",
        scalar_output.to_str().unwrap(),
    ]);
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
        "--t-degrees-of-freedom-file",
        per_gene_df.to_str().unwrap(),
        "--output",
        per_gene_output.to_str().unwrap(),
    ]);

    assert_deseq_results_table(&residual_output);
    assert_deseq_results_table(&scalar_output);
    assert_deseq_results_table(&per_gene_output);
}

#[test]
fn cli_wald_accepts_result_filter_controls() {
    let dir = temp_dir("wald-result-controls");
    let output = dir.join("wald.tsv");

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
        "10",
        "--disable-independent-filtering",
        "--independent-filtering-alpha",
        "0.05",
        "--independent-filtering-theta",
        "0,0.5,1",
        "--output",
        output.to_str().unwrap(),
    ]);

    assert_deseq_results_table(&output);
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

#[test]
fn cli_lrt_rejects_ambiguous_r_cleaned_coefficient_name_alias() {
    let dir = temp_dir("lrt-coefficient-name-ambiguous");
    let design = dir.join("design.tsv");
    let reduced_design = dir.join("reduced.tsv");
    let output = dir.join("lrt.tsv");
    write_ambiguous_coefficient_design_fixture(&design);
    write_intercept_design_fixture(&reduced_design);

    run_cli_failure_contains(
        &[
            "lrt",
            "--counts",
            reference_data_path("counts.tsv").to_str().unwrap(),
            "--design",
            design.to_str().unwrap(),
            "--reduced-design",
            reduced_design.to_str().unwrap(),
            "--fit-type",
            "mean",
            "--coefficient-name",
            "(Intercept)",
            "--output",
            output.to_str().unwrap(),
        ],
        "invalid options: coefficient '(Intercept)' resolves ambiguously after R-style cleanup",
    );
}

#[test]
fn cli_lrt_rejects_factor_level_contrast_without_sample_levels() {
    let dir = temp_dir("lrt-contrast-factor-missing-levels");
    let design = dir.join("design.tsv");
    let output = dir.join("lrt.tsv");
    write_standard_contrast_design_fixture(&design);

    run_cli_failure_contains(
        &[
            "lrt",
            "--counts",
            reference_data_path("counts.tsv").to_str().unwrap(),
            "--design",
            design.to_str().unwrap(),
            "--reduced-design",
            reference_data_path("design_reduced.tsv").to_str().unwrap(),
            "--fit-type",
            "mean",
            "--contrast-factor",
            "condition",
            "--contrast-numerator",
            "B",
            "--contrast-denominator",
            "A",
            "--output",
            output.to_str().unwrap(),
        ],
        "invalid options: factor-level results contrast requires --contrast-sample-levels",
    );
}

#[test]
fn cli_lrt_rejects_ambiguous_cleaned_factor_level_alias() {
    let dir = temp_dir("lrt-contrast-factor-ambiguous-level");
    let design = dir.join("design.tsv");
    let reduced_design = dir.join("reduced.tsv");
    let levels = dir.join("levels.tsv");
    let output = dir.join("lrt.tsv");
    write_cell_type_contrast_design_fixture(&design);
    write_intercept_design_fixture(&reduced_design);
    write_ambiguous_cell_type_sample_level_fixture(&levels);

    run_cli_failure_contains(
        &[
            "lrt",
            "--counts",
            reference_data_path("counts.tsv").to_str().unwrap(),
            "--design",
            design.to_str().unwrap(),
            "--reduced-design",
            reduced_design.to_str().unwrap(),
            "--fit-type",
            "mean",
            "--contrast-factor",
            "cell.type",
            "--contrast-numerator",
            "B.cell",
            "--contrast-denominator",
            "A.cell",
            "--contrast-sample-levels",
            levels.to_str().unwrap(),
            "--output",
            output.to_str().unwrap(),
        ],
        "denominator level 'A.cell' resolves ambiguously",
    );
}

#[test]
fn cli_lrt_rejects_sample_levels_without_factor_level_contrast() {
    let dir = temp_dir("lrt-sample-levels-without-factor");
    let levels = dir.join("levels.tsv");
    let output = dir.join("lrt.tsv");
    write_sample_level_fixture(&levels);

    run_cli_failure_contains(
        &[
            "lrt",
            "--counts",
            reference_data_path("counts.tsv").to_str().unwrap(),
            "--design",
            reference_data_path("design_full.tsv").to_str().unwrap(),
            "--reduced-design",
            reference_data_path("design_reduced.tsv").to_str().unwrap(),
            "--fit-type",
            "mean",
            "--contrast-sample-levels",
            levels.to_str().unwrap(),
            "--output",
            output.to_str().unwrap(),
        ],
        "invalid dimensions for factor-level contrast inputs",
    );
}

#[test]
fn cli_lrt_accepts_result_filter_controls() {
    let dir = temp_dir("lrt-result-controls");
    let output = dir.join("lrt.tsv");

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
        "--disable-cooks-cutoff",
        "--disable-independent-filtering",
        "--independent-filtering-theta",
        "0,0.5,1",
        "--output",
        output.to_str().unwrap(),
    ]);

    assert_deseq_results_table(&output);
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

fn assert_numeric_tsv(path: &Path, expected: &[f64]) {
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(path)
        .unwrap();
    let values = reader
        .records()
        .flat_map(|r| {
            r.unwrap()
                .iter()
                .skip(1)
                .map(|s| s.parse::<f64>().unwrap())
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    assert_eq!(values.len(), expected.len());
    for (a, b) in values.iter().zip(expected) {
        common::assert_float_close(*a, *b, 1e-12, 1e-12, "CLI numerical output");
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

#[path = "cli/contrast_options.rs"]
mod contrast_options;
#[path = "cli/normalization_options.rs"]
mod normalization_options;
#[path = "cli/numerical_assertions.rs"]
mod numerical_assertions;
