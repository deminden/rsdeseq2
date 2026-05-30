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

fn write_geometric_mean_fixture(path: &Path) {
    fs::write(
        path,
        "\
gene\tgeo_mean
gene1\t3
gene2\t4
gene3\t6
gene4\t8
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

fn write_unit_observation_weight_fixture(path: &Path) {
    fs::write(
        path,
        "\
gene\tsample1\tsample2\tsample3\tsample4
gene1\t1\t1\t1\t1
gene2\t1\t1\t1\t1
gene3\t1\t1\t1\t1
gene4\t1\t1\t1\t1
",
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

fn write_reserved_contrast_design_fixture(path: &Path) {
    fs::write(
        path,
        "\
sample\tIntercept\tcondition_if._vs_TRUE.
sample3\t1\t1
sample1\t1\t0
sample4\t1\t1
sample2\t1\t0
",
    )
    .unwrap();
}

fn write_reserved_coefficient_design_fixture(path: &Path) {
    fs::write(
        path,
        "\
sample\tIntercept\tif.
sample3\t1\t1
sample1\t1\t0
sample4\t1\t1
sample2\t1\t0
",
    )
    .unwrap();
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

fn write_reserved_sample_level_fixture(path: &Path) {
    fs::write(
        path,
        "\
sample\tcondition
sample3\tif
sample1\tTRUE
sample4\tif
sample2\tTRUE
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

fn run_cli_failure(args: &[&str]) {
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
}

fn assert_deseq_results_table(path: &Path) {
    let table = fs::read_to_string(path).unwrap();
    let mut lines = table.lines();
    let header = lines.next().unwrap();
    assert!(
        header == "gene\tbaseMean\tlog2FoldChange\tlfcSE\tstat\tpvalue\tpadj\tdispersion\tconverged"
            || header
                == "gene\tbaseMean\tlog2FoldChange\tlfcSE\tstat\tpvalue\tpadj\tdispersion\tconverged\tfiltered",
        "unexpected result header: {header}"
    );

    let rows = lines.collect::<Vec<_>>();
    assert_eq!(rows.len(), 4);
    assert!(rows.iter().any(|row| row.starts_with("gene1\t")));
    assert!(rows.iter().any(|row| row.starts_with("gene4\t")));
}

fn assert_base_mean_table(path: &Path) {
    let table = fs::read_to_string(path).unwrap();
    let mut lines = table.lines();
    assert_eq!(lines.next().unwrap(), "gene\tbase_mean");
    assert_eq!(lines.count(), 4);
}

fn assert_normalized_counts_table(path: &Path) {
    let table = fs::read_to_string(path).unwrap();
    let mut lines = table.lines();
    assert_eq!(
        lines.next().unwrap(),
        "gene\tsample1\tsample2\tsample3\tsample4"
    );
    assert_eq!(lines.count(), 4);
}

fn assert_size_factor_table(path: &Path) {
    let table = fs::read_to_string(path).unwrap();
    let mut lines = table.lines();
    assert_eq!(lines.next().unwrap(), "sample\tsize_factor");
    assert_eq!(lines.count(), 4);
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
fn cli_size_factors_accepts_control_genes() {
    let dir = temp_dir("size-factors-control");
    let output = dir.join("size_factors.tsv");

    run_cli(&[
        "size-factors",
        "--counts",
        reference_data_path("counts.tsv").to_str().unwrap(),
        "--control-genes",
        "0,2",
        "--output",
        output.to_str().unwrap(),
    ]);

    assert_size_factor_table(&output);
}

#[test]
fn cli_size_factors_accepts_geometric_means() {
    let dir = temp_dir("size-factors-geo");
    let geometric_means = dir.join("geometric_means.tsv");
    let output = dir.join("size_factors.tsv");
    write_geometric_mean_fixture(&geometric_means);

    run_cli(&[
        "size-factors",
        "--counts",
        reference_data_path("counts.tsv").to_str().unwrap(),
        "--geometric-means",
        geometric_means.to_str().unwrap(),
        "--output",
        output.to_str().unwrap(),
    ]);

    assert_size_factor_table(&output);
}

#[test]
fn cli_base_mean_accepts_normalization_factors() {
    let dir = temp_dir("base-mean-nf");
    let output = dir.join("base_mean.tsv");

    run_cli(&[
        "base-mean",
        "--counts",
        reference_data_path("counts.tsv").to_str().unwrap(),
        "--normalization-factors",
        reference_data_path("normalization_factors.tsv")
            .to_str()
            .unwrap(),
        "--output",
        output.to_str().unwrap(),
    ]);

    assert_base_mean_table(&output);
}

#[test]
fn cli_base_mean_accepts_size_factors() {
    let dir = temp_dir("base-mean-sf");
    let size_factors = dir.join("size_factors.tsv");
    let output = dir.join("base_mean.tsv");
    write_size_factor_fixture(&size_factors);

    run_cli(&[
        "base-mean",
        "--counts",
        reference_data_path("counts.tsv").to_str().unwrap(),
        "--size-factors",
        size_factors.to_str().unwrap(),
        "--output",
        output.to_str().unwrap(),
    ]);

    assert_base_mean_table(&output);
}

#[test]
fn cli_base_mean_accepts_observation_weights() {
    let dir = temp_dir("base-mean-weights");
    let output = dir.join("base_mean.tsv");

    run_cli(&[
        "base-mean",
        "--counts",
        reference_data_path("counts.tsv").to_str().unwrap(),
        "--observation-weights",
        reference_data_path("observation_weights.tsv")
            .to_str()
            .unwrap(),
        "--output",
        output.to_str().unwrap(),
    ]);

    assert_base_mean_table(&output);
}

#[test]
fn cli_base_mean_accepts_control_genes() {
    let dir = temp_dir("base-mean-control");
    let output = dir.join("base_mean.tsv");

    run_cli(&[
        "base-mean",
        "--counts",
        reference_data_path("counts.tsv").to_str().unwrap(),
        "--control-genes",
        "0,2",
        "--output",
        output.to_str().unwrap(),
    ]);

    assert_base_mean_table(&output);
}

#[test]
fn cli_base_mean_accepts_geometric_means() {
    let dir = temp_dir("base-mean-geo");
    let geometric_means = dir.join("geometric_means.tsv");
    let output = dir.join("base_mean.tsv");
    write_geometric_mean_fixture(&geometric_means);

    run_cli(&[
        "base-mean",
        "--counts",
        reference_data_path("counts.tsv").to_str().unwrap(),
        "--geometric-means",
        geometric_means.to_str().unwrap(),
        "--output",
        output.to_str().unwrap(),
    ]);

    assert_base_mean_table(&output);
}

#[test]
fn cli_normalized_counts_accepts_size_factors() {
    let dir = temp_dir("normalized-sf");
    let size_factors = dir.join("size_factors.tsv");
    let output = dir.join("normalized.tsv");
    write_size_factor_fixture(&size_factors);

    run_cli(&[
        "normalized-counts",
        "--counts",
        reference_data_path("counts.tsv").to_str().unwrap(),
        "--size-factors",
        size_factors.to_str().unwrap(),
        "--output",
        output.to_str().unwrap(),
    ]);

    assert_normalized_counts_table(&output);
}

#[test]
fn cli_normalized_counts_accepts_normalization_factors() {
    let dir = temp_dir("normalized-nf");
    let output = dir.join("normalized.tsv");

    run_cli(&[
        "normalized-counts",
        "--counts",
        reference_data_path("counts.tsv").to_str().unwrap(),
        "--normalization-factors",
        reference_data_path("normalization_factors.tsv")
            .to_str()
            .unwrap(),
        "--output",
        output.to_str().unwrap(),
    ]);

    assert_normalized_counts_table(&output);
}

#[test]
fn cli_normalized_counts_accepts_control_genes() {
    let dir = temp_dir("normalized-control");
    let output = dir.join("normalized.tsv");

    run_cli(&[
        "normalized-counts",
        "--counts",
        reference_data_path("counts.tsv").to_str().unwrap(),
        "--control-genes",
        "0,2",
        "--output",
        output.to_str().unwrap(),
    ]);

    assert_normalized_counts_table(&output);
}

#[test]
fn cli_normalized_counts_accepts_geometric_means() {
    let dir = temp_dir("normalized-geo");
    let geometric_means = dir.join("geometric_means.tsv");
    let output = dir.join("normalized.tsv");
    write_geometric_mean_fixture(&geometric_means);

    run_cli(&[
        "normalized-counts",
        "--counts",
        reference_data_path("counts.tsv").to_str().unwrap(),
        "--geometric-means",
        geometric_means.to_str().unwrap(),
        "--output",
        output.to_str().unwrap(),
    ]);

    assert_normalized_counts_table(&output);
}

#[test]
fn cli_vst_runs_blind_mean_fit() {
    let dir = temp_dir("vst-blind");
    let output = dir.join("vst.tsv");

    run_cli(&[
        "vst",
        "--counts",
        reference_data_path("counts.tsv").to_str().unwrap(),
        "--fit-type",
        "mean",
        "--output",
        output.to_str().unwrap(),
    ]);

    assert_matrix_table(&output);
}

#[test]
fn cli_vst_accepts_control_genes() {
    let dir = temp_dir("vst-control");
    let output = dir.join("vst.tsv");

    run_cli(&[
        "vst",
        "--counts",
        reference_data_path("counts.tsv").to_str().unwrap(),
        "--fit-type",
        "mean",
        "--control-genes",
        "0,2",
        "--output",
        output.to_str().unwrap(),
    ]);

    assert_matrix_table(&output);
}

#[test]
fn cli_vst_accepts_geometric_means() {
    let dir = temp_dir("vst-geo");
    let geometric_means = dir.join("geometric_means.tsv");
    let output = dir.join("vst.tsv");
    write_geometric_mean_fixture(&geometric_means);

    run_cli(&[
        "vst",
        "--counts",
        reference_data_path("counts.tsv").to_str().unwrap(),
        "--fit-type",
        "mean",
        "--geometric-means",
        geometric_means.to_str().unwrap(),
        "--output",
        output.to_str().unwrap(),
    ]);

    assert_matrix_table(&output);
}

#[test]
fn cli_vst_runs_design_aware_mean_fit() {
    let dir = temp_dir("vst-design");
    let output = dir.join("vst.tsv");

    run_cli(&[
        "vst",
        "--counts",
        reference_data_path("counts.tsv").to_str().unwrap(),
        "--design",
        reference_data_path("design_full.tsv").to_str().unwrap(),
        "--blind=false",
        "--fit-type",
        "mean",
        "--output",
        output.to_str().unwrap(),
    ]);

    assert_matrix_table(&output);
}

#[test]
fn cli_rlog_runs_blind_mean_fit() {
    let dir = temp_dir("rlog-blind");
    let output = dir.join("rlog.tsv");

    run_cli(&[
        "rlog",
        "--counts",
        reference_data_path("counts.tsv").to_str().unwrap(),
        "--fit-type",
        "mean",
        "--output",
        output.to_str().unwrap(),
    ]);

    assert_matrix_table(&output);
}

#[test]
fn cli_rlog_runs_design_aware_mean_fit() {
    let dir = temp_dir("rlog-design");
    let output = dir.join("rlog.tsv");

    run_cli(&[
        "rlog",
        "--counts",
        reference_data_path("counts.tsv").to_str().unwrap(),
        "--design",
        reference_data_path("design_full.tsv").to_str().unwrap(),
        "--blind=false",
        "--fit-type",
        "mean",
        "--output",
        output.to_str().unwrap(),
    ]);

    assert_matrix_table(&output);
}

#[test]
fn cli_rlog_accepts_normalization_factors() {
    let dir = temp_dir("rlog-nf");
    let output = dir.join("rlog.tsv");

    run_cli(&[
        "rlog",
        "--counts",
        reference_data_path("counts.tsv").to_str().unwrap(),
        "--normalization-factors",
        reference_data_path("normalization_factors.tsv")
            .to_str()
            .unwrap(),
        "--fit-type",
        "mean",
        "--output",
        output.to_str().unwrap(),
    ]);

    assert_matrix_table(&output);
}

#[test]
fn cli_rlog_accepts_size_factors() {
    let dir = temp_dir("rlog-sf");
    let size_factors = dir.join("size_factors.tsv");
    let output = dir.join("rlog.tsv");
    write_size_factor_fixture(&size_factors);

    run_cli(&[
        "rlog",
        "--counts",
        reference_data_path("counts.tsv").to_str().unwrap(),
        "--size-factors",
        size_factors.to_str().unwrap(),
        "--fit-type",
        "mean",
        "--output",
        output.to_str().unwrap(),
    ]);

    assert_matrix_table(&output);
}

#[test]
fn cli_rlog_accepts_control_genes() {
    let dir = temp_dir("rlog-control");
    let output = dir.join("rlog.tsv");

    run_cli(&[
        "rlog",
        "--counts",
        reference_data_path("counts.tsv").to_str().unwrap(),
        "--control-genes",
        "0,2",
        "--fit-type",
        "mean",
        "--output",
        output.to_str().unwrap(),
    ]);

    assert_matrix_table(&output);
}

#[test]
fn cli_rlog_accepts_geometric_means() {
    let dir = temp_dir("rlog-geo");
    let geometric_means = dir.join("geometric_means.tsv");
    let output = dir.join("rlog.tsv");
    write_geometric_mean_fixture(&geometric_means);

    run_cli(&[
        "rlog",
        "--counts",
        reference_data_path("counts.tsv").to_str().unwrap(),
        "--geometric-means",
        geometric_means.to_str().unwrap(),
        "--fit-type",
        "mean",
        "--output",
        output.to_str().unwrap(),
    ]);

    assert_matrix_table(&output);
}

#[test]
fn cli_rlog_accepts_observation_weights() {
    let dir = temp_dir("rlog-weights");
    let weights = dir.join("observation_weights.tsv");
    let output = dir.join("rlog.tsv");
    write_unit_observation_weight_fixture(&weights);

    run_cli(&[
        "rlog",
        "--counts",
        reference_data_path("counts.tsv").to_str().unwrap(),
        "--observation-weights",
        weights.to_str().unwrap(),
        "--fit-type",
        "mean",
        "--output",
        output.to_str().unwrap(),
    ]);

    assert_matrix_table(&output);
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
fn cli_rlog_frozen_intercepts_require_prior_variance() {
    let dir = temp_dir("rlog-frozen-missing-prior");
    let frozen_intercept = dir.join("frozen_intercept.tsv");
    let output = dir.join("rlog.tsv");
    write_frozen_intercept_fixture(&frozen_intercept);

    run_cli_failure(&[
        "rlog",
        "--counts",
        reference_data_path("counts.tsv").to_str().unwrap(),
        "--frozen-intercept",
        frozen_intercept.to_str().unwrap(),
        "--fit-type",
        "mean",
        "--output",
        output.to_str().unwrap(),
    ]);
}

#[test]
fn cli_rlog_prior_variance_requires_frozen_intercepts() {
    let dir = temp_dir("rlog-prior-without-frozen");
    let output = dir.join("rlog.tsv");

    run_cli_failure(&[
        "rlog",
        "--counts",
        reference_data_path("counts.tsv").to_str().unwrap(),
        "--rlog-prior-variance",
        "2.5",
        "--fit-type",
        "mean",
        "--output",
        output.to_str().unwrap(),
    ]);
}

#[test]
fn cli_rlog_requires_design_when_not_blind() {
    let dir = temp_dir("rlog-missing-design");
    let output = dir.join("rlog.tsv");

    run_cli_failure(&[
        "rlog",
        "--counts",
        reference_data_path("counts.tsv").to_str().unwrap(),
        "--blind=false",
        "--fit-type",
        "mean",
        "--output",
        output.to_str().unwrap(),
    ]);
}

#[test]
fn cli_wald_accepts_control_genes() {
    let dir = temp_dir("wald-control");
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
        "--control-genes",
        "0,1,2,3",
        "--output",
        output.to_str().unwrap(),
    ]);

    assert_deseq_results_table(&output);
}

#[test]
fn cli_wald_accepts_geometric_means() {
    let dir = temp_dir("wald-geo");
    let geometric_means = dir.join("geometric_means.tsv");
    let output = dir.join("wald.tsv");
    write_geometric_mean_fixture(&geometric_means);

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
        "--geometric-means",
        geometric_means.to_str().unwrap(),
        "--output",
        output.to_str().unwrap(),
    ]);

    assert_deseq_results_table(&output);
}

#[test]
fn cli_wald_writes_deseq_results_table() {
    let dir = temp_dir("wald");
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
        "--output",
        output.to_str().unwrap(),
    ]);

    assert_deseq_results_table(&output);
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

    run_cli_failure(&[
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
    ]);
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
fn cli_wald_accepts_numeric_contrast() {
    let dir = temp_dir("wald-contrast");
    let output = dir.join("wald.tsv");

    run_cli(&[
        "wald",
        "--counts",
        reference_data_path("counts.tsv").to_str().unwrap(),
        "--design",
        reference_data_path("design_full.tsv").to_str().unwrap(),
        "--fit-type",
        "mean",
        "--contrast",
        "0,1",
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

    run_cli_failure(&[
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

    run_cli_failure(&[
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
    ]);
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

    run_cli_failure(&[
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
}

#[test]
fn cli_wald_rejects_incomplete_expanded_beta_prior_inputs() {
    let dir = temp_dir("wald-beta-prior-expanded-incomplete");
    let standard = dir.join("standard.tsv");
    let expanded = dir.join("expanded.tsv");
    let output = dir.join("wald.tsv");
    write_standard_contrast_design_fixture(&standard);
    write_expanded_factor_design_fixture(&expanded);

    run_cli_failure(&[
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
    ]);
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

    run_cli_failure(&[
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
    ]);
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

    run_cli_failure(&[
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
    ]);
}

#[test]
fn cli_wald_rejects_mixed_beta_prior_matrix_and_factor_inputs() {
    let dir = temp_dir("wald-beta-prior-mixed-inputs");
    let standard = dir.join("standard.tsv");
    let expanded = dir.join("expanded.tsv");
    let levels = dir.join("levels.tsv");
    let dispersions = dir.join("dispersions.tsv");
    let base_mean = dir.join("base_mean.tsv");
    let disp_fit = dir.join("disp_fit.tsv");
    let output = dir.join("wald.tsv");
    write_standard_contrast_design_fixture(&standard);
    write_expanded_factor_design_fixture(&expanded);
    write_sample_level_fixture(&levels);
    write_gene_numeric_fixture(&dispersions, "dispersion", &[0.08, 0.12, 0.05, 0.2]);
    write_gene_numeric_fixture(&base_mean, "baseMean", &[16.5, 3.0, 97.5, 7.5]);
    write_gene_numeric_fixture(&disp_fit, "dispFit", &[0.08, 0.12, 0.05, 0.2]);

    run_cli_failure(&[
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
}

#[test]
fn cli_wald_rejects_mixed_beta_prior_matrix_and_additive_inputs() {
    let dir = temp_dir("wald-beta-prior-matrix-additive-inputs");
    let standard = dir.join("standard.tsv");
    let expanded = dir.join("expanded.tsv");
    let condition_levels = dir.join("condition.tsv");
    let batch_levels = dir.join("batch.tsv");
    let dispersions = dir.join("dispersions.tsv");
    let base_mean = dir.join("base_mean.tsv");
    let disp_fit = dir.join("disp_fit.tsv");
    let output = dir.join("wald.tsv");
    write_additive_contrast_design_fixture(&standard);
    write_expanded_factor_design_fixture(&expanded);
    write_sample_level_fixture(&condition_levels);
    write_batch_sample_level_fixture(&batch_levels);
    write_gene_numeric_fixture(&dispersions, "dispersion", &[0.08, 0.12, 0.05, 0.2]);
    write_gene_numeric_fixture(&base_mean, "baseMean", &[16.5, 3.0, 97.5, 7.5]);
    write_gene_numeric_fixture(&disp_fit, "dispFit", &[0.08, 0.12, 0.05, 0.2]);
    let sample_levels = format!("{},{}", condition_levels.display(), batch_levels.display());

    run_cli_failure(&[
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
}

#[test]
fn cli_wald_rejects_mixed_beta_prior_factor_and_additive_inputs() {
    let dir = temp_dir("wald-beta-prior-factor-additive-inputs");
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

    run_cli_failure(&[
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
        condition_levels.to_str().unwrap(),
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

    run_cli_failure(&[
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
    ]);
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

    run_cli_failure(&[
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
    ]);
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

    run_cli_failure(&[
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
    ]);
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

    run_cli_failure(&[
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
    ]);
}

#[test]
fn cli_wald_accepts_named_contrast() {
    let dir = temp_dir("wald-contrast-name");
    let output = dir.join("wald.tsv");

    run_cli(&[
        "wald",
        "--counts",
        reference_data_path("counts.tsv").to_str().unwrap(),
        "--design",
        reference_data_path("design_full.tsv").to_str().unwrap(),
        "--fit-type",
        "mean",
        "--contrast-name",
        "conditionB",
        "--output",
        output.to_str().unwrap(),
    ]);

    assert_deseq_results_table(&output);
}

#[test]
fn cli_wald_accepts_r_cleaned_named_contrast_alias() {
    let dir = temp_dir("wald-contrast-name-cleaned");
    let design = dir.join("design.tsv");
    let output = dir.join("wald.tsv");
    write_reserved_coefficient_design_fixture(&design);

    run_cli(&[
        "wald",
        "--counts",
        reference_data_path("counts.tsv").to_str().unwrap(),
        "--design",
        design.to_str().unwrap(),
        "--fit-type",
        "mean",
        "--contrast-name",
        "if",
        "--output",
        output.to_str().unwrap(),
    ]);

    assert_deseq_results_table(&output);
}

#[test]
fn cli_wald_accepts_coefficient_name() {
    let dir = temp_dir("wald-coefficient-name");
    let output = dir.join("wald.tsv");

    run_cli(&[
        "wald",
        "--counts",
        reference_data_path("counts.tsv").to_str().unwrap(),
        "--design",
        reference_data_path("design_full.tsv").to_str().unwrap(),
        "--fit-type",
        "mean",
        "--coefficient-name",
        "conditionB",
        "--output",
        output.to_str().unwrap(),
    ]);

    assert_deseq_results_table(&output);
}

#[test]
fn cli_wald_accepts_r_cleaned_coefficient_name_alias() {
    let dir = temp_dir("wald-coefficient-name-cleaned");
    let design = dir.join("design.tsv");
    let output = dir.join("wald.tsv");
    write_reserved_coefficient_design_fixture(&design);

    run_cli(&[
        "wald",
        "--counts",
        reference_data_path("counts.tsv").to_str().unwrap(),
        "--design",
        design.to_str().unwrap(),
        "--fit-type",
        "mean",
        "--coefficient-name",
        "if",
        "--output",
        output.to_str().unwrap(),
    ]);

    assert_deseq_results_table(&output);
}

#[test]
fn cli_wald_rejects_ambiguous_r_cleaned_coefficient_name_alias() {
    let dir = temp_dir("wald-coefficient-name-ambiguous");
    let design = dir.join("design.tsv");
    let output = dir.join("wald.tsv");
    write_ambiguous_coefficient_design_fixture(&design);

    run_cli_failure(&[
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
    ]);
}

#[test]
fn cli_wald_accepts_list_contrast() {
    let dir = temp_dir("wald-contrast-list");
    let output = dir.join("wald.tsv");

    run_cli(&[
        "wald",
        "--counts",
        reference_data_path("counts.tsv").to_str().unwrap(),
        "--design",
        reference_data_path("design_full.tsv").to_str().unwrap(),
        "--fit-type",
        "mean",
        "--contrast-positive",
        "conditionB",
        "--contrast-negative",
        "(Intercept)",
        "--contrast-positive-weight",
        "1",
        "--contrast-negative-weight=-0.5",
        "--output",
        output.to_str().unwrap(),
    ]);

    assert_deseq_results_table(&output);
}

#[test]
fn cli_wald_accepts_factor_level_contrast() {
    let dir = temp_dir("wald-contrast-factor");
    let design = dir.join("design.tsv");
    let levels = dir.join("levels.tsv");
    let output = dir.join("wald.tsv");
    write_standard_contrast_design_fixture(&design);
    write_sample_level_fixture(&levels);

    run_cli(&[
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
        "--contrast-sample-levels",
        levels.to_str().unwrap(),
        "--output",
        output.to_str().unwrap(),
    ]);

    assert_deseq_results_table(&output);
}

#[test]
fn cli_wald_accepts_factor_level_contrast_with_r_cleaned_level_names() {
    let dir = temp_dir("wald-contrast-factor-cleaned");
    let design = dir.join("design.tsv");
    let levels = dir.join("levels.tsv");
    let output = dir.join("wald.tsv");
    write_reserved_contrast_design_fixture(&design);
    write_reserved_sample_level_fixture(&levels);

    run_cli(&[
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
        "if",
        "--contrast-denominator",
        "TRUE",
        "--contrast-sample-levels",
        levels.to_str().unwrap(),
        "--output",
        output.to_str().unwrap(),
    ]);

    assert_deseq_results_table(&output);
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
fn cli_wald_accepts_size_factors() {
    let dir = temp_dir("wald-sf");
    let size_factors = dir.join("size_factors.tsv");
    let output = dir.join("wald.tsv");
    write_size_factor_fixture(&size_factors);

    run_cli(&[
        "wald",
        "--counts",
        reference_data_path("counts.tsv").to_str().unwrap(),
        "--design",
        reference_data_path("design_full.tsv").to_str().unwrap(),
        "--size-factors",
        size_factors.to_str().unwrap(),
        "--fit-type",
        "mean",
        "--coefficient",
        "1",
        "--output",
        output.to_str().unwrap(),
    ]);

    assert_deseq_results_table(&output);
}

#[test]
fn cli_wald_accepts_observation_weights() {
    let dir = temp_dir("wald-weights");
    let output = dir.join("wald.tsv");

    run_cli(&[
        "wald",
        "--counts",
        reference_data_path("counts.tsv").to_str().unwrap(),
        "--design",
        reference_data_path("design_full.tsv").to_str().unwrap(),
        "--observation-weights",
        reference_data_path("observation_weights.tsv")
            .to_str()
            .unwrap(),
        "--fit-type",
        "mean",
        "--coefficient",
        "1",
        "--output",
        output.to_str().unwrap(),
    ]);

    assert_deseq_results_table(&output);
}

#[test]
fn cli_lrt_accepts_control_genes() {
    let dir = temp_dir("lrt-control");
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
        "--control-genes",
        "0,1,2,3",
        "--output",
        output.to_str().unwrap(),
    ]);

    assert_deseq_results_table(&output);
}

#[test]
fn cli_lrt_accepts_geometric_means() {
    let dir = temp_dir("lrt-geo");
    let geometric_means = dir.join("geometric_means.tsv");
    let output = dir.join("lrt.tsv");
    write_geometric_mean_fixture(&geometric_means);

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
        "--geometric-means",
        geometric_means.to_str().unwrap(),
        "--output",
        output.to_str().unwrap(),
    ]);

    assert_deseq_results_table(&output);
}

#[test]
fn cli_lrt_writes_deseq_results_table() {
    let dir = temp_dir("lrt");
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
fn cli_lrt_accepts_coefficient_name() {
    let dir = temp_dir("lrt-coefficient-name");
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
        "--coefficient-name",
        "conditionB",
        "--output",
        output.to_str().unwrap(),
    ]);

    assert_deseq_results_table(&output);
}

#[test]
fn cli_lrt_accepts_r_cleaned_coefficient_name_alias() {
    let dir = temp_dir("lrt-coefficient-name-cleaned");
    let design = dir.join("design.tsv");
    let reduced_design = dir.join("reduced.tsv");
    let output = dir.join("lrt.tsv");
    write_reserved_coefficient_design_fixture(&design);
    write_intercept_design_fixture(&reduced_design);

    run_cli(&[
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
        "if",
        "--output",
        output.to_str().unwrap(),
    ]);

    assert_deseq_results_table(&output);
}

#[test]
fn cli_lrt_rejects_ambiguous_r_cleaned_coefficient_name_alias() {
    let dir = temp_dir("lrt-coefficient-name-ambiguous");
    let design = dir.join("design.tsv");
    let reduced_design = dir.join("reduced.tsv");
    let output = dir.join("lrt.tsv");
    write_ambiguous_coefficient_design_fixture(&design);
    write_intercept_design_fixture(&reduced_design);

    run_cli_failure(&[
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
    ]);
}

#[test]
fn cli_lrt_accepts_numeric_contrast() {
    let dir = temp_dir("lrt-contrast");
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
        "--contrast",
        "0,1",
        "--output",
        output.to_str().unwrap(),
    ]);

    assert_deseq_results_table(&output);
}

#[test]
fn cli_lrt_accepts_named_contrast() {
    let dir = temp_dir("lrt-contrast-name");
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
        "--contrast-name",
        "conditionB",
        "--output",
        output.to_str().unwrap(),
    ]);

    assert_deseq_results_table(&output);
}

#[test]
fn cli_lrt_accepts_r_cleaned_named_contrast_alias() {
    let dir = temp_dir("lrt-contrast-name-cleaned");
    let design = dir.join("design.tsv");
    let reduced_design = dir.join("reduced.tsv");
    let output = dir.join("lrt.tsv");
    write_reserved_coefficient_design_fixture(&design);
    write_intercept_design_fixture(&reduced_design);

    run_cli(&[
        "lrt",
        "--counts",
        reference_data_path("counts.tsv").to_str().unwrap(),
        "--design",
        design.to_str().unwrap(),
        "--reduced-design",
        reduced_design.to_str().unwrap(),
        "--fit-type",
        "mean",
        "--contrast-name",
        "if",
        "--output",
        output.to_str().unwrap(),
    ]);

    assert_deseq_results_table(&output);
}

#[test]
fn cli_lrt_accepts_list_contrast() {
    let dir = temp_dir("lrt-contrast-list");
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
        "--contrast-positive",
        "conditionB",
        "--contrast-negative",
        "(Intercept)",
        "--contrast-positive-weight",
        "1",
        "--contrast-negative-weight=-0.5",
        "--output",
        output.to_str().unwrap(),
    ]);

    assert_deseq_results_table(&output);
}

#[test]
fn cli_lrt_accepts_factor_level_contrast() {
    let dir = temp_dir("lrt-contrast-factor");
    let design = dir.join("design.tsv");
    let levels = dir.join("levels.tsv");
    let output = dir.join("lrt.tsv");
    write_standard_contrast_design_fixture(&design);
    write_sample_level_fixture(&levels);

    run_cli(&[
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
        "--contrast-sample-levels",
        levels.to_str().unwrap(),
        "--output",
        output.to_str().unwrap(),
    ]);

    assert_deseq_results_table(&output);
}

#[test]
fn cli_lrt_accepts_factor_level_contrast_with_r_cleaned_level_names() {
    let dir = temp_dir("lrt-contrast-factor-cleaned");
    let design = dir.join("design.tsv");
    let reduced_design = dir.join("reduced.tsv");
    let levels = dir.join("levels.tsv");
    let output = dir.join("lrt.tsv");
    write_reserved_contrast_design_fixture(&design);
    write_intercept_design_fixture(&reduced_design);
    write_reserved_sample_level_fixture(&levels);

    run_cli(&[
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
        "condition",
        "--contrast-numerator",
        "if",
        "--contrast-denominator",
        "TRUE",
        "--contrast-sample-levels",
        levels.to_str().unwrap(),
        "--output",
        output.to_str().unwrap(),
    ]);

    assert_deseq_results_table(&output);
}

#[test]
fn cli_lrt_rejects_sample_levels_without_factor_level_contrast() {
    let dir = temp_dir("lrt-sample-levels-without-factor");
    let levels = dir.join("levels.tsv");
    let output = dir.join("lrt.tsv");
    write_sample_level_fixture(&levels);

    run_cli_failure(&[
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
    ]);
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
