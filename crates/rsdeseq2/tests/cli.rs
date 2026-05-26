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
sample2\t1
sample3\t1
sample4\t1
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

fn run_cli(args: &[&str]) {
    let status = Command::new(env!("CARGO_BIN_EXE_rsdeseq2"))
        .args(args)
        .status()
        .unwrap();
    assert!(status.success(), "CLI exited with status {status}");
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
