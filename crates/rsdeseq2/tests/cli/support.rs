use super::*;

static TEMP_COUNTER: AtomicUsize = AtomicUsize::new(0);

pub(super) fn temp_dir(label: &str) -> PathBuf {
    let id = TEMP_COUNTER.fetch_add(1, Ordering::Relaxed);
    let path =
        std::env::temp_dir().join(format!("rsdeseq2-cli-{label}-{}-{id}", std::process::id()));
    fs::create_dir_all(&path).unwrap();
    path
}

pub(super) fn reference_data_path(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/deseq2_reference")
        .join(name)
}

pub(super) fn write_size_factor_fixture(path: &Path) {
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

pub(super) fn write_frozen_intercept_fixture(path: &Path) {
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

pub(super) fn write_observation_weight_fixture(path: &Path) {
    let text = fs::read_to_string(reference_data_path("observation_weights.tsv")).unwrap();
    fs::write(
        path,
        text.replace("gene2\t1\t1\t0\t0", "gene2\t1\t1\t0.25\t0.5"),
    )
    .unwrap();
}

pub(super) fn write_standard_contrast_design_fixture(path: &Path) {
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

pub(super) fn write_additive_contrast_design_fixture(path: &Path) {
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

pub(super) fn write_mismatched_standard_contrast_design_fixture(path: &Path) {
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

pub(super) fn write_expanded_factor_design_fixture(path: &Path) {
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

pub(super) fn write_gene_numeric_fixture(path: &Path, column: &str, values: &[f64]) {
    let genes = ["gene1", "gene2", "gene3", "gene4"];
    let mut text = format!("gene\t{column}\n");
    for (gene, value) in genes.iter().zip(values.iter()) {
        text.push_str(&format!("{gene}\t{value}\n"));
    }
    fs::write(path, text).unwrap();
}

pub(super) fn write_ambiguous_coefficient_design_fixture(path: &Path) {
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

pub(super) fn write_intercept_design_fixture(path: &Path) {
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

pub(super) fn write_sample_level_fixture(path: &Path) {
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

pub(super) fn write_batch_sample_level_fixture(path: &Path) {
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

pub(super) fn write_cell_type_contrast_design_fixture(path: &Path) {
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

pub(super) fn write_ambiguous_cell_type_sample_level_fixture(path: &Path) {
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

pub(super) fn run_cli(args: &[&str]) {
    let status = Command::new(env!("CARGO_BIN_EXE_rsdeseq2"))
        .args(args)
        .status()
        .unwrap();
    assert!(status.success(), "CLI exited with status {status}");
}

pub(super) fn run_cli_failure_contains(args: &[&str], needle: &str) {
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

pub(super) fn assert_deseq_results_table(path: &Path) {
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

pub(super) fn assert_tsv_starts_with(path: &Path, header: &str) {
    let table = fs::read_to_string(path).unwrap();
    assert!(
        table.starts_with(header),
        "unexpected TSV header for {}:\n{table}",
        path.display()
    );
}

pub(super) fn assert_tsv_contains(path: &Path, needle: &str) {
    let table = fs::read_to_string(path).unwrap();
    assert!(
        table.contains(needle),
        "expected {} to contain {needle:?}, got:\n{table}",
        path.display()
    );
}

pub(super) fn assert_matrix_table(path: &Path) {
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

pub(super) fn assert_numeric_tsv(path: &Path, expected: &[f64]) {
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
