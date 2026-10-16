//! CLI routing checks use direct library inputs as their oracle. Independent R
//! references test the numerical implementation; these cases test flag wiring.
use super::*;

pub(super) fn invoke(label: &str, args: &[String]) {
    eprintln!("CLI routing case: {label}");
    let output = Command::new(env!("CARGO_BIN_EXE_rsdeseq2"))
        .args(args)
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{label}: CLI failed: {args:?}\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
}

pub(super) fn read_export(path: &Path) -> (Vec<String>, Vec<Vec<String>>) {
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(path)
        .unwrap();
    let headers = reader
        .headers()
        .unwrap()
        .iter()
        .map(str::to_owned)
        .collect();
    let rows = reader
        .records()
        .map(|row| row.unwrap().iter().map(str::to_owned).collect())
        .collect();
    (headers, rows)
}

pub(super) fn number(text: &str) -> Option<f64> {
    if text == "NA" {
        return None;
    }
    let value: f64 = text.parse().expect("numeric CLI cell or NA");
    assert!(value.is_finite(), "non-finite CLI cell: {text}");
    Some(value)
}

pub(super) fn assert_result_export(path: &Path, expected: &DeseqResults, label: &str) {
    let (headers, rows) = read_export(path);
    assert_eq!(
        headers,
        std::iter::once("gene")
            .chain(expected.column_names())
            .collect::<Vec<_>>(),
        "{label}: columns"
    );
    assert_eq!(rows.len(), expected.rows.len(), "{label}: row count");
    for (actual, expected) in rows.iter().zip(&expected.rows) {
        assert_eq!(
            actual[0],
            expected.gene.as_deref().unwrap(),
            "{label}: gene order"
        );
        for (column, cell) in headers.iter().zip(actual).skip(1) {
            let numeric = match column.as_str() {
                "baseMean" => Some(expected.base_mean),
                "log2FoldChange" => expected.log2_fold_change,
                "lfcSE" => expected.lfc_se,
                "stat" => expected.stat,
                "pvalue" => expected.pvalue,
                "padj" => expected.padj,
                "dispersion" => expected.dispersion,
                "maxCooks" => expected.max_cooks,
                logical => {
                    let value = match logical {
                        "converged" => expected.converged,
                        "cooksOutlier" => expected.cooks_outlier,
                        "filtered" => expected.filtered,
                        _ => panic!("unhandled result column {logical}"),
                    };
                    assert_eq!(
                        cell,
                        match value {
                            Some(true) => "TRUE",
                            Some(false) => "FALSE",
                            None => "NA",
                        },
                        "{label}: {logical}"
                    );
                    continue;
                }
            };
            common::assert_option_close(
                number(cell),
                numeric,
                1e-12,
                1e-12,
                &format!("{label}: {} {column}", actual[0]),
            );
        }
    }
}

pub(super) fn result_values(results: &DeseqResults) -> Vec<Option<f64>> {
    results
        .rows
        .iter()
        .flat_map(|r| {
            [
                Some(r.base_mean),
                r.log2_fold_change,
                r.lfc_se,
                r.stat,
                r.pvalue,
                r.padj,
                r.dispersion,
            ]
        })
        .collect()
}

pub(super) fn numerically_distinct(a: &[Option<f64>], b: &[Option<f64>]) -> bool {
    assert_eq!(a.len(), b.len());
    a.iter().zip(b).any(|(a, b)| match (a, b) {
        (Some(a), Some(b)) => (a - b).abs() > 1e-6,
        (None, None) => false,
        _ => true,
    })
}

pub(super) fn test_builder() -> DeseqBuilder {
    DeseqBuilder::new()
        .automatic_dispersion_means()
        .fit_type(FitType::Mean)
        .disable_cooks_cutoff()
        .disable_independent_filtering()
}

pub(super) fn basic_args(command: &str, counts: &Path, output: &Path) -> Vec<String> {
    [
        command,
        "--counts",
        counts.to_str().unwrap(),
        "--output",
        output.to_str().unwrap(),
    ]
    .map(str::to_owned)
    .to_vec()
}

pub(super) fn model_args(args: &mut Vec<String>, command: &str, design: &Path, reduced: &Path) {
    args.extend(["--fit-type", "mean"].map(str::to_owned));
    if matches!(command, "wald" | "lrt") {
        args.extend(
            [
                "--design",
                design.to_str().unwrap(),
                "--disable-cooks-cutoff",
                "--disable-independent-filtering",
            ]
            .map(str::to_owned),
        );
        if command == "lrt" {
            args.extend(["--reduced-design", reduced.to_str().unwrap()].map(str::to_owned));
        }
    }
}

pub(super) fn assert_matrix_export(
    path: &Path,
    counts: &CountMatrix,
    expected: &RowMajorMatrix<f64>,
    label: &str,
) {
    assert_eq!(expected.n_rows(), counts.n_genes(), "{label}: oracle rows");
    assert_eq!(
        expected.n_cols(),
        counts.n_samples(),
        "{label}: oracle columns"
    );
    let (headers, rows) = read_export(path);
    assert_eq!(
        headers,
        std::iter::once("gene")
            .chain(counts.sample_names().unwrap().iter().map(String::as_str))
            .collect::<Vec<_>>(),
        "{label}: columns"
    );
    assert_eq!(rows.len(), counts.n_genes(), "{label}: genes");
    for (i, row) in rows.iter().enumerate() {
        assert_eq!(row.len(), headers.len(), "{label}: row width");
        assert_eq!(
            &row[0],
            &counts.gene_names().unwrap()[i],
            "{label}: gene order"
        );
        for (actual, expected) in row.iter().skip(1).zip(expected.row(i).unwrap()) {
            common::assert_option_close(number(actual), Some(*expected), 1e-12, 1e-12, label);
        }
    }
}
