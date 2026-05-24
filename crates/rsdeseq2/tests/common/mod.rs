#![allow(dead_code)]

use csv::ReaderBuilder;
use rsdeseq2::prelude::*;
use std::collections::HashMap;
use std::path::{Path, PathBuf};

pub type TsvRow = HashMap<String, String>;

pub fn reference_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/deseq2_reference")
}

pub fn optional_reference_file(name: &str) -> Option<PathBuf> {
    let path = reference_dir().join(name);
    if path.exists() {
        Some(path)
    } else {
        eprintln!(
            "skipping optional DESeq2 golden test; generate {} with scripts/generate_deseq2_references.R",
            path.display()
        );
        None
    }
}

pub fn read_optional_tsv(name: &str) -> Option<Vec<TsvRow>> {
    optional_reference_file(name).map(|path| read_tsv(&path))
}

pub fn read_tsv(path: &Path) -> Vec<TsvRow> {
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(path)
        .unwrap_or_else(|error| panic!("failed to open {}: {error}", path.display()));
    let headers = reader
        .headers()
        .unwrap_or_else(|error| panic!("failed to read headers from {}: {error}", path.display()))
        .clone();
    reader
        .records()
        .map(|record| {
            let record = record.unwrap_or_else(|error| {
                panic!("failed to read row from {}: {error}", path.display())
            });
            headers
                .iter()
                .zip(record.iter())
                .map(|(key, value)| (key.to_string(), value.to_string()))
                .collect::<TsvRow>()
        })
        .collect()
}

pub fn reference_counts() -> CountMatrix {
    CountMatrix::from_row_major_u32_with_names(
        4,
        4,
        vec![
            10, 12, 20, 24, //
            0, 0, 5, 7, //
            100, 80, 90, 120, //
            3, 6, 9, 12,
        ],
        Some(reference_gene_names()),
        Some(reference_sample_names()),
    )
    .unwrap()
}

pub fn reference_full_design() -> DesignMatrix {
    DesignMatrix::from_row_major(
        4,
        2,
        vec![
            1.0, 0.0, //
            1.0, 0.0, //
            1.0, 1.0, //
            1.0, 1.0,
        ],
        Some(vec!["(Intercept)".to_string(), "conditionB".to_string()]),
    )
    .unwrap()
}

pub fn reference_reduced_design() -> DesignMatrix {
    DesignMatrix::from_row_major(
        4,
        1,
        vec![1.0, 1.0, 1.0, 1.0],
        Some(vec!["(Intercept)".to_string()]),
    )
    .unwrap()
}

pub fn reference_gene_names() -> Vec<String> {
    (1..=4).map(|idx| format!("gene{idx}")).collect()
}

pub fn reference_sample_names() -> Vec<String> {
    (1..=4).map(|idx| format!("sample{idx}")).collect()
}

pub fn read_size_factors(name: &str) -> Option<Vec<f64>> {
    let rows = read_optional_tsv(name)?;
    Some(
        rows.iter()
            .map(|row| parse_required_f64(row, "size_factor"))
            .collect(),
    )
}

pub fn read_fixed_dispersions() -> Option<Vec<f64>> {
    let rows = read_optional_tsv("fixed_dispersions.tsv")?;
    Some(
        rows.iter()
            .map(|row| parse_required_f64(row, "dispersion"))
            .collect(),
    )
}

pub fn read_reference_matrix(name: &str) -> Option<RowMajorMatrix<f64>> {
    let rows = read_optional_tsv(name)?;
    let samples = reference_sample_names();
    let values = rows
        .iter()
        .flat_map(|row| samples.iter().map(|sample| parse_required_f64(row, sample)))
        .collect::<Vec<_>>();
    Some(RowMajorMatrix::from_row_major(rows.len(), samples.len(), values).unwrap())
}

pub fn parse_required_f64(row: &TsvRow, column: &str) -> f64 {
    let value = row
        .get(column)
        .unwrap_or_else(|| panic!("missing column {column} in reference row"));
    value
        .parse::<f64>()
        .unwrap_or_else(|error| panic!("invalid float in column {column}: {value}: {error}"))
}

pub fn parse_optional_f64(row: &TsvRow, column: &str) -> Option<f64> {
    let value = row
        .get(column)
        .unwrap_or_else(|| panic!("missing column {column} in reference row"));
    match value.as_str() {
        "" | "NA" | "NaN" => None,
        _ => {
            Some(value.parse::<f64>().unwrap_or_else(|error| {
                panic!("invalid float in column {column}: {value}: {error}")
            }))
        }
    }
}

pub fn parse_required_usize(row: &TsvRow, column: &str) -> usize {
    let value = row
        .get(column)
        .unwrap_or_else(|| panic!("missing column {column} in reference row"));
    value
        .parse::<usize>()
        .unwrap_or_else(|error| panic!("invalid integer in column {column}: {value}: {error}"))
}

pub fn parse_required_bool(row: &TsvRow, column: &str) -> bool {
    match row
        .get(column)
        .unwrap_or_else(|| panic!("missing column {column} in reference row"))
        .as_str()
    {
        "TRUE" | "true" => true,
        "FALSE" | "false" => false,
        value => panic!("invalid bool in column {column}: {value}"),
    }
}

pub fn assert_float_close(actual: f64, expected: f64, atol: f64, rtol: f64, label: &str) {
    if expected.is_nan() {
        assert!(actual.is_nan(), "{label}: expected NaN, got {actual}");
        return;
    }
    let diff = (actual - expected).abs();
    let allowed = atol + rtol * expected.abs().max(1.0);
    assert!(
        diff <= allowed,
        "{label}: actual={actual}, expected={expected}, diff={diff}, allowed={allowed}"
    );
}

pub fn assert_option_close(
    actual: Option<f64>,
    expected: Option<f64>,
    atol: f64,
    rtol: f64,
    label: &str,
) {
    match (actual, expected) {
        (Some(actual), Some(expected)) => assert_float_close(actual, expected, atol, rtol, label),
        (None, None) => {}
        _ => panic!("{label}: actual={actual:?}, expected={expected:?}"),
    }
}
