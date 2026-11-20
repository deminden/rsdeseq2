//! Shared test inputs, fixture readers and numerical assertions.
pub(crate) mod builtin_checks;
mod formula_reference;
pub(crate) mod posterior;
pub(crate) mod regression_checks;

pub(crate) use formula_reference::formula_reference;

use csv::ReaderBuilder;
use rsdeseq2::errors::NameScope;
use rsdeseq2::prelude::*;
use std::collections::HashMap;
use std::path::{Path, PathBuf};

pub type TsvRow = HashMap<String, String>;

#[track_caller]
pub(crate) fn assert_ambiguous_name(
    error: DeseqError,
    expected_domain: ErrorKind,
    expected_scope: NameScope,
    expected_name: &str,
    expected_matches: &[&str],
) {
    assert_eq!(error.kind(), expected_domain);
    std::assert_matches!(
        error,
        DeseqError::AmbiguousName { scope, name, matches, .. }
            if scope == expected_scope && name == expected_name && matches == expected_matches
    );
}

#[track_caller]
pub(crate) fn assert_name_not_found(
    error: DeseqError,
    expected_domain: ErrorKind,
    expected_scope: NameScope,
    expected_name: &str,
) {
    assert_eq!(error.kind(), expected_domain);
    std::assert_matches!(
        error,
        DeseqError::NameNotFound { scope, name, .. }
            if scope == expected_scope && name == expected_name
    );
}

#[track_caller]
pub(crate) fn assert_ordered_factor_contrast_error(error: DeseqError, factor: &str) {
    let expected_operation = format!("factor-level contrast for ordered factor '{factor}'");
    std::assert_matches!(
        error,
        DeseqError::MissingPrerequisite { domain: ErrorKind::Options, operation, prerequisite }
            if operation == expected_operation
                && prerequisite == "a named or numeric contrast on its polynomial columns"
    );
}

macro_rules! json_fixture {
    ($file:literal) => {
        json_fixture!($file, serde_json::Value)
    };
    ($file:literal, $value:ty) => {
        serde_json::from_str::<$value>(include_str!(concat!("data/", $file)))
            .unwrap_or_else(|error| panic!("invalid fixture {}: {error}", $file))
    };
}
pub(crate) use json_fixture;

pub struct TsvTable {
    pub headers: Vec<String>,
    pub rows: Vec<Vec<String>>,
}

impl TsvTable {
    pub fn read(path: &Path) -> Self {
        let (headers, records) = tsv_records(path);
        let headers = headers.iter().map(str::to_owned).collect();
        let rows = records
            .map(|record| record.iter().map(str::to_owned).collect())
            .collect();
        Self { headers, rows }
    }

    pub fn column(&self, name: &str) -> usize {
        self.headers
            .iter()
            .position(|header| header == name)
            .unwrap_or_else(|| panic!("missing {name:?} in {:?}", self.headers))
    }
}

fn tsv_records(
    path: &Path,
) -> (
    csv::StringRecord,
    impl Iterator<Item = csv::StringRecord> + '_,
) {
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(path)
        .unwrap_or_else(|error| panic!("failed to open {}: {error}", path.display()));
    let headers = reader
        .headers()
        .unwrap_or_else(|error| panic!("failed to read headers from {}: {error}", path.display()))
        .clone();
    assert!(!headers.is_empty(), "empty header: {}", path.display());
    let unique = headers.iter().collect::<std::collections::HashSet<_>>();
    assert_eq!(
        unique.len(),
        headers.len(),
        "duplicate reference columns: {}",
        path.display()
    );
    let records = reader.into_records().map(move |record| {
        record.unwrap_or_else(|error| panic!("failed to read row from {}: {error}", path.display()))
    });
    (headers, records)
}

pub fn reference_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/deseq2_reference")
}

pub fn reference_file(name: &str) -> PathBuf {
    let path = reference_dir().join(name);
    assert!(
        path.is_file(),
        "required DESeq2 reference missing: {}; regenerate with scripts/references/generate_deseq2_references.R",
        path.display()
    );
    path
}

pub fn read_reference_tsv(name: &str) -> Vec<TsvRow> {
    let rows = read_tsv(&reference_file(name));
    assert!(!rows.is_empty(), "empty DESeq2 reference: {name}");
    rows
}

pub fn read_tsv(path: &Path) -> Vec<TsvRow> {
    let (headers, records) = tsv_records(path);
    records
        .map(|record| {
            headers
                .iter()
                .zip(record.iter())
                .map(|(key, value)| (key.to_string(), value.to_string()))
                .collect()
        })
        .collect()
}

pub fn two_group_design(samples: usize, coefficient: Option<&str>) -> DesignMatrix {
    let values = (0..samples)
        .flat_map(|sample| [1.0, if sample < samples / 2 { 0.0 } else { 1.0 }])
        .collect();
    let names = coefficient.map(|name| vec!["Intercept".to_string(), name.to_string()]);
    DesignMatrix::from_row_major(samples, 2, values, names).unwrap()
}

pub fn decimal_numbers(values: &[String]) -> Vec<f64> {
    values.iter().map(|value| value.parse().unwrap()).collect()
}

pub fn json_prefixed_bits(values: &serde_json::Value) -> Vec<f64> {
    values
        .as_array()
        .unwrap()
        .iter()
        .map(|value| {
            let bits = value.as_str().unwrap().strip_prefix("bits:").unwrap();
            f64::from_bits(u64::from_str_radix(bits, 16).unwrap())
        })
        .collect()
}

pub fn optional_bits(values: &[Option<f64>]) -> Vec<Option<u64>> {
    values.iter().map(|value| value.map(f64::to_bits)).collect()
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

pub fn read_size_factors(name: &str) -> Vec<f64> {
    let rows = read_reference_tsv(name);

    rows.iter()
        .map(|row| parse_required_f64(row, "size_factor"))
        .collect()
}

pub fn read_fixed_dispersions() -> Vec<f64> {
    let rows = read_reference_tsv("fixed_dispersions.tsv");

    rows.iter()
        .map(|row| parse_required_f64(row, "dispersion"))
        .collect()
}

pub fn read_reference_matrix(name: &str) -> RowMajorMatrix<f64> {
    let rows = read_reference_tsv(name);
    let samples = reference_sample_names();
    let values = rows
        .iter()
        .flat_map(|row| samples.iter().map(|sample| parse_required_f64(row, sample)))
        .collect::<Vec<_>>();
    RowMajorMatrix::from_row_major(rows.len(), samples.len(), values).unwrap()
}

pub fn parse_required_f64(row: &TsvRow, column: &str) -> f64 {
    let value = row
        .get(column)
        .unwrap_or_else(|| panic!("missing column {column} in reference row"));
    let number = value
        .parse::<f64>()
        .unwrap_or_else(|error| panic!("invalid float in column {column}: {value}: {error}"));
    assert!(
        number.is_finite(),
        "non-finite reference value in {column}: {value}"
    );
    number
}

pub fn parse_optional_f64(row: &TsvRow, column: &str) -> Option<f64> {
    let value = row
        .get(column)
        .unwrap_or_else(|| panic!("missing column {column} in reference row"));
    match value.as_str() {
        "" | "NA" | "NaN" => None,
        _ => Some(parse_required_f64(row, column)),
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
    assert!(
        actual.is_finite() && expected.is_finite(),
        "{label}: non-finite comparison {actual}, {expected}"
    );
    let diff = (actual - expected).abs();
    let allowed = atol + rtol * expected.abs().max(1.0);
    assert!(
        diff <= allowed,
        "{label}: actual={actual}, expected={expected}, diff={diff}, allowed={allowed}"
    );
}

/// Compare a scalar to a required reference value or its explicit missing row.
pub fn assert_reference_value(
    actual: f64,
    expected: Option<f64>,
    atol: f64,
    rtol: f64,
    label: &str,
) {
    match expected {
        Some(expected) => assert_float_close(actual, expected, atol, rtol, label),
        None => assert!(actual.is_nan(), "{label}: expected missing, got {actual}"),
    }
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

/// Log10 gap between valid probabilities, with zero requiring an exact match.
pub fn probability_log10_error(actual: f64, expected: f64) -> f64 {
    match (actual, expected) {
        (0.0, 0.0) => 0.0,
        (a, e) if a > 0.0 && a <= 1.0 && e > 0.0 && e <= 1.0 => (a.log10() - e.log10()).abs(),
        _ => f64::INFINITY,
    }
}

/// Genes whose significance at `alpha` differs between `actual` and
/// `expected` adjusted p-values; a missing value counts as not significant.
/// Parity tests treat any flip as a failure.
pub fn significance_flips(actual: &[Option<f64>], expected: &[Option<f64>], alpha: f64) -> usize {
    assert_eq!(
        actual.len(),
        expected.len(),
        "padj vectors differ in length"
    );
    actual
        .iter()
        .zip(expected)
        .filter(|(a, e)| a.is_some_and(|v| v < alpha) != e.is_some_and(|v| v < alpha))
        .count()
}

/// Assert `actual` equals `expected` to an absolute 1e-12 (or relative
/// machine epsilon), treating an expected NaN as "must be NaN".
pub fn assert_float_close_or_nan(actual: f64, expected: f64, label: &str) {
    if expected.is_nan() {
        assert!(actual.is_nan(), "{label}: expected NaN, got {actual}");
    } else {
        approx::assert_relative_eq!(actual, expected, epsilon = 1e-12);
    }
}

/// [`assert_float_close_or_nan`] element by element, after a length check.
pub fn assert_slice_close_or_nan(actual: &[f64], expected: &[f64], label: &str) {
    assert_eq!(actual.len(), expected.len(), "{label}: length mismatch");
    for (index, (actual, expected)) in actual.iter().zip(expected).enumerate() {
        assert_float_close_or_nan(*actual, *expected, &format!("{label}[{index}]"));
    }
}

/// [`assert_slice_close_or_nan`] over two matrices of the same shape.
pub fn assert_matrix_close_or_nan(
    actual: &RowMajorMatrix<f64>,
    expected: &RowMajorMatrix<f64>,
    label: &str,
) {
    assert_eq!(actual.n_rows(), expected.n_rows(), "{label}: row mismatch");
    assert_eq!(
        actual.n_cols(),
        expected.n_cols(),
        "{label}: column mismatch"
    );
    assert_slice_close_or_nan(actual.as_slice(), expected.as_slice(), label);
}

/// Parse a number as R writes it: `bits:<hex>` for exact bits, `NA` and
/// `NaN` as NaN, `Inf` and `-Inf`, otherwise decimal.
pub fn r_number(value: &str) -> f64 {
    if let Some(bits) = value.strip_prefix("bits:") {
        return f64::from_bits(u64::from_str_radix(bits, 16).unwrap());
    }
    match value {
        "NA" | "NaN" => f64::NAN,
        "Inf" => f64::INFINITY,
        "-Inf" => f64::NEG_INFINITY,
        _ => value.parse().unwrap(),
    }
}

/// [`r_number`] over a list.
pub fn r_numbers(values: &[String]) -> Vec<f64> {
    values.iter().map(|value| r_number(value)).collect()
}

/// The 16-digit lower-case hex of a float's bits, as the JSON fixtures store it.
pub fn hex_bits(value: f64) -> String {
    format!("{:016x}", value.to_bits())
}
