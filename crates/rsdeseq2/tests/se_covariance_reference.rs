mod common;

use std::collections::HashMap;
use std::path::Path;

use common::*;
use csv::ReaderBuilder;
use rsdeseq2::prelude::*;

#[test]
fn hard_high_dispersion_fixed_dispersions_match_deseq2_se_covariance() {
    let fixture = reference_dir();
    let (genes, samples, counts_values) =
        read_count_fixture(&fixture.join("hard_se_covariance_counts.tsv"));
    let (coefficient_names, design_values) =
        read_design_fixture(&fixture.join("hard_se_covariance_design.tsv"));
    let size_factors =
        read_size_factor_fixture(&fixture.join("hard_se_covariance_size_factors.tsv"));
    let dispersions = read_gene_value_fixture(
        &fixture.join("hard_se_covariance_dispersions.tsv"),
        "dispersion",
    );
    let expected_beta = read_matrix_fixture(&fixture.join("hard_se_covariance_beta.tsv"));
    let expected_beta_se = read_matrix_fixture(&fixture.join("hard_se_covariance_beta_se.tsv"));
    let expected_mu = read_matrix_fixture(&fixture.join("hard_se_covariance_mu.tsv"));
    let expected_hat = read_matrix_fixture(&fixture.join("hard_se_covariance_hat.tsv"));

    let counts = CountMatrix::from_row_major_u32_with_names(
        genes.len(),
        samples.len(),
        counts_values,
        Some(genes.clone()),
        Some(samples.clone()),
    )
    .unwrap();
    let design = DesignMatrix::from_row_major(
        samples.len(),
        coefficient_names.len(),
        design_values,
        Some(coefficient_names.clone()),
    )
    .unwrap();
    let normalization = RowMajorMatrix::from_row_major(
        genes.len(),
        samples.len(),
        genes
            .iter()
            .flat_map(|_| samples.iter().map(|sample| size_factors[sample]))
            .collect(),
    )
    .unwrap();

    let fit = fit_fixed_dispersion_irls_with_normalization_factors(
        &counts,
        &design,
        &normalization,
        &genes
            .iter()
            .map(|gene| dispersions[gene])
            .collect::<Vec<_>>(),
        IrlsOptions {
            use_optim: false,
            solver: IrlsSolver::NormalEquations,
            ..IrlsOptions::default()
        },
    )
    .unwrap();

    assert_matrix_close(
        "hard SE covariance beta",
        &fit.beta,
        &expected_beta,
        1.0e-12,
        1.0e-12,
    );
    assert_matrix_close(
        "hard SE covariance betaSE",
        &fit.beta_se,
        &expected_beta_se,
        1.0e-13,
        1.0e-13,
    );
    assert_matrix_close(
        "hard SE covariance mu",
        &fit.mu,
        &expected_mu,
        1.0e-6,
        1.0e-13,
    );
    assert_matrix_close(
        "hard SE covariance hat",
        &fit.hat_diagonal,
        &expected_hat,
        1.0e-14,
        1.0e-13,
    );
}

struct MatrixFixture {
    values: RowMajorMatrix<f64>,
}

fn read_tsv_records(path: &Path) -> (Vec<String>, Vec<Vec<String>>) {
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(path)
        .unwrap_or_else(|error| panic!("failed to open {}: {error}", path.display()));
    let header = reader
        .headers()
        .unwrap_or_else(|error| panic!("failed to read headers from {}: {error}", path.display()))
        .iter()
        .map(str::to_string)
        .collect::<Vec<_>>();
    let rows = reader
        .records()
        .map(|record| {
            record
                .unwrap_or_else(|error| {
                    panic!("failed to read row from {}: {error}", path.display())
                })
                .iter()
                .map(str::to_string)
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    (header, rows)
}

fn read_count_fixture(path: &Path) -> (Vec<String>, Vec<String>, Vec<u32>) {
    let (header, rows) = read_tsv_records(path);
    let samples = header[1..].to_vec();
    let mut genes = Vec::with_capacity(rows.len());
    let mut values = Vec::with_capacity(rows.len() * samples.len());
    for row in rows {
        genes.push(row[0].clone());
        values.extend(row[1..].iter().map(|value| {
            value
                .parse::<u32>()
                .unwrap_or_else(|error| panic!("invalid count {value}: {error}"))
        }));
    }
    (genes, samples, values)
}

fn read_design_fixture(path: &Path) -> (Vec<String>, Vec<f64>) {
    let (header, rows) = read_tsv_records(path);
    let coefficient_names = header[1..].to_vec();
    let values = rows
        .iter()
        .flat_map(|row| row[1..].iter().map(|value| parse_fixture_f64(value)))
        .collect::<Vec<_>>();
    (coefficient_names, values)
}

fn read_size_factor_fixture(path: &Path) -> HashMap<String, f64> {
    let (_, rows) = read_tsv_records(path);
    rows.into_iter()
        .map(|row| (row[0].clone(), parse_fixture_f64(&row[1])))
        .collect()
}

fn read_gene_value_fixture(path: &Path, column: &str) -> HashMap<String, f64> {
    let (header, rows) = read_tsv_records(path);
    let column_idx = header
        .iter()
        .position(|name| name == column)
        .unwrap_or_else(|| panic!("missing column {column} in {}", path.display()));
    rows.into_iter()
        .map(|row| (row[0].clone(), parse_fixture_f64(&row[column_idx])))
        .collect()
}

fn read_matrix_fixture(path: &Path) -> MatrixFixture {
    let (header, rows) = read_tsv_records(path);
    let n_cols = header.len() - 1;
    let values = rows
        .iter()
        .flat_map(|row| row[1..].iter().map(|value| parse_fixture_f64(value)))
        .collect::<Vec<_>>();
    MatrixFixture {
        values: RowMajorMatrix::from_row_major(rows.len(), n_cols, values).unwrap(),
    }
}

fn assert_matrix_close(
    label: &str,
    actual: &RowMajorMatrix<f64>,
    expected: &MatrixFixture,
    atol: f64,
    rtol: f64,
) {
    assert_eq!(actual.n_rows(), expected.values.n_rows(), "{label} rows");
    assert_eq!(actual.n_cols(), expected.values.n_cols(), "{label} columns");
    for row in 0..actual.n_rows() {
        let actual_row = actual.row(row).unwrap();
        let expected_row = expected.values.row(row).unwrap();
        for col in 0..actual.n_cols() {
            assert_float_close(
                actual_row[col],
                expected_row[col],
                atol,
                rtol,
                &format!("{label} row {row} col {col}"),
            );
        }
    }
}

fn parse_fixture_f64(value: &str) -> f64 {
    match value {
        "NA" | "" => f64::NAN,
        _ => value
            .parse()
            .unwrap_or_else(|error| panic!("invalid float {value}: {error}")),
    }
}
