use std::collections::HashMap;
use std::env;
use std::fs;
use std::path::{Path, PathBuf};

use rsdeseq2::prelude::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let fixture = env::args().nth(1).map(PathBuf::from).unwrap_or_else(|| {
        PathBuf::from("../../results/fixtures/se_covariance_hard_real_2026-06-05")
    });

    let (genes, samples, counts_values) = read_count_matrix(&fixture.join("selected_counts.tsv"))?;
    let counts = CountMatrix::from_row_major_u32_with_names(
        genes.len(),
        samples.len(),
        counts_values,
        Some(genes.clone()),
        Some(samples.clone()),
    )?;
    let (coefficient_names, design_values) =
        read_design_matrix(&fixture.join("design_matrix.tsv"))?;
    let design = DesignMatrix::from_row_major(
        samples.len(),
        coefficient_names.len(),
        design_values,
        Some(coefficient_names.clone()),
    )?;
    let size_factors = read_size_factors(&fixture.join("size_factors.tsv"))?;
    let normalization = RowMajorMatrix::from_row_major(
        genes.len(),
        samples.len(),
        (0..genes.len())
            .flat_map(|_| samples.iter().map(|sample| size_factors[sample]))
            .collect(),
    )?;
    let dispersions = read_gene_values(&fixture.join("selected_dispersions.tsv"), "dispersion")?;
    let expected_beta = read_named_matrix(&fixture.join("beta_no_optim.tsv"))?;
    let expected_beta_se = read_named_matrix(&fixture.join("beta_se_no_optim.tsv"))?;
    let expected_mu = read_named_matrix(&fixture.join("mu_no_optim.tsv"))?;
    let expected_hat = read_named_matrix(&fixture.join("hat_no_optim.tsv"))?;

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
    )?;

    summarize_matrix(
        "beta",
        &genes,
        &coefficient_names,
        fit.beta.as_slice(),
        &expected_beta,
    );
    summarize_matrix(
        "betaSE",
        &genes,
        &expected_beta_se.column_names,
        fit.beta_se.as_slice(),
        &expected_beta_se,
    );
    summarize_matrix("mu", &genes, &samples, fit.mu.as_slice(), &expected_mu);
    summarize_matrix(
        "hat",
        &genes,
        &samples,
        fit.hat_diagonal.as_slice(),
        &expected_hat,
    );
    Ok(())
}

struct NamedMatrix {
    column_names: Vec<String>,
    rows: HashMap<String, Vec<f64>>,
}

fn summarize_matrix(
    label: &str,
    genes: &[String],
    columns: &[String],
    actual: &[f64],
    expected: &NamedMatrix,
) {
    let mut diffs = Vec::new();
    let mut worst = ("", "", 0.0, 0.0, 0.0);
    for (gene_idx, gene) in genes.iter().enumerate() {
        let Some(expected_row) = expected.rows.get(gene) else {
            continue;
        };
        for (col_idx, column) in columns.iter().enumerate() {
            let actual_value = actual[gene_idx * columns.len() + col_idx];
            let expected_value = expected_row[col_idx];
            if !actual_value.is_finite() || !expected_value.is_finite() {
                continue;
            }
            let diff = (actual_value - expected_value).abs();
            diffs.push(diff);
            if diff > worst.2 {
                worst = (
                    gene.as_str(),
                    column.as_str(),
                    diff,
                    actual_value,
                    expected_value,
                );
            }
        }
    }
    diffs.sort_by(|left, right| left.total_cmp(right));
    let mean = diffs.iter().sum::<f64>() / diffs.len() as f64;
    let median = diffs[diffs.len() / 2];
    let p99 = diffs[((diffs.len() as f64 * 0.99).ceil() as usize).saturating_sub(1)];
    println!(
        "{label}\tmean={mean:.6e}\tmedian={median:.6e}\tp99={p99:.6e}\tmax={:.6e}\tworst_gene={}\tworst_col={}\tactual={:.16e}\texpected={:.16e}",
        worst.2, worst.0, worst.1, worst.3, worst.4
    );
}

type TsvTable = (Vec<String>, Vec<Vec<String>>);

fn read_tsv(path: &Path) -> Result<TsvTable, Box<dyn std::error::Error>> {
    let text = fs::read_to_string(path)?;
    let mut lines = text.lines();
    let header = lines
        .next()
        .ok_or_else(|| format!("empty TSV: {}", path.display()))?
        .split('\t')
        .map(str::to_string)
        .collect::<Vec<_>>();
    let rows = lines
        .filter(|line| !line.is_empty())
        .map(|line| line.split('\t').map(str::to_string).collect::<Vec<_>>())
        .collect::<Vec<_>>();
    Ok((header, rows))
}

type CountMatrixParts = (Vec<String>, Vec<String>, Vec<u32>);

fn read_count_matrix(path: &Path) -> Result<CountMatrixParts, Box<dyn std::error::Error>> {
    let (header, rows) = read_tsv(path)?;
    let samples = header[1..].to_vec();
    let mut genes = Vec::with_capacity(rows.len());
    let mut values = Vec::with_capacity(rows.len() * samples.len());
    for row in rows {
        genes.push(row[0].clone());
        for value in &row[1..] {
            values.push(value.parse()?);
        }
    }
    Ok((genes, samples, values))
}

fn read_design_matrix(path: &Path) -> Result<(Vec<String>, Vec<f64>), Box<dyn std::error::Error>> {
    let (header, rows) = read_tsv(path)?;
    let coefficient_names = header[1..].to_vec();
    let mut values = Vec::with_capacity(rows.len() * coefficient_names.len());
    for row in rows {
        for value in &row[1..] {
            values.push(value.parse()?);
        }
    }
    Ok((coefficient_names, values))
}

fn read_size_factors(path: &Path) -> Result<HashMap<String, f64>, Box<dyn std::error::Error>> {
    let (_, rows) = read_tsv(path)?;
    let mut values = HashMap::with_capacity(rows.len());
    for row in rows {
        values.insert(row[0].clone(), parse_f64(&row[1])?);
    }
    Ok(values)
}

fn read_gene_values(
    path: &Path,
    column: &str,
) -> Result<HashMap<String, f64>, Box<dyn std::error::Error>> {
    let (header, rows) = read_tsv(path)?;
    let col_idx = header
        .iter()
        .position(|name| name == column)
        .ok_or_else(|| format!("missing column {column} in {}", path.display()))?;
    let mut values = HashMap::with_capacity(rows.len());
    for row in rows {
        values.insert(row[0].clone(), parse_f64(&row[col_idx])?);
    }
    Ok(values)
}

fn read_named_matrix(path: &Path) -> Result<NamedMatrix, Box<dyn std::error::Error>> {
    let (header, rows) = read_tsv(path)?;
    let mut values = HashMap::with_capacity(rows.len());
    for row in rows {
        values.insert(
            row[0].clone(),
            row[1..]
                .iter()
                .map(|value| parse_f64(value))
                .collect::<Result<Vec<f64>, _>>()?,
        );
    }
    Ok(NamedMatrix {
        column_names: header[1..].to_vec(),
        rows: values,
    })
}

fn parse_f64(value: &str) -> Result<f64, Box<dyn std::error::Error>> {
    if value == "NA" || value.is_empty() {
        Ok(f64::NAN)
    } else {
        Ok(value.parse()?)
    }
}
