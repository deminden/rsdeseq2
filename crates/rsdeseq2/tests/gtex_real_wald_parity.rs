use rsdeseq2::prelude::*;
use std::path::PathBuf;

#[derive(Debug)]
struct Tsv {
    headers: Vec<String>,
    rows: Vec<Vec<String>>,
}

impl Tsv {
    fn read(name: &str) -> Self {
        let path = fixture_dir().join(name);
        let mut reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .from_path(&path)
            .unwrap_or_else(|error| panic!("read {}: {error}", path.display()));
        let headers = reader
            .headers()
            .unwrap_or_else(|error| panic!("read headers from {}: {error}", path.display()))
            .iter()
            .map(str::to_owned)
            .collect();
        let rows = reader
            .records()
            .map(|record| {
                record
                    .unwrap_or_else(|error| panic!("read row from {}: {error}", path.display()))
                    .iter()
                    .map(str::to_owned)
                    .collect()
            })
            .collect();
        Self { headers, rows }
    }

    fn column(&self, name: &str) -> usize {
        self.headers
            .iter()
            .position(|header| header == name)
            .unwrap_or_else(|| panic!("missing {name:?} in {:?}", self.headers))
    }
}

fn fixture_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/gtex_v8_real_wald")
}

fn parse_optional(value: &str) -> Option<f64> {
    match value {
        "NA" | "NaN" | "" => None,
        value => Some(value.parse().expect("numeric fixture value")),
    }
}

#[test]
fn gtex_v8_all_gene_counts_match_deseq2_fixed_dispersion_wald_reference() {
    let counts_table = Tsv::read("counts.tsv");
    let design_table = Tsv::read("design_matrix.tsv");
    let reference = Tsv::read("reference_results.tsv");
    let reference_size_factors = Tsv::read("reference_size_factors.tsv");

    let gene_names = counts_table
        .rows
        .iter()
        .map(|row| row[0].clone())
        .collect::<Vec<_>>();
    let sample_names = counts_table.headers[1..].to_vec();
    let n_genes = gene_names.len();
    let n_samples = sample_names.len();
    assert_eq!(n_genes, reference.rows.len(), "reference gene count");
    assert_eq!(n_samples, design_table.rows.len(), "design sample count");
    assert_eq!(
        n_samples,
        reference_size_factors.rows.len(),
        "size-factor sample count"
    );
    assert_eq!(design_table.rows[0][0], sample_names[0]);

    let mut count_values = Vec::with_capacity(n_genes * n_samples);
    for row in &counts_table.rows {
        assert_eq!(row.len(), n_samples + 1, "count row width");
        count_values.extend(row[1..].iter().map(|value| {
            value
                .parse::<u32>()
                .unwrap_or_else(|error| panic!("invalid integer count {value:?}: {error}"))
        }));
    }
    let counts = CountMatrix::from_row_major_u32_with_names(
        n_genes,
        n_samples,
        count_values,
        Some(gene_names.clone()),
        Some(sample_names.clone()),
    )
    .unwrap();

    let design_names = design_table.headers[1..].to_vec();
    let mut design_values = Vec::with_capacity(n_samples * design_names.len());
    for (index, row) in design_table.rows.iter().enumerate() {
        assert_eq!(row[0], sample_names[index], "design sample alignment");
        design_values.extend(row[1..].iter().map(|value| value.parse::<f64>().unwrap()));
    }
    let design = DesignMatrix::from_row_major(
        n_samples,
        design_names.len(),
        design_values,
        Some(design_names.clone()),
    )
    .unwrap();

    let dispersion_column = reference.column("dispersion");
    let dispersions = reference
        .rows
        .iter()
        .map(|row| parse_optional(&row[dispersion_column]).unwrap_or(1.0))
        .collect::<Vec<_>>();
    let result_coefficient = design_names.len() - 1;
    let (fit, results) = DeseqBuilder::new()
        .disable_cooks_cutoff()
        .disable_independent_filtering()
        .fit_fixed_dispersion_wald(&counts, &design, &dispersions, result_coefficient)
        .unwrap();

    assert_eq!(results.rows.len(), n_genes, "Rust result gene count");
    let sf_column = reference_size_factors.column("sizeFactor");
    let mut max_size_factor_relative_error: f64 = 0.0;
    for (index, expected_row) in reference_size_factors.rows.iter().enumerate() {
        assert_eq!(reference_size_factors.rows[index][0], sample_names[index]);
        let expected = expected_row[sf_column].parse::<f64>().unwrap();
        let actual = fit.size_factors[index];
        max_size_factor_relative_error =
            max_size_factor_relative_error.max((actual / expected - 1.0).abs());
    }

    let columns = ["baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue"];
    let indexes = columns.map(|name| reference.column(name));
    let mut max_errors = [0.0_f64; 5];
    let mut comparable_rows = 0usize;
    for (index, (actual, expected_row)) in results.rows.iter().zip(&reference.rows).enumerate() {
        assert_eq!(actual.gene.as_deref(), Some(gene_names[index].as_str()));
        let expected_values = indexes.map(|column| parse_optional(&expected_row[column]));
        for (slot, expected) in expected_values.iter().enumerate() {
            if let Some(expected) = expected {
                let actual = match slot {
                    0 => Some(actual.base_mean),
                    1 => actual.log2_fold_change,
                    2 => actual.lfc_se,
                    3 => actual.stat,
                    4 => actual.pvalue,
                    _ => unreachable!(),
                };
                let actual = actual.unwrap_or_else(|| {
                    panic!(
                        "{} unexpectedly missing {}",
                        gene_names[index], columns[slot]
                    )
                });
                max_errors[slot] = max_errors[slot].max((actual - expected).abs());
            }
        }
        if expected_values[1].is_some()
            && expected_values[2].is_some()
            && expected_values[3].is_some()
        {
            comparable_rows += 1;
        }
    }

    eprintln!(
        "GTEx v8 parity: {n_genes} genes x {n_samples} samples; {comparable_rows} finite Wald rows; max errors baseMean={:.3e}, log2FC={:.3e}, lfcSE={:.3e}, stat={:.3e}, pvalue={:.3e}; size-factor max relative error={:.3e}",
        max_errors[0],
        max_errors[1],
        max_errors[2],
        max_errors[3],
        max_errors[4],
        max_size_factor_relative_error
    );
    assert!(max_size_factor_relative_error < 1e-10);
    assert!(max_errors[0] < 1e-6);
    assert!(comparable_rows > n_genes / 2);
    assert!(
        max_errors[1] < 1e-5,
        "log2FC max abs error: {}",
        max_errors[1]
    );
    assert!(
        max_errors[2] < 1e-6,
        "lfcSE max abs error: {}",
        max_errors[2]
    );
    assert!(
        max_errors[3] < 1e-5,
        "stat max abs error: {}",
        max_errors[3]
    );
    assert!(
        max_errors[4] < 1e-6,
        "pvalue max abs error: {}",
        max_errors[4]
    );
}
