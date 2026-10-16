use rsdeseq2::prelude::*;
use std::path::PathBuf;

#[derive(Debug)]
struct Tsv {
    headers: Vec<String>,
    rows: Vec<Vec<String>>,
}

impl Tsv {
    fn read(dir: &std::path::Path, name: &str) -> Self {
        let path = dir.join(name);
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
        value => {
            let number: f64 = value.parse().expect("numeric fixture value");
            assert!(number.is_finite(), "non-finite numeric fixture: {value}");
            Some(number)
        }
    }
}

fn assert_value(actual: Option<f64>, expected: Option<f64>, tolerance: f64, label: &str) -> f64 {
    match (actual, expected) {
        (None, None) => 0.0,
        (Some(a), Some(b)) => {
            assert!(
                a.is_finite() && b.is_finite(),
                "{label}: non-finite result {a}, {b}"
            );
            let error = (a - b).abs();
            assert!(
                error < tolerance,
                "{label}: {a} vs {b}, error {error} >= {tolerance}"
            );
            error
        }
        _ => panic!("{label}: missingness mismatch {actual:?} vs {expected:?}"),
    }
}

#[test]
fn compact_gtex_wald_matches_all_six_deseq2_fields() {
    let manifest: serde_json::Value = serde_json::from_str(
        &std::fs::read_to_string(fixture_dir().join("compact/selection.json")).unwrap(),
    )
    .unwrap();
    assert_eq!(manifest["source_genes"], 56200);
    assert_eq!(manifest["selected_genes"], 512);
    assert_eq!(manifest["samples"], 61);
    check_wald(true);
}

#[test]
#[ignore = "full 56,200-gene release replay; compact panel runs in normal CI"]
fn full_gtex_wald_matches_deseq2() {
    check_wald(false);
}

#[test]
#[should_panic(expected = "non-finite result")]
fn parity_rejects_nan() {
    assert_value(Some(f64::NAN), Some(1.0), 1e-6, "regression");
}

#[test]
#[should_panic(expected = "missingness mismatch")]
fn parity_rejects_unexpected_value_for_na() {
    assert_value(Some(1.0), None, 1e-6, "regression");
}

fn check_wald(compact: bool) {
    let dir = if compact {
        fixture_dir().join("compact")
    } else {
        fixture_dir()
    };
    let counts_table = Tsv::read(&dir, "counts.tsv");
    let design_table = Tsv::read(&fixture_dir(), "design_matrix.tsv");
    let reference = Tsv::read(&dir, "reference_results.tsv");
    let reference_size_factors = Tsv::read(&fixture_dir(), "reference_size_factors.tsv");
    let gene_names = counts_table
        .rows
        .iter()
        .map(|row| row[0].clone())
        .collect::<Vec<_>>();
    let sample_names = counts_table.headers[1..].to_vec();
    let n_genes = gene_names.len();
    let n_samples = sample_names.len();
    assert_eq!(n_genes, if compact { 512 } else { 56200 });
    assert_eq!(n_samples, 61);
    assert_eq!(
        gene_names
            .iter()
            .collect::<std::collections::HashSet<_>>()
            .len(),
        n_genes
    );
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
    let sf_column = reference_size_factors.column("sizeFactor");
    let factors = reference_size_factors
        .rows
        .iter()
        .map(|r| r[sf_column].parse::<f64>().unwrap())
        .collect::<Vec<_>>();
    let builder = if compact {
        DeseqBuilder::new().size_factors(factors)
    } else {
        DeseqBuilder::new()
    };
    let (fit, results) = builder
        .disable_cooks_cutoff()
        .disable_independent_filtering()
        .fit_fixed_dispersion_wald(&counts, &design, &dispersions, result_coefficient)
        .unwrap();

    assert_eq!(results.rows.len(), n_genes, "Rust result gene count");
    let mut max_size_factor_relative_error: f64 = 0.0;
    for (index, expected_row) in reference_size_factors.rows.iter().enumerate() {
        assert_eq!(reference_size_factors.rows[index][0], sample_names[index]);
        let expected = expected_row[sf_column].parse::<f64>().unwrap();
        let actual = fit.size_factors[index];
        assert!(actual.is_finite() && expected.is_finite() && expected > 0.0);
        max_size_factor_relative_error =
            max_size_factor_relative_error.max((actual / expected - 1.0).abs());
    }

    let columns = ["baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue"];
    let indexes = columns.map(|name| reference.column(name));
    let mut max_errors = [0.0_f64; 6];
    let padj = reference_padj(&reference);
    let mut comparable_rows = 0usize;
    for (index, (actual, expected_row)) in results.rows.iter().zip(&reference.rows).enumerate() {
        assert_eq!(actual.gene.as_deref(), Some(gene_names[index].as_str()));
        let expected_values = indexes.map(|column| parse_optional(&expected_row[column]));
        assert_eq!(
            expected_row[0], gene_names[index],
            "reference gene alignment"
        );
        let actual_values = [
            Some(actual.base_mean),
            actual.log2_fold_change,
            actual.lfc_se,
            actual.stat,
            actual.pvalue,
            actual.padj,
        ];
        let expected_values_all = [
            expected_values[0],
            expected_values[1],
            expected_values[2],
            expected_values[3],
            expected_values[4],
            padj[index],
        ];
        for slot in 0..6 {
            let label = format!(
                "{} {}",
                gene_names[index],
                ["baseMean", "log2FC", "lfcSE", "stat", "pvalue", "padj"][slot]
            );
            max_errors[slot] = max_errors[slot].max(assert_value(
                actual_values[slot],
                expected_values_all[slot],
                [1e-6, 1e-5, 1e-6, 1e-5, 1e-6, 1e-6][slot],
                &label,
            ));
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
    eprintln!("padj max absolute error={:.3e}", max_errors[5]);
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

// The original full reference predates the padj column. Compute its BH oracle
// from the R p-values, independently of the candidate's p-values and BH code.
// The compact oracle is frozen directly by stats::p.adjust in R.
fn reference_padj(reference: &Tsv) -> Vec<Option<f64>> {
    if let Some(col) = reference.headers.iter().position(|s| s == "padj") {
        return reference
            .rows
            .iter()
            .map(|r| parse_optional(&r[col]))
            .collect();
    }
    let col = reference.column("pvalue");
    let mut ranked = reference
        .rows
        .iter()
        .enumerate()
        .filter_map(|(i, r)| parse_optional(&r[col]).map(|p| (i, p)))
        .collect::<Vec<_>>();
    ranked.sort_by(|a, b| a.1.total_cmp(&b.1));
    let mut out = vec![None; reference.rows.len()];
    let mut minimum = 1.0_f64;
    for rank in (0..ranked.len()).rev() {
        minimum = minimum.min(ranked[rank].1 * ranked.len() as f64 / (rank + 1) as f64);
        out[ranked[rank].0] = Some(minimum);
    }
    out
}

#[test]
fn full_gtex_normalization_matches_deseq2_without_fitting_glms() {
    let counts = rsdeseq2::io::read_count_matrix_tsv(fixture_dir().join("counts.tsv")).unwrap();
    assert_eq!((counts.n_genes(), counts.n_samples()), (56200, 61));
    let fit = DeseqBuilder::new()
        .fit_size_factors_and_base_means(&counts)
        .unwrap();
    let expected = Tsv::read(&fixture_dir(), "reference_size_factors.tsv");
    assert_eq!(fit.size_factors.len(), expected.rows.len());
    for (a, r) in fit.size_factors.iter().zip(&expected.rows) {
        let b = r[1].parse::<f64>().unwrap();
        assert!(a.is_finite() && (a / b - 1.0).abs() < 1e-10);
    }
}
