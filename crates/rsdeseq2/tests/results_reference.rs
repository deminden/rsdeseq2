mod common;

use common::*;
use rsdeseq2::prelude::*;

#[test]
fn normalization_stage_matches_optional_deseq2_reference() {
    let Some(expected_size_factors) = read_size_factors("size_factors_ratio.tsv") else {
        return;
    };
    let Some(expected_base) = read_optional_tsv("base_metadata_ratio.tsv") else {
        return;
    };

    let counts = reference_counts();
    let fit = DeseqBuilder::new()
        .fit_size_factors_and_base_means(&counts)
        .unwrap();

    assert_eq!(fit.size_factors.len(), expected_size_factors.len());
    for (sample, (actual, expected)) in fit
        .size_factors
        .iter()
        .copied()
        .zip(expected_size_factors.iter().copied())
        .enumerate()
    {
        assert_float_close(
            actual,
            expected,
            1e-12,
            1e-12,
            &format!("size factor sample {sample}"),
        );
    }

    assert_eq!(fit.base_mean.len(), expected_base.len());
    for (gene, row) in expected_base.iter().enumerate() {
        assert_eq!(
            row.get("gene").map(String::as_str),
            Some(reference_gene_names()[gene].as_str())
        );
        assert_float_close(
            fit.base_mean[gene],
            parse_required_f64(row, "baseMean"),
            1e-10,
            1e-10,
            &format!("baseMean gene {gene}"),
        );
        assert_float_close(
            fit.base_var[gene],
            parse_required_f64(row, "baseVar"),
            1e-10,
            1e-10,
            &format!("baseVar gene {gene}"),
        );
        assert_eq!(fit.all_zero[gene], parse_required_bool(row, "allZero"));
    }
}

#[test]
fn normalized_counts_match_optional_deseq2_reference() {
    let Some(expected_size_factors) = read_size_factors("size_factors_ratio.tsv") else {
        return;
    };
    let Some(expected_counts) = read_optional_tsv("normalized_counts_ratio.tsv") else {
        return;
    };

    let counts = reference_counts();
    let normalized = normalized_counts(&counts, &expected_size_factors).unwrap();
    let samples = reference_sample_names();

    assert_eq!(normalized.n_rows(), expected_counts.len());
    for (gene, row) in expected_counts.iter().enumerate() {
        assert_eq!(
            row.get("gene").map(String::as_str),
            Some(reference_gene_names()[gene].as_str())
        );
        for (sample, sample_name) in samples.iter().enumerate() {
            assert_float_close(
                normalized.row(gene).unwrap()[sample],
                parse_required_f64(row, sample_name),
                1e-10,
                1e-10,
                &format!("normalized count gene {gene} sample {sample}"),
            );
        }
    }
}

#[test]
fn normalization_factors_match_optional_deseq2_reference() {
    let Some(factor_rows) = read_optional_tsv("normalization_factors.tsv") else {
        return;
    };
    let Some(expected_counts) = read_optional_tsv("normalized_counts_nf.tsv") else {
        return;
    };
    let Some(expected_base) = read_optional_tsv("base_metadata_nf.tsv") else {
        return;
    };

    let counts = reference_counts();
    let samples = reference_sample_names();
    let factor_values = factor_rows
        .iter()
        .flat_map(|row| samples.iter().map(|sample| parse_required_f64(row, sample)))
        .collect::<Vec<_>>();
    let normalization_factors =
        RowMajorMatrix::from_row_major(factor_rows.len(), samples.len(), factor_values).unwrap();

    let normalized = normalized_counts_with_factors(&counts, &normalization_factors).unwrap();
    for (gene, row) in expected_counts.iter().enumerate() {
        assert_eq!(
            row.get("gene").map(String::as_str),
            Some(reference_gene_names()[gene].as_str())
        );
        for (sample, sample_name) in samples.iter().enumerate() {
            assert_float_close(
                normalized.row(gene).unwrap()[sample],
                parse_required_f64(row, sample_name),
                1e-10,
                1e-10,
                &format!("normalization-factor count gene {gene} sample {sample}"),
            );
        }
    }

    let fit = DeseqBuilder::new()
        .normalization_factors(normalization_factors)
        .fit_size_factors_and_base_means(&counts)
        .unwrap();
    for (gene, row) in expected_base.iter().enumerate() {
        assert_float_close(
            fit.base_mean[gene],
            parse_required_f64(row, "baseMean"),
            1e-10,
            1e-10,
            &format!("normalization-factor baseMean gene {gene}"),
        );
        assert_float_close(
            fit.base_var[gene],
            parse_required_f64(row, "baseVar"),
            1e-10,
            1e-10,
            &format!("normalization-factor baseVar gene {gene}"),
        );
        assert_eq!(fit.all_zero[gene], parse_required_bool(row, "allZero"));
    }
}

#[test]
fn weighted_base_metadata_matches_optional_deseq2_reference() {
    let Some(expected_base) = read_optional_tsv("base_metadata_weighted.tsv") else {
        return;
    };
    let Some(expected_weights) = read_reference_matrix("observation_weights_normalized.tsv") else {
        return;
    };
    let Some(observation_weights) = read_reference_matrix("observation_weights.tsv") else {
        return;
    };
    let Some(size_factors) = read_size_factors("size_factors_ratio.tsv") else {
        return;
    };

    let fit = DeseqBuilder::new()
        .size_factors(size_factors)
        .observation_weights(observation_weights)
        .fit_size_factors_and_base_means_with_design(&reference_counts(), &reference_full_design())
        .unwrap();

    assert_eq!(
        fit.weights_design_rank,
        Some(reference_full_design().n_coefficients())
    );
    assert_eq!(fit.base_mean.len(), expected_base.len());
    for (gene, row) in expected_base.iter().enumerate() {
        assert_float_close(
            fit.base_mean[gene],
            parse_required_f64(row, "baseMean"),
            1e-10,
            1e-10,
            &format!("weighted baseMean gene {gene}"),
        );
        assert_float_close(
            fit.base_var[gene],
            parse_required_f64(row, "baseVar"),
            1e-10,
            1e-10,
            &format!("weighted baseVar gene {gene}"),
        );
        assert_eq!(fit.all_zero[gene], parse_required_bool(row, "allZero"));
        assert_eq!(
            fit.weights_fail.as_ref().unwrap()[gene],
            parse_required_bool(row, "weightsFail")
        );
        for sample in 0..expected_weights.n_cols() {
            assert_float_close(
                fit.observation_weights.as_ref().unwrap().row(gene).unwrap()[sample],
                expected_weights.row(gene).unwrap()[sample],
                1e-12,
                1e-12,
                &format!("normalized observation weight gene {gene} sample {sample}"),
            );
        }
    }
}

#[test]
fn full_deseq2_results_reference_shape_is_documented() {
    let Some(rows) = read_optional_tsv("results_wald_ratio.tsv") else {
        return;
    };
    let genes = reference_gene_names();

    assert_eq!(rows.len(), genes.len());
    for (idx, row) in rows.iter().enumerate() {
        assert_eq!(
            row.get("gene").map(String::as_str),
            Some(genes[idx].as_str())
        );
        for column in [
            "baseMean",
            "log2FoldChange",
            "lfcSE",
            "stat",
            "pvalue",
            "padj",
        ] {
            assert!(
                row.contains_key(column),
                "missing DESeq2 results column {column}"
            );
        }
    }
}
