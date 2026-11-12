use super::*;

#[test]
fn native_normalization_factor_dispersion_intermediates_match_deseq2_reference() {
    let reference_rows = read_reference_tsv("native_nf_dispersion_reference.tsv");
    let mu_rows = read_reference_tsv("native_nf_mu.tsv");
    let factor_rows = read_reference_tsv("normalization_factors.tsv");

    let counts = reference_counts();
    let samples = reference_sample_names();
    let factor_values = factor_rows
        .iter()
        .flat_map(|row| samples.iter().map(|sample| parse_required_f64(row, sample)))
        .collect::<Vec<_>>();
    let normalization_factors =
        RowMajorMatrix::from_row_major(counts.n_genes(), counts.n_samples(), factor_values)
            .unwrap();
    let normalized = normalized_counts_with_factors(&counts, &normalization_factors).unwrap();
    let base_mean = base_mean(&normalized).unwrap();
    let base_var = base_variance(&normalized).unwrap();
    let all_zero = counts.all_zero_flags();
    let size_factors = vec![1.0; counts.n_samples()];
    let rough = rough_dispersion_estimates(&normalized, &reference_full_design()).unwrap();
    let moments = moments_dispersion_estimates_with_normalization_factors(
        &base_mean,
        &base_var,
        &normalization_factors,
        Some(&all_zero),
    )
    .unwrap();
    let initial = initial_dispersion_estimates(&rough, &moments, 1e-8, 10.0).unwrap();
    let output = estimate_gene_wise_dispersions_linear_mu(
        GeneWiseDispersionInput {
            counts: &counts,
            design: &reference_full_design(),
            size_factors: &size_factors,
            normalization_factors: Some(&normalization_factors),
            normalized_counts: &normalized,
            base_mean: &base_mean,
            base_var: &base_var,
            all_zero: &all_zero,
            observation_weights: None,
        },
        GeneWiseDispersionOptions {
            fit_method: GeneWiseDispersionFitMethod::Grid,
            use_cox_reid: false,
            ..GeneWiseDispersionOptions::default()
        },
    )
    .unwrap();

    assert_eq!(reference_rows.len(), counts.n_genes());
    for (gene, row) in reference_rows.iter().enumerate() {
        assert_eq!(
            row.get("gene").map(String::as_str),
            Some(reference_gene_names()[gene].as_str())
        );
        assert_float_close(
            base_mean[gene],
            parse_required_f64(row, "baseMean"),
            1e-10,
            1e-10,
            &format!("native NF baseMean gene {gene}"),
        );
        assert_float_close(
            base_var[gene],
            parse_required_f64(row, "baseVar"),
            1e-10,
            1e-10,
            &format!("native NF baseVar gene {gene}"),
        );
        assert_eq!(all_zero[gene], parse_required_bool(row, "allZero"));
        assert_float_close(
            rough[gene],
            parse_required_f64(row, "roughDisp"),
            1e-8,
            1e-8,
            &format!("native NF roughDisp gene {gene}"),
        );
        assert_float_close(
            moments[gene],
            parse_required_f64(row, "momentsDisp"),
            1e-10,
            1e-10,
            &format!("native NF momentsDisp gene {gene}"),
        );
        assert_float_close(
            initial[gene],
            parse_required_f64(row, "dispInit"),
            1e-10,
            1e-10,
            &format!("native NF dispInit gene {gene}"),
        );
    }

    for (gene, row) in mu_rows.iter().enumerate() {
        assert_eq!(
            row.get("gene").map(String::as_str),
            Some(reference_gene_names()[gene].as_str())
        );
        for (sample, sample_name) in samples.iter().enumerate() {
            assert_float_close(
                output.mu.row(gene).unwrap()[sample],
                parse_required_f64(row, sample_name),
                1e-8,
                1e-8,
                &format!("native NF mu gene {gene} sample {sample}"),
            );
        }
    }
}
