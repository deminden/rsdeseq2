use super::*;

#[test]
fn native_glm_mu_cox_reid_gene_wise_matches_deseq2_reference() {
    let reference_rows = read_reference_tsv("native_glm_mu_cr_reference.tsv");
    let mu_rows = read_reference_tsv("native_glm_mu_cr_dispersion_mu.tsv");
    let size_factors = read_size_factors("size_factors_ratio.tsv");

    let counts = reference_counts();
    let design = reference_full_design();
    let fit = DeseqBuilder::new()
        .size_factors(size_factors)
        .gene_wise_dispersion_options(GeneWiseDispersionOptions {
            niter: 2,
            ..GeneWiseDispersionOptions::default()
        })
        .fit_gene_wise_dispersions_glm_mu(&counts, &design)
        .unwrap();
    let disp_gene_est = fit.disp_gene_est.as_ref().unwrap();
    let disp_gene_iter = fit.disp_gene_iter.as_ref().unwrap();
    let mu = fit.mu.as_ref().unwrap();
    let samples = reference_sample_names();

    assert_eq!(reference_rows.len(), counts.n_genes());
    assert_eq!(mu_rows.len(), counts.n_genes());
    for (gene, row) in reference_rows.iter().enumerate() {
        assert_eq!(
            row.get("gene").map(String::as_str),
            Some(reference_gene_names()[gene].as_str())
        );
        assert_eq!(fit.all_zero[gene], parse_required_bool(row, "allZero"));
        assert_float_close(
            fit.base_mean[gene],
            parse_required_f64(row, "baseMean"),
            1e-10,
            1e-10,
            &format!("GLM-mu Cox-Reid baseMean gene {gene}"),
        );
        assert_float_close(
            fit.base_var[gene],
            parse_required_f64(row, "baseVar"),
            1e-10,
            1e-10,
            &format!("GLM-mu Cox-Reid baseVar gene {gene}"),
        );
        assert_f64_or_missing(
            disp_gene_est[gene],
            parse_optional_f64(row, "dispGeneEst"),
            1e-8,
            1e-8,
            &format!("GLM-mu Cox-Reid dispGeneEst gene {gene}"),
        );
        if !fit.all_zero[gene] {
            assert!(
                disp_gene_iter[gene] > 0,
                "GLM-mu Cox-Reid gene-wise iterations gene {gene}"
            );
            assert!(
                parse_required_f64(row, "dispGeneIter") > 0.0,
                "DESeq2 GLM-mu Cox-Reid iterations gene {gene}"
            );
        }
    }

    for (gene, row) in mu_rows.iter().enumerate() {
        for (sample, sample_name) in samples.iter().enumerate() {
            assert_f64_or_missing(
                mu.row(gene).unwrap()[sample],
                parse_optional_f64(row, sample_name),
                1e-4,
                1e-4,
                &format!("GLM-mu Cox-Reid dispersion mu gene {gene} sample {sample}"),
            );
        }
    }
}

#[test]
fn native_weighted_glm_mu_cox_reid_gene_wise_matches_deseq2_reference() {
    let reference_rows = read_reference_tsv("native_weighted_glm_mu_cr_reference.tsv");
    let mu_rows = read_reference_tsv("native_weighted_glm_mu_cr_dispersion_mu.tsv");
    let size_factors = read_size_factors("size_factors_ratio.tsv");
    let observation_weights = read_reference_matrix("observation_weights.tsv");

    let counts = reference_counts();
    let design = reference_full_design();
    let fit = DeseqBuilder::new()
        .size_factors(size_factors)
        .observation_weights(observation_weights)
        .gene_wise_dispersion_options(GeneWiseDispersionOptions {
            niter: 2,
            ..GeneWiseDispersionOptions::default()
        })
        .fit_gene_wise_dispersions_glm_mu(&counts, &design)
        .unwrap();
    let disp_gene_est = fit.disp_gene_est.as_ref().unwrap();
    let disp_gene_iter = fit.disp_gene_iter.as_ref().unwrap();
    let mu = fit.mu.as_ref().unwrap();

    let samples = reference_sample_names();

    assert_eq!(reference_rows.len(), counts.n_genes());
    assert_eq!(mu_rows.len(), counts.n_genes());
    for (gene, row) in reference_rows.iter().enumerate() {
        assert_eq!(
            row.get("gene").map(String::as_str),
            Some(reference_gene_names()[gene].as_str())
        );
        assert_eq!(fit.all_zero[gene], parse_required_bool(row, "allZero"));
        assert_eq!(
            fit.weights_fail.as_ref().unwrap()[gene],
            parse_required_bool(row, "weightsFail")
        );
        assert_float_close(
            fit.base_mean[gene],
            parse_required_f64(row, "baseMean"),
            1e-10,
            1e-10,
            &format!("weighted GLM-mu Cox-Reid baseMean gene {gene}"),
        );
        assert_float_close(
            fit.base_var[gene],
            parse_required_f64(row, "baseVar"),
            1e-10,
            1e-10,
            &format!("weighted GLM-mu Cox-Reid baseVar gene {gene}"),
        );
        assert_f64_or_missing(
            disp_gene_est[gene],
            parse_optional_f64(row, "dispGeneEst"),
            1e-8,
            1e-8,
            &format!("weighted GLM-mu Cox-Reid dispGeneEst gene {gene}"),
        );
        if !fit.all_zero[gene] {
            assert!(
                disp_gene_iter[gene] > 0,
                "weighted GLM-mu Cox-Reid gene-wise iterations gene {gene}"
            );
            assert!(
                parse_required_f64(row, "dispGeneIter") > 0.0,
                "DESeq2 weighted GLM-mu Cox-Reid iterations gene {gene}"
            );
        }
    }

    for (gene, row) in mu_rows.iter().enumerate() {
        for (sample, sample_name) in samples.iter().enumerate() {
            assert_f64_or_missing(
                mu.row(gene).unwrap()[sample],
                parse_optional_f64(row, sample_name),
                1e-4,
                1e-4,
                &format!("weighted GLM-mu Cox-Reid dispersion mu gene {gene} sample {sample}"),
            );
        }
    }
}
