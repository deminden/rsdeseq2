use super::*;

#[test]
fn native_weighted_glm_mu_wald_matches_deseq2_reference() {
    let rows = read_reference_tsv("native_weighted_glm_mu_reference.tsv");
    let wald_mu_rows = read_reference_tsv("native_weighted_glm_mu_wald_mu.tsv");
    let wald_hat_rows = read_reference_tsv("native_weighted_glm_mu_wald_hat.tsv");
    let size_factors = read_size_factors("size_factors_ratio.tsv");
    let observation_weights = read_reference_matrix("observation_weights.tsv");

    let counts = reference_counts();
    let design = reference_full_design();
    let (fit, results) = DeseqBuilder::new()
        .size_factors(size_factors)
        .observation_weights(observation_weights)
        .fit_type(FitType::Mean)
        .gene_wise_dispersion_options(GeneWiseDispersionOptions {
            use_cox_reid: false,
            niter: 2,
            ..GeneWiseDispersionOptions::default()
        })
        .disable_cooks_cutoff()
        .disable_independent_filtering()
        .fit_wald_glm_mu(&counts, &design, 1)
        .unwrap();

    let beta = fit.beta.as_ref().unwrap();
    let beta_se = fit.beta_se.as_ref().unwrap();
    let beta_converged = fit.beta_converged.as_ref().unwrap();
    let beta_iter = fit.beta_iter.as_ref().unwrap();
    let log_like = fit.log_like.as_ref().unwrap();
    let dispersion = fit.dispersion.as_ref().unwrap();
    let mu = fit.mu.as_ref().unwrap();
    let hat = fit.hat_diagonal.as_ref().unwrap();
    let wald = fit.wald.as_ref().unwrap();
    let samples = reference_sample_names();
    let diagnostics = fit.deseq2_mcols_diagnostics();
    let disp_gene_iter = diagnostics.disp_gene_iter.as_ref().unwrap();

    assert_eq!(
        diagnostics.disp_gene_iter.as_ref(),
        fit.disp_gene_iter.as_ref()
    );

    assert_eq!(rows.len(), results.rows.len());
    for (gene, row) in rows.iter().enumerate() {
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
            &format!("native weighted GLM-mu Wald baseMean gene {gene}"),
        );
        assert_f64_or_missing(
            dispersion[gene],
            parse_optional_f64(row, "dispersion"),
            &format!("native weighted GLM-mu Wald dispersion gene {gene}"),
        );
        assert_f64_or_missing(
            beta.row(gene).unwrap()[0],
            parse_optional_f64(row, "beta_intercept"),
            &format!("native weighted GLM-mu Wald beta intercept gene {gene}"),
        );
        assert_f64_or_missing(
            beta.row(gene).unwrap()[1],
            parse_optional_f64(row, "beta_conditionB"),
            &format!("native weighted GLM-mu Wald beta conditionB gene {gene}"),
        );
        assert_f64_or_missing(
            beta_se.row(gene).unwrap()[0],
            parse_optional_f64(row, "beta_se_intercept"),
            &format!("native weighted GLM-mu Wald beta SE intercept gene {gene}"),
        );
        assert_f64_or_missing(
            beta_se.row(gene).unwrap()[1],
            parse_optional_f64(row, "beta_se_conditionB"),
            &format!("native weighted GLM-mu Wald beta SE conditionB gene {gene}"),
        );
        assert_option_close(
            wald.stat[gene],
            parse_optional_f64(row, "stat_conditionB"),
            1e-4,
            1e-4,
            &format!("native weighted GLM-mu Wald stat gene {gene}"),
        );
        assert_option_close(
            wald.pvalue[gene],
            parse_optional_f64(row, "pvalue_conditionB"),
            1e-4,
            1e-4,
            &format!("native weighted GLM-mu Wald pvalue gene {gene}"),
        );
        assert_option_close(
            results.rows[gene].pvalue,
            parse_optional_f64(row, "pvalue_conditionB"),
            1e-4,
            1e-4,
            &format!("native weighted GLM-mu result pvalue gene {gene}"),
        );
        assert_option_close(
            results.rows[gene].padj,
            parse_optional_f64(row, "padj_conditionB"),
            1e-4,
            1e-4,
            &format!("native weighted GLM-mu result padj gene {gene}"),
        );
        if row.contains_key("log_like") {
            assert_f64_or_missing(
                log_like[gene],
                parse_optional_f64(row, "log_like"),
                &format!("native weighted GLM-mu log-likelihood gene {gene}"),
            );
        }
        if !fit.all_zero[gene] {
            assert!(
                disp_gene_iter[gene] > 0,
                "native weighted GLM-mu gene-wise iterations gene {gene}"
            );
            assert!(
                parse_required_f64(row, "dispGeneIter") > 0.0,
                "DESeq2 native weighted GLM-mu gene-wise iterations gene {gene}"
            );
            assert_eq!(
                beta_converged[gene],
                parse_required_bool(row, "beta_converged"),
                "native weighted GLM-mu beta convergence gene {gene}"
            );
            assert_eq!(
                beta_iter[gene],
                parse_required_f64(row, "beta_iterations") as usize,
                "native weighted GLM-mu beta iterations gene {gene}"
            );
        }
    }

    for (gene, (mu_row, hat_row)) in wald_mu_rows.iter().zip(wald_hat_rows.iter()).enumerate() {
        for (sample, sample_name) in samples.iter().enumerate() {
            assert_f64_or_missing(
                mu.row(gene).unwrap()[sample],
                parse_optional_f64(mu_row, sample_name),
                &format!("native weighted GLM-mu Wald mu gene {gene} sample {sample}"),
            );
            assert_f64_or_missing(
                hat.row(gene).unwrap()[sample],
                parse_optional_f64(hat_row, sample_name),
                &format!("native weighted GLM-mu Wald hat gene {gene} sample {sample}"),
            );
        }
    }
}

#[test]
fn native_weighted_glm_mu_local_wald_matches_deseq2_reference() {
    let rows = read_reference_tsv("native_weighted_glm_mu_local_reference.tsv");
    let wald_mu_rows = read_reference_tsv("native_weighted_glm_mu_local_wald_mu.tsv");
    let wald_hat_rows = read_reference_tsv("native_weighted_glm_mu_local_wald_hat.tsv");
    let size_factors = read_size_factors("size_factors_ratio.tsv");
    let observation_weights = read_reference_matrix("observation_weights.tsv");

    let counts = reference_counts();
    let design = reference_full_design();
    let (fit, results) = DeseqBuilder::new()
        .size_factors(size_factors)
        .observation_weights(observation_weights)
        .fit_type(FitType::Local)
        .gene_wise_dispersion_options(GeneWiseDispersionOptions {
            use_cox_reid: false,
            niter: 2,
            ..GeneWiseDispersionOptions::default()
        })
        .disable_cooks_cutoff()
        .disable_independent_filtering()
        .fit_wald_glm_mu(&counts, &design, 1)
        .unwrap();

    let beta = fit.beta.as_ref().unwrap();
    let beta_se = fit.beta_se.as_ref().unwrap();
    let beta_converged = fit.beta_converged.as_ref().unwrap();
    let beta_iter = fit.beta_iter.as_ref().unwrap();
    let log_like = fit.log_like.as_ref().unwrap();
    let dispersion = fit.dispersion.as_ref().unwrap();
    let mu = fit.mu.as_ref().unwrap();
    let hat = fit.hat_diagonal.as_ref().unwrap();
    let wald = fit.wald.as_ref().unwrap();
    let samples = reference_sample_names();
    let diagnostics = fit.deseq2_mcols_diagnostics();
    let disp_gene_iter = diagnostics.disp_gene_iter.as_ref().unwrap();

    assert_eq!(
        diagnostics.disp_gene_iter.as_ref(),
        fit.disp_gene_iter.as_ref()
    );

    assert_eq!(rows.len(), results.rows.len());
    for (gene, row) in rows.iter().enumerate() {
        let all_zero = parse_required_bool(row, "allZero");
        let weights_fail = parse_required_bool(row, "weightsFail");
        let skipped = all_zero || weights_fail;
        assert_eq!(fit.all_zero[gene], skipped);
        assert_eq!(fit.weights_fail.as_ref().unwrap()[gene], weights_fail);
        assert_float_close(
            fit.base_mean[gene],
            parse_required_f64(row, "baseMean"),
            1e-10,
            1e-10,
            &format!("native weighted GLM-mu local Wald baseMean gene {gene}"),
        );
        assert_f64_or_missing(
            dispersion[gene],
            parse_optional_f64(row, "dispersion"),
            &format!("native weighted GLM-mu local Wald dispersion gene {gene}"),
        );
        assert_f64_or_missing(
            beta.row(gene).unwrap()[0],
            parse_optional_f64(row, "beta_intercept"),
            &format!("native weighted GLM-mu local Wald beta intercept gene {gene}"),
        );
        assert_f64_or_missing(
            beta.row(gene).unwrap()[1],
            parse_optional_f64(row, "beta_conditionB"),
            &format!("native weighted GLM-mu local Wald beta conditionB gene {gene}"),
        );
        assert_f64_or_missing(
            beta_se.row(gene).unwrap()[0],
            parse_optional_f64(row, "beta_se_intercept"),
            &format!("native weighted GLM-mu local Wald beta SE intercept gene {gene}"),
        );
        assert_f64_or_missing(
            beta_se.row(gene).unwrap()[1],
            parse_optional_f64(row, "beta_se_conditionB"),
            &format!("native weighted GLM-mu local Wald beta SE conditionB gene {gene}"),
        );
        assert_option_close(
            wald.stat[gene],
            parse_optional_f64(row, "stat_conditionB"),
            1e-4,
            1e-4,
            &format!("native weighted GLM-mu local Wald stat gene {gene}"),
        );
        assert_option_close(
            wald.pvalue[gene],
            parse_optional_f64(row, "pvalue_conditionB"),
            1e-4,
            1e-4,
            &format!("native weighted GLM-mu local Wald pvalue gene {gene}"),
        );
        assert_option_close(
            results.rows[gene].pvalue,
            parse_optional_f64(row, "pvalue_conditionB"),
            1e-4,
            1e-4,
            &format!("native weighted GLM-mu local result pvalue gene {gene}"),
        );
        assert_option_close(
            results.rows[gene].padj,
            parse_optional_f64(row, "padj_conditionB"),
            1e-4,
            1e-4,
            &format!("native weighted GLM-mu local result padj gene {gene}"),
        );
        if row.contains_key("log_like") {
            assert_f64_or_missing(
                log_like[gene],
                parse_optional_f64(row, "log_like"),
                &format!("native weighted GLM-mu local log-likelihood gene {gene}"),
            );
        }
        if !skipped {
            assert!(
                disp_gene_iter[gene] > 0,
                "native weighted GLM-mu local gene-wise iterations gene {gene}"
            );
            assert!(
                parse_required_f64(row, "dispGeneIter") > 0.0,
                "DESeq2 native weighted GLM-mu local gene-wise iterations gene {gene}"
            );
            assert_eq!(
                beta_converged[gene],
                parse_required_bool(row, "beta_converged"),
                "native weighted GLM-mu local beta convergence gene {gene}"
            );
            assert_eq!(
                beta_iter[gene],
                parse_required_f64(row, "beta_iterations") as usize,
                "native weighted GLM-mu local beta iterations gene {gene}"
            );
        }
    }

    for (gene, (mu_row, hat_row)) in wald_mu_rows.iter().zip(wald_hat_rows.iter()).enumerate() {
        for (sample, sample_name) in samples.iter().enumerate() {
            assert_f64_or_missing(
                mu.row(gene).unwrap()[sample],
                parse_optional_f64(mu_row, sample_name),
                &format!("native weighted GLM-mu local Wald mu gene {gene} sample {sample}"),
            );
            assert_f64_or_missing(
                hat.row(gene).unwrap()[sample],
                parse_optional_f64(hat_row, sample_name),
                &format!("native weighted GLM-mu local Wald hat gene {gene} sample {sample}"),
            );
        }
    }
}

#[test]
fn native_weighted_glm_mu_mean_cox_reid_wald_matches_deseq2_reference() {
    let rows = read_reference_tsv("native_weighted_glm_mu_mean_cr_map_reference.tsv");
    let wald_mu_rows = read_reference_tsv("native_weighted_glm_mu_mean_cr_wald_mu.tsv");
    let wald_hat_rows = read_reference_tsv("native_weighted_glm_mu_mean_cr_wald_hat.tsv");
    let size_factors = read_size_factors("size_factors_ratio.tsv");
    let observation_weights = read_reference_matrix("observation_weights.tsv");

    let counts = reference_counts();
    let design = reference_full_design();
    let (fit, results) = DeseqBuilder::new()
        .size_factors(size_factors)
        .observation_weights(observation_weights)
        .fit_type(FitType::Mean)
        .gene_wise_dispersion_options(GeneWiseDispersionOptions {
            niter: 2,
            ..GeneWiseDispersionOptions::default()
        })
        .disable_cooks_cutoff()
        .disable_independent_filtering()
        .fit_wald_glm_mu(&counts, &design, 1)
        .unwrap();

    let beta = fit.beta.as_ref().unwrap();
    let beta_se = fit.beta_se.as_ref().unwrap();
    let beta_converged = fit.beta_converged.as_ref().unwrap();
    let beta_iter = fit.beta_iter.as_ref().unwrap();
    let log_like = fit.log_like.as_ref().unwrap();
    let dispersion = fit.dispersion.as_ref().unwrap();
    let mu = fit.mu.as_ref().unwrap();
    let hat = fit.hat_diagonal.as_ref().unwrap();
    let wald = fit.wald.as_ref().unwrap();
    let samples = reference_sample_names();
    let diagnostics = fit.deseq2_mcols_diagnostics();
    let disp_gene_iter = diagnostics.disp_gene_iter.as_ref().unwrap();

    assert_eq!(
        diagnostics.disp_gene_iter.as_ref(),
        fit.disp_gene_iter.as_ref()
    );

    assert_eq!(rows.len(), results.rows.len());
    for (gene, row) in rows.iter().enumerate() {
        let all_zero = parse_required_bool(row, "allZero");
        let weights_fail = parse_required_bool(row, "weightsFail");
        let skipped = all_zero || weights_fail;

        assert_eq!(fit.all_zero[gene], skipped);
        assert_eq!(fit.weights_fail.as_ref().unwrap()[gene], weights_fail);
        assert_float_close(
            fit.base_mean[gene],
            parse_required_f64(row, "baseMean"),
            1e-10,
            1e-10,
            &format!("native weighted GLM-mu Cox-Reid Wald baseMean gene {gene}"),
        );
        assert_f64_or_missing(
            dispersion[gene],
            parse_optional_f64(row, "dispersion"),
            &format!("native weighted GLM-mu Cox-Reid Wald dispersion gene {gene}"),
        );
        assert_f64_or_missing(
            beta.row(gene).unwrap()[0],
            parse_optional_f64(row, "beta_intercept"),
            &format!("native weighted GLM-mu Cox-Reid Wald beta intercept gene {gene}"),
        );
        assert_f64_or_missing(
            beta.row(gene).unwrap()[1],
            parse_optional_f64(row, "beta_conditionB"),
            &format!("native weighted GLM-mu Cox-Reid Wald beta conditionB gene {gene}"),
        );
        assert_f64_or_missing(
            beta_se.row(gene).unwrap()[0],
            parse_optional_f64(row, "beta_se_intercept"),
            &format!("native weighted GLM-mu Cox-Reid Wald beta SE intercept gene {gene}"),
        );
        assert_f64_or_missing(
            beta_se.row(gene).unwrap()[1],
            parse_optional_f64(row, "beta_se_conditionB"),
            &format!("native weighted GLM-mu Cox-Reid Wald beta SE conditionB gene {gene}"),
        );
        assert_option_close(
            wald.stat[gene],
            parse_optional_f64(row, "stat_conditionB"),
            1e-4,
            1e-4,
            &format!("native weighted GLM-mu Cox-Reid Wald stat gene {gene}"),
        );
        assert_option_close(
            wald.pvalue[gene],
            parse_optional_f64(row, "pvalue_conditionB"),
            1e-4,
            1e-4,
            &format!("native weighted GLM-mu Cox-Reid Wald pvalue gene {gene}"),
        );
        assert_option_close(
            results.rows[gene].pvalue,
            parse_optional_f64(row, "pvalue_conditionB"),
            1e-4,
            1e-4,
            &format!("native weighted GLM-mu Cox-Reid result pvalue gene {gene}"),
        );
        assert_option_close(
            results.rows[gene].padj,
            parse_optional_f64(row, "padj_conditionB"),
            1e-4,
            1e-4,
            &format!("native weighted GLM-mu Cox-Reid result padj gene {gene}"),
        );
        if row.contains_key("log_like") {
            assert_f64_or_missing(
                log_like[gene],
                parse_optional_f64(row, "log_like"),
                &format!("native weighted GLM-mu Cox-Reid log-likelihood gene {gene}"),
            );
        }
        if !skipped {
            assert!(
                disp_gene_iter[gene] > 0,
                "native weighted GLM-mu Cox-Reid gene-wise iterations gene {gene}"
            );
            assert_eq!(
                beta_converged[gene],
                parse_required_bool(row, "beta_converged"),
                "native weighted GLM-mu Cox-Reid beta convergence gene {gene}"
            );
            assert_eq!(
                beta_iter[gene],
                parse_required_f64(row, "beta_iterations") as usize,
                "native weighted GLM-mu Cox-Reid beta iterations gene {gene}"
            );
        }
    }

    for (gene, (mu_row, hat_row)) in wald_mu_rows.iter().zip(wald_hat_rows.iter()).enumerate() {
        for (sample, sample_name) in samples.iter().enumerate() {
            assert_f64_or_missing(
                mu.row(gene).unwrap()[sample],
                parse_optional_f64(mu_row, sample_name),
                &format!("native weighted GLM-mu Cox-Reid Wald mu gene {gene} sample {sample}"),
            );
            assert_f64_or_missing(
                hat.row(gene).unwrap()[sample],
                parse_optional_f64(hat_row, sample_name),
                &format!("native weighted GLM-mu Cox-Reid Wald hat gene {gene} sample {sample}"),
            );
        }
    }
}

#[test]
fn native_weighted_glm_mu_local_cox_reid_wald_matches_deseq2_reference() {
    let rows = read_reference_tsv("native_weighted_glm_mu_local_cr_map_reference.tsv");
    let wald_mu_rows = read_reference_tsv("native_weighted_glm_mu_local_cr_wald_mu.tsv");
    let wald_hat_rows = read_reference_tsv("native_weighted_glm_mu_local_cr_wald_hat.tsv");
    let size_factors = read_size_factors("size_factors_ratio.tsv");
    let observation_weights = read_reference_matrix("observation_weights.tsv");

    let counts = reference_counts();
    let design = reference_full_design();
    let (fit, results) = DeseqBuilder::new()
        .size_factors(size_factors)
        .observation_weights(observation_weights)
        .fit_type(FitType::Local)
        .gene_wise_dispersion_options(GeneWiseDispersionOptions {
            niter: 2,
            ..GeneWiseDispersionOptions::default()
        })
        .disable_cooks_cutoff()
        .disable_independent_filtering()
        .fit_wald_glm_mu(&counts, &design, 1)
        .unwrap();

    let beta = fit.beta.as_ref().unwrap();
    let beta_se = fit.beta_se.as_ref().unwrap();
    let beta_converged = fit.beta_converged.as_ref().unwrap();
    let beta_iter = fit.beta_iter.as_ref().unwrap();
    let log_like = fit.log_like.as_ref().unwrap();
    let dispersion = fit.dispersion.as_ref().unwrap();
    let mu = fit.mu.as_ref().unwrap();
    let hat = fit.hat_diagonal.as_ref().unwrap();
    let wald = fit.wald.as_ref().unwrap();
    let samples = reference_sample_names();
    let diagnostics = fit.deseq2_mcols_diagnostics();
    let disp_gene_iter = diagnostics.disp_gene_iter.as_ref().unwrap();

    assert_eq!(
        diagnostics.disp_gene_iter.as_ref(),
        fit.disp_gene_iter.as_ref()
    );

    assert_eq!(rows.len(), results.rows.len());
    for (gene, row) in rows.iter().enumerate() {
        let all_zero = parse_required_bool(row, "allZero");
        let weights_fail = parse_required_bool(row, "weightsFail");
        let skipped = all_zero || weights_fail;

        assert_eq!(fit.all_zero[gene], skipped);
        assert_eq!(fit.weights_fail.as_ref().unwrap()[gene], weights_fail);
        assert_float_close(
            fit.base_mean[gene],
            parse_required_f64(row, "baseMean"),
            1e-10,
            1e-10,
            &format!("native weighted GLM-mu Cox-Reid local Wald baseMean gene {gene}"),
        );
        assert_f64_or_missing(
            dispersion[gene],
            parse_optional_f64(row, "dispersion"),
            &format!("native weighted GLM-mu Cox-Reid local Wald dispersion gene {gene}"),
        );
        assert_f64_or_missing(
            beta.row(gene).unwrap()[0],
            parse_optional_f64(row, "beta_intercept"),
            &format!("native weighted GLM-mu Cox-Reid local Wald beta intercept gene {gene}"),
        );
        assert_f64_or_missing(
            beta.row(gene).unwrap()[1],
            parse_optional_f64(row, "beta_conditionB"),
            &format!("native weighted GLM-mu Cox-Reid local Wald beta conditionB gene {gene}"),
        );
        assert_f64_or_missing(
            beta_se.row(gene).unwrap()[0],
            parse_optional_f64(row, "beta_se_intercept"),
            &format!("native weighted GLM-mu Cox-Reid local Wald beta SE intercept gene {gene}"),
        );
        assert_f64_or_missing(
            beta_se.row(gene).unwrap()[1],
            parse_optional_f64(row, "beta_se_conditionB"),
            &format!("native weighted GLM-mu Cox-Reid local Wald beta SE conditionB gene {gene}"),
        );
        assert_option_close(
            wald.stat[gene],
            parse_optional_f64(row, "stat_conditionB"),
            1e-4,
            1e-4,
            &format!("native weighted GLM-mu Cox-Reid local Wald stat gene {gene}"),
        );
        assert_option_close(
            wald.pvalue[gene],
            parse_optional_f64(row, "pvalue_conditionB"),
            1e-4,
            1e-4,
            &format!("native weighted GLM-mu Cox-Reid local Wald pvalue gene {gene}"),
        );
        assert_option_close(
            results.rows[gene].pvalue,
            parse_optional_f64(row, "pvalue_conditionB"),
            1e-4,
            1e-4,
            &format!("native weighted GLM-mu Cox-Reid local result pvalue gene {gene}"),
        );
        assert_option_close(
            results.rows[gene].padj,
            parse_optional_f64(row, "padj_conditionB"),
            1e-4,
            1e-4,
            &format!("native weighted GLM-mu Cox-Reid local result padj gene {gene}"),
        );
        if row.contains_key("log_like") {
            assert_f64_or_missing(
                log_like[gene],
                parse_optional_f64(row, "log_like"),
                &format!("native weighted GLM-mu Cox-Reid local log-likelihood gene {gene}"),
            );
        }
        if !skipped {
            assert!(
                disp_gene_iter[gene] > 0,
                "native weighted GLM-mu Cox-Reid local gene-wise iterations gene {gene}"
            );
            assert_eq!(
                beta_converged[gene],
                parse_required_bool(row, "beta_converged"),
                "native weighted GLM-mu Cox-Reid local beta convergence gene {gene}"
            );
            assert_eq!(
                beta_iter[gene],
                parse_required_f64(row, "beta_iterations") as usize,
                "native weighted GLM-mu Cox-Reid local beta iterations gene {gene}"
            );
        }
    }

    for (gene, (mu_row, hat_row)) in wald_mu_rows.iter().zip(wald_hat_rows.iter()).enumerate() {
        for (sample, sample_name) in samples.iter().enumerate() {
            assert_f64_or_missing(
                mu.row(gene).unwrap()[sample],
                parse_optional_f64(mu_row, sample_name),
                &format!(
                    "native weighted GLM-mu Cox-Reid local Wald mu gene {gene} sample {sample}"
                ),
            );
            assert_f64_or_missing(
                hat.row(gene).unwrap()[sample],
                parse_optional_f64(hat_row, sample_name),
                &format!(
                    "native weighted GLM-mu Cox-Reid local Wald hat gene {gene} sample {sample}"
                ),
            );
        }
    }
}
