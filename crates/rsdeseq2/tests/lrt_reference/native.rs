use super::*;

#[test]
fn native_glm_mu_mean_lrt_matches_deseq2_reference() {
    let rows = read_reference_tsv("native_glm_mu_mean_lrt_reference.tsv");
    let size_factors = read_size_factors("size_factors_ratio.tsv");

    let counts = reference_counts();
    let full_design = reference_full_design();
    let reduced_design = reference_reduced_design();
    let (fit, results) = DeseqBuilder::new()
        .size_factors(size_factors)
        .fit_type(FitType::Mean)
        .gene_wise_dispersion_options(GeneWiseDispersionOptions {
            use_cox_reid: false,
            niter: 2,
            ..GeneWiseDispersionOptions::default()
        })
        .disable_cooks_cutoff()
        .disable_independent_filtering()
        .fit_lrt_glm_mu(&counts, &full_design, &reduced_design, 1)
        .unwrap();

    let beta = fit.beta.as_ref().unwrap();
    let beta_se = fit.beta_se.as_ref().unwrap();
    let beta_converged = fit.beta_converged.as_ref().unwrap();
    let beta_iter = fit.beta_iter.as_ref().unwrap();
    let log_like = fit.log_like.as_ref().unwrap();
    let reduced_log_like = fit.reduced_log_like.as_ref().unwrap();
    let reduced_beta_converged = fit.reduced_beta_converged.as_ref().unwrap();
    let reduced_beta_iter = fit.reduced_beta_iter.as_ref().unwrap();
    let dispersion = fit.dispersion.as_ref().unwrap();
    let lrt = fit.lrt.as_ref().unwrap();
    let diagnostics = fit.deseq2_mcols_diagnostics();
    let disp_gene_iter = diagnostics.disp_gene_iter.as_ref().unwrap();

    assert_eq!(lrt.degrees_of_freedom, 1);
    assert_eq!(rows.len(), results.rows.len());
    assert_eq!(rows.len(), beta.n_rows());
    assert_eq!(reduced_beta_converged, &lrt.reduced_converged);
    assert_eq!(reduced_beta_iter.len(), rows.len());
    assert_eq!(
        diagnostics.disp_gene_iter.as_ref(),
        fit.disp_gene_iter.as_ref()
    );
    for (gene, row) in rows.iter().enumerate() {
        assert_eq!(fit.all_zero[gene], parse_required_bool(row, "allZero"));
        assert_float_close(
            fit.base_mean[gene],
            parse_required_f64(row, "baseMean"),
            1e-10,
            1e-10,
            &format!("native GLM-mu mean LRT baseMean gene {gene}"),
        );
        assert_f64_or_missing(
            dispersion[gene],
            parse_optional_f64(row, "dispersion"),
            &format!("native GLM-mu mean LRT dispersion gene {gene}"),
        );
        assert_f64_or_missing(
            beta.row(gene).unwrap()[0],
            parse_optional_f64(row, "beta_intercept"),
            &format!("native GLM-mu mean LRT full beta intercept gene {gene}"),
        );
        assert_f64_or_missing(
            beta.row(gene).unwrap()[1],
            parse_optional_f64(row, "beta_conditionB"),
            &format!("native GLM-mu mean LRT full beta conditionB gene {gene}"),
        );
        assert_f64_or_missing(
            beta_se.row(gene).unwrap()[0],
            parse_optional_f64(row, "beta_se_intercept"),
            &format!("native GLM-mu mean LRT full beta SE intercept gene {gene}"),
        );
        assert_f64_or_missing(
            beta_se.row(gene).unwrap()[1],
            parse_optional_f64(row, "beta_se_conditionB"),
            &format!("native GLM-mu mean LRT full beta SE conditionB gene {gene}"),
        );
        assert_option_close(
            lrt.deviance[gene],
            parse_optional_f64(row, "lrt_stat"),
            1e-5,
            1e-5,
            &format!("native GLM-mu mean LRT statistic gene {gene}"),
        );
        assert_f64_or_missing(
            log_like[gene],
            parse_optional_f64(row, "log_like_full"),
            &format!("native GLM-mu mean LRT full log-likelihood gene {gene}"),
        );
        assert_f64_or_missing(
            fit.full_deviance.as_ref().unwrap()[gene],
            parse_optional_f64(row, "log_like_full").map(|log_like| -2.0 * log_like),
            &format!("native GLM-mu mean LRT full deviance gene {gene}"),
        );
        assert_f64_or_missing(
            reduced_log_like[gene],
            parse_optional_f64(row, "log_like_reduced"),
            &format!("native GLM-mu mean LRT reduced log-likelihood gene {gene}"),
        );
        assert_option_close(
            lrt.pvalue[gene],
            parse_optional_f64(row, "pvalue"),
            1e-7,
            1e-5,
            &format!("native GLM-mu mean LRT p-value gene {gene}"),
        );
        assert_option_close(
            results.rows[gene].stat,
            parse_optional_f64(row, "lrt_stat"),
            1e-5,
            1e-5,
            &format!("native GLM-mu mean LRT result statistic gene {gene}"),
        );
        assert_eq!(
            results.rows[gene].pvalue, lrt.pvalue[gene],
            "native GLM-mu mean LRT result p-value gene {gene}"
        );
        assert_option_close(
            results.rows[gene].padj,
            parse_optional_f64(row, "padj"),
            1e-7,
            1e-5,
            &format!("native GLM-mu mean LRT result padj gene {gene}"),
        );
        if fit.all_zero[gene] {
            assert_eq!(disp_gene_iter[gene], 0);
            assert_eq!(beta_iter[gene], 0);
            assert_eq!(reduced_beta_iter[gene], 0);
        } else {
            assert!(
                disp_gene_iter[gene] > 0,
                "native GLM-mu mean LRT gene-wise iterations gene {gene}"
            );
            assert_eq!(
                beta_converged[gene],
                parse_required_bool(row, "full_converged"),
                "native GLM-mu mean LRT full convergence gene {gene}"
            );
            assert_eq!(
                reduced_beta_converged[gene],
                parse_required_bool(row, "reduced_converged"),
                "native GLM-mu mean LRT reduced convergence gene {gene}"
            );
            assert_eq!(
                beta_iter[gene],
                parse_required_f64(row, "full_iterations") as usize,
                "native GLM-mu mean LRT full iterations gene {gene}"
            );
            assert_eq!(
                reduced_beta_iter[gene],
                parse_required_f64(row, "reduced_iterations") as usize,
                "native GLM-mu mean LRT reduced iterations gene {gene}"
            );
        }
    }
}

#[test]
fn native_glm_mu_local_lrt_matches_deseq2_reference() {
    let rows = read_reference_tsv("native_glm_mu_local_lrt_reference.tsv");
    let size_factors = read_size_factors("size_factors_ratio.tsv");

    let counts = reference_counts();
    let full_design = reference_full_design();
    let reduced_design = reference_reduced_design();
    let (fit, results) = DeseqBuilder::new()
        .size_factors(size_factors)
        .fit_type(FitType::Local)
        .gene_wise_dispersion_options(GeneWiseDispersionOptions {
            use_cox_reid: false,
            niter: 2,
            ..GeneWiseDispersionOptions::default()
        })
        .disable_cooks_cutoff()
        .disable_independent_filtering()
        .fit_lrt_glm_mu(&counts, &full_design, &reduced_design, 1)
        .unwrap();

    let beta = fit.beta.as_ref().unwrap();
    let beta_se = fit.beta_se.as_ref().unwrap();
    let beta_converged = fit.beta_converged.as_ref().unwrap();
    let beta_iter = fit.beta_iter.as_ref().unwrap();
    let log_like = fit.log_like.as_ref().unwrap();
    let reduced_log_like = fit.reduced_log_like.as_ref().unwrap();
    let reduced_beta_converged = fit.reduced_beta_converged.as_ref().unwrap();
    let reduced_beta_iter = fit.reduced_beta_iter.as_ref().unwrap();
    let dispersion = fit.dispersion.as_ref().unwrap();
    let lrt = fit.lrt.as_ref().unwrap();
    let diagnostics = fit.deseq2_mcols_diagnostics();
    let disp_gene_iter = diagnostics.disp_gene_iter.as_ref().unwrap();

    assert_eq!(lrt.degrees_of_freedom, 1);
    assert_eq!(rows.len(), results.rows.len());
    assert_eq!(rows.len(), beta.n_rows());
    assert_eq!(reduced_beta_converged, &lrt.reduced_converged);
    assert_eq!(reduced_beta_iter.len(), rows.len());
    assert_eq!(
        diagnostics.disp_gene_iter.as_ref(),
        fit.disp_gene_iter.as_ref()
    );
    for (gene, row) in rows.iter().enumerate() {
        assert_eq!(fit.all_zero[gene], parse_required_bool(row, "allZero"));
        assert_float_close(
            fit.base_mean[gene],
            parse_required_f64(row, "baseMean"),
            1e-10,
            1e-10,
            &format!("native GLM-mu local LRT baseMean gene {gene}"),
        );
        assert_f64_or_missing(
            dispersion[gene],
            parse_optional_f64(row, "dispersion"),
            &format!("native GLM-mu local LRT dispersion gene {gene}"),
        );
        assert_f64_or_missing(
            beta.row(gene).unwrap()[0],
            parse_optional_f64(row, "beta_intercept"),
            &format!("native GLM-mu local LRT full beta intercept gene {gene}"),
        );
        assert_f64_or_missing(
            beta.row(gene).unwrap()[1],
            parse_optional_f64(row, "beta_conditionB"),
            &format!("native GLM-mu local LRT full beta conditionB gene {gene}"),
        );
        assert_f64_or_missing(
            beta_se.row(gene).unwrap()[0],
            parse_optional_f64(row, "beta_se_intercept"),
            &format!("native GLM-mu local LRT full beta SE intercept gene {gene}"),
        );
        assert_f64_or_missing(
            beta_se.row(gene).unwrap()[1],
            parse_optional_f64(row, "beta_se_conditionB"),
            &format!("native GLM-mu local LRT full beta SE conditionB gene {gene}"),
        );
        assert_option_close(
            lrt.deviance[gene],
            parse_optional_f64(row, "lrt_stat"),
            1e-5,
            1e-5,
            &format!("native GLM-mu local LRT statistic gene {gene}"),
        );
        assert_f64_or_missing(
            log_like[gene],
            parse_optional_f64(row, "log_like_full"),
            &format!("native GLM-mu local LRT full log-likelihood gene {gene}"),
        );
        assert_f64_or_missing(
            fit.full_deviance.as_ref().unwrap()[gene],
            parse_optional_f64(row, "log_like_full").map(|log_like| -2.0 * log_like),
            &format!("native GLM-mu local LRT full deviance gene {gene}"),
        );
        assert_f64_or_missing(
            reduced_log_like[gene],
            parse_optional_f64(row, "log_like_reduced"),
            &format!("native GLM-mu local LRT reduced log-likelihood gene {gene}"),
        );
        assert_option_close(
            lrt.pvalue[gene],
            parse_optional_f64(row, "pvalue"),
            1e-7,
            1e-5,
            &format!("native GLM-mu local LRT p-value gene {gene}"),
        );
        assert_option_close(
            results.rows[gene].stat,
            parse_optional_f64(row, "lrt_stat"),
            1e-5,
            1e-5,
            &format!("native GLM-mu local LRT result statistic gene {gene}"),
        );
        assert_eq!(
            results.rows[gene].pvalue, lrt.pvalue[gene],
            "native GLM-mu local LRT result p-value gene {gene}"
        );
        assert_option_close(
            results.rows[gene].padj,
            parse_optional_f64(row, "padj"),
            1e-7,
            1e-5,
            &format!("native GLM-mu local LRT result padj gene {gene}"),
        );
        if fit.all_zero[gene] {
            assert_eq!(disp_gene_iter[gene], 0);
            assert_eq!(beta_iter[gene], 0);
            assert_eq!(reduced_beta_iter[gene], 0);
        } else {
            assert!(
                disp_gene_iter[gene] > 0,
                "native GLM-mu local LRT gene-wise iterations gene {gene}"
            );
            assert_eq!(
                beta_converged[gene],
                parse_required_bool(row, "full_converged"),
                "native GLM-mu local LRT full convergence gene {gene}"
            );
            assert_eq!(
                reduced_beta_converged[gene],
                parse_required_bool(row, "reduced_converged"),
                "native GLM-mu local LRT reduced convergence gene {gene}"
            );
            assert_eq!(
                beta_iter[gene],
                parse_required_f64(row, "full_iterations") as usize,
                "native GLM-mu local LRT full iterations gene {gene}"
            );
            assert_eq!(
                reduced_beta_iter[gene],
                parse_required_f64(row, "reduced_iterations") as usize,
                "native GLM-mu local LRT reduced iterations gene {gene}"
            );
        }
    }
}

#[test]
fn native_glm_mu_mean_cox_reid_lrt_matches_deseq2_reference() {
    let rows = read_reference_tsv("native_glm_mu_mean_cr_lrt_reference.tsv");
    let size_factors = read_size_factors("size_factors_ratio.tsv");

    let counts = reference_counts();
    let full_design = reference_full_design();
    let reduced_design = reference_reduced_design();
    let (fit, results) = DeseqBuilder::new()
        .size_factors(size_factors)
        .fit_type(FitType::Mean)
        .gene_wise_dispersion_options(GeneWiseDispersionOptions {
            niter: 2,
            ..GeneWiseDispersionOptions::default()
        })
        .disable_cooks_cutoff()
        .disable_independent_filtering()
        .fit_lrt_glm_mu(&counts, &full_design, &reduced_design, 1)
        .unwrap();

    let beta = fit.beta.as_ref().unwrap();
    let beta_se = fit.beta_se.as_ref().unwrap();
    let beta_converged = fit.beta_converged.as_ref().unwrap();
    let beta_iter = fit.beta_iter.as_ref().unwrap();
    let log_like = fit.log_like.as_ref().unwrap();
    let reduced_log_like = fit.reduced_log_like.as_ref().unwrap();
    let reduced_beta_converged = fit.reduced_beta_converged.as_ref().unwrap();
    let reduced_beta_iter = fit.reduced_beta_iter.as_ref().unwrap();
    let dispersion = fit.dispersion.as_ref().unwrap();
    let lrt = fit.lrt.as_ref().unwrap();
    let diagnostics = fit.deseq2_mcols_diagnostics();
    let disp_gene_iter = diagnostics.disp_gene_iter.as_ref().unwrap();

    assert_eq!(lrt.degrees_of_freedom, 1);
    assert_eq!(rows.len(), results.rows.len());
    assert_eq!(rows.len(), beta.n_rows());
    assert_eq!(reduced_beta_converged, &lrt.reduced_converged);
    assert_eq!(reduced_beta_iter.len(), rows.len());
    assert_eq!(
        diagnostics.disp_gene_iter.as_ref(),
        fit.disp_gene_iter.as_ref()
    );
    for (gene, row) in rows.iter().enumerate() {
        assert_eq!(fit.all_zero[gene], parse_required_bool(row, "allZero"));
        assert_float_close(
            fit.base_mean[gene],
            parse_required_f64(row, "baseMean"),
            1e-10,
            1e-10,
            &format!("native GLM-mu Cox-Reid mean LRT baseMean gene {gene}"),
        );
        assert_f64_or_missing(
            dispersion[gene],
            parse_optional_f64(row, "dispersion"),
            &format!("native GLM-mu Cox-Reid mean LRT dispersion gene {gene}"),
        );
        assert_f64_or_missing(
            beta.row(gene).unwrap()[0],
            parse_optional_f64(row, "beta_intercept"),
            &format!("native GLM-mu Cox-Reid mean LRT full beta intercept gene {gene}"),
        );
        assert_f64_or_missing(
            beta.row(gene).unwrap()[1],
            parse_optional_f64(row, "beta_conditionB"),
            &format!("native GLM-mu Cox-Reid mean LRT full beta conditionB gene {gene}"),
        );
        assert_f64_or_missing(
            beta_se.row(gene).unwrap()[0],
            parse_optional_f64(row, "beta_se_intercept"),
            &format!("native GLM-mu Cox-Reid mean LRT full beta SE intercept gene {gene}"),
        );
        assert_f64_or_missing(
            beta_se.row(gene).unwrap()[1],
            parse_optional_f64(row, "beta_se_conditionB"),
            &format!("native GLM-mu Cox-Reid mean LRT full beta SE conditionB gene {gene}"),
        );
        assert_option_close(
            lrt.deviance[gene],
            parse_optional_f64(row, "lrt_stat"),
            1e-5,
            1e-5,
            &format!("native GLM-mu Cox-Reid mean LRT statistic gene {gene}"),
        );
        assert_f64_or_missing(
            log_like[gene],
            parse_optional_f64(row, "log_like_full"),
            &format!("native GLM-mu Cox-Reid mean LRT full log-likelihood gene {gene}"),
        );
        assert_f64_or_missing(
            fit.full_deviance.as_ref().unwrap()[gene],
            parse_optional_f64(row, "log_like_full").map(|log_like| -2.0 * log_like),
            &format!("native GLM-mu Cox-Reid mean LRT full deviance gene {gene}"),
        );
        assert_f64_or_missing(
            reduced_log_like[gene],
            parse_optional_f64(row, "log_like_reduced"),
            &format!("native GLM-mu Cox-Reid mean LRT reduced log-likelihood gene {gene}"),
        );
        assert_option_close(
            lrt.pvalue[gene],
            parse_optional_f64(row, "pvalue"),
            1e-7,
            1e-5,
            &format!("native GLM-mu Cox-Reid mean LRT p-value gene {gene}"),
        );
        assert_option_close(
            results.rows[gene].stat,
            parse_optional_f64(row, "lrt_stat"),
            1e-5,
            1e-5,
            &format!("native GLM-mu Cox-Reid mean LRT result statistic gene {gene}"),
        );
        assert_eq!(
            results.rows[gene].pvalue, lrt.pvalue[gene],
            "native GLM-mu Cox-Reid mean LRT result p-value gene {gene}"
        );
        assert_option_close(
            results.rows[gene].padj,
            parse_optional_f64(row, "padj"),
            1e-7,
            1e-5,
            &format!("native GLM-mu Cox-Reid mean LRT result padj gene {gene}"),
        );
        if fit.all_zero[gene] {
            assert_eq!(disp_gene_iter[gene], 0);
            assert_eq!(beta_iter[gene], 0);
            assert_eq!(reduced_beta_iter[gene], 0);
        } else {
            assert!(
                disp_gene_iter[gene] > 0,
                "native GLM-mu Cox-Reid mean LRT gene-wise iterations gene {gene}"
            );
            assert_eq!(
                beta_converged[gene],
                parse_required_bool(row, "full_converged"),
                "native GLM-mu Cox-Reid mean LRT full convergence gene {gene}"
            );
            assert_eq!(
                reduced_beta_converged[gene],
                parse_required_bool(row, "reduced_converged"),
                "native GLM-mu Cox-Reid mean LRT reduced convergence gene {gene}"
            );
            assert_eq!(
                beta_iter[gene],
                parse_required_f64(row, "full_iterations") as usize,
                "native GLM-mu Cox-Reid mean LRT full iterations gene {gene}"
            );
            assert_eq!(
                reduced_beta_iter[gene],
                parse_required_f64(row, "reduced_iterations") as usize,
                "native GLM-mu Cox-Reid mean LRT reduced iterations gene {gene}"
            );
        }
    }
}

#[test]
fn native_glm_mu_local_cox_reid_lrt_matches_deseq2_reference() {
    let rows = read_reference_tsv("native_glm_mu_local_cr_lrt_reference.tsv");
    let size_factors = read_size_factors("size_factors_ratio.tsv");

    let counts = reference_counts();
    let full_design = reference_full_design();
    let reduced_design = reference_reduced_design();
    let (fit, results) = DeseqBuilder::new()
        .size_factors(size_factors)
        .fit_type(FitType::Local)
        .gene_wise_dispersion_options(GeneWiseDispersionOptions {
            niter: 2,
            ..GeneWiseDispersionOptions::default()
        })
        .disable_cooks_cutoff()
        .disable_independent_filtering()
        .fit_lrt_glm_mu(&counts, &full_design, &reduced_design, 1)
        .unwrap();

    let beta = fit.beta.as_ref().unwrap();
    let beta_se = fit.beta_se.as_ref().unwrap();
    let beta_converged = fit.beta_converged.as_ref().unwrap();
    let beta_iter = fit.beta_iter.as_ref().unwrap();
    let log_like = fit.log_like.as_ref().unwrap();
    let reduced_log_like = fit.reduced_log_like.as_ref().unwrap();
    let reduced_beta_converged = fit.reduced_beta_converged.as_ref().unwrap();
    let reduced_beta_iter = fit.reduced_beta_iter.as_ref().unwrap();
    let dispersion = fit.dispersion.as_ref().unwrap();
    let lrt = fit.lrt.as_ref().unwrap();
    let diagnostics = fit.deseq2_mcols_diagnostics();
    let disp_gene_iter = diagnostics.disp_gene_iter.as_ref().unwrap();

    assert_eq!(lrt.degrees_of_freedom, 1);
    assert_eq!(rows.len(), results.rows.len());
    assert_eq!(rows.len(), beta.n_rows());
    assert_eq!(reduced_beta_converged, &lrt.reduced_converged);
    assert_eq!(reduced_beta_iter.len(), rows.len());
    assert_eq!(
        diagnostics.disp_gene_iter.as_ref(),
        fit.disp_gene_iter.as_ref()
    );
    for (gene, row) in rows.iter().enumerate() {
        assert_eq!(fit.all_zero[gene], parse_required_bool(row, "allZero"));
        assert_float_close(
            fit.base_mean[gene],
            parse_required_f64(row, "baseMean"),
            1e-10,
            1e-10,
            &format!("native GLM-mu Cox-Reid local LRT baseMean gene {gene}"),
        );
        assert_f64_or_missing(
            dispersion[gene],
            parse_optional_f64(row, "dispersion"),
            &format!("native GLM-mu Cox-Reid local LRT dispersion gene {gene}"),
        );
        assert_f64_or_missing(
            beta.row(gene).unwrap()[0],
            parse_optional_f64(row, "beta_intercept"),
            &format!("native GLM-mu Cox-Reid local LRT full beta intercept gene {gene}"),
        );
        assert_f64_or_missing(
            beta.row(gene).unwrap()[1],
            parse_optional_f64(row, "beta_conditionB"),
            &format!("native GLM-mu Cox-Reid local LRT full beta conditionB gene {gene}"),
        );
        assert_f64_or_missing(
            beta_se.row(gene).unwrap()[0],
            parse_optional_f64(row, "beta_se_intercept"),
            &format!("native GLM-mu Cox-Reid local LRT full beta SE intercept gene {gene}"),
        );
        assert_f64_or_missing(
            beta_se.row(gene).unwrap()[1],
            parse_optional_f64(row, "beta_se_conditionB"),
            &format!("native GLM-mu Cox-Reid local LRT full beta SE conditionB gene {gene}"),
        );
        assert_option_close(
            lrt.deviance[gene],
            parse_optional_f64(row, "lrt_stat"),
            1e-5,
            1e-5,
            &format!("native GLM-mu Cox-Reid local LRT statistic gene {gene}"),
        );
        assert_f64_or_missing(
            log_like[gene],
            parse_optional_f64(row, "log_like_full"),
            &format!("native GLM-mu Cox-Reid local LRT full log-likelihood gene {gene}"),
        );
        assert_f64_or_missing(
            fit.full_deviance.as_ref().unwrap()[gene],
            parse_optional_f64(row, "log_like_full").map(|log_like| -2.0 * log_like),
            &format!("native GLM-mu Cox-Reid local LRT full deviance gene {gene}"),
        );
        assert_f64_or_missing(
            reduced_log_like[gene],
            parse_optional_f64(row, "log_like_reduced"),
            &format!("native GLM-mu Cox-Reid local LRT reduced log-likelihood gene {gene}"),
        );
        assert_option_close(
            lrt.pvalue[gene],
            parse_optional_f64(row, "pvalue"),
            1e-7,
            1e-5,
            &format!("native GLM-mu Cox-Reid local LRT p-value gene {gene}"),
        );
        assert_option_close(
            results.rows[gene].stat,
            parse_optional_f64(row, "lrt_stat"),
            1e-5,
            1e-5,
            &format!("native GLM-mu Cox-Reid local LRT result statistic gene {gene}"),
        );
        assert_eq!(
            results.rows[gene].pvalue, lrt.pvalue[gene],
            "native GLM-mu Cox-Reid local LRT result p-value gene {gene}"
        );
        assert_option_close(
            results.rows[gene].padj,
            parse_optional_f64(row, "padj"),
            1e-7,
            1e-5,
            &format!("native GLM-mu Cox-Reid local LRT result padj gene {gene}"),
        );
        if fit.all_zero[gene] {
            assert_eq!(disp_gene_iter[gene], 0);
            assert_eq!(beta_iter[gene], 0);
            assert_eq!(reduced_beta_iter[gene], 0);
        } else {
            assert!(
                disp_gene_iter[gene] > 0,
                "native GLM-mu Cox-Reid local LRT gene-wise iterations gene {gene}"
            );
            assert_eq!(
                beta_converged[gene],
                parse_required_bool(row, "full_converged"),
                "native GLM-mu Cox-Reid local LRT full convergence gene {gene}"
            );
            assert_eq!(
                reduced_beta_converged[gene],
                parse_required_bool(row, "reduced_converged"),
                "native GLM-mu Cox-Reid local LRT reduced convergence gene {gene}"
            );
            assert_eq!(
                beta_iter[gene],
                parse_required_f64(row, "full_iterations") as usize,
                "native GLM-mu Cox-Reid local LRT full iterations gene {gene}"
            );
            assert_eq!(
                reduced_beta_iter[gene],
                parse_required_f64(row, "reduced_iterations") as usize,
                "native GLM-mu Cox-Reid local LRT reduced iterations gene {gene}"
            );
        }
    }
}
