mod common;

use common::*;
use rsdeseq2::prelude::*;

fn assert_f64_or_missing(actual: f64, expected: Option<f64>, label: &str) {
    match expected {
        Some(expected) => assert_float_close(actual, expected, 1e-5, 1e-5, label),
        None => assert!(actual.is_nan(), "{label}: expected missing, got {actual}"),
    }
}

#[test]
fn native_glm_mu_mean_wald_matches_optional_deseq2_reference() {
    let Some(rows) = read_optional_tsv("native_glm_mu_mean_reference.tsv") else {
        return;
    };
    let Some(wald_mu_rows) = read_optional_tsv("native_glm_mu_mean_wald_mu.tsv") else {
        return;
    };
    let Some(wald_hat_rows) = read_optional_tsv("native_glm_mu_mean_wald_hat.tsv") else {
        return;
    };
    let Some(size_factors) = read_size_factors("size_factors_ratio.tsv") else {
        return;
    };

    let counts = reference_counts();
    let design = reference_full_design();
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

    assert_eq!(rows.len(), results.rows.len());
    for (gene, row) in rows.iter().enumerate() {
        assert_eq!(fit.all_zero[gene], parse_required_bool(row, "allZero"));
        assert_float_close(
            fit.base_mean[gene],
            parse_required_f64(row, "baseMean"),
            1e-10,
            1e-10,
            &format!("native GLM-mu mean Wald baseMean gene {gene}"),
        );
        assert_f64_or_missing(
            dispersion[gene],
            parse_optional_f64(row, "dispersion"),
            &format!("native GLM-mu mean Wald dispersion gene {gene}"),
        );
        assert_f64_or_missing(
            beta.row(gene).unwrap()[0],
            parse_optional_f64(row, "beta_intercept"),
            &format!("native GLM-mu mean Wald beta intercept gene {gene}"),
        );
        assert_f64_or_missing(
            beta.row(gene).unwrap()[1],
            parse_optional_f64(row, "beta_conditionB"),
            &format!("native GLM-mu mean Wald beta conditionB gene {gene}"),
        );
        assert_f64_or_missing(
            beta_se.row(gene).unwrap()[0],
            parse_optional_f64(row, "beta_se_intercept"),
            &format!("native GLM-mu mean Wald beta SE intercept gene {gene}"),
        );
        assert_f64_or_missing(
            beta_se.row(gene).unwrap()[1],
            parse_optional_f64(row, "beta_se_conditionB"),
            &format!("native GLM-mu mean Wald beta SE conditionB gene {gene}"),
        );
        assert_option_close(
            wald.stat[gene],
            parse_optional_f64(row, "stat_conditionB"),
            1e-5,
            1e-5,
            &format!("native GLM-mu mean Wald stat gene {gene}"),
        );
        assert_option_close(
            wald.pvalue[gene],
            parse_optional_f64(row, "pvalue_conditionB"),
            1e-5,
            1e-5,
            &format!("native GLM-mu mean Wald pvalue gene {gene}"),
        );
        assert_option_close(
            results.rows[gene].pvalue,
            parse_optional_f64(row, "pvalue_conditionB"),
            1e-5,
            1e-5,
            &format!("native GLM-mu mean result pvalue gene {gene}"),
        );
        assert_option_close(
            results.rows[gene].padj,
            parse_optional_f64(row, "padj_conditionB"),
            1e-5,
            1e-5,
            &format!("native GLM-mu mean result padj gene {gene}"),
        );
        assert_f64_or_missing(
            log_like[gene],
            parse_optional_f64(row, "log_like"),
            &format!("native GLM-mu mean log-likelihood gene {gene}"),
        );
        if !fit.all_zero[gene] {
            assert_eq!(
                beta_converged[gene],
                parse_required_bool(row, "beta_converged"),
                "native GLM-mu mean beta convergence gene {gene}"
            );
            assert_eq!(
                beta_iter[gene],
                parse_required_f64(row, "beta_iterations") as usize,
                "native GLM-mu mean beta iterations gene {gene}"
            );
        }
    }

    for (gene, (mu_row, hat_row)) in wald_mu_rows.iter().zip(wald_hat_rows.iter()).enumerate() {
        for (sample, sample_name) in samples.iter().enumerate() {
            assert_f64_or_missing(
                mu.row(gene).unwrap()[sample],
                parse_optional_f64(mu_row, sample_name),
                &format!("native GLM-mu mean Wald mu gene {gene} sample {sample}"),
            );
            assert_f64_or_missing(
                hat.row(gene).unwrap()[sample],
                parse_optional_f64(hat_row, sample_name),
                &format!("native GLM-mu mean Wald hat gene {gene} sample {sample}"),
            );
        }
    }
}

#[test]
fn native_glm_mu_local_wald_matches_optional_deseq2_reference() {
    let Some(rows) = read_optional_tsv("native_glm_mu_local_reference.tsv") else {
        return;
    };
    let Some(wald_mu_rows) = read_optional_tsv("native_glm_mu_local_wald_mu.tsv") else {
        return;
    };
    let Some(wald_hat_rows) = read_optional_tsv("native_glm_mu_local_wald_hat.tsv") else {
        return;
    };
    let Some(size_factors) = read_size_factors("size_factors_ratio.tsv") else {
        return;
    };

    let counts = reference_counts();
    let design = reference_full_design();
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

    assert_eq!(rows.len(), results.rows.len());
    for (gene, row) in rows.iter().enumerate() {
        assert_eq!(fit.all_zero[gene], parse_required_bool(row, "allZero"));
        assert_float_close(
            fit.base_mean[gene],
            parse_required_f64(row, "baseMean"),
            1e-10,
            1e-10,
            &format!("native GLM-mu local Wald baseMean gene {gene}"),
        );
        assert_f64_or_missing(
            dispersion[gene],
            parse_optional_f64(row, "dispersion"),
            &format!("native GLM-mu local Wald dispersion gene {gene}"),
        );
        assert_f64_or_missing(
            beta.row(gene).unwrap()[0],
            parse_optional_f64(row, "beta_intercept"),
            &format!("native GLM-mu local Wald beta intercept gene {gene}"),
        );
        assert_f64_or_missing(
            beta.row(gene).unwrap()[1],
            parse_optional_f64(row, "beta_conditionB"),
            &format!("native GLM-mu local Wald beta conditionB gene {gene}"),
        );
        assert_f64_or_missing(
            beta_se.row(gene).unwrap()[0],
            parse_optional_f64(row, "beta_se_intercept"),
            &format!("native GLM-mu local Wald beta SE intercept gene {gene}"),
        );
        assert_f64_or_missing(
            beta_se.row(gene).unwrap()[1],
            parse_optional_f64(row, "beta_se_conditionB"),
            &format!("native GLM-mu local Wald beta SE conditionB gene {gene}"),
        );
        assert_option_close(
            wald.stat[gene],
            parse_optional_f64(row, "stat_conditionB"),
            1e-5,
            1e-5,
            &format!("native GLM-mu local Wald stat gene {gene}"),
        );
        assert_option_close(
            wald.pvalue[gene],
            parse_optional_f64(row, "pvalue_conditionB"),
            1e-5,
            1e-5,
            &format!("native GLM-mu local Wald pvalue gene {gene}"),
        );
        assert_option_close(
            results.rows[gene].pvalue,
            parse_optional_f64(row, "pvalue_conditionB"),
            1e-5,
            1e-5,
            &format!("native GLM-mu local result pvalue gene {gene}"),
        );
        assert_option_close(
            results.rows[gene].padj,
            parse_optional_f64(row, "padj_conditionB"),
            1e-5,
            1e-5,
            &format!("native GLM-mu local result padj gene {gene}"),
        );
        assert_f64_or_missing(
            log_like[gene],
            parse_optional_f64(row, "log_like"),
            &format!("native GLM-mu local log-likelihood gene {gene}"),
        );
        if !fit.all_zero[gene] {
            assert_eq!(
                beta_converged[gene],
                parse_required_bool(row, "beta_converged"),
                "native GLM-mu local beta convergence gene {gene}"
            );
            assert_eq!(
                beta_iter[gene],
                parse_required_f64(row, "beta_iterations") as usize,
                "native GLM-mu local beta iterations gene {gene}"
            );
        }
    }

    for (gene, (mu_row, hat_row)) in wald_mu_rows.iter().zip(wald_hat_rows.iter()).enumerate() {
        for (sample, sample_name) in samples.iter().enumerate() {
            assert_f64_or_missing(
                mu.row(gene).unwrap()[sample],
                parse_optional_f64(mu_row, sample_name),
                &format!("native GLM-mu local Wald mu gene {gene} sample {sample}"),
            );
            assert_f64_or_missing(
                hat.row(gene).unwrap()[sample],
                parse_optional_f64(hat_row, sample_name),
                &format!("native GLM-mu local Wald hat gene {gene} sample {sample}"),
            );
        }
    }
}

#[test]
fn native_glm_mu_mean_cox_reid_wald_matches_optional_deseq2_reference() {
    let Some(rows) = read_optional_tsv("native_glm_mu_mean_cr_map_reference.tsv") else {
        return;
    };
    let Some(wald_mu_rows) = read_optional_tsv("native_glm_mu_mean_cr_wald_mu.tsv") else {
        return;
    };
    let Some(wald_hat_rows) = read_optional_tsv("native_glm_mu_mean_cr_wald_hat.tsv") else {
        return;
    };
    let Some(size_factors) = read_size_factors("size_factors_ratio.tsv") else {
        return;
    };

    let counts = reference_counts();
    let design = reference_full_design();
    let (fit, results) = DeseqBuilder::new()
        .size_factors(size_factors)
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

    assert_eq!(rows.len(), results.rows.len());
    for (gene, row) in rows.iter().enumerate() {
        assert_eq!(fit.all_zero[gene], parse_required_bool(row, "allZero"));
        assert_float_close(
            fit.base_mean[gene],
            parse_required_f64(row, "baseMean"),
            1e-10,
            1e-10,
            &format!("native GLM-mu Cox-Reid mean Wald baseMean gene {gene}"),
        );
        assert_f64_or_missing(
            dispersion[gene],
            parse_optional_f64(row, "dispersion"),
            &format!("native GLM-mu Cox-Reid mean Wald dispersion gene {gene}"),
        );
        assert_f64_or_missing(
            beta.row(gene).unwrap()[0],
            parse_optional_f64(row, "beta_intercept"),
            &format!("native GLM-mu Cox-Reid mean Wald beta intercept gene {gene}"),
        );
        assert_f64_or_missing(
            beta.row(gene).unwrap()[1],
            parse_optional_f64(row, "beta_conditionB"),
            &format!("native GLM-mu Cox-Reid mean Wald beta conditionB gene {gene}"),
        );
        assert_f64_or_missing(
            beta_se.row(gene).unwrap()[0],
            parse_optional_f64(row, "beta_se_intercept"),
            &format!("native GLM-mu Cox-Reid mean Wald beta SE intercept gene {gene}"),
        );
        assert_f64_or_missing(
            beta_se.row(gene).unwrap()[1],
            parse_optional_f64(row, "beta_se_conditionB"),
            &format!("native GLM-mu Cox-Reid mean Wald beta SE conditionB gene {gene}"),
        );
        assert_option_close(
            wald.stat[gene],
            parse_optional_f64(row, "stat_conditionB"),
            1e-5,
            1e-5,
            &format!("native GLM-mu Cox-Reid mean Wald stat gene {gene}"),
        );
        assert_option_close(
            wald.pvalue[gene],
            parse_optional_f64(row, "pvalue_conditionB"),
            1e-5,
            1e-5,
            &format!("native GLM-mu Cox-Reid mean Wald pvalue gene {gene}"),
        );
        assert_option_close(
            results.rows[gene].pvalue,
            parse_optional_f64(row, "pvalue_conditionB"),
            1e-5,
            1e-5,
            &format!("native GLM-mu Cox-Reid mean result pvalue gene {gene}"),
        );
        assert_option_close(
            results.rows[gene].padj,
            parse_optional_f64(row, "padj_conditionB"),
            1e-5,
            1e-5,
            &format!("native GLM-mu Cox-Reid mean result padj gene {gene}"),
        );
        assert_f64_or_missing(
            log_like[gene],
            parse_optional_f64(row, "log_like"),
            &format!("native GLM-mu Cox-Reid mean log-likelihood gene {gene}"),
        );
        if !fit.all_zero[gene] {
            assert_eq!(
                beta_converged[gene],
                parse_required_bool(row, "beta_converged"),
                "native GLM-mu Cox-Reid mean beta convergence gene {gene}"
            );
            assert_eq!(
                beta_iter[gene],
                parse_required_f64(row, "beta_iterations") as usize,
                "native GLM-mu Cox-Reid mean beta iterations gene {gene}"
            );
        }
    }

    for (gene, (mu_row, hat_row)) in wald_mu_rows.iter().zip(wald_hat_rows.iter()).enumerate() {
        for (sample, sample_name) in samples.iter().enumerate() {
            assert_f64_or_missing(
                mu.row(gene).unwrap()[sample],
                parse_optional_f64(mu_row, sample_name),
                &format!("native GLM-mu Cox-Reid mean Wald mu gene {gene} sample {sample}"),
            );
            assert_f64_or_missing(
                hat.row(gene).unwrap()[sample],
                parse_optional_f64(hat_row, sample_name),
                &format!("native GLM-mu Cox-Reid mean Wald hat gene {gene} sample {sample}"),
            );
        }
    }
}

#[test]
fn native_glm_mu_local_cox_reid_wald_matches_optional_deseq2_reference() {
    let Some(rows) = read_optional_tsv("native_glm_mu_local_cr_map_reference.tsv") else {
        return;
    };
    let Some(wald_mu_rows) = read_optional_tsv("native_glm_mu_local_cr_wald_mu.tsv") else {
        return;
    };
    let Some(wald_hat_rows) = read_optional_tsv("native_glm_mu_local_cr_wald_hat.tsv") else {
        return;
    };
    let Some(size_factors) = read_size_factors("size_factors_ratio.tsv") else {
        return;
    };

    let counts = reference_counts();
    let design = reference_full_design();
    let (fit, results) = DeseqBuilder::new()
        .size_factors(size_factors)
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

    assert_eq!(rows.len(), results.rows.len());
    for (gene, row) in rows.iter().enumerate() {
        assert_eq!(fit.all_zero[gene], parse_required_bool(row, "allZero"));
        assert_float_close(
            fit.base_mean[gene],
            parse_required_f64(row, "baseMean"),
            1e-10,
            1e-10,
            &format!("native GLM-mu Cox-Reid local Wald baseMean gene {gene}"),
        );
        assert_f64_or_missing(
            dispersion[gene],
            parse_optional_f64(row, "dispersion"),
            &format!("native GLM-mu Cox-Reid local Wald dispersion gene {gene}"),
        );
        assert_f64_or_missing(
            beta.row(gene).unwrap()[0],
            parse_optional_f64(row, "beta_intercept"),
            &format!("native GLM-mu Cox-Reid local Wald beta intercept gene {gene}"),
        );
        assert_f64_or_missing(
            beta.row(gene).unwrap()[1],
            parse_optional_f64(row, "beta_conditionB"),
            &format!("native GLM-mu Cox-Reid local Wald beta conditionB gene {gene}"),
        );
        assert_f64_or_missing(
            beta_se.row(gene).unwrap()[0],
            parse_optional_f64(row, "beta_se_intercept"),
            &format!("native GLM-mu Cox-Reid local Wald beta SE intercept gene {gene}"),
        );
        assert_f64_or_missing(
            beta_se.row(gene).unwrap()[1],
            parse_optional_f64(row, "beta_se_conditionB"),
            &format!("native GLM-mu Cox-Reid local Wald beta SE conditionB gene {gene}"),
        );
        assert_option_close(
            wald.stat[gene],
            parse_optional_f64(row, "stat_conditionB"),
            1e-5,
            1e-5,
            &format!("native GLM-mu Cox-Reid local Wald stat gene {gene}"),
        );
        assert_option_close(
            wald.pvalue[gene],
            parse_optional_f64(row, "pvalue_conditionB"),
            1e-5,
            1e-5,
            &format!("native GLM-mu Cox-Reid local Wald pvalue gene {gene}"),
        );
        assert_option_close(
            results.rows[gene].pvalue,
            parse_optional_f64(row, "pvalue_conditionB"),
            1e-5,
            1e-5,
            &format!("native GLM-mu Cox-Reid local result pvalue gene {gene}"),
        );
        assert_option_close(
            results.rows[gene].padj,
            parse_optional_f64(row, "padj_conditionB"),
            1e-5,
            1e-5,
            &format!("native GLM-mu Cox-Reid local result padj gene {gene}"),
        );
        assert_f64_or_missing(
            log_like[gene],
            parse_optional_f64(row, "log_like"),
            &format!("native GLM-mu Cox-Reid local log-likelihood gene {gene}"),
        );
        if !fit.all_zero[gene] {
            assert_eq!(
                beta_converged[gene],
                parse_required_bool(row, "beta_converged"),
                "native GLM-mu Cox-Reid local beta convergence gene {gene}"
            );
            assert_eq!(
                beta_iter[gene],
                parse_required_f64(row, "beta_iterations") as usize,
                "native GLM-mu Cox-Reid local beta iterations gene {gene}"
            );
        }
    }

    for (gene, (mu_row, hat_row)) in wald_mu_rows.iter().zip(wald_hat_rows.iter()).enumerate() {
        for (sample, sample_name) in samples.iter().enumerate() {
            assert_f64_or_missing(
                mu.row(gene).unwrap()[sample],
                parse_optional_f64(mu_row, sample_name),
                &format!("native GLM-mu Cox-Reid local Wald mu gene {gene} sample {sample}"),
            );
            assert_f64_or_missing(
                hat.row(gene).unwrap()[sample],
                parse_optional_f64(hat_row, sample_name),
                &format!("native GLM-mu Cox-Reid local Wald hat gene {gene} sample {sample}"),
            );
        }
    }
}

#[test]
fn fixed_dispersion_wald_matches_optional_deseq2_internal_reference() {
    let Some(rows) = read_optional_tsv("fixed_wald_reference.tsv") else {
        return;
    };
    let Some(size_factors) = read_size_factors("size_factors_ratio.tsv") else {
        return;
    };

    let dispersions = rows
        .iter()
        .map(|row| parse_required_f64(row, "dispersion"))
        .collect::<Vec<_>>();
    let (fit, results) = DeseqBuilder::new()
        .size_factors(size_factors)
        .disable_cooks_cutoff()
        .disable_independent_filtering()
        .fit_fixed_dispersion_wald(
            &reference_counts(),
            &reference_full_design(),
            &dispersions,
            1,
        )
        .unwrap();

    let beta = fit.beta.as_ref().unwrap();
    let beta_se = fit.beta_se.as_ref().unwrap();
    let log_like = fit.log_like.as_ref().unwrap();
    let wald = fit.wald.as_ref().unwrap();
    let genes = reference_gene_names();

    assert_eq!(rows.len(), results.rows.len());
    for (gene, row) in rows.iter().enumerate() {
        assert_eq!(
            row.get("gene").map(String::as_str),
            Some(genes[gene].as_str())
        );
        assert_float_close(
            beta.row(gene).unwrap()[0],
            parse_required_f64(row, "beta_intercept"),
            1e-5,
            1e-5,
            &format!("wald beta intercept gene {gene}"),
        );
        assert_float_close(
            beta.row(gene).unwrap()[1],
            parse_required_f64(row, "beta_conditionB"),
            1e-5,
            1e-5,
            &format!("wald beta conditionB gene {gene}"),
        );
        assert_float_close(
            beta_se.row(gene).unwrap()[0],
            parse_required_f64(row, "beta_se_intercept"),
            1e-5,
            1e-5,
            &format!("wald beta SE intercept gene {gene}"),
        );
        assert_float_close(
            beta_se.row(gene).unwrap()[1],
            parse_required_f64(row, "beta_se_conditionB"),
            1e-5,
            1e-5,
            &format!("wald beta SE conditionB gene {gene}"),
        );
        assert_option_close(
            wald.stat[gene],
            Some(parse_required_f64(row, "stat_conditionB")),
            1e-5,
            1e-5,
            &format!("wald stat gene {gene}"),
        );
        assert_option_close(
            wald.pvalue[gene],
            Some(parse_required_f64(row, "pvalue_conditionB")),
            1e-7,
            1e-5,
            &format!("wald p-value gene {gene}"),
        );
        if row.contains_key("log_like") {
            assert_float_close(
                log_like[gene],
                parse_required_f64(row, "log_like"),
                1e-5,
                1e-5,
                &format!("wald log-likelihood gene {gene}"),
            );
        }
        assert_option_close(
            results.rows[gene].log2_fold_change,
            Some(parse_required_f64(row, "beta_conditionB")),
            1e-5,
            1e-5,
            &format!("wald result LFC gene {gene}"),
        );
        assert_option_close(
            results.rows[gene].lfc_se,
            Some(parse_required_f64(row, "beta_se_conditionB")),
            1e-5,
            1e-5,
            &format!("wald result LFC SE gene {gene}"),
        );
    }
}

#[test]
fn fixed_dispersion_wald_t_matches_optional_deseq2_internal_reference() {
    let Some(rows) = read_optional_tsv("fixed_wald_t_reference.tsv") else {
        return;
    };
    let Some(size_factors) = read_size_factors("size_factors_ratio.tsv") else {
        return;
    };

    let dispersions = rows
        .iter()
        .map(|row| parse_required_f64(row, "dispersion"))
        .collect::<Vec<_>>();
    let (fit, results) = DeseqBuilder::new()
        .size_factors(size_factors)
        .disable_cooks_cutoff()
        .disable_independent_filtering()
        .wald_t_residual_degrees_of_freedom()
        .fit_fixed_dispersion_wald(
            &reference_counts(),
            &reference_full_design(),
            &dispersions,
            1,
        )
        .unwrap();

    let wald = fit.wald.as_ref().unwrap();
    let degrees_of_freedom = wald.degrees_of_freedom.as_ref().unwrap();

    assert_eq!(rows.len(), results.rows.len());
    for (gene, row) in rows.iter().enumerate() {
        assert_option_close(
            wald.stat[gene],
            Some(parse_required_f64(row, "stat_conditionB")),
            1e-5,
            1e-5,
            &format!("t-Wald stat gene {gene}"),
        );
        assert_option_close(
            degrees_of_freedom[gene],
            Some(parse_required_f64(row, "df_conditionB")),
            1e-12,
            1e-12,
            &format!("t-Wald df gene {gene}"),
        );
        assert_option_close(
            wald.pvalue[gene],
            Some(parse_required_f64(row, "pvalue_conditionB")),
            1e-7,
            1e-5,
            &format!("t-Wald p-value gene {gene}"),
        );
        assert_eq!(results.rows[gene].pvalue, wald.pvalue[gene]);
    }
}

#[test]
fn fixed_dispersion_force_optim_wald_matches_optional_deseq2_internal_reference() {
    let Some(rows) = read_optional_tsv("fixed_force_optim_wald_reference.tsv") else {
        return;
    };
    let Some(expected_mu) = read_optional_tsv("fixed_force_optim_mu_full.tsv") else {
        return;
    };
    let Some(expected_hat) = read_optional_tsv("fixed_force_optim_hat_full.tsv") else {
        return;
    };
    let Some(size_factors) = read_size_factors("size_factors_ratio.tsv") else {
        return;
    };

    let dispersions = rows
        .iter()
        .map(|row| parse_required_f64(row, "dispersion"))
        .collect::<Vec<_>>();
    let (fit, results) = DeseqBuilder::new()
        .size_factors(size_factors)
        .disable_cooks_cutoff()
        .disable_independent_filtering()
        .irls_options(IrlsOptions {
            use_optim: true,
            force_optim: true,
            ..IrlsOptions::default()
        })
        .fit_fixed_dispersion_wald(
            &reference_counts(),
            &reference_full_design(),
            &dispersions,
            1,
        )
        .unwrap();

    let beta = fit.beta.as_ref().unwrap();
    let beta_se = fit.beta_se.as_ref().unwrap();
    let beta_converged = fit.beta_converged.as_ref().unwrap();
    let beta_iter = fit.beta_iter.as_ref().unwrap();
    let log_like = fit.log_like.as_ref().unwrap();
    let mu = fit.mu.as_ref().unwrap();
    let hat = fit.hat_diagonal.as_ref().unwrap();
    let wald = fit.wald.as_ref().unwrap();
    let samples = reference_sample_names();

    assert_eq!(rows.len(), results.rows.len());
    assert_eq!(expected_mu.len(), results.rows.len());
    assert_eq!(expected_hat.len(), results.rows.len());
    for (gene, row) in rows.iter().enumerate() {
        assert_float_close(
            beta.row(gene).unwrap()[0],
            parse_required_f64(row, "beta_intercept"),
            1e-1,
            1e-3,
            &format!("force-optim Wald beta intercept gene {gene}"),
        );
        assert_float_close(
            beta.row(gene).unwrap()[1],
            parse_required_f64(row, "beta_conditionB"),
            1e-1,
            1e-3,
            &format!("force-optim Wald beta conditionB gene {gene}"),
        );
        assert_float_close(
            beta_se.row(gene).unwrap()[0],
            parse_required_f64(row, "beta_se_intercept"),
            2e-4,
            2e-4,
            &format!("force-optim Wald beta SE intercept gene {gene}"),
        );
        assert_float_close(
            beta_se.row(gene).unwrap()[1],
            parse_required_f64(row, "beta_se_conditionB"),
            2e-4,
            2e-4,
            &format!("force-optim Wald beta SE conditionB gene {gene}"),
        );
        assert_option_close(
            wald.stat[gene],
            Some(parse_required_f64(row, "stat_conditionB")),
            1e-1,
            1e-3,
            &format!("force-optim Wald stat gene {gene}"),
        );
        assert_option_close(
            wald.pvalue[gene],
            Some(parse_required_f64(row, "pvalue_conditionB")),
            1e-6,
            1e-5,
            &format!("force-optim Wald p-value gene {gene}"),
        );
        assert_float_close(
            log_like[gene],
            parse_required_f64(row, "log_like"),
            2e-4,
            2e-4,
            &format!("force-optim Wald log-likelihood gene {gene}"),
        );
        assert_eq!(
            beta_converged[gene],
            parse_required_bool(row, "converged"),
            "force-optim convergence gene {gene}"
        );
        assert_eq!(
            beta_iter[gene],
            parse_required_usize(row, "iterations"),
            "force-optim IRLS iterations gene {gene}"
        );

        let mu_row = &expected_mu[gene];
        let hat_row = &expected_hat[gene];
        for (sample, sample_name) in samples.iter().enumerate() {
            assert_float_close(
                mu.row(gene).unwrap()[sample],
                parse_required_f64(mu_row, sample_name),
                2e-4,
                2e-4,
                &format!("force-optim fitted mean gene {gene} sample {sample}"),
            );
            assert_float_close(
                hat.row(gene).unwrap()[sample],
                parse_required_f64(hat_row, sample_name),
                2e-4,
                2e-4,
                &format!("force-optim hat diagonal gene {gene} sample {sample}"),
            );
        }
    }
}

#[test]
fn fixed_dispersion_cooks_match_optional_deseq2_internal_reference() {
    let Some(expected_cooks) = read_optional_tsv("fixed_cooks_full.tsv") else {
        return;
    };
    let Some(size_factors) = read_size_factors("size_factors_ratio.tsv") else {
        return;
    };
    let Some(dispersions) = read_fixed_dispersions() else {
        return;
    };

    let (fit, _results) = DeseqBuilder::new()
        .size_factors(size_factors)
        .disable_cooks_cutoff()
        .disable_independent_filtering()
        .fit_fixed_dispersion_wald(
            &reference_counts(),
            &reference_full_design(),
            &dispersions,
            1,
        )
        .unwrap();
    let cooks = fit.cooks.as_ref().unwrap();
    let samples = reference_sample_names();

    assert_eq!(cooks.n_rows(), expected_cooks.len());
    for (gene, row) in expected_cooks.iter().enumerate() {
        for (sample, sample_name) in samples.iter().enumerate() {
            assert_float_close(
                cooks.row(gene).unwrap()[sample],
                parse_required_f64(row, sample_name),
                1e-4,
                1e-4,
                &format!("Cook's distance gene {gene} sample {sample}"),
            );
        }
    }
}

#[test]
fn fixed_dispersion_weighted_wald_matches_optional_deseq2_internal_reference() {
    let Some(rows) = read_optional_tsv("fixed_wald_weighted_reference.tsv") else {
        return;
    };
    let Some(size_factors) = read_size_factors("size_factors_ratio.tsv") else {
        return;
    };
    let Some(observation_weights) = read_reference_matrix("observation_weights.tsv") else {
        return;
    };

    let dispersions = rows
        .iter()
        .map(|row| parse_required_f64(row, "dispersion"))
        .collect::<Vec<_>>();
    let (fit, results) = DeseqBuilder::new()
        .size_factors(size_factors)
        .observation_weights(observation_weights)
        .disable_cooks_cutoff()
        .disable_independent_filtering()
        .fit_fixed_dispersion_wald(
            &reference_counts(),
            &reference_full_design(),
            &dispersions,
            1,
        )
        .unwrap();

    let beta = fit.beta.as_ref().unwrap();
    let beta_se = fit.beta_se.as_ref().unwrap();
    let log_like = fit.log_like.as_ref().unwrap();
    let wald = fit.wald.as_ref().unwrap();

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
            &format!("weighted Wald baseMean gene {gene}"),
        );
        assert_f64_or_missing(
            beta.row(gene).unwrap()[0],
            parse_optional_f64(row, "beta_intercept"),
            &format!("weighted Wald beta intercept gene {gene}"),
        );
        assert_f64_or_missing(
            beta.row(gene).unwrap()[1],
            parse_optional_f64(row, "beta_conditionB"),
            &format!("weighted Wald beta conditionB gene {gene}"),
        );
        assert_f64_or_missing(
            beta_se.row(gene).unwrap()[1],
            parse_optional_f64(row, "beta_se_conditionB"),
            &format!("weighted Wald beta SE conditionB gene {gene}"),
        );
        assert_option_close(
            wald.stat[gene],
            parse_optional_f64(row, "stat_conditionB"),
            1e-5,
            1e-5,
            &format!("weighted Wald stat gene {gene}"),
        );
        assert_option_close(
            wald.pvalue[gene],
            parse_optional_f64(row, "pvalue_conditionB"),
            1e-7,
            1e-5,
            &format!("weighted Wald p-value gene {gene}"),
        );
        if row.contains_key("log_like") {
            assert_f64_or_missing(
                log_like[gene],
                parse_optional_f64(row, "log_like"),
                &format!("weighted Wald log-likelihood gene {gene}"),
            );
        }
        assert_eq!(results.rows[gene].pvalue, wald.pvalue[gene]);
    }
}

#[test]
fn native_weighted_glm_mu_wald_matches_optional_deseq2_reference() {
    let Some(rows) = read_optional_tsv("native_weighted_glm_mu_reference.tsv") else {
        return;
    };
    let Some(wald_mu_rows) = read_optional_tsv("native_weighted_glm_mu_wald_mu.tsv") else {
        return;
    };
    let Some(wald_hat_rows) = read_optional_tsv("native_weighted_glm_mu_wald_hat.tsv") else {
        return;
    };
    let Some(size_factors) = read_size_factors("size_factors_ratio.tsv") else {
        return;
    };
    let Some(observation_weights) = read_reference_matrix("observation_weights.tsv") else {
        return;
    };

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
fn native_weighted_glm_mu_local_wald_matches_optional_deseq2_reference() {
    let Some(rows) = read_optional_tsv("native_weighted_glm_mu_local_reference.tsv") else {
        return;
    };
    let Some(wald_mu_rows) = read_optional_tsv("native_weighted_glm_mu_local_wald_mu.tsv") else {
        return;
    };
    let Some(wald_hat_rows) = read_optional_tsv("native_weighted_glm_mu_local_wald_hat.tsv") else {
        return;
    };
    let Some(size_factors) = read_size_factors("size_factors_ratio.tsv") else {
        return;
    };
    let Some(observation_weights) = read_reference_matrix("observation_weights.tsv") else {
        return;
    };

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
fn native_weighted_glm_mu_mean_cox_reid_wald_matches_optional_deseq2_reference() {
    let Some(rows) = read_optional_tsv("native_weighted_glm_mu_mean_cr_map_reference.tsv") else {
        return;
    };
    let Some(wald_mu_rows) = read_optional_tsv("native_weighted_glm_mu_mean_cr_wald_mu.tsv") else {
        return;
    };
    let Some(wald_hat_rows) = read_optional_tsv("native_weighted_glm_mu_mean_cr_wald_hat.tsv")
    else {
        return;
    };
    let Some(size_factors) = read_size_factors("size_factors_ratio.tsv") else {
        return;
    };
    let Some(observation_weights) = read_reference_matrix("observation_weights.tsv") else {
        return;
    };

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
fn native_weighted_glm_mu_local_cox_reid_wald_matches_optional_deseq2_reference() {
    let Some(rows) = read_optional_tsv("native_weighted_glm_mu_local_cr_map_reference.tsv") else {
        return;
    };
    let Some(wald_mu_rows) = read_optional_tsv("native_weighted_glm_mu_local_cr_wald_mu.tsv")
    else {
        return;
    };
    let Some(wald_hat_rows) = read_optional_tsv("native_weighted_glm_mu_local_cr_wald_hat.tsv")
    else {
        return;
    };
    let Some(size_factors) = read_size_factors("size_factors_ratio.tsv") else {
        return;
    };
    let Some(observation_weights) = read_reference_matrix("observation_weights.tsv") else {
        return;
    };

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
