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
fn fixed_dispersion_lrt_matches_optional_deseq2_internal_reference() {
    let Some(rows) = read_optional_tsv("fixed_lrt_reference.tsv") else {
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
        .fit_fixed_dispersion_lrt(
            &reference_counts(),
            &reference_full_design(),
            &reference_reduced_design(),
            &dispersions,
            1,
        )
        .unwrap();

    let beta = fit.beta.as_ref().unwrap();
    let beta_se = fit.beta_se.as_ref().unwrap();
    let log_like = fit.log_like.as_ref().unwrap();
    let full_deviance = fit.full_deviance.as_ref().unwrap();
    let reduced_log_like = fit.reduced_log_like.as_ref().unwrap();
    let reduced_beta_converged = fit.reduced_beta_converged.as_ref().unwrap();
    let reduced_beta_iter = fit.reduced_beta_iter.as_ref().unwrap();
    let lrt = fit.lrt.as_ref().unwrap();
    let genes = reference_gene_names();

    assert_eq!(lrt.degrees_of_freedom, 1);
    assert_eq!(rows.len(), results.rows.len());
    assert_eq!(reduced_beta_converged, &lrt.reduced_converged);
    assert_eq!(reduced_beta_iter.len(), rows.len());
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
            &format!("LRT full beta intercept gene {gene}"),
        );
        assert_float_close(
            beta.row(gene).unwrap()[1],
            parse_required_f64(row, "beta_conditionB"),
            1e-5,
            1e-5,
            &format!("LRT full beta conditionB gene {gene}"),
        );
        assert_float_close(
            beta_se.row(gene).unwrap()[1],
            parse_required_f64(row, "beta_se_conditionB"),
            1e-5,
            1e-5,
            &format!("LRT full beta SE conditionB gene {gene}"),
        );
        assert_option_close(
            lrt.deviance[gene],
            Some(parse_required_f64(row, "lrt_stat")),
            1e-5,
            1e-5,
            &format!("LRT statistic gene {gene}"),
        );
        assert_float_close(
            log_like[gene],
            parse_required_f64(row, "log_like_full"),
            1e-5,
            1e-5,
            &format!("LRT full log-likelihood gene {gene}"),
        );
        assert_float_close(
            full_deviance[gene],
            -2.0 * parse_required_f64(row, "log_like_full"),
            1e-5,
            1e-5,
            &format!("LRT full deviance gene {gene}"),
        );
        assert_float_close(
            reduced_log_like[gene],
            parse_required_f64(row, "log_like_reduced"),
            1e-5,
            1e-5,
            &format!("LRT reduced log-likelihood gene {gene}"),
        );
        assert_option_close(
            lrt.pvalue[gene],
            Some(parse_required_f64(row, "pvalue")),
            1e-7,
            1e-5,
            &format!("LRT p-value gene {gene}"),
        );
        assert_option_close(
            results.rows[gene].stat,
            Some(parse_required_f64(row, "lrt_stat")),
            1e-5,
            1e-5,
            &format!("LRT result statistic gene {gene}"),
        );
        assert_eq!(
            fit.beta_converged.as_ref().unwrap()[gene],
            parse_required_bool(row, "full_converged")
        );
        assert_eq!(
            reduced_beta_converged[gene],
            parse_required_bool(row, "reduced_converged")
        );
        assert!(reduced_beta_iter[gene] > 0);
    }
}

#[test]
fn fixed_dispersion_weighted_lrt_matches_optional_deseq2_internal_reference() {
    let Some(rows) = read_optional_tsv("fixed_lrt_weighted_reference.tsv") else {
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
        .fit_fixed_dispersion_lrt(
            &reference_counts(),
            &reference_full_design(),
            &reference_reduced_design(),
            &dispersions,
            1,
        )
        .unwrap();

    let beta = fit.beta.as_ref().unwrap();
    let beta_se = fit.beta_se.as_ref().unwrap();
    let log_like = fit.log_like.as_ref().unwrap();
    let full_deviance = fit.full_deviance.as_ref().unwrap();
    let reduced_log_like = fit.reduced_log_like.as_ref().unwrap();
    let reduced_beta_converged = fit.reduced_beta_converged.as_ref().unwrap();
    let reduced_beta_iter = fit.reduced_beta_iter.as_ref().unwrap();
    let lrt = fit.lrt.as_ref().unwrap();

    assert_eq!(lrt.degrees_of_freedom, 1);
    assert_eq!(rows.len(), results.rows.len());
    assert_eq!(reduced_beta_converged, &lrt.reduced_converged);
    assert_eq!(reduced_beta_iter.len(), rows.len());
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
            &format!("weighted LRT baseMean gene {gene}"),
        );
        assert_f64_or_missing(
            beta.row(gene).unwrap()[0],
            parse_optional_f64(row, "beta_intercept"),
            &format!("weighted LRT full beta intercept gene {gene}"),
        );
        assert_f64_or_missing(
            beta.row(gene).unwrap()[1],
            parse_optional_f64(row, "beta_conditionB"),
            &format!("weighted LRT full beta conditionB gene {gene}"),
        );
        assert_f64_or_missing(
            beta_se.row(gene).unwrap()[1],
            parse_optional_f64(row, "beta_se_conditionB"),
            &format!("weighted LRT full beta SE conditionB gene {gene}"),
        );
        assert_option_close(
            lrt.deviance[gene],
            parse_optional_f64(row, "lrt_stat"),
            1e-5,
            1e-5,
            &format!("weighted LRT statistic gene {gene}"),
        );
        assert_f64_or_missing(
            log_like[gene],
            parse_optional_f64(row, "log_like_full"),
            &format!("weighted LRT full log-likelihood gene {gene}"),
        );
        assert_f64_or_missing(
            full_deviance[gene],
            parse_optional_f64(row, "log_like_full").map(|log_like| -2.0 * log_like),
            &format!("weighted LRT full deviance gene {gene}"),
        );
        assert_f64_or_missing(
            reduced_log_like[gene],
            parse_optional_f64(row, "log_like_reduced"),
            &format!("weighted LRT reduced log-likelihood gene {gene}"),
        );
        assert_option_close(
            lrt.pvalue[gene],
            parse_optional_f64(row, "pvalue"),
            1e-7,
            1e-5,
            &format!("weighted LRT p-value gene {gene}"),
        );
        assert_eq!(
            reduced_beta_converged[gene],
            parse_required_bool(row, "reduced_converged")
        );
        if fit.all_zero[gene] || fit.weights_fail.as_ref().unwrap()[gene] {
            assert_eq!(reduced_beta_iter[gene], 0);
        } else {
            assert!(reduced_beta_iter[gene] > 0);
        }
        assert_eq!(results.rows[gene].pvalue, lrt.pvalue[gene]);
    }
}

#[test]
fn native_weighted_glm_mu_lrt_matches_optional_deseq2_reference() {
    let Some(rows) = read_optional_tsv("native_weighted_glm_mu_lrt_reference.tsv") else {
        return;
    };
    let Some(size_factors) = read_size_factors("size_factors_ratio.tsv") else {
        return;
    };
    let Some(observation_weights) = read_reference_matrix("observation_weights.tsv") else {
        return;
    };

    let counts = reference_counts();
    let full_design = reference_full_design();
    let reduced_design = reference_reduced_design();
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

    assert_eq!(lrt.degrees_of_freedom, 1);
    assert_eq!(rows.len(), results.rows.len());
    assert_eq!(rows.len(), beta.n_rows());
    assert_eq!(reduced_beta_converged, &lrt.reduced_converged);
    assert_eq!(reduced_beta_iter.len(), rows.len());
    for (gene, row) in rows.iter().enumerate() {
        let all_zero = parse_required_bool(row, "allZero");
        let weights_fail = parse_required_bool(row, "weightsFail");
        let skipped = all_zero || weights_fail;

        assert_eq!(
            fit.all_zero[gene], skipped,
            "native weighted GLM-mu LRT effective all-zero/weights-fail mask gene {gene}"
        );
        assert_eq!(
            fit.weights_fail.as_ref().unwrap()[gene],
            weights_fail,
            "native weighted GLM-mu LRT weightsFail gene {gene}"
        );
        assert_float_close(
            fit.base_mean[gene],
            parse_required_f64(row, "baseMean"),
            1e-10,
            1e-10,
            &format!("native weighted GLM-mu LRT baseMean gene {gene}"),
        );
        assert_f64_or_missing(
            dispersion[gene],
            parse_optional_f64(row, "dispersion"),
            &format!("native weighted GLM-mu LRT dispersion gene {gene}"),
        );
        assert_f64_or_missing(
            beta.row(gene).unwrap()[0],
            parse_optional_f64(row, "beta_intercept"),
            &format!("native weighted GLM-mu LRT full beta intercept gene {gene}"),
        );
        assert_f64_or_missing(
            beta.row(gene).unwrap()[1],
            parse_optional_f64(row, "beta_conditionB"),
            &format!("native weighted GLM-mu LRT full beta conditionB gene {gene}"),
        );
        assert_f64_or_missing(
            beta_se.row(gene).unwrap()[0],
            parse_optional_f64(row, "beta_se_intercept"),
            &format!("native weighted GLM-mu LRT full beta SE intercept gene {gene}"),
        );
        assert_f64_or_missing(
            beta_se.row(gene).unwrap()[1],
            parse_optional_f64(row, "beta_se_conditionB"),
            &format!("native weighted GLM-mu LRT full beta SE conditionB gene {gene}"),
        );
        assert_option_close(
            lrt.deviance[gene],
            parse_optional_f64(row, "lrt_stat"),
            1e-4,
            1e-4,
            &format!("native weighted GLM-mu LRT statistic gene {gene}"),
        );
        assert_f64_or_missing(
            log_like[gene],
            parse_optional_f64(row, "log_like_full"),
            &format!("native weighted GLM-mu LRT full log-likelihood gene {gene}"),
        );
        assert_f64_or_missing(
            fit.full_deviance.as_ref().unwrap()[gene],
            parse_optional_f64(row, "log_like_full").map(|log_like| -2.0 * log_like),
            &format!("native weighted GLM-mu LRT full deviance gene {gene}"),
        );
        assert_f64_or_missing(
            reduced_log_like[gene],
            parse_optional_f64(row, "log_like_reduced"),
            &format!("native weighted GLM-mu LRT reduced log-likelihood gene {gene}"),
        );
        assert_option_close(
            lrt.pvalue[gene],
            parse_optional_f64(row, "pvalue"),
            1e-4,
            1e-4,
            &format!("native weighted GLM-mu LRT p-value gene {gene}"),
        );
        assert_option_close(
            results.rows[gene].stat,
            parse_optional_f64(row, "lrt_stat"),
            1e-4,
            1e-4,
            &format!("native weighted GLM-mu LRT result statistic gene {gene}"),
        );
        assert_eq!(
            results.rows[gene].pvalue, lrt.pvalue[gene],
            "native weighted GLM-mu LRT result p-value gene {gene}"
        );

        if skipped {
            assert_eq!(beta_iter[gene], 0);
            assert_eq!(reduced_beta_iter[gene], 0);
        } else {
            assert_eq!(
                beta_converged[gene],
                parse_required_bool(row, "full_converged"),
                "native weighted GLM-mu LRT full convergence gene {gene}"
            );
            assert_eq!(
                reduced_beta_converged[gene],
                parse_required_bool(row, "reduced_converged"),
                "native weighted GLM-mu LRT reduced convergence gene {gene}"
            );
            assert_eq!(
                beta_iter[gene],
                parse_required_f64(row, "full_iterations") as usize,
                "native weighted GLM-mu LRT full iterations gene {gene}"
            );
            assert_eq!(
                reduced_beta_iter[gene],
                parse_required_f64(row, "reduced_iterations") as usize,
                "native weighted GLM-mu LRT reduced iterations gene {gene}"
            );
        }
    }
}
