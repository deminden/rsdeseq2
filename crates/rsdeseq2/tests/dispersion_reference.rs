mod common;

use common::*;
use rsdeseq2::prelude::*;

fn assert_f64_or_missing(actual: f64, expected: Option<f64>, atol: f64, rtol: f64, label: &str) {
    match expected {
        Some(expected) => assert_float_close(actual, expected, atol, rtol, label),
        None => assert!(actual.is_nan(), "{label}: expected missing, got {actual}"),
    }
}

#[test]
fn supplied_fixed_dispersions_from_optional_reference_are_preserved() {
    let Some(dispersions) = read_fixed_dispersions() else {
        return;
    };
    let Some(size_factors) = read_size_factors("size_factors_ratio.tsv") else {
        return;
    };

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

    let fit_dispersions = fit.dispersion.as_ref().unwrap();
    assert_eq!(fit_dispersions.len(), dispersions.len());
    for (gene, (actual, expected)) in fit_dispersions
        .iter()
        .copied()
        .zip(dispersions.iter().copied())
        .enumerate()
    {
        assert_float_close(
            actual,
            expected,
            1e-15,
            1e-15,
            &format!("supplied dispersion gene {gene}"),
        );
        assert_eq!(results.rows[gene].dispersion, Some(expected));
    }
}

#[test]
fn native_dispersion_reference_is_explicit_future_work() {
    let Some(rows) = read_optional_tsv("metadata.tsv") else {
        return;
    };
    assert!(
        rows.iter().any(|row| {
            row.get("key").map(String::as_str) == Some("fixed_glm_reference_mode")
                && row
                    .get("value")
                    .is_some_and(|value| value.contains("supplied dispersions"))
        }),
        "metadata should document that current GLM references use supplied dispersions"
    );
}

#[test]
fn native_weighted_glm_mu_map_matches_optional_deseq2_reference() {
    let Some(reference_rows) = read_optional_tsv("native_weighted_glm_mu_reference.tsv") else {
        return;
    };
    let Some(mu_rows) = read_optional_tsv("native_weighted_glm_mu_dispersion_mu.tsv") else {
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
    let fit = DeseqBuilder::new()
        .size_factors(size_factors)
        .observation_weights(observation_weights)
        .fit_type(FitType::Mean)
        .gene_wise_dispersion_options(GeneWiseDispersionOptions {
            use_cox_reid: false,
            niter: 2,
            ..GeneWiseDispersionOptions::default()
        })
        .fit_map_dispersions_glm_mu(&counts, &design)
        .unwrap();

    let disp_gene_est = fit.disp_gene_est.as_ref().unwrap();
    let disp_gene_iter = fit.disp_gene_iter.as_ref().unwrap();
    let disp_fit = fit.disp_fit.as_ref().unwrap();
    let disp_map = fit.disp_map.as_ref().unwrap();
    let dispersion = fit.dispersion.as_ref().unwrap();
    let disp_iter = fit.disp_iter.as_ref().unwrap();
    let disp_outlier = fit.disp_outlier.as_ref().unwrap();
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
            &format!("weighted GLM-mu baseMean gene {gene}"),
        );
        assert_float_close(
            fit.base_var[gene],
            parse_required_f64(row, "baseVar"),
            1e-10,
            1e-10,
            &format!("weighted GLM-mu baseVar gene {gene}"),
        );
        assert_f64_or_missing(
            disp_gene_est[gene],
            parse_optional_f64(row, "dispGeneEst"),
            1e-4,
            1e-4,
            &format!("weighted GLM-mu dispGeneEst gene {gene}"),
        );
        assert_f64_or_missing(
            disp_fit[gene],
            parse_optional_f64(row, "dispFit"),
            1e-4,
            1e-4,
            &format!("weighted GLM-mu dispFit gene {gene}"),
        );
        assert_f64_or_missing(
            disp_map[gene],
            parse_optional_f64(row, "dispMAP"),
            1e-4,
            1e-4,
            &format!("weighted GLM-mu dispMAP gene {gene}"),
        );
        assert_f64_or_missing(
            dispersion[gene],
            parse_optional_f64(row, "dispersion"),
            1e-4,
            1e-4,
            &format!("weighted GLM-mu final dispersion gene {gene}"),
        );
        if !fit.all_zero[gene] {
            // The Armijo path can take a different number of equivalent
            // line-search steps while landing on the same dispersion and mu.
            assert!(
                disp_gene_iter[gene] > 0,
                "weighted GLM-mu gene-wise iterations gene {gene}"
            );
            assert!(
                parse_required_f64(row, "dispGeneIter") > 0.0,
                "DESeq2 weighted GLM-mu gene-wise iterations gene {gene}"
            );
            assert_eq!(
                disp_iter[gene],
                parse_required_f64(row, "dispIter") as usize,
                "weighted GLM-mu MAP iterations gene {gene}"
            );
            assert_eq!(
                disp_outlier[gene],
                parse_required_bool(row, "dispOutlier"),
                "weighted GLM-mu MAP outlier gene {gene}"
            );
        }
    }

    assert_float_close(
        fit.disp_prior_var.unwrap(),
        parse_required_f64(&reference_rows[0], "dispPriorVar"),
        1e-4,
        1e-4,
        "weighted GLM-mu dispersion prior variance",
    );

    for (gene, row) in mu_rows.iter().enumerate() {
        for (sample, sample_name) in samples.iter().enumerate() {
            assert_f64_or_missing(
                mu.row(gene).unwrap()[sample],
                parse_optional_f64(row, sample_name),
                1e-4,
                1e-4,
                &format!("weighted GLM-mu dispersion mu gene {gene} sample {sample}"),
            );
        }
    }
}

#[test]
fn native_weighted_glm_mu_cox_reid_gene_wise_matches_optional_deseq2_reference() {
    let Some(reference_rows) = read_optional_tsv("native_weighted_glm_mu_cr_reference.tsv") else {
        return;
    };
    let Some(mu_rows) = read_optional_tsv("native_weighted_glm_mu_cr_dispersion_mu.tsv") else {
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

#[test]
fn native_normalization_factor_dispersion_intermediates_match_optional_deseq2_reference() {
    let Some(reference_rows) = read_optional_tsv("native_nf_dispersion_reference.tsv") else {
        return;
    };
    let Some(mu_rows) = read_optional_tsv("native_nf_mu.tsv") else {
        return;
    };
    let Some(factor_rows) = read_optional_tsv("normalization_factors.tsv") else {
        return;
    };

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

#[test]
fn parametric_dispersion_trend_matches_optional_deseq2_reference() {
    let Some(rows) = read_optional_tsv("parametric_trend_reference.tsv") else {
        return;
    };
    let means = rows
        .iter()
        .map(|row| parse_required_f64(row, "baseMean"))
        .collect::<Vec<_>>();
    let disps = rows
        .iter()
        .map(|row| parse_required_f64(row, "dispGeneEst"))
        .collect::<Vec<_>>();

    let fit = fit_parametric_dispersion_trend(
        &means,
        &disps,
        ParametricDispersionTrendOptions::default(),
    )
    .unwrap();
    let expected_asympt = parse_required_f64(&rows[0], "asymptDisp");
    let expected_extra = parse_required_f64(&rows[0], "extraPois");

    assert_float_close(
        fit.trend.asympt_disp,
        expected_asympt,
        1e-8,
        1e-8,
        "parametric trend asymptDisp",
    );
    assert_float_close(
        fit.trend.extra_pois,
        expected_extra,
        1e-8,
        1e-8,
        "parametric trend extraPois",
    );
    for (idx, row) in rows.iter().enumerate() {
        assert_eq!(fit.use_for_fit[idx], parse_required_bool(row, "useForFit"));
        assert_float_close(
            fit.disp_fit[idx],
            parse_required_f64(row, "dispFit"),
            1e-8,
            1e-8,
            &format!("parametric trend dispFit row {idx}"),
        );
    }
}

#[test]
fn parametric_dispersion_trend_predicts_optional_deseq2_reference_means() {
    let Some(fit_rows) = read_optional_tsv("parametric_trend_reference.tsv") else {
        return;
    };
    let Some(prediction_rows) = read_optional_tsv("parametric_trend_prediction_reference.tsv")
    else {
        return;
    };
    let means = fit_rows
        .iter()
        .map(|row| parse_required_f64(row, "baseMean"))
        .collect::<Vec<_>>();
    let disps = fit_rows
        .iter()
        .map(|row| parse_required_f64(row, "dispGeneEst"))
        .collect::<Vec<_>>();
    let prediction_means = prediction_rows
        .iter()
        .map(|row| parse_required_f64(row, "mean"))
        .collect::<Vec<_>>();

    let fit = fit_parametric_dispersion_trend(
        &means,
        &disps,
        ParametricDispersionTrendOptions::default(),
    )
    .unwrap();
    let predicted = fit
        .trend
        .evaluate_many_allow_missing(&prediction_means)
        .unwrap();

    for (idx, row) in prediction_rows.iter().enumerate() {
        assert_float_close(
            predicted[idx],
            parse_required_f64(row, "dispFit"),
            1e-8,
            1e-8,
            &format!("parametric trend prediction row {idx}"),
        );
    }
}

#[test]
fn mean_dispersion_trend_matches_optional_deseq2_reference() {
    let Some(rows) = read_optional_tsv("mean_trend_reference.tsv") else {
        return;
    };
    let means = rows
        .iter()
        .map(|row| parse_required_f64(row, "baseMean"))
        .collect::<Vec<_>>();
    let disps = rows
        .iter()
        .map(|row| parse_required_f64(row, "dispGeneEst"))
        .collect::<Vec<_>>();

    let fit =
        fit_mean_dispersion_trend(&means, &disps, MeanDispersionTrendOptions::default()).unwrap();
    let expected_mean = parse_required_f64(&rows[0], "meanDisp");

    assert_float_close(
        fit.trend.mean_disp,
        expected_mean,
        1e-12,
        1e-12,
        "mean trend meanDisp",
    );
    for (idx, row) in rows.iter().enumerate() {
        assert_eq!(fit.use_for_fit[idx], parse_required_bool(row, "useForFit"));
        assert_eq!(
            fit.use_for_mean[idx],
            parse_required_bool(row, "useForMean")
        );
        assert_float_close(
            fit.disp_fit[idx],
            parse_required_f64(row, "dispFit"),
            1e-12,
            1e-12,
            &format!("mean trend dispFit row {idx}"),
        );
    }
}

#[test]
fn local_dispersion_trend_matches_optional_deseq2_reference_shape() {
    let Some(rows) = read_optional_tsv("local_trend_reference.tsv") else {
        return;
    };
    let means = rows
        .iter()
        .map(|row| parse_required_f64(row, "baseMean"))
        .collect::<Vec<_>>();
    let disps = rows
        .iter()
        .map(|row| parse_required_f64(row, "dispGeneEst"))
        .collect::<Vec<_>>();

    let fit =
        fit_local_dispersion_trend(&means, &disps, LocalDispersionTrendOptions::default()).unwrap();

    assert_eq!(
        fit.used_min_disp_floor,
        parse_required_bool(&rows[0], "usedMinDispFloor")
    );
    for (idx, row) in rows.iter().enumerate() {
        assert_eq!(fit.use_for_fit[idx], parse_required_bool(row, "useForFit"));
        assert_float_close(
            fit.disp_fit[idx],
            parse_required_f64(row, "dispFit"),
            1e-3,
            1e-3,
            &format!("local trend dispFit row {idx}"),
        );
    }
}

#[test]
fn local_dispersion_trend_predicts_optional_deseq2_reference_means() {
    let Some(fit_rows) = read_optional_tsv("local_trend_reference.tsv") else {
        return;
    };
    let Some(prediction_rows) = read_optional_tsv("local_trend_prediction_reference.tsv") else {
        return;
    };
    let means = fit_rows
        .iter()
        .map(|row| parse_required_f64(row, "baseMean"))
        .collect::<Vec<_>>();
    let disps = fit_rows
        .iter()
        .map(|row| parse_required_f64(row, "dispGeneEst"))
        .collect::<Vec<_>>();
    let prediction_means = prediction_rows
        .iter()
        .map(|row| parse_required_f64(row, "mean"))
        .collect::<Vec<_>>();

    let fit =
        fit_local_dispersion_trend(&means, &disps, LocalDispersionTrendOptions::default()).unwrap();
    let predicted = fit
        .trend
        .evaluate_many_allow_missing(&prediction_means)
        .unwrap();

    for (idx, row) in prediction_rows.iter().enumerate() {
        assert_float_close(
            predicted[idx],
            parse_required_f64(row, "dispFit"),
            2e-3,
            2e-3,
            &format!("local trend prediction row {idx}"),
        );
    }
}

#[test]
fn local_dispersion_trend_floor_matches_optional_deseq2_reference() {
    let Some(rows) = read_optional_tsv("local_trend_floor_reference.tsv") else {
        return;
    };
    let means = rows
        .iter()
        .map(|row| parse_required_f64(row, "baseMean"))
        .collect::<Vec<_>>();
    let disps = rows
        .iter()
        .map(|row| parse_required_f64(row, "dispGeneEst"))
        .collect::<Vec<_>>();

    let fit =
        fit_local_dispersion_trend(&means, &disps, LocalDispersionTrendOptions::default()).unwrap();

    assert_eq!(
        fit.used_min_disp_floor,
        parse_required_bool(&rows[0], "usedMinDispFloor")
    );
    assert_eq!(fit.genes_used, 0);
    for (idx, row) in rows.iter().enumerate() {
        assert_eq!(fit.use_for_fit[idx], parse_required_bool(row, "useForFit"));
        assert_float_close(
            fit.disp_fit[idx],
            parse_required_f64(row, "dispFit"),
            1e-15,
            1e-15,
            &format!("local trend floor dispFit row {idx}"),
        );
    }
}

#[test]
fn local_dispersion_trend_mixed_threshold_matches_optional_deseq2_reference() {
    let Some(rows) = read_optional_tsv("local_trend_mixed_threshold_reference.tsv") else {
        return;
    };
    let means = rows
        .iter()
        .map(|row| parse_required_f64(row, "baseMean"))
        .collect::<Vec<_>>();
    let disps = rows
        .iter()
        .map(|row| parse_required_f64(row, "dispGeneEst"))
        .collect::<Vec<_>>();

    let fit =
        fit_local_dispersion_trend(&means, &disps, LocalDispersionTrendOptions::default()).unwrap();

    assert_eq!(
        fit.used_min_disp_floor,
        parse_required_bool(&rows[0], "usedMinDispFloor")
    );
    assert!(fit.genes_used > 0);
    assert!(fit.genes_used < rows.len());
    for (idx, row) in rows.iter().enumerate() {
        assert_eq!(fit.use_for_fit[idx], parse_required_bool(row, "useForFit"));
        assert_float_close(
            fit.disp_fit[idx],
            parse_required_f64(row, "dispFit"),
            2e-2,
            2e-2,
            &format!("local trend mixed threshold dispFit row {idx}"),
        );
    }
}

#[test]
fn dispersion_prior_variance_matches_optional_deseq2_reference() {
    let Some(rows) = read_optional_tsv("dispersion_prior_variance_reference.tsv") else {
        return;
    };
    let disp_gene_est = rows
        .iter()
        .map(|row| parse_required_f64(row, "dispGeneEst"))
        .collect::<Vec<_>>();
    let disp_fit = rows
        .iter()
        .map(|row| parse_required_f64(row, "dispFit"))
        .collect::<Vec<_>>();
    let n_samples = parse_required_f64(&rows[0], "nSamples") as usize;
    let n_coefficients = parse_required_f64(&rows[0], "nCoefficients") as usize;

    let output = estimate_dispersion_prior_variance(
        &disp_gene_est,
        &disp_fit,
        1e-8,
        n_samples,
        n_coefficients,
    )
    .unwrap();

    assert_eq!(
        output.residual_degrees_of_freedom,
        parse_required_f64(&rows[0], "residualDf") as usize
    );
    assert_float_close(
        output.var_log_disp_estimates,
        parse_required_f64(&rows[0], "varLogDispEsts"),
        1e-10,
        1e-10,
        "dispersion prior varLogDispEsts",
    );
    assert_float_close(
        output.expected_log_dispersion_variance,
        parse_required_f64(&rows[0], "expectedLogDispVariance"),
        1e-10,
        1e-10,
        "dispersion prior expected sampling variance",
    );
    assert_float_close(
        output.disp_prior_var,
        parse_required_f64(&rows[0], "dispPriorVar"),
        1e-10,
        1e-10,
        "dispersion prior variance",
    );
    for (idx, row) in rows.iter().enumerate() {
        assert_eq!(
            output.above_min_disp[idx],
            parse_required_bool(row, "aboveMinDisp")
        );
    }
}

#[test]
fn low_df_dispersion_prior_variance_matches_optional_deseq2_reference_inputs() {
    let Some(rows) = read_optional_tsv("dispersion_prior_variance_low_df_reference.tsv") else {
        return;
    };
    let disp_gene_est = rows
        .iter()
        .map(|row| parse_required_f64(row, "dispGeneEst"))
        .collect::<Vec<_>>();
    let disp_fit = rows
        .iter()
        .map(|row| parse_required_f64(row, "dispFit"))
        .collect::<Vec<_>>();
    let n_samples = parse_required_f64(&rows[0], "nSamples") as usize;
    let n_coefficients = parse_required_f64(&rows[0], "nCoefficients") as usize;

    let output = estimate_dispersion_prior_variance(
        &disp_gene_est,
        &disp_fit,
        1e-8,
        n_samples,
        n_coefficients,
    )
    .unwrap();

    assert_eq!(output.residual_degrees_of_freedom, 2);
    for (idx, row) in rows.iter().enumerate() {
        assert_eq!(
            output.above_min_disp[idx],
            parse_required_bool(row, "aboveMinDisp")
        );
    }
    assert_float_close(
        output.var_log_disp_estimates,
        parse_required_f64(&rows[0], "varLogDispEsts"),
        1e-12,
        1e-12,
        "low-df dispersion prior varLogDispEsts",
    );
    assert_float_close(
        output.expected_log_dispersion_variance,
        parse_required_f64(&rows[0], "expectedLogDispVariance"),
        1e-12,
        1e-12,
        "low-df dispersion prior expected sampling variance",
    );
    let deseq2_prior_var = parse_required_f64(&rows[0], "dispPriorVar");
    assert!(deseq2_prior_var >= 0.25);
    assert!(deseq2_prior_var <= 8.0);
    assert!(output.disp_prior_var >= 0.25);
    assert!(output.disp_prior_var <= 8.0);
}

#[test]
fn map_dispersion_matches_optional_deseq2_internal_reference() {
    let Some(rows) = read_optional_tsv("map_dispersion_reference.tsv") else {
        return;
    };
    let row = &rows[0];
    let counts = CountMatrix::from_row_major_u32(1, 4, vec![10, 30, 10, 30]).unwrap();
    let mu = RowMajorMatrix::from_row_major(1, 4, vec![20.0, 20.0, 20.0, 20.0]).unwrap();

    let output = estimate_map_dispersions(
        MapDispersionInput {
            counts: &counts,
            design: &reference_full_design(),
            mu: &mu,
            disp_gene_est: &[parse_required_f64(row, "dispGeneEst")],
            disp_fit: &[parse_required_f64(row, "dispFit")],
            all_zero: &[false],
            observation_weights: None,
            disp_prior_var: parse_required_f64(row, "dispPriorVar"),
            var_log_disp_estimates: parse_required_f64(row, "varLogDispEsts"),
        },
        MapDispersionOptions {
            use_cox_reid: parse_required_bool(row, "useCR"),
            ..MapDispersionOptions::default()
        },
    )
    .unwrap();

    assert_float_close(
        output.disp_init[0],
        parse_required_f64(row, "dispInit"),
        1e-12,
        1e-12,
        "MAP dispInit",
    );
    assert_float_close(
        output.disp_map[0],
        parse_required_f64(row, "dispMAP"),
        1e-6,
        1e-6,
        "MAP dispMAP",
    );
    assert_float_close(
        output.dispersion[0],
        parse_required_f64(row, "dispersion"),
        1e-6,
        1e-6,
        "MAP final dispersion",
    );
    assert_eq!(
        output.disp_iter[0],
        parse_required_f64(row, "dispIter") as usize
    );
    assert_eq!(
        output.disp_outlier[0],
        parse_required_bool(row, "dispOutlier")
    );
    assert_eq!(output.converged[0], parse_required_bool(row, "converged"));
}
