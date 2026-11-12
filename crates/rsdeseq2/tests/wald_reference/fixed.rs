use super::*;

#[test]
fn fixed_dispersion_wald_matches_deseq2_internal_reference() {
    let rows = read_reference_tsv("fixed_wald_reference.tsv");
    let size_factors = read_size_factors("size_factors_ratio.tsv");

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
fn fixed_dispersion_wald_t_matches_deseq2_internal_reference() {
    let rows = read_reference_tsv("fixed_wald_t_reference.tsv");
    let size_factors = read_size_factors("size_factors_ratio.tsv");

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

    let request = ResultRequest::new(
        ResultsContrast::Numeric(vec![0.0, 1.0]),
        WaldTestOptions::t_residual_degrees_of_freedom(),
        CooksCutoff::Disabled,
        IndependentFilteringOptions {
            enabled: false,
            ..Default::default()
        },
    )
    .unwrap();
    let selected = fit
        .results_for_request(&reference_counts(), &request, None)
        .unwrap();
    let mut inference_only = fit.clone();
    inference_only.mu = None;
    inference_only.hat_diagonal = None;
    inference_only.log_like = None;
    inference_only.beta_optim_start = None;
    assert_eq!(
        selected,
        inference_only
            .results_for_request(&reference_counts(), &request, None)
            .unwrap()
    );
    let masked_request = ResultRequest::new(
        ResultsContrast::Numeric(vec![0.0, 1.0]),
        WaldTestOptions::t_residual_degrees_of_freedom(),
        CooksCutoff::Threshold(0.0),
        IndependentFilteringOptions::default(),
    )
    .unwrap();
    inference_only.max_cooks = Some(vec![Some(1.0); rows.len()]);
    let masked = inference_only
        .results_for_request(&reference_counts(), &masked_request, None)
        .unwrap();
    assert!(
        masked
            .rows
            .iter()
            .all(|row| row.pvalue.is_none() && row.padj.is_none())
    );
    assert!(
        ResultRequest::new(
            ResultsContrast::Numeric(vec![0.0, 1.0]),
            WaldTestOptions::normal(),
            CooksCutoff::Threshold(-1.0),
            IndependentFilteringOptions::default()
        )
        .is_err()
    );
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
        assert_option_close(
            selected.rows[gene].pvalue,
            Some(parse_required_f64(row, "pvalue_conditionB")),
            1e-7,
            1e-5,
            &format!("stored request t-Wald p-value gene {gene}"),
        );
    }
}

#[test]
fn fixed_dispersion_force_optim_wald_matches_deseq2_internal_reference() {
    let rows = read_reference_tsv("fixed_force_optim_wald_reference.tsv");
    let expected_mu = read_reference_tsv("fixed_force_optim_mu_full.tsv");
    let expected_hat = read_reference_tsv("fixed_force_optim_hat_full.tsv");
    let size_factors = read_size_factors("size_factors_ratio.tsv");

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
fn fixed_dispersion_cooks_match_deseq2_internal_reference() {
    let expected_cooks = read_reference_tsv("fixed_cooks_full.tsv");
    let size_factors = read_size_factors("size_factors_ratio.tsv");
    let dispersions = read_fixed_dispersions();

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
fn fixed_dispersion_weighted_wald_matches_deseq2_internal_reference() {
    let rows = read_reference_tsv("fixed_wald_weighted_reference.tsv");
    let size_factors = read_size_factors("size_factors_ratio.tsv");
    let observation_weights = read_reference_matrix("observation_weights.tsv");

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
