use super::*;

#[test]
fn native_glm_mu_lrt_preserves_diagnostics_and_all_zero_rows() {
    let counts = counts_with_zero_row();
    let full = full_design();
    let reduced = reduced_design();

    let (fit, results) = native_lrt_builder()
        .fit_lrt_glm_mu(&counts, &full, &reduced, 1)
        .unwrap();

    assert_eq!(fit.design.as_ref().unwrap(), &full);
    assert_eq!(fit.reduced_design.as_ref().unwrap(), &reduced);
    assert_eq!(
        fit.all_zero,
        vec![true, false, false, false, false, false, false, false]
    );
    assert!(fit.disp_prior_var.unwrap().is_finite());
    assert_eq!(fit.disp_gene_est.as_ref().unwrap().len(), counts.n_genes());
    assert_eq!(fit.disp_gene_iter.as_ref().unwrap().len(), counts.n_genes());
    assert_eq!(fit.disp_fit.as_ref().unwrap().len(), counts.n_genes());
    assert_eq!(fit.disp_map.as_ref().unwrap().len(), counts.n_genes());
    assert_eq!(fit.dispersion.as_ref().unwrap().len(), counts.n_genes());
    assert_eq!(fit.disp_gene_iter.as_ref().unwrap()[0], 0);
    assert!(
        fit.disp_gene_iter.as_ref().unwrap()[1..]
            .iter()
            .all(|iterations| *iterations > 0)
    );

    assert_eq!(fit.beta.as_ref().unwrap().n_cols(), full.n_coefficients());
    assert_eq!(
        fit.reduced_log_like.as_ref().unwrap().len(),
        counts.n_genes()
    );
    assert_eq!(fit.reduced_mu.as_ref().unwrap().n_rows(), counts.n_genes());
    assert_eq!(
        fit.reduced_mu.as_ref().unwrap().n_cols(),
        counts.n_samples()
    );
    assert_eq!(
        fit.reduced_hat_diagonal.as_ref().unwrap().n_rows(),
        counts.n_genes()
    );
    assert_eq!(
        fit.reduced_hat_diagonal.as_ref().unwrap().n_cols(),
        counts.n_samples()
    );
    assert!(fit.reduced_mu.as_ref().unwrap().row(0).unwrap()[0].is_nan());
    assert!(fit.reduced_hat_diagonal.as_ref().unwrap().row(0).unwrap()[0].is_nan());
    assert!(fit.reduced_mu.as_ref().unwrap().row(1).unwrap()[0].is_finite());
    assert!(fit.reduced_hat_diagonal.as_ref().unwrap().row(1).unwrap()[0].is_finite());
    assert_eq!(
        fit.reduced_beta_converged.as_deref(),
        Some(&[false, true, true, true, true, true, true, true][..])
    );
    assert_eq!(fit.lrt.as_ref().unwrap().degrees_of_freedom, 1);
    assert_eq!(fit.lrt.as_ref().unwrap().deviance[0], None);
    assert_eq!(fit.lrt.as_ref().unwrap().pvalue[0], None);
    assert_eq!(fit.cooks.as_ref().unwrap().n_rows(), counts.n_genes());

    assert_eq!(results.rows.len(), counts.n_genes());
    assert_eq!(results.rows[0].gene.as_deref(), Some("zero"));
    assert_eq!(results.rows[0].pvalue, None);
    assert_eq!(results.rows[0].padj, None);
    assert_eq!(results.rows[0].converged, None);
    assert!(results.rows[1].stat.unwrap().is_finite());
    assert!(results.rows[1].pvalue.unwrap().is_finite());
    assert_eq!(results.rows[1].pvalue, fit.lrt.as_ref().unwrap().pvalue[1]);
}

#[test]
fn native_glm_mu_lrt_matches_fixed_pipeline_when_reusing_final_dispersions() {
    let counts = counts_with_zero_row();
    let full = full_design();
    let reduced = reduced_design();
    let builder = native_lrt_builder();

    let (native_fit, native_results) = builder.fit_lrt_glm_mu(&counts, &full, &reduced, 1).unwrap();
    let final_dispersions = native_fit.dispersion.as_ref().unwrap().clone();
    let (fixed_fit, fixed_results) = builder
        .fit_fixed_dispersion_lrt(&counts, &full, &reduced, &final_dispersions, 1)
        .unwrap();

    for gene in 0..counts.n_genes() {
        let native_disp = native_fit.dispersion.as_ref().unwrap()[gene];
        let fixed_disp = fixed_fit.dispersion.as_ref().unwrap()[gene];
        if native_disp.is_nan() {
            assert!(fixed_disp.is_nan());
        } else {
            assert_relative_eq!(native_disp, fixed_disp, epsilon = 1e-12);
        }
        assert_eq!(
            native_fit.lrt.as_ref().unwrap().deviance[gene],
            fixed_fit.lrt.as_ref().unwrap().deviance[gene]
        );
        assert_eq!(
            native_results.rows[gene].pvalue,
            fixed_results.rows[gene].pvalue
        );
        assert_eq!(
            native_results.rows[gene].padj,
            fixed_results.rows[gene].padj
        );
        assert_eq!(
            native_results.rows[gene].log2_fold_change,
            fixed_results.rows[gene].log2_fold_change
        );
        for sample in 0..counts.n_samples() {
            let native_mu = native_fit.reduced_mu.as_ref().unwrap().row(gene).unwrap()[sample];
            let fixed_mu = fixed_fit.reduced_mu.as_ref().unwrap().row(gene).unwrap()[sample];
            if native_mu.is_nan() {
                assert!(fixed_mu.is_nan());
            } else {
                assert_relative_eq!(native_mu, fixed_mu, epsilon = 1e-12);
            }
            let native_hat = native_fit
                .reduced_hat_diagonal
                .as_ref()
                .unwrap()
                .row(gene)
                .unwrap()[sample];
            let fixed_hat = fixed_fit
                .reduced_hat_diagonal
                .as_ref()
                .unwrap()
                .row(gene)
                .unwrap()[sample];
            if native_hat.is_nan() {
                assert!(fixed_hat.is_nan());
            } else {
                assert_relative_eq!(native_hat, fixed_hat, epsilon = 1e-12);
            }
        }
    }
}

#[test]
fn top_level_lrt_fit_only_helpers_match_result_helpers() {
    let counts = counts_with_zero_row();
    let full = full_design();
    let reduced = reduced_design();
    let builder = native_lrt_builder();
    let spec = ContrastSpec::coefficient_name("condition_B_vs_A");
    let levels = ["A", "A", "A", "A", "B", "B", "B", "B"]
        .into_iter()
        .map(String::from)
        .collect::<Vec<_>>();
    let formula_builder = builder.clone().model_frame(FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "condition".to_string(),
            sample_levels: levels.clone(),
            levels: None,
            reference: None,
        }],
        numeric_covariates: Vec::new(),
    });

    let default_fit = builder.fit_lrt(&counts, &full, &reduced).unwrap();
    let (default_result_fit, _default_results) = builder
        .fit_lrt_with_results(&counts, &full, &reduced)
        .unwrap();
    assert_lrt_fit_state_matches(&default_fit, &default_result_fit, "fit-only default LRT");
    let formula_fit = formula_builder
        .fit_lrt_formula(&counts, "~ condition", "~ 1")
        .unwrap();
    let (formula_result_fit, _formula_results) = formula_builder
        .fit_lrt_formula_with_results(&counts, "~ condition", "~ 1")
        .unwrap();
    assert_lrt_fit_state_matches(&formula_fit, &default_fit, "fit-only formula LRT");
    assert_lrt_fit_state_matches(
        &formula_fit,
        &formula_result_fit,
        "fit-only formula LRT result companion",
    );

    let named_fit = builder
        .fit_lrt_name(&counts, &full, &reduced, "condition_B_vs_A")
        .unwrap();
    let (named_result_fit, _named_results) = builder
        .fit_lrt_with_results_name(&counts, &full, &reduced, "condition_B_vs_A")
        .unwrap();
    assert_lrt_fit_state_matches(&named_fit, &named_result_fit, "fit-only named LRT");
    let formula_named_fit = formula_builder
        .fit_lrt_formula_name(&counts, "~ condition", "~ 1", "condition_B_vs_A")
        .unwrap();
    let (formula_named_result_fit, _formula_named_results) = formula_builder
        .fit_lrt_formula_with_results_name(&counts, "~ condition", "~ 1", "condition_B_vs_A")
        .unwrap();
    assert_lrt_fit_state_matches(&formula_named_fit, &named_fit, "fit-only formula named LRT");
    assert_lrt_fit_state_matches(
        &formula_named_fit,
        &formula_named_result_fit,
        "fit-only formula named LRT result companion",
    );

    let contrast_fit = builder
        .fit_lrt_contrast(&counts, &full, &reduced, &[0.0, 1.0])
        .unwrap();
    let (contrast_result_fit, _contrast_results) = builder
        .fit_lrt_with_results_contrast(&counts, &full, &reduced, &[0.0, 1.0])
        .unwrap();
    assert_lrt_fit_state_matches(
        &contrast_fit,
        &contrast_result_fit,
        "fit-only numeric LRT contrast",
    );

    let spec_fit = builder
        .fit_lrt_contrast_spec(&counts, &full, &reduced, &spec)
        .unwrap();
    let (spec_result_fit, _spec_results) = builder
        .fit_lrt_with_results_contrast_spec(&counts, &full, &reduced, &spec)
        .unwrap();
    assert_lrt_fit_state_matches(&spec_fit, &spec_result_fit, "fit-only named LRT contrast");

    let factor_contrast = FactorLevelContrast {
        factor: "condition",
        numerator: "B",
        denominator: "A",
        reference: None,
        sample_levels: &levels,
    };
    let factor_fit = builder
        .fit_lrt_factor_level_contrast(&counts, &full, &reduced, factor_contrast)
        .unwrap();
    let (factor_result_fit, _factor_results) = builder
        .fit_lrt_with_results_factor_level_contrast(&counts, &full, &reduced, factor_contrast)
        .unwrap();
    assert_lrt_fit_state_matches(
        &factor_fit,
        &factor_result_fit,
        "fit-only factor-level LRT contrast",
    );
}

#[test]
fn top_level_fit_lrt_runs_default_automatic_lrt_for_last_coefficient() {
    let counts = counts_with_zero_row();
    let full = full_design();
    let reduced = reduced_design();
    let builder = native_lrt_builder().automatic_dispersion_means();

    let top_level_fit = builder.fit_lrt(&counts, &full, &reduced).unwrap();
    let (lrt_fit, _results) = builder.fit_lrt_glm_mu(&counts, &full, &reduced, 1).unwrap();

    assert_eq!(top_level_fit.counts_summary, lrt_fit.counts_summary);
    assert_eq!(top_level_fit.design, lrt_fit.design);
    assert_eq!(top_level_fit.reduced_design, lrt_fit.reduced_design);
    assert_eq!(top_level_fit.all_zero, lrt_fit.all_zero);
    assert_eq!(
        top_level_fit.reduced_beta_converged,
        lrt_fit.reduced_beta_converged
    );
    assert_lrt_fit_state_matches(&top_level_fit, &lrt_fit, "top-level LRT");
    assert_eq!(top_level_fit.lrt, lrt_fit.lrt);
    assert!(top_level_fit.dispersion.is_some());
    assert!(matches!(
        top_level_fit.dispersion_trend.as_ref(),
        Some(DispersionTrendFit::Mean(_))
    ));
    assert!(top_level_fit.beta.is_some());
    assert!(top_level_fit.lrt.is_some());

    let transformed = top_level_fit.vst(&counts).unwrap();
    assert_eq!(transformed.n_rows(), counts.n_genes());
    assert_eq!(transformed.n_cols(), counts.n_samples());
    assert_eq!(
        transformed,
        top_level_fit
            .variance_stabilizing_transform(&counts)
            .unwrap()
    );
    assert!(transformed.as_slice().iter().all(|value| value.is_finite()));
}

#[test]
fn top_level_fit_lrt_with_results_returns_default_lrt_results() {
    let counts = counts_with_zero_row();
    let full = full_design();
    let reduced = reduced_design();
    let builder = native_lrt_builder().automatic_dispersion_means();

    let (top_level_fit, top_level_results) = builder
        .fit_lrt_with_results(&counts, &full, &reduced)
        .unwrap();
    let (lrt_fit, lrt_results) = builder.fit_lrt_glm_mu(&counts, &full, &reduced, 1).unwrap();

    assert_eq!(top_level_fit.lrt, lrt_fit.lrt);
    assert_lrt_fit_state_matches(&top_level_fit, &lrt_fit, "top-level LRT results");
    assert_eq!(top_level_results, lrt_results);
    assert_eq!(top_level_results.metadata.test_type, Some(TestType::Lrt));
    assert_eq!(
        top_level_results.metadata.result_name.as_deref(),
        Some("condition_B_vs_A")
    );
}

#[test]
fn top_level_fit_lrt_with_results_accepts_coefficient_name() {
    let counts = counts_with_zero_row();
    let full = full_design();
    let reduced = reduced_design();
    let builder = native_lrt_builder().automatic_dispersion_means();

    let (named_fit, named_results) = builder
        .fit_lrt_with_results_name(&counts, &full, &reduced, "condition_B_vs_A")
        .unwrap();
    let (indexed_fit, indexed_results) =
        builder.fit_lrt_glm_mu(&counts, &full, &reduced, 1).unwrap();

    assert_eq!(named_fit.lrt, indexed_fit.lrt);
    assert_lrt_fit_state_matches(&named_fit, &indexed_fit, "named LRT");
    assert_eq!(named_results, indexed_results);
    assert!(
        builder
            .fit_lrt_with_results_name(&counts, &full, &reduced, "missing")
            .is_err()
    );
}

#[test]
fn top_level_fit_lrt_with_results_accepts_r_cleaned_coefficient_name_alias() {
    let counts = counts_with_zero_row();
    let full = full_design_with_r_cleaned_coefficient();
    let reduced = reduced_design();
    let builder = native_lrt_builder().automatic_dispersion_means();

    let (named_fit, named_results) = builder
        .fit_lrt_with_results_name(&counts, &full, &reduced, "condition-B 1")
        .unwrap();
    let (indexed_fit, indexed_results) =
        builder.fit_lrt_glm_mu(&counts, &full, &reduced, 1).unwrap();

    assert_eq!(named_fit.lrt, indexed_fit.lrt);
    assert_lrt_fit_state_matches(&named_fit, &indexed_fit, "R-cleaned named LRT");
    assert_eq!(named_results, indexed_results);
}

#[test]
fn native_linear_mu_lrt_matches_fixed_pipeline_when_reusing_final_dispersions() {
    let counts = counts_with_zero_row();
    let full = full_design();
    let reduced = reduced_design();
    let builder = native_lrt_builder();

    let (native_fit, native_results) = builder
        .fit_lrt_linear_mu(&counts, &full, &reduced, 1)
        .unwrap();
    let final_dispersions = native_fit.dispersion.as_ref().unwrap().clone();
    let (fixed_fit, fixed_results) = builder
        .fit_fixed_dispersion_lrt(&counts, &full, &reduced, &final_dispersions, 1)
        .unwrap();

    assert!(native_fit.disp_prior_var.unwrap().is_finite());
    assert_eq!(native_fit.lrt.as_ref().unwrap().degrees_of_freedom, 1);
    assert_eq!(native_results.rows[0].pvalue, None);

    for gene in 0..counts.n_genes() {
        assert_eq!(
            native_fit.lrt.as_ref().unwrap().pvalue[gene],
            fixed_fit.lrt.as_ref().unwrap().pvalue[gene]
        );
        assert_eq!(
            native_results.rows[gene].pvalue,
            fixed_results.rows[gene].pvalue
        );
        assert_eq!(
            native_results.rows[gene].stat,
            fixed_results.rows[gene].stat
        );
    }
}

#[test]
fn native_linear_mu_local_lrt_matches_fixed_pipeline_when_reusing_final_dispersions() {
    let counts = counts_with_zero_row();
    let full = full_design();
    let reduced = reduced_design();
    let builder = native_lrt_builder().fit_type(FitType::Local);

    let (native_fit, native_results) = builder
        .fit_lrt_linear_mu(&counts, &full, &reduced, 1)
        .unwrap();
    let final_dispersions = native_fit.dispersion.as_ref().unwrap().clone();
    let (fixed_fit, fixed_results) = builder
        .fit_fixed_dispersion_lrt(&counts, &full, &reduced, &final_dispersions, 1)
        .unwrap();

    assert!(native_fit.disp_prior_var.unwrap().is_finite());
    assert_eq!(
        native_fit.disp_fit.as_ref().unwrap().len(),
        counts.n_genes()
    );
    assert_eq!(native_fit.lrt.as_ref().unwrap().degrees_of_freedom, 1);
    assert_eq!(native_results.rows[0].pvalue, None);

    for gene in 0..counts.n_genes() {
        assert_eq!(
            native_fit.lrt.as_ref().unwrap().pvalue[gene],
            fixed_fit.lrt.as_ref().unwrap().pvalue[gene]
        );
        assert_eq!(
            native_results.rows[gene].pvalue,
            fixed_results.rows[gene].pvalue
        );
        assert_eq!(
            native_results.rows[gene].stat,
            fixed_results.rows[gene].stat
        );
    }
}

#[test]
fn native_glm_mu_local_lrt_matches_fixed_pipeline_when_reusing_final_dispersions() {
    let counts = counts_with_zero_row();
    let full = full_design();
    let reduced = reduced_design();
    let builder = native_lrt_builder().fit_type(FitType::Local);

    let (native_fit, native_results) = builder.fit_lrt_glm_mu(&counts, &full, &reduced, 1).unwrap();
    let final_dispersions = native_fit.dispersion.as_ref().unwrap().clone();
    let (fixed_fit, fixed_results) = builder
        .fit_fixed_dispersion_lrt(&counts, &full, &reduced, &final_dispersions, 1)
        .unwrap();

    assert!(native_fit.disp_prior_var.unwrap().is_finite());
    assert_eq!(
        native_fit.disp_fit.as_ref().unwrap().len(),
        counts.n_genes()
    );
    assert_eq!(native_fit.lrt.as_ref().unwrap().degrees_of_freedom, 1);
    assert_eq!(native_results.rows[0].pvalue, None);

    for gene in 0..counts.n_genes() {
        let native_disp = native_fit.dispersion.as_ref().unwrap()[gene];
        let fixed_disp = fixed_fit.dispersion.as_ref().unwrap()[gene];
        if native_disp.is_nan() {
            assert!(fixed_disp.is_nan());
        } else {
            assert_relative_eq!(native_disp, fixed_disp, epsilon = 1e-12);
        }
        assert_eq!(
            native_fit.lrt.as_ref().unwrap().pvalue[gene],
            fixed_fit.lrt.as_ref().unwrap().pvalue[gene]
        );
        assert_eq!(
            native_results.rows[gene].pvalue,
            fixed_results.rows[gene].pvalue
        );
        assert_eq!(
            native_results.rows[gene].stat,
            fixed_results.rows[gene].stat
        );
    }
}
