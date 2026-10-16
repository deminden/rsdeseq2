use super::*;

#[test]
fn native_linear_mu_parametric_wald_preserves_dispersion_intermediates() {
    let counts = native_wald_counts_with_zero_row();
    let design = two_group_design();

    let (fit, results) = native_wald_builder()
        .fit_wald_linear_mu_parametric(&counts, &design, 1)
        .unwrap();

    assert_eq!(fit.design.as_ref().unwrap(), &design);
    assert_eq!(results.rows.len(), counts.n_genes());
    assert_eq!(
        fit.all_zero,
        vec![true, false, false, false, false, false, false, false]
    );

    let disp_gene_est = fit.disp_gene_est.as_ref().unwrap();
    let disp_gene_iter = fit.disp_gene_iter.as_ref().unwrap();
    let disp_fit = fit.disp_fit.as_ref().unwrap();
    let disp_map = fit.disp_map.as_ref().unwrap();
    let dispersion = fit.dispersion.as_ref().unwrap();
    let disp_iter = fit.disp_iter.as_ref().unwrap();
    let disp_outlier = fit.disp_outlier.as_ref().unwrap();
    let disp_converged = fit.dispersion_converged.as_ref().unwrap();

    assert_eq!(disp_gene_est.len(), counts.n_genes());
    assert_eq!(disp_gene_iter.len(), counts.n_genes());
    assert_eq!(disp_fit.len(), counts.n_genes());
    assert_eq!(disp_map.len(), counts.n_genes());
    assert_eq!(dispersion.len(), counts.n_genes());
    assert_eq!(disp_iter.len(), counts.n_genes());
    assert_eq!(disp_outlier.len(), counts.n_genes());
    assert_eq!(disp_converged.len(), counts.n_genes());
    assert!(fit.disp_prior_var.unwrap().is_finite());

    assert!(disp_gene_est[0].is_nan());
    assert_eq!(disp_gene_iter[0], 0);
    assert!(disp_fit[0].is_nan());
    assert!(disp_map[0].is_nan());
    assert!(dispersion[0].is_nan());
    assert_eq!(disp_iter[0], 0);
    assert!(!disp_outlier[0]);
    assert!(!disp_converged[0]);

    for gene in 1..counts.n_genes() {
        assert!(disp_gene_est[gene].is_finite());
        assert!(disp_gene_iter[gene] > 0);
        assert!(disp_fit[gene].is_finite());
        assert!(disp_map[gene].is_finite());
        assert!(dispersion[gene].is_finite());
        assert!(dispersion[gene] > 0.0);
    }

    assert_eq!(fit.beta.as_ref().unwrap().n_cols(), design.n_coefficients());
    assert_eq!(fit.beta_se.as_ref().unwrap().n_rows(), counts.n_genes());
    assert_eq!(
        fit.beta_covariance.as_ref().unwrap().n_rows(),
        counts.n_genes()
    );
    assert_eq!(
        fit.beta_covariance.as_ref().unwrap().n_cols(),
        design.n_coefficients() * design.n_coefficients()
    );
    assert!(fit.beta_covariance.as_ref().unwrap().row(0).unwrap()[0].is_nan());
    assert_eq!(fit.mu.as_ref().unwrap().n_cols(), counts.n_samples());
    assert_eq!(
        fit.hat_diagonal.as_ref().unwrap().n_cols(),
        counts.n_samples()
    );
    assert_eq!(fit.cooks.as_ref().unwrap().n_rows(), counts.n_genes());
    assert_eq!(fit.max_cooks.as_ref().unwrap().len(), counts.n_genes());
    assert_eq!(fit.wald.as_ref().unwrap().pvalue.len(), counts.n_genes());

    assert_eq!(results.rows[0].gene.as_deref(), Some("zero"));
    assert_eq!(results.rows[0].log2_fold_change, None);
    assert_eq!(results.rows[0].pvalue, None);
    assert_eq!(results.rows[0].padj, None);
    assert_eq!(results.rows[0].dispersion, None);
    assert_eq!(results.rows[0].converged, None);

    assert_eq!(results.rows[1].gene.as_deref(), Some("up"));
    assert!(results.rows[1].log2_fold_change.unwrap().is_finite());
    assert!(results.rows[1].lfc_se.unwrap().is_finite());
    assert!(results.rows[1].stat.unwrap().is_finite());
    assert!(results.rows[1].pvalue.unwrap().is_finite());
    assert!(results.rows[1].padj.unwrap().is_finite());
    assert_eq!(results.rows[1].dispersion, Some(dispersion[1]));
    assert_eq!(results.rows[1].pvalue, fit.wald.as_ref().unwrap().pvalue[1]);
}

#[test]
fn top_level_fit_runs_default_automatic_wald_for_last_coefficient() {
    let counts = native_wald_counts_with_zero_row();
    let design = two_group_design();
    let builder = glm_mu_native_wald_builder().automatic_dispersion_means();

    let top_level_fit = builder.fit(&counts, &design).unwrap();
    let (wald_fit, _results) = builder.fit_wald_glm_mu(&counts, &design, 1).unwrap();

    assert_eq!(top_level_fit.counts_summary, wald_fit.counts_summary);
    assert_eq!(top_level_fit.design, wald_fit.design);
    assert_eq!(top_level_fit.all_zero, wald_fit.all_zero);
    assert_eq!(top_level_fit.beta_converged, wald_fit.beta_converged);
    assert_wald_fit_state_matches(&top_level_fit, &wald_fit, "top-level Wald");
    assert_eq!(top_level_fit.wald, wald_fit.wald);
    assert!(top_level_fit.dispersion.is_some());
    assert!(top_level_fit.beta.is_some());
    assert!(top_level_fit.wald.is_some());
}

#[test]
fn top_level_fit_with_results_returns_default_wald_results() {
    let counts = native_wald_counts_with_zero_row();
    let design = two_group_design();
    let builder = glm_mu_native_wald_builder().automatic_dispersion_means();

    let (top_level_fit, top_level_results) = builder.fit_with_results(&counts, &design).unwrap();
    let (wald_fit, wald_results) = builder.fit_wald_glm_mu(&counts, &design, 1).unwrap();

    assert_eq!(top_level_fit.wald, wald_fit.wald);
    assert_wald_fit_state_matches(&top_level_fit, &wald_fit, "top-level Wald results");
    assert_eq!(top_level_results, wald_results);
    assert_eq!(top_level_results.metadata.test_type, Some(TestType::Wald));
    assert_eq!(
        top_level_results.metadata.result_name.as_deref(),
        Some("condition_B_vs_A")
    );
}

#[test]
fn top_level_fit_with_results_accepts_coefficient_name() {
    let counts = native_wald_counts_with_zero_row();
    let design = two_group_design();
    let builder = glm_mu_native_wald_builder().automatic_dispersion_means();

    let (named_fit, named_results) = builder
        .fit_with_results_name(&counts, &design, "condition_B_vs_A")
        .unwrap();
    let (indexed_fit, indexed_results) = builder.fit_wald_glm_mu(&counts, &design, 1).unwrap();

    assert_eq!(named_fit.wald, indexed_fit.wald);
    assert_wald_fit_state_matches(&named_fit, &indexed_fit, "named top-level Wald");
    assert_eq!(named_results, indexed_results);
    assert!(
        builder
            .fit_with_results_name(&counts, &design, "missing")
            .is_err()
    );
}

#[test]
fn top_level_fit_with_results_accepts_r_cleaned_coefficient_name_alias() {
    let counts = native_wald_counts_with_zero_row();
    let design = two_group_design_with_r_cleaned_coefficient();
    let builder = glm_mu_native_wald_builder().automatic_dispersion_means();

    let (named_fit, named_results) = builder
        .fit_with_results_name(&counts, &design, "condition-B 1")
        .unwrap();
    let (indexed_fit, indexed_results) = builder.fit_wald_glm_mu(&counts, &design, 1).unwrap();

    assert_eq!(named_fit.wald, indexed_fit.wald);
    assert_wald_fit_state_matches(&named_fit, &indexed_fit, "R-cleaned named top-level Wald");
    assert_eq!(named_results, indexed_results);
}

#[test]
fn top_level_formula_with_results_accepts_r_cleaned_quoted_coefficient_alias() {
    let counts = native_wald_counts_with_zero_row();
    let levels = [
        "T cell", "T cell", "T cell", "T cell", "B cell", "B cell", "B cell", "B cell",
    ]
    .into_iter()
    .map(String::from)
    .collect::<Vec<_>>();
    let model_frame = FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "cell type".to_string(),
            sample_levels: levels,
            levels: Some(vec!["T cell".to_string(), "B cell".to_string()]),
            reference: None,
        }],
        numeric_covariates: Vec::new(),
    };
    let builder = glm_mu_native_wald_builder()
        .disable_independent_filtering()
        .model_frame(model_frame);

    let (cleaned_fit, cleaned_results) = builder
        .clone()
        .fit_formula_with_results_name(&counts, "~ `cell type`", "cell.type_B.cell_vs_T.cell")
        .unwrap();
    let (raw_fit, raw_results) = builder
        .fit_formula_with_results_name(&counts, "~ `cell type`", "cell type_B cell_vs_T cell")
        .unwrap();

    assert_wald_fit_state_matches(
        &cleaned_fit,
        &raw_fit,
        "quoted formula R-cleaned named top-level Wald",
    );
    assert_eq!(cleaned_results, raw_results);
}

#[test]
fn fixed_wald_supports_stored_model_frame_dot_power_design() {
    let counts = native_wald_counts_with_zero_row();
    let condition = [
        "T cell", "T cell", "T cell", "T cell", "B cell", "B cell", "B cell", "B cell",
    ]
    .into_iter()
    .map(String::from)
    .collect::<Vec<_>>();
    let batch = ["X 0", "Y-1", "X 0", "Y-1", "X 0", "Y-1", "X 0", "Y-1"]
        .into_iter()
        .map(String::from)
        .collect::<Vec<_>>();
    let model_frame = FormulaModelFrame {
        factors: vec![
            FormulaFactorColumn {
                name: "cell type".to_string(),
                sample_levels: condition,
                levels: Some(vec!["T cell".to_string(), "B cell".to_string()]),
                reference: None,
            },
            FormulaFactorColumn {
                name: "batch name".to_string(),
                sample_levels: batch,
                levels: Some(vec!["X 0".to_string(), "Y-1".to_string()]),
                reference: None,
            },
        ],
        numeric_covariates: vec![FormulaNumericColumn {
            name: "dose value".to_string(),
            values: vec![0.0, 1.0, 2.0, 3.0, 0.5, 1.5, 2.5, 3.5],
        }],
    };
    let builder = glm_mu_native_wald_builder()
        .disable_independent_filtering()
        .model_frame(model_frame.clone());

    let dot_design = builder.expanded_formula_design("~ .^2").unwrap();
    let explicit_design = builder
        .expanded_formula_design(
            "~ `cell type` + `batch name` + `dose value` + `cell type`:`batch name` + `cell type`:`dose value` + `batch name`:`dose value`",
        )
        .unwrap();
    assert_eq!(dot_design, explicit_design);

    let coefficient = resolve_coefficient_index(
        &dot_design.standard_design,
        "cell.type_B.cell_vs_T.cell:batch.name_Y.1_vs_X.0",
    )
    .unwrap();
    let dispersions = vec![0.2; counts.n_genes()];
    let (dot_fit, dot_results) = builder
        .clone()
        .fit_fixed_dispersion_wald(
            &counts,
            &dot_design.standard_design,
            &dispersions,
            coefficient,
        )
        .unwrap();
    let (explicit_fit, explicit_results) = builder
        .fit_fixed_dispersion_wald(
            &counts,
            &explicit_design.standard_design,
            &dispersions,
            coefficient,
        )
        .unwrap();

    assert_wald_fit_state_matches(
        &dot_fit,
        &explicit_fit,
        "stored model-frame dot-power fixed Wald",
    );
    assert_eq!(dot_results, explicit_results);
    assert_eq!(dot_fit.current_model_frame(), Some(&model_frame));
}

#[test]
fn top_level_formula_fit_state_retains_model_frame_metadata() {
    let counts = native_wald_counts_with_zero_row();
    let levels = [
        "T cell", "T cell", "T cell", "T cell", "B cell", "B cell", "B cell", "B cell",
    ]
    .into_iter()
    .map(String::from)
    .collect::<Vec<_>>();
    let model_frame = FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "cell type".to_string(),
            sample_levels: levels,
            levels: Some(vec!["T cell".to_string(), "B cell".to_string()]),
            reference: Some("T cell".to_string()),
        }],
        numeric_covariates: vec![FormulaNumericColumn {
            name: "dose value".to_string(),
            values: vec![0.5, 1.5, 0.25, 1.25, 0.75, 1.75, 0.4, 1.4],
        }],
    };
    let builder = glm_mu_native_wald_builder()
        .disable_independent_filtering()
        .model_frame(model_frame.clone());
    let options = CooksReplacementOptions::new(f64::MAX);

    let (fit, _results) = builder
        .clone()
        .fit_formula_with_results(&counts, "~ `cell type`")
        .unwrap();
    let replacement = builder
        .fit_formula_with_results_with_cooks_replacement(&counts, "~ `cell type`", &options)
        .unwrap();

    assert_eq!(fit.current_model_frame(), Some(&model_frame));
    assert_eq!(
        replacement.original_fit.current_model_frame(),
        Some(&model_frame)
    );
    assert_eq!(
        replacement
            .refit_fit
            .as_ref()
            .map(DeseqFit::current_model_frame),
        None
    );

    let transformed = builder
        .clone()
        .fit_formula_with_results(&counts, "~ relevel(`cell type`, ref='B cell')")
        .unwrap()
        .0;
    let transformed_frame = transformed.current_model_frame().unwrap();
    assert_eq!(
        transformed_frame
            .factors
            .iter()
            .map(|factor| factor.name.as_str())
            .collect::<Vec<_>>(),
        vec!["cell type", "relevel(cell type, ref = \"B cell\")"]
    );
    assert_eq!(
        transformed_frame.factors[1].reference.as_deref(),
        Some("B cell")
    );

    let transformed_coefficient = "relevel(cell type, ref = \"B cell\")_T cell_vs_B cell";
    let (named_transformed, _named_results) = builder
        .clone()
        .fit_formula_with_results_name(
            &counts,
            "~ relevel(`cell type`, ref='B cell')",
            transformed_coefficient,
        )
        .unwrap();
    let fit_only_transformed = builder
        .fit_formula_name(
            &counts,
            "~ relevel(`cell type`, ref='B cell')",
            transformed_coefficient,
        )
        .unwrap();
    let replacement_transformed = builder
        .fit_formula_with_results_name_with_cooks_replacement(
            &counts,
            "~ relevel(`cell type`, ref='B cell')",
            "relevel.cell.type..ref....B.cell.._T.cell_vs_B.cell",
            &options,
        )
        .unwrap();

    assert_eq!(
        named_transformed.current_model_frame(),
        transformed.current_model_frame()
    );
    assert_eq!(
        fit_only_transformed.current_model_frame(),
        transformed.current_model_frame()
    );
    assert_eq!(
        replacement_transformed.original_fit.current_model_frame(),
        transformed.current_model_frame()
    );

    let numeric_formula = "~ `cell type` + as.numeric(`dose value`)";
    let (numeric_fit, _numeric_results) = builder
        .fit_formula_with_results_name(&counts, numeric_formula, "dose value_as.numeric")
        .expect("formula with derived numeric metadata should fit");
    let numeric_frame = numeric_fit.current_model_frame().unwrap();
    assert_eq!(
        numeric_frame
            .numeric_covariates
            .iter()
            .map(|covariate| covariate.name.as_str())
            .collect::<Vec<_>>(),
        vec!["dose value", "dose value_as.numeric"]
    );
    assert_eq!(
        numeric_frame
            .resolved_numeric_covariate_by_alias("dose.value_as.numeric")
            .unwrap()
            .unwrap(),
        ResolvedFormulaNumericCovariate {
            name: "dose value_as.numeric",
            values: &[0.5, 1.5, 0.25, 1.25, 0.75, 1.75, 0.4, 1.4],
        }
    );

    let factor_code_formula = "~ as.double(`cell type`)";
    let (factor_code_fit, _factor_code_results) = builder
        .clone()
        .fit_formula_with_results_name(&counts, factor_code_formula, "cell type_as.double")
        .expect("formula with factor-code numeric metadata should fit");
    let factor_code_replacement = builder
        .fit_formula_with_results_name_with_cooks_replacement(
            &counts,
            factor_code_formula,
            "cell.type_as.double",
            &options,
        )
        .expect("factor-code formula replacement should resolve cleaned coefficient alias");
    let factor_code_frame = factor_code_fit.current_model_frame().unwrap();
    assert_eq!(
        factor_code_frame
            .numeric_covariates
            .iter()
            .map(|covariate| covariate.name.as_str())
            .collect::<Vec<_>>(),
        vec!["dose value", "cell type_as.double"]
    );
    assert_eq!(
        factor_code_frame
            .resolved_numeric_covariate_by_alias("cell.type_as.double")
            .unwrap()
            .unwrap(),
        ResolvedFormulaNumericCovariate {
            name: "cell type_as.double",
            values: &[1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0],
        }
    );
    assert_eq!(factor_code_replacement.refit_plan.n_refit, 0);
    assert_eq!(
        factor_code_replacement.original_fit.current_model_frame(),
        factor_code_fit.current_model_frame()
    );
    assert_eq!(
        factor_code_replacement
            .results
            .metadata
            .result_name
            .as_deref(),
        Some("cell type_as.double")
    );
}

#[test]
fn top_level_formula_list_contrast_accepts_r_cleaned_quoted_coefficient_alias() {
    let counts = native_wald_counts_with_zero_row();
    let levels = [
        "T cell", "T cell", "T cell", "T cell", "B cell", "B cell", "B cell", "B cell",
    ]
    .into_iter()
    .map(String::from)
    .collect::<Vec<_>>();
    let model_frame = FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "cell type".to_string(),
            sample_levels: levels,
            levels: Some(vec!["T cell".to_string(), "B cell".to_string()]),
            reference: None,
        }],
        numeric_covariates: Vec::new(),
    };
    let builder = glm_mu_native_wald_builder()
        .disable_independent_filtering()
        .model_frame(model_frame);

    let (cleaned_fit, cleaned_results) = builder
        .clone()
        .fit_formula_with_results_contrast_request(
            &counts,
            "~ `cell type`",
            &ResultsContrast::list(vec!["cell.type_B.cell_vs_T.cell".into()], Vec::new()),
        )
        .unwrap();
    let (raw_fit, raw_results) = builder
        .fit_formula_with_results_contrast_request(
            &counts,
            "~ `cell type`",
            &ResultsContrast::list(vec!["cell type_B cell_vs_T cell".into()], Vec::new()),
        )
        .unwrap();

    assert_wald_fit_state_matches(
        &cleaned_fit,
        &raw_fit,
        "quoted formula R-cleaned list top-level Wald",
    );
    assert_eq!(cleaned_results.rows, raw_results.rows);
    assert_eq!(
        cleaned_results.metadata.comparison.as_deref(),
        Some("coefficient list contrast: cell.type_B.cell_vs_T.cell effect")
    );
    assert_eq!(
        raw_results.metadata.comparison.as_deref(),
        Some("coefficient list contrast: cell type_B cell_vs_T cell effect")
    );
}

#[test]
fn top_level_formula_list_contrast_accepts_cleaned_positive_and_negative_aliases() {
    let counts = native_wald_counts_with_zero_row();
    let levels = [
        "T cell", "T cell", "T cell", "B cell", "B cell", "B cell", "NK cell", "NK cell",
    ]
    .into_iter()
    .map(String::from)
    .collect::<Vec<_>>();
    let model_frame = FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "cell type".to_string(),
            sample_levels: levels,
            levels: Some(vec![
                "T cell".to_string(),
                "B cell".to_string(),
                "NK cell".to_string(),
            ]),
            reference: None,
        }],
        numeric_covariates: Vec::new(),
    };
    let builder = glm_mu_native_wald_builder()
        .disable_independent_filtering()
        .model_frame(model_frame);

    let (cleaned_fit, cleaned_results) = builder
        .clone()
        .fit_formula_with_results_contrast_request(
            &counts,
            "~ `cell type`",
            &ResultsContrast::list(
                vec!["cell.type_NK.cell_vs_T.cell".into()],
                vec!["cell.type_B.cell_vs_T.cell".into()],
            ),
        )
        .unwrap();
    let (raw_fit, raw_results) = builder
        .fit_formula_with_results_contrast_request(
            &counts,
            "~ `cell type`",
            &ResultsContrast::list(
                vec!["cell type_NK cell_vs_T cell".into()],
                vec!["cell type_B cell_vs_T cell".into()],
            ),
        )
        .unwrap();

    assert_wald_fit_state_matches(
        &cleaned_fit,
        &raw_fit,
        "three-level formula R-cleaned list Wald",
    );
    assert_eq!(cleaned_results.rows, raw_results.rows);
    assert_eq!(
        cleaned_results.metadata.comparison.as_deref(),
        Some(
            "coefficient list contrast: cell.type_NK.cell_vs_T.cell vs cell.type_B.cell_vs_T.cell"
        )
    );
    assert_eq!(
        raw_results.metadata.comparison.as_deref(),
        Some(
            "coefficient list contrast: cell type_NK cell_vs_T cell vs cell type_B cell_vs_T cell"
        )
    );
}
