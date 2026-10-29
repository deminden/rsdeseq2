use super::*;

#[test]
fn native_glm_mu_lrt_cooks_replacement_refit_merges_refit_rows() {
    let counts = counts_with_zero_row();
    let full = full_design();
    let reduced = reduced_design();
    let builder = native_lrt_builder();

    let output = builder
        .fit_lrt_glm_mu_with_cooks_replacement(
            &counts,
            &full,
            &reduced,
            1,
            &CooksReplacementOptions {
                trim: 0.2,
                cooks_cutoff: 0.0,
                min_replicates: 3,
                which_samples: Some(vec![true, false, false, false, false, false, false, false]),
            },
        )
        .unwrap();

    assert!(output.original_fit.lrt.is_some());
    assert!(output.refit_plan.n_refit > 0);
    assert!(output.refit_plan.should_refit);
    assert!(!output.refit_plan.refit_rows.is_empty());
    assert!(output.refit_fit.as_ref().unwrap().lrt.is_some());
    assert!(output.refit_results.is_some());
    assert_ne!(
        output.refit_plan.replacement.replaced_counts.as_slice(),
        counts.as_slice()
    );
    assert_eq!(
        output.refit_plan.replacement.replaceable_samples,
        vec![true, false, false, false, false, false, false, false]
    );

    let refit_results = output.refit_results.as_ref().unwrap();
    for gene in output.refit_plan.refit_rows.iter().copied() {
        assert_eq!(
            output.results.rows[gene].log2_fold_change,
            refit_results.rows[gene].log2_fold_change
        );
        assert_eq!(
            output.results.rows[gene].stat,
            refit_results.rows[gene].stat
        );
        assert_eq!(
            output.results.rows[gene].pvalue,
            refit_results.rows[gene].pvalue
        );
        assert_eq!(
            output.results.rows[gene].max_cooks,
            output.refit_plan.post_refit_max_cooks[gene]
        );
        assert_relative_eq!(
            output.results.rows[gene].base_mean,
            output.refit_plan.replaced_base_mean[gene],
            epsilon = 1e-12
        );
    }
    for gene in output.refit_plan.new_all_zero_rows.iter().copied() {
        assert_eq!(output.results.rows[gene].pvalue, Some(1.0));
        assert_eq!(output.results.rows[gene].stat, Some(0.0));
        assert_eq!(output.results.rows[gene].log2_fold_change, Some(0.0));
        assert_eq!(output.results.rows[gene].lfc_se, Some(0.0));
        assert_eq!(output.results.rows[gene].dispersion, None);
    }
}

#[test]
fn native_glm_mu_lrt_contrast_cooks_replacement_refit_merges_refit_rows() {
    let counts = counts_with_zero_row();
    let full = full_design();
    let reduced = reduced_design();
    let builder = native_lrt_builder();
    let options = CooksReplacementOptions {
        trim: 0.2,
        cooks_cutoff: 0.0,
        min_replicates: 3,
        which_samples: Some(vec![true, false, false, false, false, false, false, false]),
    };

    let output = builder
        .fit_lrt_glm_mu_contrast_with_cooks_replacement(
            &counts,
            &full,
            &reduced,
            &[0.0, 1.0],
            &options,
        )
        .unwrap();
    let named_output = builder
        .fit_lrt_glm_mu_contrast_spec_with_cooks_replacement(
            &counts,
            &full,
            &reduced,
            &ContrastSpec::coefficient_name("condition_B_vs_A"),
            &options,
        )
        .unwrap();

    assert!(output.original_fit.lrt.is_some());
    assert!(output.refit_plan.should_refit);
    assert!(output.refit_fit.as_ref().unwrap().lrt.is_some());
    let refit_results = output.refit_results.as_ref().unwrap();
    for gene in output.refit_plan.refit_rows.iter().copied() {
        assert_eq!(
            output.results.rows[gene].log2_fold_change,
            refit_results.rows[gene].log2_fold_change
        );
        assert_eq!(
            output.results.rows[gene].stat,
            refit_results.rows[gene].stat
        );
        assert_eq!(
            output.results.rows[gene].pvalue,
            refit_results.rows[gene].pvalue
        );
    }
    assert_eq!(named_output.results.rows, output.results.rows);
    assert_eq!(
        named_output.results.metadata.result_name.as_deref(),
        Some("condition_B_vs_A")
    );
    assert_eq!(
        named_output
            .refit_results
            .as_ref()
            .unwrap()
            .metadata
            .comparison
            .as_deref(),
        Some("coefficient condition_B_vs_A")
    );
}

#[test]
fn native_glm_mu_lrt_factor_level_contrast_cooks_replacement_preserves_metadata() {
    let counts = counts_with_zero_row();
    let full = full_design();
    let reduced = reduced_design();
    let builder = native_lrt_builder();
    let levels = ["A", "A", "A", "A", "B", "B", "B", "B"]
        .into_iter()
        .map(String::from)
        .collect::<Vec<_>>();

    let (expected_fit, expected_results) = builder
        .clone()
        .fit_lrt_glm_mu_contrast(&counts, &full, &reduced, &[0.0, 1.0])
        .unwrap();
    let output = builder
        .fit_lrt_glm_mu_factor_level_contrast_with_cooks_replacement(
            &counts,
            &full,
            &reduced,
            FactorLevelContrast::new("condition", "B", "A", &levels),
            &CooksReplacementOptions {
                trim: 0.2,
                cooks_cutoff: f64::MAX,
                min_replicates: 3,
                which_samples: None,
            },
        )
        .unwrap();

    assert_eq!(output.refit_plan.n_refit, 0);
    assert!(output.refit_fit.is_none());
    assert!(output.refit_results.is_none());
    assert_lrt_fit_state_matches(
        &output.original_fit,
        &expected_fit,
        "native GLM-mu LRT factor-level contrast no-refit original",
    );
    assert_eq!(output.original_results.rows, expected_results.rows);
    assert_eq!(output.results.rows, expected_results.rows);
    assert_eq!(
        output.original_results.metadata.result_name.as_deref(),
        Some("condition_B_vs_A")
    );
    assert_eq!(
        output.original_results.metadata.comparison.as_deref(),
        Some("factor-level contrast: condition B vs A")
    );
    assert_eq!(
        output.results.metadata.result_name.as_deref(),
        Some("condition_B_vs_A")
    );
    assert_eq!(
        output.results.metadata.comparison.as_deref(),
        Some("factor-level contrast: condition B vs A")
    );
}

#[test]
fn native_lrt_replacement_numeric_contrasts_do_not_infer_low_count_gate_from_multifactor_frame() {
    let counts = CountMatrix::from_row_major_u32_with_names(
        3,
        6,
        vec![
            1, 20, 21, 20, 20, 20, //
            10, 11, 12, 15, 16, 17, //
            30, 32, 31, 50, 52, 51,
        ],
        Some(vec!["low_count".into(), "stable".into(), "up".into()]),
        None,
    )
    .unwrap();
    let full = DesignMatrix::from_row_major(
        6,
        3,
        vec![
            1.0, 0.0, 0.0, //
            1.0, 0.0, 1.0, //
            1.0, 0.0, 0.0, //
            1.0, 1.0, 1.0, //
            1.0, 1.0, 0.0, //
            1.0, 1.0, 1.0,
        ],
        Some(vec![
            "Intercept".into(),
            "condition_B_vs_A".into(),
            "batch_Y_vs_X".into(),
        ]),
    )
    .unwrap();
    let reduced = DesignMatrix::intercept_only(counts.n_samples()).unwrap();
    let model_frame = FormulaModelFrame {
        factors: vec![
            FormulaFactorColumn {
                name: "condition".to_string(),
                sample_levels: ["A", "A", "A", "B", "B", "B"]
                    .into_iter()
                    .map(String::from)
                    .collect(),
                levels: Some(vec!["A".to_string(), "B".to_string()]),
                reference: None,
            },
            FormulaFactorColumn {
                name: "batch".to_string(),
                sample_levels: ["X", "Y", "X", "Y", "X", "Y"]
                    .into_iter()
                    .map(String::from)
                    .collect(),
                levels: Some(vec!["X".to_string(), "Y".to_string()]),
                reference: None,
            },
        ],
        numeric_covariates: Vec::new(),
    };
    let options = CooksReplacementOptions::new(f64::MAX);
    let builder = DeseqBuilder::new()
        .fit_type(FitType::Mean)
        .size_factors(vec![1.0; 6])
        .gene_wise_dispersion_options(GeneWiseDispersionOptions {
            use_cox_reid: false,
            fit_method: GeneWiseDispersionFitMethod::Grid,
            niter: 2,
            ..GeneWiseDispersionOptions::default()
        })
        .cooks_cutoff_threshold(0.0)
        .disable_independent_filtering();

    let generic = builder
        .clone()
        .fit_lrt_with_results_contrast_with_cooks_replacement(
            &counts,
            &full,
            &reduced,
            &[0.0, 1.0, 0.0],
            &options,
        )
        .unwrap();
    let formula = builder
        .model_frame(model_frame)
        .fit_lrt_with_results_contrast_with_cooks_replacement(
            &counts,
            &full,
            &reduced,
            &[0.0, 1.0, 0.0],
            &options,
        )
        .unwrap();

    assert_eq!(generic.results.rows[0].gene.as_deref(), Some("low_count"));
    assert_eq!(formula.results.rows, generic.results.rows);
    assert_eq!(formula.results.metadata, generic.results.metadata);
}

#[test]
fn top_level_fit_lrt_with_results_cooks_replacement_runs_default_automatic_lrt() {
    let counts = counts_with_zero_row();
    let full = full_design();
    let reduced = reduced_design();
    let levels = ["A", "A", "A", "A", "B", "B", "B", "B"]
        .into_iter()
        .map(String::from)
        .collect::<Vec<_>>();
    let model_frame = FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "condition".to_string(),
            sample_levels: levels,
            levels: None,
            reference: None,
        }],
        numeric_covariates: vec![FormulaNumericColumn {
            name: "exposure".to_string(),
            values: vec![0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7],
        }],
    };
    let builder = native_lrt_builder()
        .model_frame(model_frame)
        .automatic_dispersion_means();
    let replacement_options = CooksReplacementOptions {
        trim: 0.2,
        cooks_cutoff: 0.0,
        min_replicates: 3,
        which_samples: Some(vec![true, false, false, false, false, false, false, false]),
    };

    let top_level_output = builder
        .fit_lrt_with_results_with_cooks_replacement(&counts, &full, &reduced, &replacement_options)
        .unwrap();
    let formula_output = builder
        .fit_lrt_formula_with_results_with_cooks_replacement(
            &counts,
            "~ condition",
            "~ 1",
            &replacement_options,
        )
        .unwrap();
    let lrt_output = builder
        .fit_lrt_glm_mu_with_cooks_replacement(&counts, &full, &reduced, 1, &replacement_options)
        .unwrap();
    let offset_values = builder.current_model_frame().unwrap().numeric_covariates[0]
        .values
        .to_vec();
    let mut normalization_values = Vec::with_capacity(counts.n_genes() * counts.n_samples());
    for _gene in 0..counts.n_genes() {
        for offset in &offset_values {
            normalization_values.push(offset.exp());
        }
    }
    let normalization_factors =
        RowMajorMatrix::from_row_major(counts.n_genes(), counts.n_samples(), normalization_values)
            .unwrap();
    let explicit_offset_output = builder
        .clone()
        .normalization_factors(normalization_factors.clone())
        .fit_lrt_formula_with_results_with_cooks_replacement(
            &counts,
            "~ condition",
            "~ 1",
            &replacement_options,
        )
        .unwrap();
    let formula_offset_output = builder
        .fit_lrt_formula_with_results_with_cooks_replacement(
            &counts,
            "~ condition + offset(exposure)",
            "~ 1",
            &replacement_options,
        )
        .unwrap();

    assert_eq!(top_level_output.refit_plan, lrt_output.refit_plan);
    assert_eq!(formula_output.refit_plan, lrt_output.refit_plan);
    assert_eq!(
        formula_offset_output.refit_plan,
        explicit_offset_output.refit_plan
    );
    assert_eq!(
        formula_offset_output.original_fit.normalization_factors,
        Some(normalization_factors)
    );
    assert_lrt_fit_state_matches(
        &formula_offset_output.original_fit,
        &explicit_offset_output.original_fit,
        "formula LRT replacement offset original",
    );
    assert_eq!(
        formula_offset_output.original_results,
        explicit_offset_output.original_results
    );
    assert_eq!(
        formula_offset_output.refit_results,
        explicit_offset_output.refit_results
    );
    assert_lrt_fit_state_matches(
        formula_offset_output.refit_fit.as_ref().unwrap(),
        explicit_offset_output.refit_fit.as_ref().unwrap(),
        "formula LRT replacement offset refit",
    );
    assert_eq!(
        formula_offset_output.results,
        explicit_offset_output.results
    );
    assert_eq!(
        top_level_output.original_results,
        lrt_output.original_results
    );
    assert_eq!(formula_output.original_results, lrt_output.original_results);
    assert_lrt_fit_state_matches(
        &top_level_output.original_fit,
        &lrt_output.original_fit,
        "top-level LRT replacement original",
    );
    assert_lrt_fit_state_matches(
        &formula_output.original_fit,
        &lrt_output.original_fit,
        "formula LRT replacement original",
    );
    assert_eq!(top_level_output.refit_results, lrt_output.refit_results);
    assert_eq!(formula_output.refit_results, lrt_output.refit_results);
    assert_lrt_fit_state_matches(
        top_level_output.refit_fit.as_ref().unwrap(),
        lrt_output.refit_fit.as_ref().unwrap(),
        "top-level LRT replacement refit",
    );
    assert_lrt_fit_state_matches(
        formula_output.refit_fit.as_ref().unwrap(),
        lrt_output.refit_fit.as_ref().unwrap(),
        "formula LRT replacement refit",
    );
    assert_eq!(top_level_output.results, lrt_output.results);
    assert_eq!(formula_output.results, lrt_output.results);
    assert!(top_level_output.refit_plan.should_refit);
    assert!(formula_output.refit_plan.should_refit);
    assert!(top_level_output.refit_fit.as_ref().unwrap().lrt.is_some());
    assert!(formula_output.refit_fit.as_ref().unwrap().lrt.is_some());
    assert_eq!(
        top_level_output.results.metadata.test_type,
        Some(TestType::Lrt)
    );
    assert_eq!(
        top_level_output.results.metadata.result_name.as_deref(),
        Some("condition_B_vs_A")
    );
}

#[test]
fn top_level_fit_lrt_with_results_cooks_replacement_accepts_coefficient_name() {
    let counts = counts_with_zero_row();
    let full = full_design();
    let reduced = reduced_design();
    let levels = ["A", "A", "A", "A", "B", "B", "B", "B"]
        .into_iter()
        .map(String::from)
        .collect::<Vec<_>>();
    let model_frame = FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "condition".to_string(),
            sample_levels: levels,
            levels: None,
            reference: None,
        }],
        numeric_covariates: Vec::new(),
    };
    let builder = native_lrt_builder()
        .model_frame(model_frame)
        .automatic_dispersion_means();
    let replacement_options = CooksReplacementOptions {
        trim: 0.2,
        cooks_cutoff: 0.0,
        min_replicates: 3,
        which_samples: Some(vec![true, false, false, false, false, false, false, false]),
    };

    let named_output = builder
        .fit_lrt_with_results_name_with_cooks_replacement(
            &counts,
            &full,
            &reduced,
            "condition_B_vs_A",
            &replacement_options,
        )
        .unwrap();
    let indexed_output = builder
        .fit_lrt_glm_mu_with_cooks_replacement(&counts, &full, &reduced, 1, &replacement_options)
        .unwrap();
    let formula_named_output = builder
        .fit_lrt_formula_with_results_name_with_cooks_replacement(
            &counts,
            "~ condition",
            "~ 1",
            "condition_B_vs_A",
            &replacement_options,
        )
        .unwrap();

    assert_eq!(named_output.refit_plan, indexed_output.refit_plan);
    assert_eq!(formula_named_output.refit_plan, indexed_output.refit_plan);
    assert_eq!(named_output.results, indexed_output.results);
    assert_eq!(formula_named_output.results, indexed_output.results);
    assert_lrt_fit_state_matches(
        &named_output.original_fit,
        &indexed_output.original_fit,
        "named LRT replacement original",
    );
    assert_lrt_fit_state_matches(
        &formula_named_output.original_fit,
        &indexed_output.original_fit,
        "formula named LRT replacement original",
    );
}

#[test]
fn top_level_fit_lrt_with_results_cooks_replacement_validates_reduced_design() {
    let counts = counts_with_zero_row();
    let full = full_design();
    let invalid_reduced = full_design();
    let err = native_lrt_builder()
        .fit_lrt_with_results_with_cooks_replacement(
            &counts,
            &full,
            &invalid_reduced,
            &CooksReplacementOptions::new(f64::MAX),
        )
        .unwrap_err();

    assert!(matches!(err, DeseqError::InvalidDimensions { .. }));
}

#[test]
fn native_glm_mu_lrt_cooks_replacement_refit_skips_when_no_rows_are_marked() {
    let counts = counts_with_zero_row();
    let full = full_design();
    let reduced = reduced_design();
    let builder = native_lrt_builder();

    let (expected_fit, expected_results) = builder
        .clone()
        .fit_lrt_glm_mu(&counts, &full, &reduced, 1)
        .unwrap();
    let output = builder
        .fit_lrt_glm_mu_with_cooks_replacement(
            &counts,
            &full,
            &reduced,
            1,
            &CooksReplacementOptions {
                trim: 0.2,
                cooks_cutoff: f64::MAX,
                min_replicates: 3,
                which_samples: None,
            },
        )
        .unwrap();

    assert_eq!(output.refit_plan.n_refit, 0);
    assert!(!output.refit_plan.should_refit);
    assert!(output.refit_fit.is_none());
    assert!(output.refit_results.is_none());
    assert_lrt_fit_state_matches(
        &output.original_fit,
        &expected_fit,
        "native GLM-mu LRT no-refit original",
    );
    assert_eq!(output.original_results, expected_results);
    assert_eq!(output.results, expected_results);
    assert_eq!(output.results, output.original_results);
    assert_eq!(
        output.refit_plan.replacement.replaced_counts.as_slice(),
        counts.as_slice()
    );
}
