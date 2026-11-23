use super::*;
use crate::support::run::{run_fit, run_results, run_wald_replacement};
use rsdeseq2::prelude::{Analysis, MeanModel, ResultColumn};

#[test]
fn top_level_fit_with_results_cooks_replacement_runs_default_automatic_wald() {
    let counts = native_wald_counts_with_zero_row();
    let design = two_group_design();
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
    let builder = glm_mu_native_wald_builder()
        .model_frame(model_frame)
        .automatic_dispersion_means();
    let replacement_options = CooksReplacementOptions {
        trim: 0.2,
        cooks_cutoff: 0.0,
        min_replicates: 3,
        which_samples: Some(vec![true, false, false, false, false, false, false, false]),
    };

    let top_level_output = run_wald_replacement(
        builder.run(&Analysis::new(&counts, &design).cooks_replacement(&replacement_options)),
    )
    .unwrap();
    let wald_output = run_wald_replacement(
        builder.run(
            &Analysis::new(&counts, &design)
                .means(MeanModel::GlmMu)
                .result(ResultColumn::Coefficient(1))
                .cooks_replacement(&replacement_options),
        ),
    )
    .unwrap();
    let formula_output =
        run_wald_replacement(builder.run(
            &Analysis::formula(&counts, "~ condition").cooks_replacement(&replacement_options),
        ))
        .unwrap();
    let transformed_formula_output = run_wald_replacement(
        builder.run(
            &Analysis::formula(&counts, "~ factor(condition, levels=c('A','B'))")
                .cooks_replacement(&replacement_options),
        ),
    )
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
    let explicit_offset_output = run_wald_replacement(
        builder
            .clone()
            .normalization_factors(normalization_factors.clone())
            .run(
                &Analysis::formula(&counts, "~ condition").cooks_replacement(&replacement_options),
            ),
    )
    .unwrap();
    let formula_offset_output = run_wald_replacement(
        builder.run(
            &Analysis::formula(&counts, "~ condition + offset(exposure)")
                .cooks_replacement(&replacement_options),
        ),
    )
    .unwrap();

    assert_eq!(top_level_output.refit_plan, wald_output.refit_plan);
    assert_eq!(formula_output.refit_plan, wald_output.refit_plan);
    assert_eq!(
        transformed_formula_output.refit_plan,
        wald_output.refit_plan
    );
    assert_eq!(
        formula_offset_output.refit_plan,
        explicit_offset_output.refit_plan
    );
    assert_eq!(
        formula_offset_output.original_fit.normalization_factors,
        Some(normalization_factors.clone())
    );
    assert_wald_fit_state_matches(
        &formula_offset_output.original_fit,
        &explicit_offset_output.original_fit,
        "formula top-level Wald replacement offset original",
    );
    assert_eq!(
        formula_offset_output.original_results,
        explicit_offset_output.original_results
    );
    assert_eq!(
        formula_offset_output.refit_results,
        explicit_offset_output.refit_results
    );
    assert_wald_fit_state_matches(
        formula_offset_output.refit_fit.as_ref().unwrap(),
        explicit_offset_output.refit_fit.as_ref().unwrap(),
        "formula top-level Wald replacement offset refit",
    );
    assert_eq!(
        formula_offset_output.results,
        explicit_offset_output.results
    );
    assert_eq!(
        top_level_output.original_results,
        wald_output.original_results
    );
    assert_eq!(
        formula_output.original_results,
        wald_output.original_results
    );
    assert_wald_fit_state_matches(
        &top_level_output.original_fit,
        &wald_output.original_fit,
        "top-level Wald replacement original",
    );
    assert_wald_fit_state_matches(
        &formula_output.original_fit,
        &wald_output.original_fit,
        "formula top-level Wald replacement original",
    );
    assert_eq!(top_level_output.refit_results, wald_output.refit_results);
    assert_eq!(formula_output.refit_results, wald_output.refit_results);
    assert_wald_fit_state_matches(
        top_level_output.refit_fit.as_ref().unwrap(),
        wald_output.refit_fit.as_ref().unwrap(),
        "top-level Wald replacement refit",
    );
    assert_wald_fit_state_matches(
        formula_output.refit_fit.as_ref().unwrap(),
        wald_output.refit_fit.as_ref().unwrap(),
        "formula top-level Wald replacement refit",
    );
    let transformed_model_frame = transformed_formula_output
        .original_fit
        .current_model_frame()
        .unwrap();
    assert_eq!(
        transformed_model_frame
            .factors
            .iter()
            .map(|factor| factor.name.as_str())
            .collect::<Vec<_>>(),
        vec!["condition", "factor(condition, levels = c(\"A\", \"B\"))"]
    );
    assert_eq!(
        transformed_formula_output
            .refit_fit
            .as_ref()
            .and_then(DeseqFit::current_model_frame),
        Some(transformed_model_frame)
    );
    assert_eq!(top_level_output.results, wald_output.results);
    assert_eq!(formula_output.results, wald_output.results);
    assert_eq!(
        transformed_formula_output.results.rows,
        wald_output.results.rows
    );
    assert!(top_level_output.refit_plan.should_refit);
    assert_eq!(
        top_level_output.results.metadata.result_name.as_deref(),
        Some("condition_B_vs_A")
    );
}

#[test]
fn top_level_fit_with_results_cooks_replacement_accepts_coefficient_name() {
    let counts = native_wald_counts_with_zero_row();
    let design = two_group_design();
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
    let builder = glm_mu_native_wald_builder()
        .model_frame(model_frame)
        .automatic_dispersion_means();
    let replacement_options = CooksReplacementOptions {
        trim: 0.2,
        cooks_cutoff: 0.0,
        min_replicates: 3,
        which_samples: Some(vec![true, false, false, false, false, false, false, false]),
    };

    let named_output = run_wald_replacement(
        builder.run(
            &Analysis::new(&counts, &design)
                .result(ResultColumn::Name("condition_B_vs_A"))
                .cooks_replacement(&replacement_options),
        ),
    )
    .unwrap();
    let indexed_output = run_wald_replacement(
        builder.run(
            &Analysis::new(&counts, &design)
                .means(MeanModel::GlmMu)
                .result(ResultColumn::Coefficient(1))
                .cooks_replacement(&replacement_options),
        ),
    )
    .unwrap();
    let formula_named_output = run_wald_replacement(
        builder.run(
            &Analysis::formula(&counts, "~ condition")
                .result(ResultColumn::Name("condition_B_vs_A"))
                .cooks_replacement(&replacement_options),
        ),
    )
    .unwrap();

    assert_eq!(named_output.refit_plan, indexed_output.refit_plan);
    assert_eq!(formula_named_output.refit_plan, indexed_output.refit_plan);
    assert_eq!(named_output.results, indexed_output.results);
    assert_eq!(formula_named_output.results, indexed_output.results);
    assert_wald_fit_state_matches(
        &named_output.original_fit,
        &indexed_output.original_fit,
        "named top-level Wald replacement original",
    );
    assert_wald_fit_state_matches(
        &formula_named_output.original_fit,
        &indexed_output.original_fit,
        "formula named top-level Wald replacement original",
    );
}

#[test]
fn top_level_fit_with_results_contrast_runs_default_automatic_wald_contrast() {
    let counts = native_wald_counts_with_zero_row();
    let design = two_group_design();
    let builder = glm_mu_native_wald_builder().automatic_dispersion_means();
    let contrast = [0.0, 1.0];

    let (top_level_fit, top_level_results) = run_results(
        builder.run(&Analysis::new(&counts, &design).result(ResultColumn::Numeric(&contrast))),
    )
    .unwrap();
    let (contrast_fit, contrast_results) = run_results(
        builder.run(
            &Analysis::new(&counts, &design)
                .means(MeanModel::GlmMu)
                .result(ResultColumn::Numeric(&contrast)),
        ),
    )
    .unwrap();

    assert_eq!(top_level_fit.wald, contrast_fit.wald);
    assert_wald_fit_state_matches(
        &top_level_fit,
        &contrast_fit,
        "top-level Wald contrast results",
    );
    assert_eq!(top_level_results, contrast_results);
    assert_eq!(top_level_results.metadata.test_type, Some(TestType::Wald));
    assert_eq!(
        top_level_results.metadata.result_name.as_deref(),
        Some("contrast")
    );
}

#[test]
fn top_level_fit_contrast_returns_default_automatic_wald_contrast_fit() {
    let counts = native_wald_counts_with_zero_row();
    let design = two_group_design();
    let builder = glm_mu_native_wald_builder().automatic_dispersion_means();
    let contrast = [0.0, 1.0];

    let top_level_fit = run_fit(
        builder.run(&Analysis::new(&counts, &design).result(ResultColumn::Numeric(&contrast))),
    )
    .unwrap();
    let (contrast_fit, _contrast_results) = run_results(
        builder.run(
            &Analysis::new(&counts, &design)
                .means(MeanModel::GlmMu)
                .result(ResultColumn::Numeric(&contrast)),
        ),
    )
    .unwrap();

    assert_eq!(top_level_fit.counts_summary, contrast_fit.counts_summary);
    assert_eq!(top_level_fit.design, contrast_fit.design);
    assert_eq!(top_level_fit.all_zero, contrast_fit.all_zero);
    assert_eq!(top_level_fit.beta_converged, contrast_fit.beta_converged);
    assert_wald_fit_state_matches(&top_level_fit, &contrast_fit, "top-level Wald contrast");
    assert_eq!(top_level_fit.wald, contrast_fit.wald);
    assert!(top_level_fit.dispersion.is_some());
    assert!(top_level_fit.beta.is_some());
    assert!(top_level_fit.wald.is_some());
}

#[test]
fn top_level_fit_with_results_contrast_cooks_replacement_runs_default_automatic_wald_contrast() {
    let counts = native_wald_counts_with_zero_row();
    let design = two_group_design();
    let builder = glm_mu_native_wald_builder().automatic_dispersion_means();
    let contrast = [0.0, 1.0];
    let replacement_options = CooksReplacementOptions {
        trim: 0.2,
        cooks_cutoff: 0.0,
        min_replicates: 3,
        which_samples: Some(vec![true, false, false, false, false, false, false, false]),
    };

    let top_level_output = run_wald_replacement(
        builder.run(
            &Analysis::new(&counts, &design)
                .result(ResultColumn::Numeric(&contrast))
                .cooks_replacement(&replacement_options),
        ),
    )
    .unwrap();
    let contrast_output = run_wald_replacement(
        builder.run(
            &Analysis::new(&counts, &design)
                .means(MeanModel::GlmMu)
                .result(ResultColumn::Numeric(&contrast))
                .cooks_replacement(&replacement_options),
        ),
    )
    .unwrap();

    assert_eq!(top_level_output.refit_plan, contrast_output.refit_plan);
    assert_eq!(
        top_level_output.original_results,
        contrast_output.original_results
    );
    assert_wald_fit_state_matches(
        &top_level_output.original_fit,
        &contrast_output.original_fit,
        "top-level Wald contrast replacement original",
    );
    assert_eq!(
        top_level_output.refit_results,
        contrast_output.refit_results
    );
    assert_wald_fit_state_matches(
        top_level_output.refit_fit.as_ref().unwrap(),
        contrast_output.refit_fit.as_ref().unwrap(),
        "top-level Wald contrast replacement refit",
    );
    assert_eq!(top_level_output.results, contrast_output.results);
    assert!(top_level_output.refit_plan.should_refit);
    assert_eq!(
        top_level_output.results.metadata.result_name.as_deref(),
        Some("contrast")
    );
}

#[test]
fn top_level_fit_with_results_contrast_specs_set_metadata() {
    let counts = native_wald_counts_with_zero_row();
    let design = two_group_design();
    let builder = glm_mu_native_wald_builder();
    let levels = ["A", "A", "A", "A", "B", "B", "B", "B"]
        .into_iter()
        .map(String::from)
        .collect::<Vec<_>>();

    let (primitive_fit, _primitive_results) = run_results(
        builder
            .clone()
            .run(&Analysis::new(&counts, &design).result(ResultColumn::Numeric(&[0.0, 1.0]))),
    )
    .unwrap();
    let (named_fit, named_results) =
        run_results(builder.clone().run(&Analysis::new(&counts, &design).result(
            ResultColumn::Spec(&ContrastSpec::coefficient_name("condition_B_vs_A")),
        )))
        .unwrap();
    let (factor_fit, factor_results) =
        run_results(builder.run(&Analysis::new(&counts, &design).result(
            ResultColumn::FactorLevel(FactorLevelContrast::new("condition", "B", "A", &levels)),
        )))
        .unwrap();

    assert_wald_fit_state_matches(
        &named_fit,
        &primitive_fit,
        "top-level named Wald contrast results",
    );
    assert_wald_fit_state_matches(
        &factor_fit,
        &primitive_fit,
        "top-level factor Wald contrast results",
    );
    assert_eq!(
        named_results.metadata.result_name.as_deref(),
        Some("condition_B_vs_A")
    );
    assert_eq!(
        named_results.metadata.comparison.as_deref(),
        Some("coefficient condition_B_vs_A")
    );
    assert_eq!(
        factor_results.metadata.result_name.as_deref(),
        Some("condition_B_vs_A")
    );
    assert_eq!(
        factor_results.metadata.comparison.as_deref(),
        Some("factor-level contrast: condition B vs A")
    );
}
