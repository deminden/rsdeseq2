use super::*;

#[test]
fn top_level_lrt_formula_list_contrast_cooks_replacement_accepts_r_cleaned_quoted_coefficient_alias()
 {
    let counts = counts_with_zero_row();
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
    let builder = native_lrt_builder().model_frame(model_frame);
    let options = CooksReplacementOptions::new(f64::MAX);

    let cleaned_output = builder
        .clone()
        .fit_lrt_formula_with_results_contrast_request_with_cooks_replacement(
            &counts,
            "~ `cell type`",
            "~ 1",
            &ResultsContrast::list(vec!["cell.type_B.cell_vs_T.cell".into()], Vec::new()),
            &options,
        )
        .unwrap();
    let raw_output = builder
        .fit_lrt_formula_with_results_contrast_request_with_cooks_replacement(
            &counts,
            "~ `cell type`",
            "~ 1",
            &ResultsContrast::list(vec!["cell type_B cell_vs_T cell".into()], Vec::new()),
            &options,
        )
        .unwrap();

    assert_lrt_fit_state_matches(
        &cleaned_output.original_fit,
        &raw_output.original_fit,
        "quoted formula R-cleaned list replacement LRT original",
    );
    assert_eq!(cleaned_output.refit_plan, raw_output.refit_plan);
    assert_eq!(cleaned_output.results.rows, raw_output.results.rows);
    assert_eq!(
        cleaned_output.results.metadata.comparison.as_deref(),
        Some("coefficient list contrast: cell.type_B.cell_vs_T.cell effect")
    );
    assert_eq!(
        raw_output.results.metadata.comparison.as_deref(),
        Some("coefficient list contrast: cell type_B cell_vs_T cell effect")
    );
}

#[test]
fn top_level_lrt_formula_list_contrast_cooks_replacement_accepts_cleaned_numeric_transform_alias() {
    let counts = counts_with_zero_row();
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
        numeric_covariates: Vec::new(),
    };
    let builder = native_lrt_builder().model_frame(model_frame);
    let options = CooksReplacementOptions::new(f64::MAX);

    let cleaned_output = builder
        .clone()
        .fit_lrt_formula_with_results_contrast_request_with_cooks_replacement(
            &counts,
            "~ as.double(`cell type`)",
            "~ 1",
            &ResultsContrast::list(vec!["cell.type_as.double".into()], Vec::new()),
            &options,
        )
        .unwrap();
    let raw_output = builder
        .fit_lrt_formula_with_results_contrast_request_with_cooks_replacement(
            &counts,
            "~ as.double(`cell type`)",
            "~ 1",
            &ResultsContrast::list(vec!["cell type_as.double".into()], Vec::new()),
            &options,
        )
        .unwrap();

    assert_lrt_fit_state_matches(
        &cleaned_output.original_fit,
        &raw_output.original_fit,
        "numeric-transform formula R-cleaned list replacement LRT original",
    );
    assert_eq!(cleaned_output.refit_plan, raw_output.refit_plan);
    assert_eq!(cleaned_output.results.rows, raw_output.results.rows);
    assert_eq!(
        cleaned_output
            .original_fit
            .current_model_frame()
            .unwrap()
            .resolved_numeric_covariate_by_alias("cell.type_as.double")
            .unwrap()
            .unwrap()
            .name,
        "cell type_as.double"
    );
    assert_eq!(
        cleaned_output.results.metadata.comparison.as_deref(),
        Some("coefficient list contrast: cell.type_as.double effect")
    );
    assert_eq!(
        raw_output.results.metadata.comparison.as_deref(),
        Some("coefficient list contrast: cell type_as.double effect")
    );
}

#[test]
fn top_level_lrt_formula_numeric_contrast_cooks_replacement_retains_numeric_transform_metadata() {
    let counts = counts_with_zero_row();
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
        numeric_covariates: Vec::new(),
    };
    let builder = native_lrt_builder().model_frame(model_frame);
    let options = CooksReplacementOptions::new(f64::MAX);

    let numeric_output = builder
        .fit_lrt_formula_with_results_contrast_request_with_cooks_replacement(
            &counts,
            "~ as.double(`cell type`)",
            "~ 1",
            &ResultsContrast::numeric(vec![0.0, 1.0]),
            &options,
        )
        .unwrap();

    assert_eq!(numeric_output.refit_plan.n_refit, 0);
    assert_eq!(
        numeric_output
            .original_fit
            .current_model_frame()
            .unwrap()
            .resolved_numeric_covariate_by_alias("cell.type_as.double")
            .unwrap()
            .unwrap()
            .name,
        "cell type_as.double"
    );
    assert_eq!(
        numeric_output.results.metadata.comparison.as_deref(),
        Some("primitive numeric contrast")
    );
}

#[test]
fn top_level_lrt_formula_list_contrast_cooks_replacement_accepts_cleaned_positive_and_negative_aliases()
 {
    let counts = counts_with_zero_row();
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
    let builder = native_lrt_builder().model_frame(model_frame);
    let options = CooksReplacementOptions::new(f64::MAX);

    let cleaned_output = builder
        .clone()
        .fit_lrt_formula_with_results_contrast_request_with_cooks_replacement(
            &counts,
            "~ `cell type`",
            "~ 1",
            &ResultsContrast::list(
                vec!["cell.type_NK.cell_vs_T.cell".into()],
                vec!["cell.type_B.cell_vs_T.cell".into()],
            ),
            &options,
        )
        .unwrap();
    let raw_output = builder
        .fit_lrt_formula_with_results_contrast_request_with_cooks_replacement(
            &counts,
            "~ `cell type`",
            "~ 1",
            &ResultsContrast::list(
                vec!["cell type_NK cell_vs_T cell".into()],
                vec!["cell type_B cell_vs_T cell".into()],
            ),
            &options,
        )
        .unwrap();

    assert_lrt_fit_state_matches(
        &cleaned_output.original_fit,
        &raw_output.original_fit,
        "three-level formula R-cleaned list replacement LRT original",
    );
    assert_eq!(cleaned_output.refit_plan, raw_output.refit_plan);
    assert_eq!(cleaned_output.results.rows, raw_output.results.rows);
    assert_eq!(
        cleaned_output.results.metadata.comparison.as_deref(),
        Some(
            "coefficient list contrast: cell.type_NK.cell_vs_T.cell vs cell.type_B.cell_vs_T.cell"
        )
    );
    assert_eq!(
        raw_output.results.metadata.comparison.as_deref(),
        Some(
            "coefficient list contrast: cell type_NK cell_vs_T cell vs cell type_B cell_vs_T cell"
        )
    );
}
