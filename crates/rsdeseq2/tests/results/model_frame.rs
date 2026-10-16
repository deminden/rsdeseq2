use super::*;

#[test]
fn fit_expanded_formula_model_frame_beta_prior_wald_normalized_matches_formula_workflow() {
    let counts =
        CountMatrix::from_row_major_u32(2, 4, vec![10, 12, 120, 24, 30, 33, 45, 54]).unwrap();
    let condition = vec![
        "A".to_string(),
        "A".to_string(),
        "B".to_string(),
        "B".to_string(),
    ];
    let batch = vec![
        "X".to_string(),
        "Y".to_string(),
        "X".to_string(),
        "Y".to_string(),
    ];
    let exposure_offset = [0.0_f64, 0.1, 0.2, 0.3];
    let factors = [
        ExpandedFactorSpec {
            factor: "condition",
            sample_levels: &condition,
            reference: "A",
            levels: None,
        },
        ExpandedFactorSpec {
            factor: "batch",
            sample_levels: &batch,
            reference: "X",
            levels: None,
        },
    ];
    let numeric_covariates = [ExpandedNumericSpec {
        name: "exposure_offset",
        values: &exposure_offset,
    }];
    let model_frame = FormulaModelFrame {
        factors: vec![
            FormulaFactorColumn {
                name: "condition".to_string(),
                sample_levels: condition.clone(),
                levels: None,
                reference: None,
            },
            FormulaFactorColumn {
                name: "batch".to_string(),
                sample_levels: batch.clone(),
                levels: None,
                reference: None,
            },
        ],
        numeric_covariates: vec![FormulaNumericColumn {
            name: "exposure_offset".to_string(),
            values: exposure_offset.to_vec(),
        }],
    };
    let normalization_factors = RowMajorMatrix::from_row_major(
        2,
        4,
        vec![
            1.0, 1.0, 1.0, 1.0, //
            1.0, 1.0, 1.0, 1.0,
        ],
    )
    .unwrap();
    let weights = RowMajorMatrix::from_row_major(
        2,
        4,
        vec![
            1.0, 0.9, 1.0, 0.8, //
            1.0, 1.0, 0.95, 0.9,
        ],
    )
    .unwrap();
    let dispersions = [0.05, 0.08];
    let base_mean = [41.5, 40.5];
    let disp_fit = [0.05, 0.08];
    let names = vec!["gene_a".to_string(), "gene_b".to_string()];
    let options = BetaPriorRefitOptions {
        fit_options: IrlsOptions::default(),
        variance_options: BetaPriorVarianceOptions {
            method: BetaPriorVarianceMethod::Quantile,
            upper_quantile: 0.5,
            ..BetaPriorVarianceOptions::default()
        },
    };
    let formula = "~ condition * batch + offset(exposure_offset)";

    let borrowed =
        fit_expanded_formula_beta_prior_wald_results_with_normalization_factors_and_weights(
            ExpandedFormulaBetaPriorWaldNormalizedResultsInput {
                counts: &counts,
                formula,
                factors: &factors,
                numeric_covariates: &numeric_covariates,
                normalization_factors: &normalization_factors,
                weights: Some(&weights),
                dispersions: &dispersions,
                base_mean: &base_mean,
                disp_fit: &disp_fit,
                gene_names: Some(&names),
                options: options.clone(),
            },
            3,
        )
        .unwrap();
    let model_frame_result =
        fit_expanded_formula_model_frame_beta_prior_wald_results_with_normalization_factors_and_weights(
            ExpandedFormulaModelFrameBetaPriorWaldNormalizedResultsInput {
                counts: &counts,
                formula,
                model_frame: &model_frame,
                normalization_factors: &normalization_factors,
                weights: Some(&weights),
                dispersions: &dispersions,
                base_mean: &base_mean,
                disp_fit: &disp_fit,
                gene_names: Some(&names),
                options: options.clone(),
            },
            3,
        )
        .unwrap();
    assert_eq!(model_frame_result, borrowed);

    let contrast = [0.0, 0.0, 0.0, 1.0];
    let borrowed_contrast =
        fit_expanded_formula_beta_prior_wald_contrast_results_with_normalization_factors_and_weights(
            ExpandedFormulaBetaPriorWaldNormalizedResultsInput {
                counts: &counts,
                formula,
                factors: &factors,
                numeric_covariates: &numeric_covariates,
                normalization_factors: &normalization_factors,
                weights: Some(&weights),
                dispersions: &dispersions,
                base_mean: &base_mean,
                disp_fit: &disp_fit,
                gene_names: Some(&names),
                options: options.clone(),
            },
            &contrast,
        )
        .unwrap();
    let model_frame_contrast =
        fit_expanded_formula_model_frame_beta_prior_wald_contrast_results_with_normalization_factors_and_weights(
            ExpandedFormulaModelFrameBetaPriorWaldNormalizedResultsInput {
                counts: &counts,
                formula,
                model_frame: &model_frame,
                normalization_factors: &normalization_factors,
                weights: Some(&weights),
                dispersions: &dispersions,
                base_mean: &base_mean,
                disp_fit: &disp_fit,
                gene_names: Some(&names),
                options: options.clone(),
            },
            &contrast,
        )
        .unwrap();
    assert_eq!(model_frame_contrast, borrowed_contrast);

    let replacement_options = CooksReplacementOptions {
        trim: 0.2,
        cooks_cutoff: 0.0,
        min_replicates: 3,
        which_samples: Some(vec![false, false, true, false]),
    };
    let borrowed_replacement =
        fit_expanded_formula_beta_prior_wald_results_with_normalization_factors_and_weights_and_cooks_replacement(
            ExpandedFormulaBetaPriorWaldNormalizedResultsInput {
                counts: &counts,
                formula,
                factors: &factors,
                numeric_covariates: &numeric_covariates,
                normalization_factors: &normalization_factors,
                weights: Some(&weights),
                dispersions: &dispersions,
                base_mean: &base_mean,
                disp_fit: &disp_fit,
                gene_names: Some(&names),
                options: options.clone(),
            },
            3,
            &replacement_options,
        )
        .unwrap();
    let model_frame_replacement =
        fit_expanded_formula_model_frame_beta_prior_wald_results_with_normalization_factors_and_weights_and_cooks_replacement(
            ExpandedFormulaModelFrameBetaPriorWaldNormalizedResultsInput {
                counts: &counts,
                formula,
                model_frame: &model_frame,
                normalization_factors: &normalization_factors,
                weights: Some(&weights),
                dispersions: &dispersions,
                base_mean: &base_mean,
                disp_fit: &disp_fit,
                gene_names: Some(&names),
                options: options.clone(),
            },
            3,
            &replacement_options,
        )
        .unwrap();
    assert_eq!(model_frame_replacement, borrowed_replacement);

    let borrowed_contrast_replacement =
        fit_expanded_formula_beta_prior_wald_contrast_results_with_normalization_factors_and_weights_and_cooks_replacement(
            ExpandedFormulaBetaPriorWaldNormalizedResultsInput {
                counts: &counts,
                formula,
                factors: &factors,
                numeric_covariates: &numeric_covariates,
                normalization_factors: &normalization_factors,
                weights: Some(&weights),
                dispersions: &dispersions,
                base_mean: &base_mean,
                disp_fit: &disp_fit,
                gene_names: Some(&names),
                options: options.clone(),
            },
            &contrast,
            &replacement_options,
        )
        .unwrap();
    let model_frame_contrast_replacement =
        fit_expanded_formula_model_frame_beta_prior_wald_contrast_results_with_normalization_factors_and_weights_and_cooks_replacement(
            ExpandedFormulaModelFrameBetaPriorWaldNormalizedResultsInput {
                counts: &counts,
                formula,
                model_frame: &model_frame,
                normalization_factors: &normalization_factors,
                weights: Some(&weights),
                dispersions: &dispersions,
                base_mean: &base_mean,
                disp_fit: &disp_fit,
                gene_names: Some(&names),
                options,
            },
            &contrast,
            &replacement_options,
        )
        .unwrap();
    assert_eq!(
        model_frame_contrast_replacement,
        borrowed_contrast_replacement
    );
}
