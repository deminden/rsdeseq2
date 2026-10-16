use super::*;

#[test]
fn fit_expanded_formula_beta_prior_wald_results_matches_additive_workflow() {
    let counts =
        CountMatrix::from_row_major_u32(2, 4, vec![10, 12, 20, 24, 30, 33, 45, 54]).unwrap();
    let condition = vec![
        "A".to_string(),
        "A".to_string(),
        "B".to_string(),
        "B".to_string(),
    ];
    let dose = [0.0, 1.0, 0.0, 1.0];
    let time = [1.0, 1.0, 2.0, 2.0];
    let factors = [ExpandedFactorSpec {
        factor: "condition",
        sample_levels: &condition,
        reference: "A",
        levels: None,
    }];
    let numeric_covariates = [
        ExpandedNumericSpec {
            name: "dose",
            values: &dose,
        },
        ExpandedNumericSpec {
            name: "time",
            values: &time,
        },
    ];
    let factor_numeric = [ExpandedFactorNumericInteractionSpec {
        factor: "condition",
        numeric: "dose",
    }];
    let numeric_interactions = [ExpandedNumericInteractionSpec {
        left_numeric: "dose",
        right_numeric: "time",
    }];
    let size_factors = [1.0, 1.0, 1.0, 1.0];
    let dispersions = [0.05, 0.08];
    let base_mean = [16.5, 40.5];
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

    let formula = fit_expanded_formula_beta_prior_wald_results(
        ExpandedFormulaBetaPriorWaldResultsInput {
            counts: &counts,
            formula: "~ condition + dose + time + condition:dose + dose:time",
            factors: &factors,
            numeric_covariates: &numeric_covariates,
            size_factors: &size_factors,
            weights: None,
            dispersions: &dispersions,
            base_mean: &base_mean,
            disp_fit: &disp_fit,
            gene_names: Some(&names),
            options: options.clone(),
        },
        4,
    )
    .unwrap();

    let additive = fit_expanded_additive_beta_prior_wald_results(
        ExpandedAdditiveBetaPriorWaldResultsInput {
            counts: &counts,
            factors: &factors,
            numeric_covariates: &numeric_covariates,
            interactions: &[],
            factor_numeric_interactions: &factor_numeric,
            numeric_interactions: &numeric_interactions,
            size_factors: &size_factors,
            weights: None,
            dispersions: &dispersions,
            base_mean: &base_mean,
            disp_fit: &disp_fit,
            gene_names: Some(&names),
            options: options.clone(),
        },
        4,
    )
    .unwrap();

    assert_eq!(formula.design, additive.design);
    assert_eq!(formula.fit, additive.fit);
    assert_eq!(formula.results, additive.results);
    assert_eq!(
        formula.results.metadata.result_name.as_deref(),
        Some("condition_B_vs_A:dose")
    );

    let formula_contrast = fit_expanded_formula_beta_prior_wald_contrast_results(
        ExpandedFormulaBetaPriorWaldResultsInput {
            counts: &counts,
            formula: "~ condition + dose + time + condition:dose + dose:time",
            factors: &factors,
            numeric_covariates: &numeric_covariates,
            size_factors: &size_factors,
            weights: None,
            dispersions: &dispersions,
            base_mean: &base_mean,
            disp_fit: &disp_fit,
            gene_names: Some(&names),
            options,
        },
        &[0.0, 0.0, 0.0, 0.0, 1.0, 0.0],
    )
    .unwrap();
    let direct_contrast_results = build_wald_contrast_results_from_expanded_beta_prior_fit(
        &base_mean,
        &formula.fit,
        &[0.0, 0.0, 0.0, 0.0, 1.0, 0.0],
        Some(&names),
        Some(&dispersions),
    )
    .unwrap();

    assert_eq!(formula_contrast.fit, formula.fit);
    assert_eq!(formula_contrast.results, direct_contrast_results);
}

#[test]
fn fit_expanded_formula_beta_prior_wald_replacement_matches_additive_workflow() {
    let counts =
        CountMatrix::from_row_major_u32(2, 4, vec![10, 12, 120, 24, 30, 33, 45, 54]).unwrap();
    let condition = vec![
        "A".to_string(),
        "A".to_string(),
        "B".to_string(),
        "B".to_string(),
    ];
    let dose = [0.0, 1.0, 0.0, 1.0];
    let time = [1.0, 1.0, 2.0, 2.0];
    let factors = [ExpandedFactorSpec {
        factor: "condition",
        sample_levels: &condition,
        reference: "A",
        levels: None,
    }];
    let numeric_covariates = [
        ExpandedNumericSpec {
            name: "dose",
            values: &dose,
        },
        ExpandedNumericSpec {
            name: "time",
            values: &time,
        },
    ];
    let factor_numeric = [ExpandedFactorNumericInteractionSpec {
        factor: "condition",
        numeric: "dose",
    }];
    let numeric_interactions = [ExpandedNumericInteractionSpec {
        left_numeric: "dose",
        right_numeric: "time",
    }];
    let size_factors = [1.0, 1.0, 1.0, 1.0];
    let dispersions = [0.05, 0.08];
    let base_mean = [41.5, 40.5];
    let disp_fit = [0.05, 0.08];
    let options = BetaPriorRefitOptions {
        fit_options: IrlsOptions::default(),
        variance_options: BetaPriorVarianceOptions {
            method: BetaPriorVarianceMethod::Quantile,
            upper_quantile: 0.5,
            ..BetaPriorVarianceOptions::default()
        },
    };
    let replacement_options = CooksReplacementOptions {
        trim: 0.2,
        cooks_cutoff: 0.0,
        min_replicates: 3,
        which_samples: Some(vec![false, false, true, false]),
    };

    let formula = fit_expanded_formula_beta_prior_wald_results_with_cooks_replacement(
        ExpandedFormulaBetaPriorWaldResultsInput {
            counts: &counts,
            formula: "~ condition + dose + time + condition:dose + dose:time",
            factors: &factors,
            numeric_covariates: &numeric_covariates,
            size_factors: &size_factors,
            weights: None,
            dispersions: &dispersions,
            base_mean: &base_mean,
            disp_fit: &disp_fit,
            gene_names: None,
            options: options.clone(),
        },
        4,
        &replacement_options,
    )
    .unwrap();

    let additive = fit_expanded_additive_beta_prior_wald_results_with_cooks_replacement(
        ExpandedAdditiveBetaPriorWaldResultsInput {
            counts: &counts,
            factors: &factors,
            numeric_covariates: &numeric_covariates,
            interactions: &[],
            factor_numeric_interactions: &factor_numeric,
            numeric_interactions: &numeric_interactions,
            size_factors: &size_factors,
            weights: None,
            dispersions: &dispersions,
            base_mean: &base_mean,
            disp_fit: &disp_fit,
            gene_names: None,
            options: options.clone(),
        },
        4,
        &replacement_options,
    )
    .unwrap();

    assert_eq!(formula.design, additive.design);
    assert_eq!(formula.replacement, additive.replacement);

    let formula_contrast =
        fit_expanded_formula_beta_prior_wald_contrast_results_with_cooks_replacement(
            ExpandedFormulaBetaPriorWaldResultsInput {
                counts: &counts,
                formula: "~ condition + dose + time + condition:dose + dose:time",
                factors: &factors,
                numeric_covariates: &numeric_covariates,
                size_factors: &size_factors,
                weights: None,
                dispersions: &dispersions,
                base_mean: &base_mean,
                disp_fit: &disp_fit,
                gene_names: None,
                options,
            },
            &[0.0, 0.0, 0.0, 0.0, 1.0, 0.0],
            &replacement_options,
        )
        .unwrap();
    assert_eq!(
        formula_contrast
            .replacement
            .refit_plan
            .replacement
            .replaced_counts,
        counts
    );
    assert_eq!(
        formula_contrast
            .replacement
            .results
            .metadata
            .result_name
            .as_deref(),
        Some("contrast")
    );
}

#[test]
fn fit_expanded_formula_beta_prior_wald_results_applies_formula_offsets() {
    let counts =
        CountMatrix::from_row_major_u32(2, 4, vec![10, 12, 20, 24, 30, 33, 45, 54]).unwrap();
    let condition = vec![
        "A".to_string(),
        "A".to_string(),
        "B".to_string(),
        "B".to_string(),
    ];
    let dose = [0.0_f64, 0.05, 0.1, 0.15];
    let exposure_offset = [0.0_f64, 0.1, 0.2, 0.3];
    let factors = [ExpandedFactorSpec {
        factor: "condition",
        sample_levels: &condition,
        reference: "A",
        levels: None,
    }];
    let numeric_covariates = [
        ExpandedNumericSpec {
            name: "dose",
            values: &dose,
        },
        ExpandedNumericSpec {
            name: "exposure_offset",
            values: &exposure_offset,
        },
    ];
    let size_factors = [1.0, 1.1, 0.9, 1.2];
    let dispersions = [0.05, 0.08];
    let base_mean = [16.5, 40.5];
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

    let formula = fit_expanded_formula_beta_prior_wald_results(
        ExpandedFormulaBetaPriorWaldResultsInput {
            counts: &counts,
            formula: "~ condition + offset(I(exposure_offset + dose))",
            factors: &factors,
            numeric_covariates: &numeric_covariates,
            size_factors: &size_factors,
            weights: None,
            dispersions: &dispersions,
            base_mean: &base_mean,
            disp_fit: &disp_fit,
            gene_names: Some(&names),
            options: options.clone(),
        },
        1,
    )
    .unwrap();

    let design = expanded_formula_design("~ condition", &factors, &[]).unwrap();
    let mut normalization_values = Vec::new();
    for _ in 0..counts.n_genes() {
        for (sample, size_factor) in size_factors.iter().copied().enumerate() {
            normalization_values.push(size_factor * (exposure_offset[sample] + dose[sample]).exp());
        }
    }
    let normalization_factors =
        RowMajorMatrix::from_row_major(counts.n_genes(), counts.n_samples(), normalization_values)
            .unwrap();
    let direct = fit_expanded_beta_prior_wald_results_with_normalization_factors_and_weights(
        ExpandedBetaPriorWaldNormalizedResultsInput {
            counts: &counts,
            design: ExpandedModelBetaPriorDesignInput {
                expanded_design: &design.expanded_design,
                standard_design: &design.standard_design,
                coefficient_groups: &design.coefficient_groups,
            },
            normalization_factors: &normalization_factors,
            weights: None,
            dispersions: &dispersions,
            base_mean: &base_mean,
            disp_fit: &disp_fit,
            gene_names: Some(&names),
            options: options.clone(),
        },
        1,
    )
    .unwrap();

    assert_eq!(formula.design, design);
    assert_eq!(formula.fit, direct.fit);
    assert_eq!(formula.results, direct.results);

    let formula_contrast = fit_expanded_formula_beta_prior_wald_contrast_results(
        ExpandedFormulaBetaPriorWaldResultsInput {
            counts: &counts,
            formula: "~ condition + offset(I(exposure_offset + dose))",
            factors: &factors,
            numeric_covariates: &numeric_covariates,
            size_factors: &size_factors,
            weights: None,
            dispersions: &dispersions,
            base_mean: &base_mean,
            disp_fit: &disp_fit,
            gene_names: Some(&names),
            options,
        },
        &[0.0, 1.0],
    )
    .unwrap();
    let direct_contrast =
        fit_expanded_beta_prior_wald_contrast_results_with_normalization_factors_and_weights(
            ExpandedBetaPriorWaldNormalizedResultsInput {
                counts: &counts,
                design: ExpandedModelBetaPriorDesignInput {
                    expanded_design: &design.expanded_design,
                    standard_design: &design.standard_design,
                    coefficient_groups: &design.coefficient_groups,
                },
                normalization_factors: &normalization_factors,
                weights: None,
                dispersions: &dispersions,
                base_mean: &base_mean,
                disp_fit: &disp_fit,
                gene_names: Some(&names),
                options: BetaPriorRefitOptions {
                    fit_options: IrlsOptions::default(),
                    variance_options: BetaPriorVarianceOptions {
                        method: BetaPriorVarianceMethod::Quantile,
                        upper_quantile: 0.5,
                        ..BetaPriorVarianceOptions::default()
                    },
                },
            },
            &[0.0, 1.0],
        )
        .unwrap();

    assert_eq!(formula_contrast.fit, direct_contrast.fit);
    assert_eq!(formula_contrast.results, direct_contrast.results);
}

#[test]
fn fit_expanded_formula_beta_prior_wald_replacement_applies_formula_offsets() {
    let counts =
        CountMatrix::from_row_major_u32(2, 4, vec![10, 12, 120, 24, 30, 33, 45, 54]).unwrap();
    let condition = vec![
        "A".to_string(),
        "A".to_string(),
        "B".to_string(),
        "B".to_string(),
    ];
    let dose = [0.0_f64, 0.05, 0.1, 0.15];
    let exposure_offset = [0.0_f64, 0.1, 0.2, 0.3];
    let factors = [ExpandedFactorSpec {
        factor: "condition",
        sample_levels: &condition,
        reference: "A",
        levels: None,
    }];
    let numeric_covariates = [
        ExpandedNumericSpec {
            name: "dose",
            values: &dose,
        },
        ExpandedNumericSpec {
            name: "exposure_offset",
            values: &exposure_offset,
        },
    ];
    let size_factors = [1.0, 1.1, 0.9, 1.2];
    let dispersions = [0.05, 0.08];
    let base_mean = [41.5, 40.5];
    let disp_fit = [0.05, 0.08];
    let options = BetaPriorRefitOptions {
        fit_options: IrlsOptions::default(),
        variance_options: BetaPriorVarianceOptions {
            method: BetaPriorVarianceMethod::Quantile,
            upper_quantile: 0.5,
            ..BetaPriorVarianceOptions::default()
        },
    };
    let replacement_options = CooksReplacementOptions {
        trim: 0.2,
        cooks_cutoff: 0.0,
        min_replicates: 3,
        which_samples: Some(vec![false, false, true, false]),
    };

    let formula = fit_expanded_formula_beta_prior_wald_results_with_cooks_replacement(
        ExpandedFormulaBetaPriorWaldResultsInput {
            counts: &counts,
            formula: "~ condition + offset(I(exposure_offset + dose))",
            factors: &factors,
            numeric_covariates: &numeric_covariates,
            size_factors: &size_factors,
            weights: None,
            dispersions: &dispersions,
            base_mean: &base_mean,
            disp_fit: &disp_fit,
            gene_names: None,
            options: options.clone(),
        },
        1,
        &replacement_options,
    )
    .unwrap();

    let design = expanded_formula_design("~ condition", &factors, &[]).unwrap();
    let mut normalization_values = Vec::new();
    for _ in 0..counts.n_genes() {
        for (sample, size_factor) in size_factors.iter().copied().enumerate() {
            normalization_values.push(size_factor * (exposure_offset[sample] + dose[sample]).exp());
        }
    }
    let normalization_factors =
        RowMajorMatrix::from_row_major(counts.n_genes(), counts.n_samples(), normalization_values)
            .unwrap();
    let direct =
        fit_expanded_beta_prior_wald_results_with_normalization_factors_and_weights_and_cooks_replacement(
            ExpandedBetaPriorWaldNormalizedResultsInput {
                counts: &counts,
                design: ExpandedModelBetaPriorDesignInput {
                    expanded_design: &design.expanded_design,
                    standard_design: &design.standard_design,
                    coefficient_groups: &design.coefficient_groups,
                },
                normalization_factors: &normalization_factors,
                weights: None,
                dispersions: &dispersions,
                base_mean: &base_mean,
                disp_fit: &disp_fit,
                gene_names: None,
                options: options.clone(),
            },
            1,
            &replacement_options,
        )
        .unwrap();

    assert_eq!(formula.design, design);
    assert_eq!(formula.replacement, direct);

    let formula_contrast =
        fit_expanded_formula_beta_prior_wald_contrast_results_with_cooks_replacement(
            ExpandedFormulaBetaPriorWaldResultsInput {
                counts: &counts,
                formula: "~ condition + offset(I(exposure_offset + dose))",
                factors: &factors,
                numeric_covariates: &numeric_covariates,
                size_factors: &size_factors,
                weights: None,
                dispersions: &dispersions,
                base_mean: &base_mean,
                disp_fit: &disp_fit,
                gene_names: None,
                options,
            },
            &[0.0, 1.0],
            &replacement_options,
        )
        .unwrap();
    assert!(formula_contrast.replacement.refit.is_some());
    assert_eq!(
        formula_contrast
            .replacement
            .results
            .metadata
            .result_name
            .as_deref(),
        Some("contrast")
    );
}
