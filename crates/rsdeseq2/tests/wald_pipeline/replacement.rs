use super::*;

#[test]
fn fixed_dispersion_wald_cooks_replacement_refits_marked_rows() {
    let counts = CountMatrix::from_row_major_u32(
        2,
        8,
        vec![
            0, 20, 1, 19, 2, 18, 3, 17, //
            12, 28, 10, 30, 14, 26, 11, 29,
        ],
    )
    .unwrap();
    let design = DesignMatrix::from_row_major(
        8,
        2,
        vec![
            1.0, 0.0, //
            1.0, 0.0, //
            1.0, 0.0, //
            1.0, 0.0, //
            1.0, 1.0, //
            1.0, 1.0, //
            1.0, 1.0, //
            1.0, 1.0,
        ],
        Some(vec!["Intercept".into(), "condition_B_vs_A".into()]),
    )
    .unwrap();
    let options = CooksReplacementOptions {
        trim: 0.2,
        cooks_cutoff: 0.0,
        min_replicates: 3,
        which_samples: Some(vec![true, false, false, false, false, false, false, false]),
    };

    let output = DeseqBuilder::new()
        .size_factors(vec![1.0; 8])
        .disable_independent_filtering()
        .fit_fixed_dispersion_wald_with_cooks_replacement(
            &counts,
            &design,
            &[0.1, 0.1],
            1,
            &options,
        )
        .unwrap();

    assert!(output.refit_plan.n_refit > 0);
    assert!(output.refit_plan.should_refit);
    assert!(output.refit_fit.is_some());
    assert!(output.refit_results.is_some());
    assert_ne!(
        output.refit_plan.replacement.replaced_counts.as_slice(),
        counts.as_slice()
    );

    let refit_results = output.refit_results.as_ref().unwrap();
    for gene in output.refit_plan.refit_rows.iter().copied() {
        assert_eq!(
            output.results.rows[gene].log2_fold_change,
            refit_results.rows[gene].log2_fold_change
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
}

#[test]
fn fixed_dispersion_wald_contrast_cooks_replacement_refits_marked_rows() {
    let counts = CountMatrix::from_row_major_u32(
        2,
        8,
        vec![
            0, 20, 1, 19, 2, 18, 3, 17, //
            12, 28, 10, 30, 14, 26, 11, 29,
        ],
    )
    .unwrap();
    let design = DesignMatrix::from_row_major(
        8,
        2,
        vec![
            1.0, 0.0, //
            1.0, 0.0, //
            1.0, 0.0, //
            1.0, 0.0, //
            1.0, 1.0, //
            1.0, 1.0, //
            1.0, 1.0, //
            1.0, 1.0,
        ],
        Some(vec!["Intercept".into(), "condition_B_vs_A".into()]),
    )
    .unwrap();
    let options = CooksReplacementOptions {
        trim: 0.2,
        cooks_cutoff: 0.0,
        min_replicates: 3,
        which_samples: Some(vec![true, false, false, false, false, false, false, false]),
    };

    let output = DeseqBuilder::new()
        .size_factors(vec![1.0; 8])
        .disable_independent_filtering()
        .fit_fixed_dispersion_wald_contrast_with_cooks_replacement(
            &counts,
            &design,
            &[0.1, 0.1],
            &[0.0, 1.0],
            &options,
        )
        .unwrap();

    assert!(output.refit_plan.n_refit > 0);
    assert!(output.refit_plan.should_refit);
    assert!(output.refit_fit.is_some());
    assert!(output.refit_results.is_some());

    let refit_results = output.refit_results.as_ref().unwrap();
    for gene in output.refit_plan.refit_rows.iter().copied() {
        assert_eq!(
            output.results.rows[gene].log2_fold_change,
            refit_results.rows[gene].log2_fold_change
        );
        assert_eq!(
            output.results.rows[gene].pvalue,
            refit_results.rows[gene].pvalue
        );
        assert_eq!(
            output.results.rows[gene].max_cooks,
            output.refit_plan.post_refit_max_cooks[gene]
        );
    }
}

#[test]
fn fixed_dispersion_wald_contrast_spec_cooks_replacement_preserves_metadata() {
    let counts = CountMatrix::from_row_major_u32(
        2,
        8,
        vec![
            0, 20, 1, 19, 2, 18, 3, 17, //
            12, 28, 10, 30, 14, 26, 11, 29,
        ],
    )
    .unwrap();
    let design = DesignMatrix::from_row_major(
        8,
        2,
        vec![
            1.0, 0.0, //
            1.0, 0.0, //
            1.0, 0.0, //
            1.0, 0.0, //
            1.0, 1.0, //
            1.0, 1.0, //
            1.0, 1.0, //
            1.0, 1.0,
        ],
        Some(vec!["Intercept".into(), "condition_B_vs_A".into()]),
    )
    .unwrap();
    let options = CooksReplacementOptions {
        trim: 0.2,
        cooks_cutoff: 0.0,
        min_replicates: 3,
        which_samples: Some(vec![true, false, false, false, false, false, false, false]),
    };

    let output = DeseqBuilder::new()
        .size_factors(vec![1.0; 8])
        .disable_independent_filtering()
        .fit_fixed_dispersion_wald_contrast_spec_with_cooks_replacement(
            &counts,
            &design,
            &[0.1, 0.1],
            &ContrastSpec::coefficient_name("condition_B_vs_A"),
            &options,
        )
        .unwrap();

    assert!(output.refit_plan.n_refit > 0);
    assert_eq!(
        output.results.metadata.result_name.as_deref(),
        Some("condition_B_vs_A")
    );
    assert_eq!(
        output.results.metadata.comparison.as_deref(),
        Some("coefficient condition_B_vs_A")
    );
    assert_eq!(
        output.results.metadata.contrast.as_deref(),
        Some(&[0.0, 1.0][..])
    );
    assert_eq!(
        output.original_results.metadata.contrast.as_deref(),
        Some(&[0.0, 1.0][..])
    );
    assert_eq!(
        output
            .refit_results
            .as_ref()
            .unwrap()
            .metadata
            .contrast
            .as_deref(),
        Some(&[0.0, 1.0][..])
    );
}

#[test]
fn fixed_dispersion_wald_factor_level_cooks_replacement_preserves_metadata() {
    let counts = CountMatrix::from_row_major_u32(
        2,
        8,
        vec![
            0, 20, 1, 19, 2, 18, 3, 17, //
            12, 28, 10, 30, 14, 26, 11, 29,
        ],
    )
    .unwrap();
    let design = DesignMatrix::from_row_major(
        8,
        2,
        vec![
            1.0, 0.0, //
            1.0, 0.0, //
            1.0, 0.0, //
            1.0, 0.0, //
            1.0, 1.0, //
            1.0, 1.0, //
            1.0, 1.0, //
            1.0, 1.0,
        ],
        Some(vec!["Intercept".into(), "condition_B_vs_A".into()]),
    )
    .unwrap();
    let levels = ["A", "A", "A", "A", "B", "B", "B", "B"]
        .into_iter()
        .map(String::from)
        .collect::<Vec<_>>();
    let options = CooksReplacementOptions {
        trim: 0.2,
        cooks_cutoff: 0.0,
        min_replicates: 3,
        which_samples: Some(vec![true, false, false, false, false, false, false, false]),
    };

    let output = DeseqBuilder::new()
        .size_factors(vec![1.0; 8])
        .disable_independent_filtering()
        .fit_fixed_dispersion_wald_factor_level_contrast_with_cooks_replacement(
            &counts,
            &design,
            &[0.1, 0.1],
            FactorLevelContrast::new("condition", "B", "A", &levels),
            &options,
        )
        .unwrap();
    let request = DeseqBuilder::new()
        .size_factors(vec![1.0; 8])
        .disable_independent_filtering()
        .fit_fixed_dispersion_wald_results_contrast_with_cooks_replacement(
            &counts,
            &design,
            &[0.1, 0.1],
            &ResultsContrast::character("condition", "B", "A"),
            Some(&levels),
            &options,
        )
        .unwrap();

    assert!(output.refit_plan.n_refit > 0);
    assert_eq!(request.results, output.results);
    assert_eq!(
        output.results.metadata.result_name.as_deref(),
        Some("condition_B_vs_A")
    );
    assert_eq!(
        output.results.metadata.comparison.as_deref(),
        Some("factor-level contrast: condition B vs A")
    );
    assert_eq!(
        output.results.metadata.contrast.as_deref(),
        Some(&[0.0, 1.0][..])
    );
    assert_eq!(
        output.original_results.metadata.contrast.as_deref(),
        Some(&[0.0, 1.0][..])
    );
    assert_eq!(
        output
            .refit_results
            .as_ref()
            .unwrap()
            .metadata
            .contrast
            .as_deref(),
        Some(&[0.0, 1.0][..])
    );

    let model_frame = FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "condition".to_string(),
            sample_levels: levels,
            levels: None,
            reference: None,
        }],
        numeric_covariates: Vec::new(),
    };
    let model_frame_request = DeseqBuilder::new()
        .size_factors(vec![1.0; 8])
        .disable_independent_filtering()
        .fit_fixed_dispersion_wald_results_contrast_from_model_frame_with_cooks_replacement(
            &counts,
            &design,
            &[0.1, 0.1],
            &ResultsContrast::character("condition", "B", "A"),
            &model_frame,
            &options,
        )
        .unwrap();
    assert_eq!(model_frame_request.results, output.results);
    assert_eq!(
        model_frame_request.original_fit.current_model_frame(),
        Some(&model_frame)
    );
}

#[test]
fn fixed_dispersion_wald_replacement_skips_when_no_rows_are_marked() {
    let counts = CountMatrix::from_row_major_u32(
        2,
        8,
        vec![
            0, 20, 1, 19, 2, 18, 3, 17, //
            12, 28, 10, 30, 14, 26, 11, 29,
        ],
    )
    .unwrap();
    let design = DesignMatrix::from_row_major(
        8,
        2,
        vec![
            1.0, 0.0, //
            1.0, 0.0, //
            1.0, 0.0, //
            1.0, 0.0, //
            1.0, 1.0, //
            1.0, 1.0, //
            1.0, 1.0, //
            1.0, 1.0,
        ],
        Some(vec!["Intercept".into(), "condition_B_vs_A".into()]),
    )
    .unwrap();

    let output = DeseqBuilder::new()
        .size_factors(vec![1.0; 8])
        .disable_independent_filtering()
        .fit_fixed_dispersion_wald_with_cooks_replacement(
            &counts,
            &design,
            &[0.1, 0.1],
            1,
            &CooksReplacementOptions::new(f64::MAX),
        )
        .unwrap();

    assert_eq!(output.refit_plan.n_refit, 0);
    assert!(!output.refit_plan.should_refit);
    assert!(output.refit_fit.is_none());
    assert!(output.refit_results.is_none());
    assert_eq!(
        output.results.rows.len(),
        output.original_results.rows.len()
    );
    for (final_row, original_row) in output
        .results
        .rows
        .iter()
        .zip(&output.original_results.rows)
    {
        assert_eq!(final_row.log2_fold_change, original_row.log2_fold_change);
        assert_eq!(final_row.pvalue, original_row.pvalue);
        assert_eq!(final_row.max_cooks, original_row.max_cooks);
        assert_eq!(final_row.cooks_outlier, Some(false));
    }
    assert_eq!(
        output.refit_plan.replacement.replaced_counts.as_slice(),
        counts.as_slice()
    );
}
