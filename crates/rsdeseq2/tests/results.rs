use approx::assert_relative_eq;
use rsdeseq2::prelude::*;

fn toy_fit(beta: Vec<f64>, beta_se: Vec<f64>, beta_converged: Vec<bool>) -> NbinomGlmFit {
    let n_genes = beta_converged.len();
    let n_samples = 2;
    NbinomGlmFit {
        log_like: vec![0.0; n_genes],
        beta_converged,
        beta: RowMajorMatrix::from_row_major(n_genes, 1, beta).unwrap(),
        beta_se: RowMajorMatrix::from_row_major(n_genes, 1, beta_se).unwrap(),
        beta_covariance: None,
        mu: RowMajorMatrix::from_row_major(n_genes, n_samples, vec![1.0; n_genes * n_samples])
            .unwrap(),
        beta_iter: vec![1; n_genes],
        model_matrix: DesignMatrix::from_row_major(
            n_samples,
            1,
            vec![1.0, 1.0],
            Some(vec!["Intercept".to_string()]),
        )
        .unwrap(),
        n_terms: 1,
        hat_diagonal: RowMajorMatrix::from_row_major(
            n_genes,
            n_samples,
            vec![0.5; n_genes * n_samples],
        )
        .unwrap(),
    }
}

#[test]
fn result_column_schema_matches_current_deseq2_shape() {
    assert_eq!(
        deseq2_result_core_column_names(),
        &[
            "baseMean",
            "log2FoldChange",
            "lfcSE",
            "stat",
            "pvalue",
            "padj"
        ]
    );
    assert_eq!(
        rsdeseq2_result_diagnostic_column_names(),
        &[
            "dispersion",
            "converged",
            "maxCooks",
            "cooksOutlier",
            "filtered"
        ]
    );
}

#[test]
fn result_table_metadata_exposes_description_label_precedence() {
    let wald = DeseqResultsTableMetadata {
        test_type: Some(TestType::Wald),
        result_name: Some("condition_B_vs_A".to_string()),
        comparison: Some("coefficient condition_B_vs_A".to_string()),
        ..DeseqResultsTableMetadata::default()
    };
    assert_eq!(
        wald.effect_description_label(),
        Some("coefficient condition_B_vs_A")
    );
    assert_eq!(
        wald.test_description_label(),
        Some("coefficient condition_B_vs_A")
    );

    let lrt = DeseqResultsTableMetadata {
        test_type: Some(TestType::Lrt),
        result_name: Some("condition_B_vs_A".to_string()),
        comparison: Some("full model versus reduced model".to_string()),
        ..DeseqResultsTableMetadata::default()
    };
    assert_eq!(lrt.effect_description_label(), Some("condition_B_vs_A"));
    assert_eq!(
        lrt.test_description_label(),
        Some("full model versus reduced model")
    );

    assert_eq!(
        lrt.scalar_metadata(),
        vec![
            DeseqResultsTableMetadataEntry {
                name: "testType".to_string(),
                value: "LRT".to_string(),
            },
            DeseqResultsTableMetadataEntry {
                name: "resultName".to_string(),
                value: "condition_B_vs_A".to_string(),
            },
            DeseqResultsTableMetadataEntry {
                name: "comparison".to_string(),
                value: "full model versus reduced model".to_string(),
            },
            DeseqResultsTableMetadataEntry {
                name: "lfcThreshold".to_string(),
                value: "0".to_string(),
            },
            DeseqResultsTableMetadataEntry {
                name: "pAdjustMethod".to_string(),
                value: "BH".to_string(),
            },
        ]
    );
}

#[test]
fn build_wald_results_populates_deseq2_shaped_columns() {
    let fit = toy_fit(vec![2.0, 1.0], vec![0.5, 1.0], vec![true, false]);
    let names = vec!["gene_a".to_string(), "gene_b".to_string()];
    let dispersions = vec![0.1, 0.2];
    let results =
        build_wald_results(&[10.0, 20.0], &fit, 0, Some(&names), Some(&dispersions)).unwrap();

    assert_eq!(results.len(), 2);
    assert!(!results.is_empty());
    assert_eq!(results.rows[0].gene.as_deref(), Some("gene_a"));
    assert_relative_eq!(results.rows[0].base_mean, 10.0, epsilon = 1e-12);
    assert_relative_eq!(
        results.rows[0].log2_fold_change.unwrap(),
        2.0,
        epsilon = 1e-12
    );
    assert_relative_eq!(results.rows[0].lfc_se.unwrap(), 0.5, epsilon = 1e-12);
    assert_relative_eq!(results.rows[0].stat.unwrap(), 4.0, epsilon = 1e-12);
    assert!(results.rows[0].pvalue.unwrap() < results.rows[1].pvalue.unwrap());
    assert!(results.rows[0].padj.unwrap() <= results.rows[1].padj.unwrap());
    assert_eq!(results.rows[0].dispersion, Some(0.1));
    assert_eq!(results.rows[1].converged, Some(false));
    assert_eq!(
        results.column_names(),
        vec![
            "baseMean",
            "log2FoldChange",
            "lfcSE",
            "stat",
            "pvalue",
            "padj",
            "dispersion",
            "converged"
        ]
    );
    assert_eq!(results.metadata.test_type, Some(TestType::Wald));
    assert_eq!(results.metadata.result_name.as_deref(), Some("Intercept"));
    assert_eq!(results.metadata.lfc_threshold, 0.0);
    assert_eq!(results.metadata.p_adjust_method, "BH");
    let threshold_metadata = results.clone().with_wald_test_options(
        &WaldTestOptions::normal().with_lfc_threshold(1.5, WaldAlternative::Less),
    );
    assert_eq!(threshold_metadata.metadata.lfc_threshold, 1.5);
    assert_eq!(
        threshold_metadata.metadata.alt_hypothesis.as_deref(),
        Some("less")
    );

    let metadata = results.deseq2_metadata();
    assert_eq!(metadata.table, results.metadata);
    assert!(metadata.independent_filtering.is_none());
    assert_eq!(metadata.columns[0].name, "baseMean");
    assert_eq!(metadata.columns[0].column_type, "results");
    assert_eq!(
        metadata.columns[0].description,
        "mean of normalized counts for all samples"
    );
    assert_eq!(
        metadata.columns[1].description,
        "log2 fold change (MLE): Intercept"
    );
    assert_eq!(metadata.columns[2].description, "standard error: Intercept");
    assert_eq!(metadata.columns[3].description, "Wald statistic: Intercept");
    assert_eq!(
        metadata.columns[4].description,
        "Wald test p-value: Intercept"
    );
    assert_eq!(metadata.columns[5].description, "BH adjusted p-values");
    assert_eq!(metadata.columns[6].name, "dispersion");
    assert_eq!(metadata.columns[6].column_type, "diagnostic");
}

#[test]
fn result_data_frame_assembles_typed_columns_and_metadata() {
    let fit = toy_fit(vec![2.0, f64::NAN], vec![0.5, 1.0], vec![true, false]);
    let names = vec!["gene_a".to_string(), "gene_b".to_string()];
    let results =
        build_wald_results(&[10.0, 20.0], &fit, 0, Some(&names), Some(&[0.1, 0.2])).unwrap();

    let frame = results.data_frame();

    assert_eq!(
        frame.row_names,
        vec![Some("gene_a".to_string()), Some("gene_b".to_string())]
    );
    assert_eq!(frame.metadata, results.deseq2_metadata());
    assert_eq!(
        frame
            .columns
            .iter()
            .map(|column| column.metadata.name.as_str())
            .collect::<Vec<_>>(),
        vec![
            "baseMean",
            "log2FoldChange",
            "lfcSE",
            "stat",
            "pvalue",
            "padj",
            "dispersion",
            "converged"
        ]
    );

    let base_mean = frame
        .columns
        .iter()
        .find(|column| column.metadata.name == "baseMean")
        .unwrap();
    assert_eq!(base_mean.values.len(), 2);
    assert_eq!(
        base_mean.values.as_numeric().unwrap(),
        &[Some(10.0), Some(20.0)]
    );
    assert!(base_mean.values.as_logical().is_none());

    let lfc = frame
        .columns
        .iter()
        .find(|column| column.metadata.name == "log2FoldChange")
        .unwrap();
    assert_eq!(lfc.values.as_numeric().unwrap(), &[Some(2.0), None]);

    let converged = frame
        .columns
        .iter()
        .find(|column| column.metadata.name == "converged")
        .unwrap();
    assert_eq!(
        converged.values.as_logical().unwrap(),
        &[Some(true), Some(false)]
    );
    assert!(converged.values.as_numeric().is_none());
}

#[test]
fn original_result_frame_keeps_optional_diagnostics_ordered_and_typed() {
    let results = DeseqResults {
        rows: vec![
            DeseqResultRow {
                gene: Some("gene_a".to_string()),
                base_mean: 10.0,
                log2_fold_change: Some(1.0),
                lfc_se: Some(0.5),
                stat: Some(2.0),
                pvalue: Some(0.04),
                padj: Some(0.08),
                dispersion: None,
                converged: Some(true),
                max_cooks: Some(1.5),
                cooks_outlier: Some(false),
                filtered: None,
            },
            DeseqResultRow {
                gene: Some("gene_b".to_string()),
                base_mean: 20.0,
                log2_fold_change: None,
                lfc_se: None,
                stat: None,
                pvalue: None,
                padj: None,
                dispersion: Some(0.2),
                converged: None,
                max_cooks: None,
                cooks_outlier: None,
                filtered: Some(true),
            },
        ],
        metadata: DeseqResultsTableMetadata::default(),
        independent_filtering: None,
    };

    assert_eq!(
        results.column_names(),
        vec![
            "baseMean",
            "log2FoldChange",
            "lfcSE",
            "stat",
            "pvalue",
            "padj",
            "dispersion",
            "converged",
            "maxCooks",
            "cooksOutlier",
            "filtered",
        ]
    );

    let frame = results.data_frame();
    let max_cooks = frame
        .columns
        .iter()
        .find(|column| column.metadata.name == "maxCooks")
        .unwrap();
    assert_eq!(max_cooks.values.as_numeric().unwrap(), &[Some(1.5), None]);
    assert!(max_cooks.values.as_logical().is_none());

    let filtered = frame
        .columns
        .iter()
        .find(|column| column.metadata.name == "filtered")
        .unwrap();
    assert_eq!(filtered.values.as_logical().unwrap(), &[None, Some(true)]);
    assert!(filtered.values.as_numeric().is_none());
}

#[test]
fn build_wald_results_preserves_missing_pvalues_and_padj() {
    let fit = toy_fit(vec![2.0, 1.0], vec![0.0, 1.0], vec![true, true]);
    let results = build_wald_results(&[10.0, 20.0], &fit, 0, None, None).unwrap();

    assert_eq!(results.rows[0].pvalue, None);
    assert_eq!(results.rows[0].padj, None);
    assert_eq!(results.rows[0].lfc_se, Some(0.0));
    assert!(results.rows[1].pvalue.is_some());
    assert!(results.rows[1].padj.is_some());
}

#[test]
fn build_wald_results_omits_optional_metadata_when_absent() {
    let fit = toy_fit(vec![0.0], vec![1.0], vec![true]);
    let results = build_wald_results(&[0.5], &fit, 0, None, None).unwrap();
    assert_eq!(results.rows[0].gene, None);
    assert_eq!(results.rows[0].dispersion, None);
    assert_eq!(results.rows[0].converged, Some(true));
    assert_relative_eq!(results.rows[0].pvalue.unwrap(), 1.0, epsilon = 1e-15);
}

#[test]
fn build_wald_contrast_results_uses_contrast_columns() {
    let fit = toy_fit(vec![1.0, 2.0], vec![0.5, 0.25], vec![true, false]);
    let contrast = WaldContrastOutput {
        log2_fold_change: vec![Some(2.0), Some(-1.0)],
        lfc_se: vec![Some(0.8), Some(0.5)],
        wald: WaldOutput {
            stat: vec![Some(2.5), Some(-2.0)],
            pvalue: vec![
                Some(two_sided_normal_pvalue(2.5)),
                Some(two_sided_normal_pvalue(-2.0)),
            ],
            degrees_of_freedom: None,
        },
    };
    let names = vec!["gene_a".to_string(), "gene_b".to_string()];
    let results = build_wald_contrast_results(
        &[10.0, 20.0],
        &fit,
        &contrast,
        Some(&names),
        Some(&[0.1, 0.2]),
    )
    .unwrap();

    assert_eq!(results.rows[0].gene.as_deref(), Some("gene_a"));
    assert_eq!(results.rows[0].log2_fold_change, Some(2.0));
    assert_eq!(results.rows[0].lfc_se, Some(0.8));
    assert_eq!(results.rows[0].stat, Some(2.5));
    assert_eq!(results.rows[0].dispersion, Some(0.1));
    assert_eq!(results.rows[1].converged, Some(false));
    assert!(results.rows[0].padj.is_some());
    assert!(results.rows[1].padj.is_some());
    assert_eq!(results.metadata.test_type, Some(TestType::Wald));
    assert_eq!(results.metadata.result_name.as_deref(), Some("contrast"));
    assert_eq!(
        results.metadata.comparison.as_deref(),
        Some("primitive numeric contrast")
    );
    let metadata = results.deseq2_metadata();
    assert_eq!(
        metadata.columns[1].description,
        "log2 fold change (MLE): primitive numeric contrast"
    );
    assert_eq!(
        metadata.columns[3].description,
        "Wald statistic: primitive numeric contrast"
    );
    assert_eq!(
        metadata.columns[4].description,
        "Wald test p-value: primitive numeric contrast"
    );
}

#[test]
fn build_wald_contrast_results_validates_dimensions() {
    let fit = toy_fit(vec![1.0, 2.0], vec![0.5, 0.25], vec![true, true]);
    let bad = WaldContrastOutput {
        log2_fold_change: vec![Some(2.0)],
        lfc_se: vec![Some(0.8), Some(0.5)],
        wald: WaldOutput {
            stat: vec![Some(2.5), Some(-2.0)],
            pvalue: vec![Some(0.01), Some(0.02)],
            degrees_of_freedom: None,
        },
    };

    assert!(build_wald_contrast_results(&[10.0, 20.0], &fit, &bad, None, None).is_err());
}

#[test]
fn build_wald_contrast_results_rejects_invalid_optional_outputs() {
    let fit = toy_fit(vec![1.0], vec![0.5], vec![true]);
    let invalid_lfc = WaldContrastOutput {
        log2_fold_change: vec![Some(f64::NAN)],
        lfc_se: vec![Some(0.5)],
        wald: WaldOutput {
            stat: vec![Some(2.0)],
            pvalue: vec![Some(0.05)],
            degrees_of_freedom: None,
        },
    };
    assert!(build_wald_contrast_results(&[10.0], &fit, &invalid_lfc, None, None).is_err());

    let invalid_lfc_se = WaldContrastOutput {
        log2_fold_change: vec![Some(1.0)],
        lfc_se: vec![Some(f64::INFINITY)],
        wald: WaldOutput {
            stat: vec![Some(2.0)],
            pvalue: vec![Some(0.05)],
            degrees_of_freedom: None,
        },
    };
    assert!(build_wald_contrast_results(&[10.0], &fit, &invalid_lfc_se, None, None).is_err());

    let invalid_pvalue = WaldContrastOutput {
        log2_fold_change: vec![Some(1.0)],
        lfc_se: vec![Some(0.5)],
        wald: WaldOutput {
            stat: vec![Some(2.0)],
            pvalue: vec![Some(1.2)],
            degrees_of_freedom: None,
        },
    };
    assert!(build_wald_contrast_results(&[10.0], &fit, &invalid_pvalue, None, None).is_err());

    let invalid_df = WaldContrastOutput {
        log2_fold_change: vec![Some(1.0)],
        lfc_se: vec![Some(0.5)],
        wald: WaldOutput {
            stat: vec![Some(2.0)],
            pvalue: vec![Some(0.05)],
            degrees_of_freedom: Some(vec![Some(0.0)]),
        },
    };
    assert!(build_wald_contrast_results(&[10.0], &fit, &invalid_df, None, None).is_err());
}

#[test]
fn build_wald_results_validates_dimensions() {
    let fit = toy_fit(vec![1.0, 2.0], vec![1.0, 1.0], vec![true, true]);
    assert!(build_wald_results(&[1.0], &fit, 0, None, None).is_err());
    assert!(build_wald_results(&[1.0, f64::NAN], &fit, 0, None, None).is_err());
    assert!(build_wald_results(&[1.0, -1.0], &fit, 0, None, None).is_err());

    let bad_names = vec!["gene_a".to_string()];
    assert!(build_wald_results(&[1.0, 2.0], &fit, 0, Some(&bad_names), None).is_err());

    assert!(build_wald_results(&[1.0, 2.0], &fit, 0, None, Some(&[0.1])).is_err());
    assert!(build_wald_results(&[1.0, 2.0], &fit, 1, None, None).is_err());
}

#[test]
fn build_wald_results_from_expanded_model_fit_reports_collapsed_coefficients() {
    let expanded_design = DesignMatrix::from_row_major(
        2,
        4,
        vec![
            1.0, 1.0, 0.0, 0.0, //
            1.0, 0.0, 1.0, 0.0,
        ],
        Some(vec![
            "Intercept".into(),
            "condition_A".into(),
            "condition_B".into(),
            "condition_C".into(),
        ]),
    )
    .unwrap();
    let standard_design = DesignMatrix::from_row_major(
        2,
        2,
        vec![
            1.0, 0.0, //
            1.0, 1.0,
        ],
        Some(vec!["Intercept".into(), "condition_B_vs_A".into()]),
    )
    .unwrap();
    let expanded_fit = NbinomGlmFit {
        log_like: vec![-10.0, -20.0],
        beta_converged: vec![true, true],
        beta: RowMajorMatrix::from_row_major(
            2,
            4,
            vec![
                4.0, 1.0, 3.0, 5.0, //
                6.0, -2.0, 2.0, 4.0,
            ],
        )
        .unwrap(),
        beta_se: RowMajorMatrix::from_row_major(2, 4, vec![1.0; 8]).unwrap(),
        beta_covariance: Some(
            RowMajorMatrix::from_row_major(
                2,
                16,
                vec![
                    4.0, 1.0, 2.0, 3.0, //
                    1.0, 9.0, 4.0, 5.0, //
                    2.0, 4.0, 16.0, 6.0, //
                    3.0, 5.0, 6.0, 25.0, //
                    1.0, 0.0, 0.0, 0.0, //
                    0.0, 4.0, 1.0, 1.0, //
                    0.0, 1.0, 9.0, 1.0, //
                    0.0, 1.0, 1.0, 16.0,
                ],
            )
            .unwrap(),
        ),
        mu: RowMajorMatrix::from_row_major(2, 2, vec![10.0, 20.0, 5.0, 6.0]).unwrap(),
        beta_iter: vec![5, 6],
        model_matrix: expanded_design,
        n_terms: 4,
        hat_diagonal: RowMajorMatrix::from_row_major(2, 2, vec![0.1; 4]).unwrap(),
    };
    let names = vec!["gene_a".to_string(), "gene_b".to_string()];

    let results = build_wald_results_from_expanded_model_fit(
        &[10.0, 20.0],
        &expanded_fit,
        &standard_design,
        &[vec![0], vec![1, 2]],
        1,
        Some(&names),
        Some(&[0.1, 0.2]),
    )
    .unwrap();

    assert_eq!(results.rows[0].gene.as_deref(), Some("gene_a"));
    assert_eq!(results.rows[0].log2_fold_change, Some(2.0));
    assert_relative_eq!(
        results.rows[0].lfc_se.unwrap(),
        8.25_f64.sqrt(),
        epsilon = 1e-12
    );
    assert_relative_eq!(
        results.rows[0].stat.unwrap(),
        2.0 / 8.25_f64.sqrt(),
        epsilon = 1e-12
    );
    assert_eq!(results.rows[0].dispersion, Some(0.1));
    assert_eq!(
        results.metadata.result_name.as_deref(),
        Some("condition_B_vs_A")
    );
}

#[test]
fn build_wald_contrast_results_from_expanded_model_fit_reports_collapsed_contrast() {
    let expanded_design = DesignMatrix::from_row_major(
        2,
        3,
        vec![
            1.0, 1.0, 0.0, //
            1.0, 0.0, 1.0,
        ],
        Some(vec![
            "Intercept".into(),
            "condition_A".into(),
            "condition_B".into(),
        ]),
    )
    .unwrap();
    let standard_design = DesignMatrix::from_row_major(
        2,
        2,
        vec![
            1.0, 0.0, //
            1.0, 1.0,
        ],
        Some(vec!["Intercept".into(), "condition_B_vs_A".into()]),
    )
    .unwrap();
    let expanded_fit = NbinomGlmFit {
        log_like: vec![-10.0, -20.0],
        beta_converged: vec![true, true],
        beta: RowMajorMatrix::from_row_major(
            2,
            3,
            vec![
                4.0, 1.0, 3.0, //
                6.0, -2.0, 2.0,
            ],
        )
        .unwrap(),
        beta_se: RowMajorMatrix::from_row_major(2, 3, vec![1.0; 6]).unwrap(),
        beta_covariance: Some(
            RowMajorMatrix::from_row_major(
                2,
                9,
                vec![
                    4.0, 1.0, 2.0, //
                    1.0, 9.0, 4.0, //
                    2.0, 4.0, 16.0, //
                    1.0, 0.0, 0.0, //
                    0.0, 4.0, 1.0, //
                    0.0, 1.0, 9.0,
                ],
            )
            .unwrap(),
        ),
        mu: RowMajorMatrix::from_row_major(2, 2, vec![10.0, 20.0, 5.0, 6.0]).unwrap(),
        beta_iter: vec![5, 6],
        model_matrix: expanded_design,
        n_terms: 3,
        hat_diagonal: RowMajorMatrix::from_row_major(2, 2, vec![0.1; 4]).unwrap(),
    };
    let names = vec!["gene_a".to_string(), "gene_b".to_string()];

    let results = build_wald_contrast_results_from_expanded_model_fit(
        &[10.0, 20.0],
        &expanded_fit,
        &standard_design,
        &[vec![0], vec![1, 2]],
        &[0.0, 1.0],
        Some(&names),
        Some(&[0.1, 0.2]),
    )
    .unwrap();

    assert_eq!(results.rows[0].gene.as_deref(), Some("gene_a"));
    assert_eq!(results.rows[0].log2_fold_change, Some(2.0));
    assert_relative_eq!(
        results.rows[0].lfc_se.unwrap(),
        8.25_f64.sqrt(),
        epsilon = 1e-12
    );
    assert_relative_eq!(
        results.rows[0].stat.unwrap(),
        2.0 / 8.25_f64.sqrt(),
        epsilon = 1e-12
    );
    assert_eq!(results.rows[0].dispersion, Some(0.1));
    assert_eq!(results.metadata.test_type, Some(TestType::Wald));
    assert_eq!(results.metadata.result_name.as_deref(), Some("contrast"));
    assert_eq!(
        results.metadata.comparison.as_deref(),
        Some("primitive numeric contrast")
    );
}

#[test]
fn build_wald_results_from_expanded_beta_prior_fit_uses_collapsed_prior_fit() {
    let standard_design = DesignMatrix::from_row_major(
        2,
        2,
        vec![
            1.0, 0.0, //
            1.0, 1.0,
        ],
        Some(vec!["Intercept".into(), "condition_B_vs_A".into()]),
    )
    .unwrap();
    let prior_fit = NbinomGlmFit {
        log_like: vec![-10.0, -20.0],
        beta_converged: vec![true, false],
        beta: RowMajorMatrix::from_row_major(2, 2, vec![4.0, 2.0, 6.0, -1.0]).unwrap(),
        beta_se: RowMajorMatrix::from_row_major(2, 2, vec![0.5, 0.25, 0.75, 0.5]).unwrap(),
        beta_covariance: Some(
            RowMajorMatrix::from_row_major(
                2,
                4,
                vec![
                    0.25, 0.0, //
                    0.0, 0.0625, //
                    0.5625, 0.0, //
                    0.0, 0.25,
                ],
            )
            .unwrap(),
        ),
        mu: RowMajorMatrix::from_row_major(2, 2, vec![10.0, 20.0, 5.0, 6.0]).unwrap(),
        beta_iter: vec![5, 6],
        model_matrix: standard_design,
        n_terms: 2,
        hat_diagonal: RowMajorMatrix::from_row_major(2, 2, vec![0.1; 4]).unwrap(),
    };
    let fit = ExpandedModelBetaPriorGlmFit {
        expanded_mle_fit: prior_fit.clone(),
        expanded_prior_fit: prior_fit.clone(),
        prior_fit: prior_fit.clone(),
        beta_prior_variance: vec![1e6, 1.0],
    };
    let names = vec!["gene_a".to_string(), "gene_b".to_string()];

    let coefficient_results = build_wald_results_from_expanded_beta_prior_fit(
        &[10.0, 20.0],
        &fit,
        1,
        Some(&names),
        Some(&[0.1, 0.2]),
    )
    .unwrap();
    let direct_results = build_wald_results(
        &[10.0, 20.0],
        &prior_fit,
        1,
        Some(&names),
        Some(&[0.1, 0.2]),
    )
    .unwrap();

    assert_eq!(coefficient_results, direct_results);
    assert_eq!(
        coefficient_results.metadata.result_name.as_deref(),
        Some("condition_B_vs_A")
    );

    let contrast_results = build_wald_contrast_results_from_expanded_beta_prior_fit(
        &[10.0, 20.0],
        &fit,
        &[0.0, 1.0],
        Some(&names),
        Some(&[0.1, 0.2]),
    )
    .unwrap();
    let direct_contrast = wald_test_contrast(&prior_fit, &[0.0, 1.0]).unwrap();
    let direct_contrast_results = build_wald_contrast_results(
        &[10.0, 20.0],
        &prior_fit,
        &direct_contrast,
        Some(&names),
        Some(&[0.1, 0.2]),
    )
    .unwrap();

    assert_eq!(contrast_results, direct_contrast_results);
    assert_eq!(
        contrast_results.metadata.result_name.as_deref(),
        Some("contrast")
    );
}

#[test]
fn fit_expanded_beta_prior_wald_results_runs_fit_and_result_workflow() {
    let counts =
        CountMatrix::from_row_major_u32(2, 4, vec![10, 12, 20, 24, 30, 33, 45, 54]).unwrap();
    let expanded_design = DesignMatrix::from_row_major(
        4,
        3,
        vec![
            1.0, 0.0, 0.0, //
            1.0, 0.0, 1.0, //
            1.0, 1.0, 0.0, //
            1.0, 1.0, 1.0,
        ],
        Some(vec![
            "Intercept".into(),
            "condition_B".into(),
            "batch_Y".into(),
        ]),
    )
    .unwrap();
    let standard_design = DesignMatrix::from_row_major(
        4,
        2,
        vec![
            1.0, 0.0, //
            1.0, 1.0, //
            1.0, 1.0, //
            1.0, 2.0,
        ],
        Some(vec!["Intercept".into(), "condition_or_batch".into()]),
    )
    .unwrap();
    let size_factors = [1.0, 1.0, 1.0, 1.0];
    let dispersions = [0.05, 0.08];
    let base_mean = [16.5, 40.5];
    let disp_fit = [0.05, 0.08];
    let groups = [vec![0], vec![1, 2]];
    let names = vec!["gene_a".to_string(), "gene_b".to_string()];
    let options = BetaPriorRefitOptions {
        fit_options: IrlsOptions::default(),
        variance_options: BetaPriorVarianceOptions {
            method: BetaPriorVarianceMethod::Quantile,
            upper_quantile: 0.5,
            ..BetaPriorVarianceOptions::default()
        },
    };
    let design = ExpandedModelBetaPriorDesignInput {
        expanded_design: &expanded_design,
        standard_design: &standard_design,
        coefficient_groups: &groups,
    };
    let input = ExpandedBetaPriorWaldResultsInput {
        counts: &counts,
        design,
        size_factors: &size_factors,
        weights: None,
        dispersions: &dispersions,
        base_mean: &base_mean,
        disp_fit: &disp_fit,
        gene_names: Some(&names),
        options,
    };

    let coefficient_workflow = fit_expanded_beta_prior_wald_results(input.clone(), 1).unwrap();
    let direct_fit = fit_expanded_glms_with_estimated_beta_prior_variance(
        &counts,
        design,
        &size_factors,
        &dispersions,
        &base_mean,
        &disp_fit,
        input.options.clone(),
    )
    .unwrap();
    let direct_results = build_wald_results_from_expanded_beta_prior_fit(
        &base_mean,
        &direct_fit,
        1,
        Some(&names),
        Some(&dispersions),
    )
    .unwrap();

    assert_eq!(coefficient_workflow.fit, direct_fit);
    assert_eq!(coefficient_workflow.results, direct_results);
    assert_eq!(
        coefficient_workflow.results.metadata.result_name.as_deref(),
        Some("condition_or_batch")
    );

    let contrast_workflow =
        fit_expanded_beta_prior_wald_contrast_results(input, &[0.0, 1.0]).unwrap();
    let direct_contrast_results = build_wald_contrast_results_from_expanded_beta_prior_fit(
        &base_mean,
        &coefficient_workflow.fit,
        &[0.0, 1.0],
        Some(&names),
        Some(&dispersions),
    )
    .unwrap();

    assert_eq!(contrast_workflow.fit, coefficient_workflow.fit);
    assert_eq!(contrast_workflow.results, direct_contrast_results);
    assert_eq!(
        contrast_workflow.results.metadata.result_name.as_deref(),
        Some("contrast")
    );
}

#[test]
fn fit_expanded_beta_prior_wald_results_accepts_offsets_and_weights() {
    let counts =
        CountMatrix::from_row_major_u32(2, 4, vec![10, 12, 20, 24, 30, 33, 45, 54]).unwrap();
    let expanded_design = DesignMatrix::from_row_major(
        4,
        3,
        vec![
            1.0, 0.0, 0.0, //
            1.0, 0.0, 1.0, //
            1.0, 1.0, 0.0, //
            1.0, 1.0, 1.0,
        ],
        Some(vec![
            "Intercept".into(),
            "condition_B".into(),
            "batch_Y".into(),
        ]),
    )
    .unwrap();
    let standard_design = DesignMatrix::from_row_major(
        4,
        2,
        vec![
            1.0, 0.0, //
            1.0, 1.0, //
            1.0, 1.0, //
            1.0, 2.0,
        ],
        Some(vec!["Intercept".into(), "condition_or_batch".into()]),
    )
    .unwrap();
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
    let base_mean = [16.5, 40.5];
    let disp_fit = [0.05, 0.08];
    let groups = [vec![0], vec![1, 2]];
    let names = vec!["gene_a".to_string(), "gene_b".to_string()];
    let options = BetaPriorRefitOptions {
        fit_options: IrlsOptions::default(),
        variance_options: BetaPriorVarianceOptions {
            method: BetaPriorVarianceMethod::Quantile,
            upper_quantile: 0.5,
            ..BetaPriorVarianceOptions::default()
        },
    };
    let design = ExpandedModelBetaPriorDesignInput {
        expanded_design: &expanded_design,
        standard_design: &standard_design,
        coefficient_groups: &groups,
    };
    let input = ExpandedBetaPriorWaldNormalizedResultsInput {
        counts: &counts,
        design,
        normalization_factors: &normalization_factors,
        weights: Some(&weights),
        dispersions: &dispersions,
        base_mean: &base_mean,
        disp_fit: &disp_fit,
        gene_names: Some(&names),
        options,
    };

    let coefficient_workflow =
        fit_expanded_beta_prior_wald_results_with_normalization_factors_and_weights(
            input.clone(),
            1,
        )
        .unwrap();
    let direct_fit =
        fit_expanded_glms_with_estimated_beta_prior_variance_and_normalization_factors_and_weights(
            &counts,
            design,
            BetaPriorNormalizationFactorWeightInput {
                normalization_factors: &normalization_factors,
                weights: Some(&weights),
            },
            &dispersions,
            &base_mean,
            &disp_fit,
            input.options.clone(),
        )
        .unwrap();
    let direct_results = build_wald_results_from_expanded_beta_prior_fit(
        &base_mean,
        &direct_fit,
        1,
        Some(&names),
        Some(&dispersions),
    )
    .unwrap();

    assert_eq!(coefficient_workflow.fit, direct_fit);
    assert_eq!(coefficient_workflow.results, direct_results);

    let contrast_workflow =
        fit_expanded_beta_prior_wald_contrast_results_with_normalization_factors_and_weights(
            input,
            &[0.0, 1.0],
        )
        .unwrap();
    let direct_contrast_results = build_wald_contrast_results_from_expanded_beta_prior_fit(
        &base_mean,
        &coefficient_workflow.fit,
        &[0.0, 1.0],
        Some(&names),
        Some(&dispersions),
    )
    .unwrap();

    assert_eq!(contrast_workflow.fit, coefficient_workflow.fit);
    assert_eq!(contrast_workflow.results, direct_contrast_results);
}

#[test]
fn fit_expanded_factor_beta_prior_wald_results_builds_design_and_matches_direct_workflow() {
    let counts =
        CountMatrix::from_row_major_u32(2, 4, vec![10, 12, 20, 24, 30, 33, 45, 54]).unwrap();
    let sample_levels = vec![
        "A".to_string(),
        "A".to_string(),
        "B".to_string(),
        "B".to_string(),
    ];
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
    let direct_design = expanded_factor_design("condition", &sample_levels, "A").unwrap();
    let direct_design_input = ExpandedModelBetaPriorDesignInput {
        expanded_design: &direct_design.expanded_design,
        standard_design: &direct_design.standard_design,
        coefficient_groups: &direct_design.coefficient_groups,
    };
    let direct_input = ExpandedBetaPriorWaldResultsInput {
        counts: &counts,
        design: direct_design_input,
        size_factors: &size_factors,
        weights: None,
        dispersions: &dispersions,
        base_mean: &base_mean,
        disp_fit: &disp_fit,
        gene_names: Some(&names),
        options: options.clone(),
    };
    let direct = fit_expanded_beta_prior_wald_results(direct_input, 1).unwrap();

    let factor = fit_expanded_factor_beta_prior_wald_results(
        ExpandedFactorBetaPriorWaldResultsInput {
            counts: &counts,
            factor: "condition",
            sample_levels: &sample_levels,
            reference: "A",
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

    assert_eq!(factor.design, direct_design);
    assert_eq!(factor.fit, direct.fit);
    assert_eq!(factor.results, direct.results);
    assert_eq!(
        factor.results.metadata.result_name.as_deref(),
        Some("condition_B_vs_A")
    );

    let factor_contrast = fit_expanded_factor_beta_prior_wald_contrast_results(
        ExpandedFactorBetaPriorWaldResultsInput {
            counts: &counts,
            factor: "condition",
            sample_levels: &sample_levels,
            reference: "A",
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
    let direct_contrast_results = build_wald_contrast_results_from_expanded_beta_prior_fit(
        &base_mean,
        &factor.fit,
        &[0.0, 1.0],
        Some(&names),
        Some(&dispersions),
    )
    .unwrap();

    assert_eq!(factor_contrast.fit, factor.fit);
    assert_eq!(factor_contrast.results, direct_contrast_results);
}

#[test]
fn fit_expanded_factor_beta_prior_wald_results_accepts_normalization_factors_and_weights() {
    let counts =
        CountMatrix::from_row_major_u32(2, 4, vec![10, 12, 20, 24, 30, 33, 45, 54]).unwrap();
    let sample_levels = vec![
        "A".to_string(),
        "A".to_string(),
        "B".to_string(),
        "B".to_string(),
    ];
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
    let direct_design = expanded_factor_design("condition", &sample_levels, "A").unwrap();
    let direct_design_input = ExpandedModelBetaPriorDesignInput {
        expanded_design: &direct_design.expanded_design,
        standard_design: &direct_design.standard_design,
        coefficient_groups: &direct_design.coefficient_groups,
    };
    let direct_input = ExpandedBetaPriorWaldNormalizedResultsInput {
        counts: &counts,
        design: direct_design_input,
        normalization_factors: &normalization_factors,
        weights: Some(&weights),
        dispersions: &dispersions,
        base_mean: &base_mean,
        disp_fit: &disp_fit,
        gene_names: Some(&names),
        options: options.clone(),
    };
    let direct = fit_expanded_beta_prior_wald_results_with_normalization_factors_and_weights(
        direct_input,
        1,
    )
    .unwrap();

    let factor =
        fit_expanded_factor_beta_prior_wald_results_with_normalization_factors_and_weights(
            ExpandedFactorBetaPriorWaldNormalizedResultsInput {
                counts: &counts,
                factor: "condition",
                sample_levels: &sample_levels,
                reference: "A",
                normalization_factors: &normalization_factors,
                weights: Some(&weights),
                dispersions: &dispersions,
                base_mean: &base_mean,
                disp_fit: &disp_fit,
                gene_names: Some(&names),
                options,
            },
            1,
        )
        .unwrap();

    assert_eq!(factor.design, direct_design);
    assert_eq!(factor.fit, direct.fit);
    assert_eq!(factor.results, direct.results);

    let factor_contrast =
        fit_expanded_factor_beta_prior_wald_contrast_results_with_normalization_factors_and_weights(
            ExpandedFactorBetaPriorWaldNormalizedResultsInput {
                counts: &counts,
                factor: "condition",
                sample_levels: &sample_levels,
                reference: "A",
                normalization_factors: &normalization_factors,
                weights: Some(&weights),
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
    let direct_contrast_results = build_wald_contrast_results_from_expanded_beta_prior_fit(
        &base_mean,
        &factor.fit,
        &[0.0, 1.0],
        Some(&names),
        Some(&dispersions),
    )
    .unwrap();

    assert_eq!(factor_contrast.fit, factor.fit);
    assert_eq!(factor_contrast.results, direct_contrast_results);
}

#[test]
fn fit_expanded_additive_beta_prior_wald_results_builds_design_and_matches_direct_workflow() {
    let counts =
        CountMatrix::from_row_major_u32(2, 4, vec![10, 12, 20, 24, 30, 33, 45, 54]).unwrap();
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
    let factors = [
        ExpandedFactorSpec {
            factor: "condition",
            sample_levels: &condition,
            reference: "A",
        },
        ExpandedFactorSpec {
            factor: "batch",
            sample_levels: &batch,
            reference: "X",
        },
    ];
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
    let direct_design = expanded_additive_factor_design(&factors).unwrap();
    let direct_design_input = ExpandedModelBetaPriorDesignInput {
        expanded_design: &direct_design.expanded_design,
        standard_design: &direct_design.standard_design,
        coefficient_groups: &direct_design.coefficient_groups,
    };
    let direct_input = ExpandedBetaPriorWaldResultsInput {
        counts: &counts,
        design: direct_design_input,
        size_factors: &size_factors,
        weights: None,
        dispersions: &dispersions,
        base_mean: &base_mean,
        disp_fit: &disp_fit,
        gene_names: Some(&names),
        options: options.clone(),
    };
    let direct = fit_expanded_beta_prior_wald_results(direct_input, 1).unwrap();

    let additive = fit_expanded_additive_beta_prior_wald_results(
        ExpandedAdditiveBetaPriorWaldResultsInput {
            counts: &counts,
            factors: &factors,
            numeric_covariates: &[],
            interactions: &[],
            factor_numeric_interactions: &[],
            numeric_interactions: &[],
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

    assert_eq!(additive.design, direct_design);
    assert_eq!(additive.fit, direct.fit);
    assert_eq!(additive.results, direct.results);
    assert_eq!(
        additive.results.metadata.result_name.as_deref(),
        Some("condition_B_vs_A")
    );

    let contrast = fit_expanded_additive_beta_prior_wald_contrast_results(
        ExpandedAdditiveBetaPriorWaldResultsInput {
            counts: &counts,
            factors: &factors,
            numeric_covariates: &[],
            interactions: &[],
            factor_numeric_interactions: &[],
            numeric_interactions: &[],
            size_factors: &size_factors,
            weights: None,
            dispersions: &dispersions,
            base_mean: &base_mean,
            disp_fit: &disp_fit,
            gene_names: Some(&names),
            options,
        },
        &[0.0, 1.0, -1.0],
    )
    .unwrap();
    let direct_contrast_results = build_wald_contrast_results_from_expanded_beta_prior_fit(
        &base_mean,
        &additive.fit,
        &[0.0, 1.0, -1.0],
        Some(&names),
        Some(&dispersions),
    )
    .unwrap();

    assert_eq!(contrast.fit, additive.fit);
    assert_eq!(contrast.results, direct_contrast_results);
}

#[test]
fn fit_expanded_additive_beta_prior_wald_results_accepts_normalization_factors_and_weights() {
    let counts =
        CountMatrix::from_row_major_u32(2, 4, vec![10, 12, 20, 24, 30, 33, 45, 54]).unwrap();
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
    let factors = [
        ExpandedFactorSpec {
            factor: "condition",
            sample_levels: &condition,
            reference: "A",
        },
        ExpandedFactorSpec {
            factor: "batch",
            sample_levels: &batch,
            reference: "X",
        },
    ];
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
    let direct_design = expanded_additive_factor_design(&factors).unwrap();
    let direct_design_input = ExpandedModelBetaPriorDesignInput {
        expanded_design: &direct_design.expanded_design,
        standard_design: &direct_design.standard_design,
        coefficient_groups: &direct_design.coefficient_groups,
    };
    let direct_input = ExpandedBetaPriorWaldNormalizedResultsInput {
        counts: &counts,
        design: direct_design_input,
        normalization_factors: &normalization_factors,
        weights: Some(&weights),
        dispersions: &dispersions,
        base_mean: &base_mean,
        disp_fit: &disp_fit,
        gene_names: Some(&names),
        options: options.clone(),
    };
    let direct = fit_expanded_beta_prior_wald_results_with_normalization_factors_and_weights(
        direct_input,
        2,
    )
    .unwrap();

    let additive =
        fit_expanded_additive_beta_prior_wald_results_with_normalization_factors_and_weights(
            ExpandedAdditiveBetaPriorWaldNormalizedResultsInput {
                counts: &counts,
                factors: &factors,
                numeric_covariates: &[],
                interactions: &[],
                factor_numeric_interactions: &[],
                numeric_interactions: &[],
                normalization_factors: &normalization_factors,
                weights: Some(&weights),
                dispersions: &dispersions,
                base_mean: &base_mean,
                disp_fit: &disp_fit,
                gene_names: Some(&names),
                options,
            },
            2,
        )
        .unwrap();

    assert_eq!(additive.design, direct_design);
    assert_eq!(additive.fit, direct.fit);
    assert_eq!(additive.results, direct.results);

    let contrast =
        fit_expanded_additive_beta_prior_wald_contrast_results_with_normalization_factors_and_weights(
            ExpandedAdditiveBetaPriorWaldNormalizedResultsInput {
                counts: &counts,
                factors: &factors,
                numeric_covariates: &[],
                interactions: &[],
                factor_numeric_interactions: &[],
                numeric_interactions: &[],
                normalization_factors: &normalization_factors,
                weights: Some(&weights),
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
            &[0.0, 1.0, -1.0],
        )
        .unwrap();
    let direct_contrast_results = build_wald_contrast_results_from_expanded_beta_prior_fit(
        &base_mean,
        &additive.fit,
        &[0.0, 1.0, -1.0],
        Some(&names),
        Some(&dispersions),
    )
    .unwrap();

    assert_eq!(contrast.fit, additive.fit);
    assert_eq!(contrast.results, direct_contrast_results);
}

#[test]
fn fit_expanded_additive_beta_prior_wald_results_accepts_numeric_covariates() {
    let counts =
        CountMatrix::from_row_major_u32(2, 4, vec![10, 12, 20, 24, 30, 33, 45, 54]).unwrap();
    let condition = vec![
        "A".to_string(),
        "A".to_string(),
        "B".to_string(),
        "B".to_string(),
    ];
    let dose = [0.0, 1.0, 0.0, 1.0];
    let factors = [ExpandedFactorSpec {
        factor: "condition",
        sample_levels: &condition,
        reference: "A",
    }];
    let numeric_covariates = [ExpandedNumericSpec {
        name: "dose",
        values: &dose,
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
    let direct_design = expanded_additive_design(&factors, &numeric_covariates).unwrap();
    let direct_design_input = ExpandedModelBetaPriorDesignInput {
        expanded_design: &direct_design.expanded_design,
        standard_design: &direct_design.standard_design,
        coefficient_groups: &direct_design.coefficient_groups,
    };
    let direct_input = ExpandedBetaPriorWaldResultsInput {
        counts: &counts,
        design: direct_design_input,
        size_factors: &size_factors,
        weights: None,
        dispersions: &dispersions,
        base_mean: &base_mean,
        disp_fit: &disp_fit,
        gene_names: Some(&names),
        options: options.clone(),
    };
    let direct = fit_expanded_beta_prior_wald_results(direct_input, 2).unwrap();

    let additive = fit_expanded_additive_beta_prior_wald_results(
        ExpandedAdditiveBetaPriorWaldResultsInput {
            counts: &counts,
            factors: &factors,
            numeric_covariates: &numeric_covariates,
            interactions: &[],
            factor_numeric_interactions: &[],
            numeric_interactions: &[],
            size_factors: &size_factors,
            weights: None,
            dispersions: &dispersions,
            base_mean: &base_mean,
            disp_fit: &disp_fit,
            gene_names: Some(&names),
            options,
        },
        2,
    )
    .unwrap();

    assert_eq!(additive.design, direct_design);
    assert_eq!(additive.fit, direct.fit);
    assert_eq!(additive.results, direct.results);
    assert_eq!(
        additive.results.metadata.result_name.as_deref(),
        Some("dose")
    );
}

#[test]
fn fit_expanded_additive_beta_prior_wald_results_accepts_factor_interactions() {
    let counts =
        CountMatrix::from_row_major_u32(2, 4, vec![10, 12, 20, 24, 30, 33, 45, 54]).unwrap();
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
    let factors = [
        ExpandedFactorSpec {
            factor: "condition",
            sample_levels: &condition,
            reference: "A",
        },
        ExpandedFactorSpec {
            factor: "batch",
            sample_levels: &batch,
            reference: "X",
        },
    ];
    let interactions = [ExpandedFactorInteractionSpec {
        left_factor: "condition",
        right_factor: "batch",
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
    let direct_design =
        expanded_additive_design_with_interactions(&factors, &[], &interactions).unwrap();
    let direct_design_input = ExpandedModelBetaPriorDesignInput {
        expanded_design: &direct_design.expanded_design,
        standard_design: &direct_design.standard_design,
        coefficient_groups: &direct_design.coefficient_groups,
    };
    let direct_input = ExpandedBetaPriorWaldResultsInput {
        counts: &counts,
        design: direct_design_input,
        size_factors: &size_factors,
        weights: None,
        dispersions: &dispersions,
        base_mean: &base_mean,
        disp_fit: &disp_fit,
        gene_names: Some(&names),
        options: options.clone(),
    };
    let direct = fit_expanded_beta_prior_wald_results(direct_input, 3).unwrap();

    let additive = fit_expanded_additive_beta_prior_wald_results(
        ExpandedAdditiveBetaPriorWaldResultsInput {
            counts: &counts,
            factors: &factors,
            numeric_covariates: &[],
            interactions: &interactions,
            factor_numeric_interactions: &[],
            numeric_interactions: &[],
            size_factors: &size_factors,
            weights: None,
            dispersions: &dispersions,
            base_mean: &base_mean,
            disp_fit: &disp_fit,
            gene_names: Some(&names),
            options,
        },
        3,
    )
    .unwrap();

    assert_eq!(additive.design, direct_design);
    assert_eq!(additive.fit, direct.fit);
    assert_eq!(additive.results, direct.results);
    assert_eq!(
        additive.results.metadata.result_name.as_deref(),
        Some("condition_B_vs_A:batch_Y_vs_X")
    );
}

#[test]
fn fit_expanded_additive_beta_prior_wald_results_accepts_numeric_interactions() {
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
    let direct_design = expanded_additive_design_with_all_interactions(
        &factors,
        &numeric_covariates,
        &[],
        &factor_numeric,
        &numeric_interactions,
    )
    .unwrap();
    let direct_design_input = ExpandedModelBetaPriorDesignInput {
        expanded_design: &direct_design.expanded_design,
        standard_design: &direct_design.standard_design,
        coefficient_groups: &direct_design.coefficient_groups,
    };
    let direct_input = ExpandedBetaPriorWaldResultsInput {
        counts: &counts,
        design: direct_design_input,
        size_factors: &size_factors,
        weights: None,
        dispersions: &dispersions,
        base_mean: &base_mean,
        disp_fit: &disp_fit,
        gene_names: Some(&names),
        options: options.clone(),
    };
    let direct = fit_expanded_beta_prior_wald_results(direct_input, 4).unwrap();

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
            options,
        },
        4,
    )
    .unwrap();

    assert_eq!(additive.design, direct_design);
    assert_eq!(additive.fit, direct.fit);
    assert_eq!(additive.results, direct.results);
    assert_eq!(
        additive.results.metadata.result_name.as_deref(),
        Some("condition_B_vs_A:dose")
    );
}

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
fn fit_expanded_formula_beta_prior_wald_results_accepts_normalization_factors_and_weights() {
    let counts =
        CountMatrix::from_row_major_u32(2, 4, vec![10, 12, 20, 24, 30, 33, 45, 54]).unwrap();
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
    let factors = [
        ExpandedFactorSpec {
            factor: "condition",
            sample_levels: &condition,
            reference: "A",
        },
        ExpandedFactorSpec {
            factor: "batch",
            sample_levels: &batch,
            reference: "X",
        },
    ];
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

    let formula =
        fit_expanded_formula_beta_prior_wald_results_with_normalization_factors_and_weights(
            ExpandedFormulaBetaPriorWaldNormalizedResultsInput {
                counts: &counts,
                formula: "~ condition * batch",
                factors: &factors,
                numeric_covariates: &[],
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

    let additive =
        fit_expanded_additive_beta_prior_wald_results_with_normalization_factors_and_weights(
            ExpandedAdditiveBetaPriorWaldNormalizedResultsInput {
                counts: &counts,
                factors: &factors,
                numeric_covariates: &[],
                interactions: &[ExpandedFactorInteractionSpec {
                    left_factor: "condition",
                    right_factor: "batch",
                }],
                factor_numeric_interactions: &[],
                numeric_interactions: &[],
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

    assert_eq!(formula.design, additive.design);
    assert_eq!(formula.fit, additive.fit);
    assert_eq!(formula.results, additive.results);

    let formula_contrast =
        fit_expanded_formula_beta_prior_wald_contrast_results_with_normalization_factors_and_weights(
            ExpandedFormulaBetaPriorWaldNormalizedResultsInput {
                counts: &counts,
                formula: "~ condition * batch",
                factors: &factors,
                numeric_covariates: &[],
                normalization_factors: &normalization_factors,
                weights: Some(&weights),
                dispersions: &dispersions,
                base_mean: &base_mean,
                disp_fit: &disp_fit,
                gene_names: Some(&names),
                options,
            },
            &[0.0, 0.0, 0.0, 1.0],
        )
        .unwrap();
    let direct_contrast_results = build_wald_contrast_results_from_expanded_beta_prior_fit(
        &base_mean,
        &formula.fit,
        &[0.0, 0.0, 0.0, 1.0],
        Some(&names),
        Some(&dispersions),
    )
    .unwrap();

    assert_eq!(formula_contrast.fit, formula.fit);
    assert_eq!(formula_contrast.results, direct_contrast_results);
}

#[test]
fn build_lrt_results_uses_full_model_beta_and_lrt_pvalue() {
    let fit = toy_fit(vec![1.0, 2.0], vec![0.5, 0.25], vec![true, true]);
    let lrt = LrtOutput {
        deviance: vec![Some(4.0), Some(1.0)],
        pvalue: vec![Some(0.04550026389635853), Some(0.31731050786291415)],
        degrees_of_freedom: 1,
        reduced_converged: vec![true, false],
    };
    let names = vec!["gene_a".to_string(), "gene_b".to_string()];

    let results = build_lrt_results(
        &[10.0, 20.0],
        &fit,
        &lrt,
        0,
        Some(&names),
        Some(&[0.1, 0.2]),
    )
    .unwrap();

    assert_eq!(results.rows[0].gene.as_deref(), Some("gene_a"));
    assert_eq!(results.rows[0].log2_fold_change, Some(1.0));
    assert_eq!(results.rows[0].lfc_se, Some(0.5));
    assert_eq!(results.rows[0].stat, Some(4.0));
    assert_eq!(results.rows[0].pvalue, Some(0.04550026389635853));
    assert_eq!(results.rows[0].dispersion, Some(0.1));
    assert!(results.rows[0].padj.unwrap() <= results.rows[1].padj.unwrap());
    assert_eq!(results.metadata.test_type, Some(TestType::Lrt));
    assert_eq!(results.metadata.result_name.as_deref(), Some("Intercept"));
    assert_eq!(
        results.metadata.comparison.as_deref(),
        Some("full model versus reduced model")
    );
    let metadata = results.deseq2_metadata();
    assert_eq!(
        metadata.columns[1].description,
        "log2 fold change (MLE): Intercept"
    );
    assert_eq!(
        metadata.columns[3].description,
        "LRT statistic: full model versus reduced model"
    );
    assert_eq!(
        metadata.columns[4].description,
        "LRT p-value: full model versus reduced model"
    );
}

#[test]
fn build_lrt_results_validates_dimensions() {
    let fit = toy_fit(vec![1.0, 2.0], vec![0.5, 0.25], vec![true, true]);
    let lrt = LrtOutput {
        deviance: vec![Some(4.0)],
        pvalue: vec![Some(0.04550026389635853), Some(0.31731050786291415)],
        degrees_of_freedom: 1,
        reduced_converged: vec![true, true],
    };
    assert!(build_lrt_results(&[10.0, 20.0], &fit, &lrt, 0, None, None).is_err());
    assert!(build_lrt_results(&[10.0, 20.0], &fit, &lrt, 1, None, None).is_err());
}

#[test]
fn build_lrt_results_rejects_invalid_optional_outputs() {
    let fit = toy_fit(vec![1.0], vec![0.5], vec![true]);
    let invalid_pvalue = LrtOutput {
        deviance: vec![Some(4.0)],
        pvalue: vec![Some(f64::NAN)],
        degrees_of_freedom: 1,
        reduced_converged: vec![true],
    };
    assert!(build_lrt_results(&[10.0], &fit, &invalid_pvalue, 0, None, None).is_err());

    let invalid_deviance = LrtOutput {
        deviance: vec![Some(f64::INFINITY)],
        pvalue: vec![Some(0.05)],
        degrees_of_freedom: 1,
        reduced_converged: vec![true],
    };
    assert!(build_lrt_results(&[10.0], &fit, &invalid_deviance, 0, None, None).is_err());

    let missing_reduced_flag = LrtOutput {
        deviance: vec![Some(4.0)],
        pvalue: vec![Some(0.05)],
        degrees_of_freedom: 1,
        reduced_converged: Vec::new(),
    };
    assert!(build_lrt_results(&[10.0], &fit, &missing_reduced_flag, 0, None, None).is_err());
}

#[test]
fn build_lrt_contrast_results_uses_contrast_effect_but_lrt_test_columns() {
    let fit = toy_fit(vec![1.0, 2.0], vec![0.5, 0.25], vec![true, false]);
    let lrt = LrtOutput {
        deviance: vec![Some(4.0), Some(1.0)],
        pvalue: vec![Some(0.04550026389635853), Some(0.31731050786291415)],
        degrees_of_freedom: 1,
        reduced_converged: vec![true, false],
    };
    let contrast = WaldContrastOutput {
        log2_fold_change: vec![Some(2.0), Some(-1.0)],
        lfc_se: vec![Some(0.8), Some(0.5)],
        wald: WaldOutput {
            stat: vec![Some(2.5), Some(-2.0)],
            pvalue: vec![
                Some(two_sided_normal_pvalue(2.5)),
                Some(two_sided_normal_pvalue(-2.0)),
            ],
            degrees_of_freedom: None,
        },
    };
    let names = vec!["gene_a".to_string(), "gene_b".to_string()];

    let results = build_lrt_contrast_results(
        &[10.0, 20.0],
        &fit,
        &lrt,
        &contrast,
        Some(&names),
        Some(&[0.1, 0.2]),
    )
    .unwrap();

    assert_eq!(results.rows[0].gene.as_deref(), Some("gene_a"));
    assert_eq!(results.rows[0].log2_fold_change, Some(2.0));
    assert_eq!(results.rows[0].lfc_se, Some(0.8));
    assert_eq!(results.rows[0].stat, Some(4.0));
    assert_eq!(results.rows[0].pvalue, Some(0.04550026389635853));
    assert_eq!(results.rows[0].dispersion, Some(0.1));
    assert_eq!(results.rows[1].log2_fold_change, Some(-1.0));
    assert_eq!(results.rows[1].stat, Some(1.0));
    assert!(results.rows[0].padj.unwrap() <= results.rows[1].padj.unwrap());
    assert_eq!(results.rows[1].converged, Some(false));
    assert_eq!(results.metadata.test_type, Some(TestType::Lrt));
    assert_eq!(results.metadata.result_name.as_deref(), Some("contrast"));
    assert_eq!(
        results.metadata.comparison.as_deref(),
        Some("primitive numeric contrast")
    );

    let metadata = results.deseq2_metadata();
    assert_eq!(
        metadata.columns[1].description,
        "log2 fold change (MLE): contrast"
    );
    assert_eq!(metadata.columns[2].description, "standard error: contrast");
    assert_eq!(
        metadata.columns[3].description,
        "LRT statistic: primitive numeric contrast"
    );
    assert_eq!(
        metadata.columns[4].description,
        "LRT p-value: primitive numeric contrast"
    );
}

#[test]
fn build_lrt_contrast_results_validates_lrt_and_contrast_outputs() {
    let fit = toy_fit(vec![1.0, 2.0], vec![0.5, 0.25], vec![true, true]);
    let valid_lrt = LrtOutput {
        deviance: vec![Some(4.0), Some(1.0)],
        pvalue: vec![Some(0.04550026389635853), Some(0.31731050786291415)],
        degrees_of_freedom: 1,
        reduced_converged: vec![true, true],
    };
    let valid_contrast = WaldContrastOutput {
        log2_fold_change: vec![Some(2.0), Some(-1.0)],
        lfc_se: vec![Some(0.8), Some(0.5)],
        wald: WaldOutput {
            stat: vec![Some(2.5), Some(-2.0)],
            pvalue: vec![Some(0.01), Some(0.02)],
            degrees_of_freedom: None,
        },
    };

    let bad_lrt = LrtOutput {
        deviance: vec![Some(f64::NAN), Some(1.0)],
        ..valid_lrt.clone()
    };
    assert!(
        build_lrt_contrast_results(&[10.0, 20.0], &fit, &bad_lrt, &valid_contrast, None, None,)
            .is_err()
    );

    let bad_contrast = WaldContrastOutput {
        lfc_se: vec![Some(0.8), Some(f64::INFINITY)],
        ..valid_contrast
    };
    assert!(
        build_lrt_contrast_results(&[10.0, 20.0], &fit, &valid_lrt, &bad_contrast, None, None,)
            .is_err()
    );
}

#[test]
fn default_cooks_cutoff_matches_deseq2_f_distribution_shape() {
    let cutoff = default_cooks_cutoff(3, 1).unwrap().unwrap();
    assert!(cutoff > 90.0);
    assert!(cutoff < 110.0);
    assert_eq!(default_cooks_cutoff(2, 2).unwrap(), None);
    assert!(default_cooks_cutoff(3, 0).is_err());
}

#[test]
fn resolve_cooks_cutoff_handles_disabled_threshold_and_invalid_values() {
    assert_eq!(
        resolve_cooks_cutoff(CooksCutoff::Disabled, 3, 1).unwrap(),
        None
    );
    assert_eq!(
        resolve_cooks_cutoff(CooksCutoff::Threshold(0.5), 3, 1).unwrap(),
        Some(0.5)
    );
    assert!(resolve_cooks_cutoff(CooksCutoff::Threshold(f64::NAN), 3, 1).is_err());
}

#[test]
fn apply_cooks_cutoff_masks_outlier_pvalues_and_recomputes_padj() {
    let fit = toy_fit(
        vec![0.0, 2.0, 1.0],
        vec![1.0, 0.5, 1.0],
        vec![true, true, true],
    );
    let mut results = build_wald_results(&[1.0, 2.0, 3.0], &fit, 0, None, None).unwrap();
    results.rows[0].max_cooks = Some(0.0);
    results.rows[1].max_cooks = Some(10.0);
    results.rows[2].max_cooks = None;

    apply_cooks_cutoff(&mut results, Some(5.0)).unwrap();

    assert!(results.column_names().contains(&"cooksOutlier"));
    assert_eq!(results.rows[0].cooks_outlier, Some(false));
    assert!(results.rows[0].pvalue.is_some());
    assert!(results.rows[0].padj.is_some());
    assert_eq!(results.rows[1].cooks_outlier, Some(true));
    assert_eq!(results.rows[1].pvalue, None);
    assert_eq!(results.rows[1].padj, None);
    assert_eq!(results.rows[2].cooks_outlier, None);
    assert!(results.rows[2].pvalue.is_some());
    assert!(results.rows[2].padj.is_some());
}

#[test]
fn recompute_padj_rejects_invalid_mutated_pvalues() {
    let fit = toy_fit(vec![0.0], vec![1.0], vec![true]);
    let mut results = build_wald_results(&[1.0], &fit, 0, None, None).unwrap();
    results.rows[0].pvalue = Some(1.2);

    assert!(recompute_padj(&mut results).is_err());
}

#[test]
fn apply_cooks_cutoff_rejects_nonfinite_max_cooks() {
    let fit = toy_fit(vec![0.0], vec![1.0], vec![true]);
    let mut results = build_wald_results(&[1.0], &fit, 0, None, None).unwrap();
    results.rows[0].max_cooks = Some(f64::INFINITY);

    assert!(apply_cooks_cutoff(&mut results, Some(5.0)).is_err());
}

#[test]
fn apply_cooks_cutoff_low_count_heuristic_spares_rows_with_three_larger_counts() {
    let fit = toy_fit(vec![0.0], vec![1.0], vec![true]);
    let mut results = build_wald_results(&[1.0], &fit, 0, None, None).unwrap();
    results.rows[0].max_cooks = Some(10.0);
    let original_pvalue = results.rows[0].pvalue;
    let counts = CountMatrix::from_row_major_u32(1, 4, vec![1, 5, 6, 7]).unwrap();
    let cooks = RowMajorMatrix::from_row_major(1, 4, vec![10.0, 0.1, 0.2, 0.3]).unwrap();

    apply_cooks_cutoff_with_low_count_heuristic(&mut results, Some(5.0), &counts, &cooks).unwrap();

    assert_eq!(results.rows[0].cooks_outlier, Some(false));
    assert_eq!(results.rows[0].pvalue, original_pvalue);
    assert!(results.rows[0].padj.is_some());
}

#[test]
fn apply_cooks_cutoff_low_count_heuristic_masks_when_outlier_count_is_not_low() {
    let fit = toy_fit(vec![0.0], vec![1.0], vec![true]);
    let mut results = build_wald_results(&[1.0], &fit, 0, None, None).unwrap();
    results.rows[0].max_cooks = Some(10.0);
    let counts = CountMatrix::from_row_major_u32(1, 4, vec![9, 5, 6, 7]).unwrap();
    let cooks = RowMajorMatrix::from_row_major(1, 4, vec![10.0, 0.1, 0.2, 0.3]).unwrap();

    apply_cooks_cutoff_with_low_count_heuristic(&mut results, Some(5.0), &counts, &cooks).unwrap();

    assert_eq!(results.rows[0].cooks_outlier, Some(true));
    assert_eq!(results.rows[0].pvalue, None);
    assert_eq!(results.rows[0].padj, None);
}

#[test]
fn apply_cooks_cutoff_low_count_heuristic_validates_inputs() {
    let fit = toy_fit(vec![0.0], vec![1.0], vec![true]);
    let mut results = build_wald_results(&[1.0], &fit, 0, None, None).unwrap();
    results.rows[0].max_cooks = Some(10.0);
    let counts = CountMatrix::from_row_major_u32(1, 4, vec![1, 5, 6, 7]).unwrap();
    let bad_cooks = RowMajorMatrix::from_row_major(1, 3, vec![10.0, 0.1, 0.2]).unwrap();

    assert!(apply_cooks_cutoff_with_low_count_heuristic(
        &mut results,
        Some(5.0),
        &counts,
        &bad_cooks,
    )
    .is_err());
    assert!(apply_cooks_cutoff_with_low_count_heuristic(
        &mut results,
        Some(f64::NAN),
        &counts,
        &RowMajorMatrix::from_row_major(1, 4, vec![10.0, 0.1, 0.2, 0.3]).unwrap(),
    )
    .is_err());
    results.rows[0].max_cooks = Some(f64::NAN);
    assert!(apply_cooks_cutoff_with_low_count_heuristic(
        &mut results,
        Some(5.0),
        &counts,
        &RowMajorMatrix::from_row_major(1, 4, vec![10.0, 0.1, 0.2, 0.3]).unwrap(),
    )
    .is_err());
}

#[test]
fn apply_cooks_cutoff_none_leaves_results_unchanged() {
    let fit = toy_fit(vec![2.0], vec![0.5], vec![true]);
    let mut results = build_wald_results(&[1.0], &fit, 0, None, None).unwrap();
    results.rows[0].max_cooks = Some(10.0);
    let before = results.clone();

    apply_cooks_cutoff(&mut results, None).unwrap();

    assert_eq!(results, before);
}
