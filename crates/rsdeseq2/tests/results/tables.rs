use super::*;

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
fn result_table_metadata_preserves_resolved_numeric_contrast() {
    let mut results = DeseqResults::default();
    results.set_resolved_contrast_metadata(
        "condition_B_vs_A",
        "factor-level contrast: condition B vs A",
        &[0.0, 1.0, -1.0],
    );

    assert_eq!(
        results.metadata.contrast.as_deref(),
        Some(&[0.0, 1.0, -1.0][..])
    );
    assert_eq!(
        results.metadata.scalar_metadata(),
        vec![
            DeseqResultsTableMetadataEntry {
                name: "resultName".to_string(),
                value: "condition_B_vs_A".to_string(),
            },
            DeseqResultsTableMetadataEntry {
                name: "comparison".to_string(),
                value: "factor-level contrast: condition B vs A".to_string(),
            },
            DeseqResultsTableMetadataEntry {
                name: "contrast".to_string(),
                value: "0,1,-1".to_string(),
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
