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
fn build_wald_results_validates_dimensions() {
    let fit = toy_fit(vec![1.0, 2.0], vec![1.0, 1.0], vec![true, true]);
    assert!(build_wald_results(&[1.0], &fit, 0, None, None).is_err());

    let bad_names = vec!["gene_a".to_string()];
    assert!(build_wald_results(&[1.0, 2.0], &fit, 0, Some(&bad_names), None).is_err());

    assert!(build_wald_results(&[1.0, 2.0], &fit, 0, None, Some(&[0.1])).is_err());
    assert!(build_wald_results(&[1.0, 2.0], &fit, 1, None, None).is_err());
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
fn default_cooks_cutoff_matches_deseq2_f_distribution_shape() {
    let cutoff = default_cooks_cutoff(3, 1).unwrap().unwrap();
    assert!(cutoff > 90.0);
    assert!(cutoff < 110.0);
    assert_eq!(default_cooks_cutoff(2, 2).unwrap(), None);
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
