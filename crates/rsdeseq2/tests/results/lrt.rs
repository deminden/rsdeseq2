use super::*;

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
