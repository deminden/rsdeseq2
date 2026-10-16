use super::*;

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

    assert!(
        apply_cooks_cutoff_with_low_count_heuristic(&mut results, Some(5.0), &counts, &bad_cooks,)
            .is_err()
    );
    assert!(
        apply_cooks_cutoff_with_low_count_heuristic(
            &mut results,
            Some(f64::NAN),
            &counts,
            &RowMajorMatrix::from_row_major(1, 4, vec![10.0, 0.1, 0.2, 0.3]).unwrap(),
        )
        .is_err()
    );
    results.rows[0].max_cooks = Some(f64::NAN);
    assert!(
        apply_cooks_cutoff_with_low_count_heuristic(
            &mut results,
            Some(5.0),
            &counts,
            &RowMajorMatrix::from_row_major(1, 4, vec![10.0, 0.1, 0.2, 0.3]).unwrap(),
        )
        .is_err()
    );
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
