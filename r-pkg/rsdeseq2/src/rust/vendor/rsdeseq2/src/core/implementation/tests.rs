#[cfg(test)]
mod tests {
    use super::{
        apply_cooks_cutoff_for_factor_level_metadata,
        factor_level_contrast_is_single_two_level_condition, full_deviance_from_log_like,
        CountMatrix,
    };
    use crate::contrasts::FactorLevelContrast;
    use crate::errors::DeseqError;
    use crate::matrix::RowMajorMatrix;
    use crate::results::{DeseqResultRow, DeseqResults};
    use std::assert_matches;

    #[test]
    fn count_matrix_rejects_bad_length() {
        let err = CountMatrix::from_row_major_u32(2, 3, vec![1, 2]).unwrap_err();
        assert!(err.to_string().contains("invalid dimensions"));
    }

    #[test]
    fn count_matrix_row_access() {
        let counts = CountMatrix::from_row_major_u32(2, 3, vec![1, 2, 3, 4, 5, 6]).unwrap();
        assert_eq!(counts.row(1).unwrap(), &[4, 5, 6]);
    }

    #[test]
    fn count_matrix_accepts_u64_when_values_fit() {
        let counts = CountMatrix::from_row_major_u64(1, 3, vec![1, 2, 3]).unwrap();
        assert_eq!(counts.as_slice(), &[1, 2, 3]);
    }

    #[test]
    fn count_matrix_index_spans_are_copy_and_reusable() {
        let counts = CountMatrix::from_row_major_u32(2, 3, vec![1, 2, 3, 4, 5, 6]).unwrap();
        let genes = counts.gene_indices();
        let samples = counts.sample_indices();

        let first_genes = genes.into_iter().collect::<Vec<_>>();
        let second_genes = genes.into_iter().collect::<Vec<_>>();
        let first_samples = samples.into_iter().collect::<Vec<_>>();
        let second_samples = samples.into_iter().collect::<Vec<_>>();

        assert_eq!(first_genes, vec![0, 1]);
        assert_eq!(second_genes, first_genes);
        assert_eq!(first_samples, vec![0, 1, 2]);
        assert_eq!(second_samples, first_samples);
    }

    #[test]
    fn count_matrix_gene_rows_accept_legacy_and_new_ranges() {
        let counts = CountMatrix::from_row_major_u32(3, 2, vec![1, 2, 3, 4, 5, 6]).unwrap();

        assert_eq!(counts.gene_rows(1..3).unwrap(), &[3, 4, 5, 6]);
        assert_eq!(
            counts
                .gene_rows(core::range::Range { start: 0, end: 2 })
                .unwrap(),
            &[1, 2, 3, 4]
        );
        assert_eq!(
            counts
                .gene_rows(core::range::RangeFrom { start: 2 })
                .unwrap(),
            &[5, 6]
        );
    }

    #[test]
    fn count_matrix_gene_rows_reports_range_errors_with_debuggable_match() {
        let counts = CountMatrix::from_row_major_u32(2, 2, vec![1, 2, 3, 4]).unwrap();
        let start = 2;
        let end = 1;

        assert_matches!(
            counts.gene_rows(start..end).unwrap_err(),
            DeseqError::InvalidDimensions { .. }
        );
        assert_matches!(
            counts
                .gene_rows(core::range::Range { start: 0, end: 3 })
                .unwrap_err(),
            DeseqError::InvalidDimensions { .. }
        );
    }

    #[test]
    fn all_zero_gene_detection() {
        let counts = CountMatrix::from_row_major_u32(2, 3, vec![0, 0, 0, 4, 5, 6]).unwrap();
        assert!(counts.is_all_zero_gene(0).unwrap());
        assert!(!counts.is_all_zero_gene(1).unwrap());
    }

    #[test]
    fn full_deviance_from_log_like_masks_nonfinite_values() {
        assert_eq!(full_deviance_from_log_like(-2.0), 4.0);
        assert!(full_deviance_from_log_like(f64::NAN).is_nan());
        assert!(full_deviance_from_log_like(f64::MAX).is_nan());
        assert!(full_deviance_from_log_like(-f64::MAX).is_nan());
    }

    #[test]
    fn factor_level_low_count_heuristic_requires_exact_two_levels() {
        let two_levels = ["A", "A", "B", "B"]
            .into_iter()
            .map(String::from)
            .collect::<Vec<_>>();
        assert!(factor_level_contrast_is_single_two_level_condition(
            FactorLevelContrast::new("condition", "B", "A", &two_levels)
        ));

        let missing_numerator = ["A", "A", "A"]
            .into_iter()
            .map(String::from)
            .collect::<Vec<_>>();
        assert!(!factor_level_contrast_is_single_two_level_condition(
            FactorLevelContrast::new("condition", "B", "A", &missing_numerator)
        ));

        let extra_level = ["A", "B", "C"]
            .into_iter()
            .map(String::from)
            .collect::<Vec<_>>();
        assert!(!factor_level_contrast_is_single_two_level_condition(
            FactorLevelContrast::new("condition", "B", "A", &extra_level)
        ));
    }

    #[test]
    fn factor_level_cooks_filter_uses_low_count_heuristic_only_for_two_levels() {
        let counts = CountMatrix::from_row_major_u32(1, 4, vec![1, 10, 11, 12]).unwrap();
        let cooks = RowMajorMatrix::from_row_major(1, 4, vec![10.0, 0.1, 0.2, 0.3]).unwrap();
        let two_levels = ["A", "A", "B", "B"]
            .into_iter()
            .map(String::from)
            .collect::<Vec<_>>();
        let mut two_level_results = one_cooks_result();
        apply_cooks_cutoff_for_factor_level_metadata(
            &mut two_level_results,
            Some(5.0),
            &counts,
            &cooks,
            FactorLevelContrast::new("condition", "B", "A", &two_levels),
        )
        .unwrap();
        assert_eq!(two_level_results.rows[0].pvalue, Some(0.01));
        assert_eq!(two_level_results.rows[0].cooks_outlier, Some(false));

        let extra_level = ["A", "A", "B", "C"]
            .into_iter()
            .map(String::from)
            .collect::<Vec<_>>();
        let mut extra_level_results = one_cooks_result();
        apply_cooks_cutoff_for_factor_level_metadata(
            &mut extra_level_results,
            Some(5.0),
            &counts,
            &cooks,
            FactorLevelContrast::new("condition", "B", "A", &extra_level),
        )
        .unwrap();
        assert_eq!(extra_level_results.rows[0].pvalue, None);
        assert_eq!(extra_level_results.rows[0].cooks_outlier, Some(true));
    }

    fn one_cooks_result() -> DeseqResults {
        DeseqResults {
            rows: vec![DeseqResultRow {
                gene: Some("gene1".to_string()),
                base_mean: 1.0,
                log2_fold_change: Some(0.0),
                lfc_se: Some(1.0),
                stat: Some(0.0),
                pvalue: Some(0.01),
                padj: Some(0.01),
                dispersion: Some(0.1),
                converged: Some(true),
                max_cooks: Some(10.0),
                cooks_outlier: None,
                filtered: None,
            }],
            ..DeseqResults::default()
        }
    }
}
