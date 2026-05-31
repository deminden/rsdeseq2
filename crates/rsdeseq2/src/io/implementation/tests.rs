#[cfg(test)]
mod tests {
    use super::*;
    use crate::diagnostics::Deseq2McolsDiagnostics;
    use crate::independent_filtering::IndependentFilteringOutput;
    use crate::results::{DeseqResultRow, DeseqResults};
    use std::fs;
    use std::path::PathBuf;

    #[test]
    fn parse_count_field_accepts_integer_scientific_notation() {
        let path = Path::new("counts.tsv");
        assert_eq!(parse_count_field("1e+05", path).unwrap(), 100_000);
        assert_eq!(parse_count_field("42.0", path).unwrap(), 42);
    }

    #[test]
    fn parse_count_field_rejects_non_integer_numeric_values() {
        let path = Path::new("counts.tsv");
        assert!(matches!(
            parse_count_field("1.5", path),
            Err(DeseqError::InvalidCounts { .. })
        ));
        assert!(matches!(
            parse_count_field("-1", path),
            Err(DeseqError::InvalidCounts { .. })
        ));
    }

    #[test]
    fn read_design_matrix_tsv_reads_numeric_matrix_and_names() {
        let path = unique_test_path("design.tsv");
        fs::write(
            &path,
            concat!(
                "sample\tIntercept\tcondition_B_vs_A\n",
                "s1\t1\t0\n",
                "s2\t1\t0\n",
                "s3\t1\t1\n",
            ),
        )
        .unwrap();

        let design = read_design_matrix_tsv(&path).unwrap();
        let _ = fs::remove_file(&path);

        assert_eq!(design.n_samples(), 3);
        assert_eq!(design.n_coefficients(), 2);
        assert_eq!(
            design.coefficient_names().unwrap(),
            &["Intercept".to_string(), "condition_B_vs_A".to_string()]
        );
        assert_eq!(design.matrix().as_slice(), &[1.0, 0.0, 1.0, 0.0, 1.0, 1.0]);
    }

    #[test]
    fn read_labeled_design_matrix_tsv_preserves_sample_names() {
        let path = unique_test_path("labeled_design.tsv");
        fs::write(
            &path,
            concat!(
                "sample\tIntercept\tcondition_B_vs_A\n",
                "s1\t1\t0\n",
                "s2\t1\t0\n",
                "s3\t1\t1\n",
            ),
        )
        .unwrap();

        let labeled = read_labeled_design_matrix_tsv(&path).unwrap();
        let _ = fs::remove_file(&path);

        assert_eq!(
            labeled.sample_names,
            vec!["s1".to_string(), "s2".to_string(), "s3".to_string()]
        );
        assert_eq!(
            labeled.design.matrix().as_slice(),
            &[1.0, 0.0, 1.0, 0.0, 1.0, 1.0]
        );
    }

    #[test]
    fn align_design_matrix_to_samples_uses_count_sample_order() {
        let design = DesignMatrix::from_row_major(
            3,
            2,
            vec![1.0, 1.0, 1.0, 0.0, 1.0, 2.0],
            Some(vec!["Intercept".to_string(), "condition".to_string()]),
        )
        .unwrap();
        let labeled = LabeledDesignMatrix {
            design,
            sample_names: vec!["s3".to_string(), "s1".to_string(), "s2".to_string()],
        };
        let samples = vec!["s1".to_string(), "s2".to_string(), "s3".to_string()];

        let aligned = align_design_matrix_to_samples(labeled, &samples).unwrap();

        assert_eq!(aligned.matrix().as_slice(), &[1.0, 0.0, 1.0, 2.0, 1.0, 1.0]);
        assert_eq!(
            aligned.coefficient_names().unwrap(),
            &["Intercept".to_string(), "condition".to_string()]
        );
    }

    #[test]
    fn align_design_matrix_to_samples_rejects_missing_and_duplicate_samples() {
        let samples = vec!["s1".to_string(), "s2".to_string()];
        let missing = LabeledDesignMatrix {
            design: DesignMatrix::from_row_major(2, 1, vec![1.0, 1.0], Some(vec!["x".to_string()]))
                .unwrap(),
            sample_names: vec!["s1".to_string(), "s3".to_string()],
        };
        assert!(align_design_matrix_to_samples(missing, &samples).is_err());

        let duplicated = LabeledDesignMatrix {
            design: DesignMatrix::from_row_major(2, 1, vec![1.0, 1.0], Some(vec!["x".to_string()]))
                .unwrap(),
            sample_names: vec!["s1".to_string(), "s1".to_string()],
        };
        assert!(align_design_matrix_to_samples(duplicated, &samples).is_err());
    }

    #[test]
    fn read_design_matrix_tsv_validates_shape_and_values() {
        let bad_value = unique_test_path("bad_design_value.tsv");
        fs::write(
            &bad_value,
            concat!("sample\tIntercept\tcondition_B_vs_A\n", "s1\t1\tNA\n",),
        )
        .unwrap();
        assert!(matches!(
            read_design_matrix_tsv(&bad_value),
            Err(DeseqError::ParseFloat { .. })
        ));
        let _ = fs::remove_file(&bad_value);

        let bad_shape = unique_test_path("bad_design_shape.tsv");
        fs::write(
            &bad_shape,
            concat!("sample\tIntercept\tcondition_B_vs_A\n", "s1\t1\n",),
        )
        .unwrap();
        assert!(read_design_matrix_tsv(&bad_shape).is_err());
        let _ = fs::remove_file(&bad_shape);
    }

    #[test]
    fn read_normalization_factors_tsv_reads_positive_matrix() {
        let path = unique_test_path("read_normalization_factors.tsv");
        fs::write(
            &path,
            concat!(
                "gene\tsample_1\tsample_2\n",
                "gene_a\t1\t2\n",
                "gene_b\t0.5\t4\n",
            ),
        )
        .unwrap();

        let factors = read_normalization_factors_tsv(&path).unwrap();
        let _ = fs::remove_file(&path);

        assert_eq!(factors.n_rows(), 2);
        assert_eq!(factors.n_cols(), 2);
        assert_eq!(factors.as_slice(), &[1.0, 2.0, 0.5, 4.0]);
    }

    #[test]
    fn read_labeled_normalization_factors_tsv_preserves_names() {
        let path = unique_test_path("read_labeled_normalization_factors.tsv");
        fs::write(
            &path,
            concat!(
                "gene\tsample_1\tsample_2\n",
                "gene_a\t1\t2\n",
                "gene_b\t0.5\t4\n",
            ),
        )
        .unwrap();

        let factors = read_labeled_normalization_factors_tsv(&path).unwrap();
        let _ = fs::remove_file(&path);

        assert_eq!(
            factors.gene_names,
            vec!["gene_a".to_string(), "gene_b".to_string()]
        );
        assert_eq!(
            factors.sample_names,
            vec!["sample_1".to_string(), "sample_2".to_string()]
        );
        assert_eq!(factors.matrix.as_slice(), &[1.0, 2.0, 0.5, 4.0]);
    }

    #[test]
    fn align_labeled_assay_matrix_to_counts_uses_gene_and_sample_order() {
        let counts = CountMatrix::from_row_major_u32_with_names(
            2,
            2,
            vec![1, 2, 3, 4],
            Some(vec!["gene_a".to_string(), "gene_b".to_string()]),
            Some(vec!["sample_1".to_string(), "sample_2".to_string()]),
        )
        .unwrap();
        let labeled = LabeledAssayMatrix {
            matrix: RowMajorMatrix::from_row_major(2, 2, vec![4.0, 3.0, 2.0, 1.0]).unwrap(),
            gene_names: vec!["gene_b".to_string(), "gene_a".to_string()],
            sample_names: vec!["sample_2".to_string(), "sample_1".to_string()],
        };

        let aligned =
            align_labeled_assay_matrix_to_counts(labeled, &counts, "test matrix").unwrap();

        assert_eq!(aligned.as_slice(), &[1.0, 2.0, 3.0, 4.0]);
    }

    #[test]
    fn align_labeled_assay_matrix_to_counts_rejects_missing_and_duplicate_names() {
        let counts = CountMatrix::from_row_major_u32_with_names(
            1,
            2,
            vec![1, 2],
            Some(vec!["gene_a".to_string()]),
            Some(vec!["sample_1".to_string(), "sample_2".to_string()]),
        )
        .unwrap();
        let missing = LabeledAssayMatrix {
            matrix: RowMajorMatrix::from_row_major(1, 2, vec![1.0, 2.0]).unwrap(),
            gene_names: vec!["gene_b".to_string()],
            sample_names: vec!["sample_1".to_string(), "sample_2".to_string()],
        };
        assert!(align_labeled_assay_matrix_to_counts(missing, &counts, "test matrix").is_err());

        let duplicated = LabeledAssayMatrix {
            matrix: RowMajorMatrix::from_row_major(1, 2, vec![1.0, 2.0]).unwrap(),
            gene_names: vec!["gene_a".to_string()],
            sample_names: vec!["sample_1".to_string(), "sample_1".to_string()],
        };
        assert!(align_labeled_assay_matrix_to_counts(duplicated, &counts, "test matrix").is_err());
    }

    #[test]
    fn read_normalization_factors_tsv_validates_shape_and_values() {
        let bad_value = unique_test_path("bad_normalization_factor_value.tsv");
        fs::write(
            &bad_value,
            concat!("gene\tsample_1\tsample_2\n", "gene_a\t1\t0\n",),
        )
        .unwrap();
        assert!(matches!(
            read_normalization_factors_tsv(&bad_value),
            Err(DeseqError::InvalidSizeFactors { .. })
        ));
        let _ = fs::remove_file(&bad_value);

        let bad_shape = unique_test_path("bad_normalization_factor_shape.tsv");
        fs::write(
            &bad_shape,
            concat!("gene\tsample_1\tsample_2\n", "gene_a\t1\n",),
        )
        .unwrap();
        assert!(read_normalization_factors_tsv(&bad_shape).is_err());
        let _ = fs::remove_file(&bad_shape);
    }

    #[test]
    fn read_observation_weights_tsv_reads_nonnegative_matrix() {
        let path = unique_test_path("read_observation_weights.tsv");
        fs::write(
            &path,
            concat!(
                "gene\tsample_1\tsample_2\n",
                "gene_a\t1\t0\n",
                "gene_b\t0.5\t4\n",
            ),
        )
        .unwrap();

        let weights = read_observation_weights_tsv(&path).unwrap();
        let _ = fs::remove_file(&path);

        assert_eq!(weights.n_rows(), 2);
        assert_eq!(weights.n_cols(), 2);
        assert_eq!(weights.as_slice(), &[1.0, 0.0, 0.5, 4.0]);
    }

    #[test]
    fn read_labeled_observation_weights_tsv_preserves_names() {
        let path = unique_test_path("read_labeled_observation_weights.tsv");
        fs::write(
            &path,
            concat!(
                "gene\tsample_1\tsample_2\n",
                "gene_a\t1\t0\n",
                "gene_b\t0.5\t4\n",
            ),
        )
        .unwrap();

        let weights = read_labeled_observation_weights_tsv(&path).unwrap();
        let _ = fs::remove_file(&path);

        assert_eq!(
            weights.gene_names,
            vec!["gene_a".to_string(), "gene_b".to_string()]
        );
        assert_eq!(
            weights.sample_names,
            vec!["sample_1".to_string(), "sample_2".to_string()]
        );
        assert_eq!(weights.matrix.as_slice(), &[1.0, 0.0, 0.5, 4.0]);
    }

    #[test]
    fn read_observation_weights_tsv_validates_shape_and_values() {
        let bad_value = unique_test_path("bad_observation_weight_value.tsv");
        fs::write(
            &bad_value,
            concat!("gene\tsample_1\tsample_2\n", "gene_a\t1\t-0.1\n",),
        )
        .unwrap();
        assert!(matches!(
            read_observation_weights_tsv(&bad_value),
            Err(DeseqError::NonFiniteValue { .. })
        ));
        let _ = fs::remove_file(&bad_value);

        let bad_shape = unique_test_path("bad_observation_weight_shape.tsv");
        fs::write(
            &bad_shape,
            concat!("gene\tsample_1\tsample_2\n", "gene_a\t1\n",),
        )
        .unwrap();
        assert!(read_observation_weights_tsv(&bad_shape).is_err());
        let _ = fs::remove_file(&bad_shape);
    }

    #[test]
    fn read_size_factors_tsv_reads_positive_values() {
        let path = unique_test_path("read_size_factors.tsv");
        fs::write(
            &path,
            concat!(
                "sample\tsize_factor\n",
                "sample_1\t1\n",
                "sample_2\t0.5\n",
                "sample_3\t2\n",
            ),
        )
        .unwrap();

        let size_factors = read_size_factors_tsv(&path).unwrap();
        let _ = fs::remove_file(&path);

        assert_eq!(size_factors, vec![1.0, 0.5, 2.0]);
    }

    #[test]
    fn read_labeled_size_factors_tsv_preserves_sample_names() {
        let path = unique_test_path("read_labeled_size_factors.tsv");
        fs::write(
            &path,
            concat!("sample\tsize_factor\n", "sample_1\t1\n", "sample_2\t0.5\n",),
        )
        .unwrap();

        let size_factors = read_labeled_size_factors_tsv(&path).unwrap();
        let _ = fs::remove_file(&path);

        assert_eq!(
            size_factors,
            vec![
                SampleNumericValue {
                    sample: "sample_1".to_string(),
                    value: 1.0
                },
                SampleNumericValue {
                    sample: "sample_2".to_string(),
                    value: 0.5
                },
            ]
        );
    }

    #[test]
    fn align_sample_numeric_values_to_samples_uses_count_sample_order() {
        let values = vec![
            SampleNumericValue {
                sample: "sample_3".to_string(),
                value: 3.0,
            },
            SampleNumericValue {
                sample: "sample_1".to_string(),
                value: 1.0,
            },
            SampleNumericValue {
                sample: "sample_2".to_string(),
                value: 2.0,
            },
        ];
        let samples = vec![
            "sample_1".to_string(),
            "sample_2".to_string(),
            "sample_3".to_string(),
        ];

        assert_eq!(
            align_sample_numeric_values_to_samples(&values, &samples, "size-factor").unwrap(),
            vec![1.0, 2.0, 3.0]
        );
    }

    #[test]
    fn align_sample_numeric_values_to_samples_rejects_missing_and_duplicate_samples() {
        let samples = vec!["sample_1".to_string(), "sample_2".to_string()];
        let missing = vec![
            SampleNumericValue {
                sample: "sample_1".to_string(),
                value: 1.0,
            },
            SampleNumericValue {
                sample: "sample_3".to_string(),
                value: 3.0,
            },
        ];
        assert!(align_sample_numeric_values_to_samples(&missing, &samples, "size-factor").is_err());

        let duplicated = vec![
            SampleNumericValue {
                sample: "sample_1".to_string(),
                value: 1.0,
            },
            SampleNumericValue {
                sample: "sample_1".to_string(),
                value: 2.0,
            },
        ];
        assert!(
            align_sample_numeric_values_to_samples(&duplicated, &samples, "size-factor").is_err()
        );
    }

    #[test]
    fn read_size_factors_tsv_validates_shape_and_values() {
        let bad_value = unique_test_path("bad_size_factor_value.tsv");
        fs::write(
            &bad_value,
            concat!("sample\tsize_factor\n", "sample_1\t0\n",),
        )
        .unwrap();
        assert!(matches!(
            read_size_factors_tsv(&bad_value),
            Err(DeseqError::InvalidSizeFactors { .. })
        ));
        let _ = fs::remove_file(&bad_value);

        let bad_shape = unique_test_path("bad_size_factor_shape.tsv");
        fs::write(
            &bad_shape,
            concat!("sample\tsize_factor\n", "sample_1\t1\textra\n",),
        )
        .unwrap();
        assert!(read_size_factors_tsv(&bad_shape).is_err());
        let _ = fs::remove_file(&bad_shape);
    }

    #[test]
    fn read_geometric_means_tsv_reads_nonnegative_values() {
        let path = unique_test_path("read_geometric_means.tsv");
        fs::write(
            &path,
            concat!(
                "gene\tgeo_mean\n",
                "gene_1\t1\n",
                "gene_2\t0\n",
                "gene_3\t2.5\n",
            ),
        )
        .unwrap();

        let geometric_means = read_geometric_means_tsv(&path).unwrap();
        let _ = fs::remove_file(&path);

        assert_eq!(geometric_means, vec![1.0, 0.0, 2.5]);
    }

    #[test]
    fn read_labeled_geometric_means_tsv_preserves_gene_names() {
        let path = unique_test_path("read_labeled_geometric_means.tsv");
        fs::write(
            &path,
            concat!("gene\tgeo_mean\n", "gene_1\t1\n", "gene_2\t0\n",),
        )
        .unwrap();

        let geometric_means = read_labeled_geometric_means_tsv(&path).unwrap();
        let _ = fs::remove_file(&path);

        assert_eq!(
            geometric_means,
            vec![
                GeneNumericValue {
                    gene: "gene_1".to_string(),
                    value: 1.0
                },
                GeneNumericValue {
                    gene: "gene_2".to_string(),
                    value: 0.0
                },
            ]
        );
    }

    #[test]
    fn align_gene_numeric_values_to_genes_uses_count_gene_order() {
        let values = vec![
            GeneNumericValue {
                gene: "gene_3".to_string(),
                value: 3.0,
            },
            GeneNumericValue {
                gene: "gene_1".to_string(),
                value: 1.0,
            },
            GeneNumericValue {
                gene: "gene_2".to_string(),
                value: 2.0,
            },
        ];
        let genes = vec![
            "gene_1".to_string(),
            "gene_2".to_string(),
            "gene_3".to_string(),
        ];

        assert_eq!(
            align_gene_numeric_values_to_genes(&values, &genes, "geometric-mean").unwrap(),
            vec![1.0, 2.0, 3.0]
        );
    }

    #[test]
    fn align_gene_numeric_values_to_genes_rejects_missing_and_duplicate_genes() {
        let genes = vec!["gene_1".to_string(), "gene_2".to_string()];
        let missing = vec![
            GeneNumericValue {
                gene: "gene_1".to_string(),
                value: 1.0,
            },
            GeneNumericValue {
                gene: "gene_3".to_string(),
                value: 3.0,
            },
        ];
        assert!(align_gene_numeric_values_to_genes(&missing, &genes, "geometric-mean").is_err());

        let duplicated = vec![
            GeneNumericValue {
                gene: "gene_1".to_string(),
                value: 1.0,
            },
            GeneNumericValue {
                gene: "gene_1".to_string(),
                value: 2.0,
            },
        ];
        assert!(align_gene_numeric_values_to_genes(&duplicated, &genes, "geometric-mean").is_err());
    }

    #[test]
    fn read_geometric_means_tsv_validates_shape_and_values() {
        let bad_value = unique_test_path("bad_geometric_mean_value.tsv");
        fs::write(&bad_value, concat!("gene\tgeo_mean\n", "gene_1\t-1\n")).unwrap();
        assert!(matches!(
            read_geometric_means_tsv(&bad_value),
            Err(DeseqError::InvalidSizeFactors { .. })
        ));
        let _ = fs::remove_file(&bad_value);

        let bad_shape = unique_test_path("bad_geometric_mean_shape.tsv");
        fs::write(
            &bad_shape,
            concat!("gene\tgeo_mean\n", "gene_1\t1\textra\n"),
        )
        .unwrap();
        assert!(read_geometric_means_tsv(&bad_shape).is_err());
        let _ = fs::remove_file(&bad_shape);
    }

    #[test]
    fn read_wald_t_degrees_of_freedom_tsv_reads_finite_values() {
        let path = unique_test_path("read_wald_t_df.tsv");
        fs::write(&path, concat!("gene\tdf\n", "gene_1\t4\n", "gene_2\t2.5\n")).unwrap();

        let degrees_of_freedom = read_wald_t_degrees_of_freedom_tsv(&path).unwrap();
        let _ = fs::remove_file(&path);

        assert_eq!(degrees_of_freedom, vec![4.0, 2.5]);
    }

    #[test]
    fn read_labeled_wald_t_degrees_of_freedom_tsv_preserves_gene_names() {
        let path = unique_test_path("read_labeled_wald_t_df.tsv");
        fs::write(&path, concat!("gene\tdf\n", "gene_1\t4\n", "gene_2\t2.5\n")).unwrap();

        let degrees_of_freedom = read_labeled_wald_t_degrees_of_freedom_tsv(&path).unwrap();
        let _ = fs::remove_file(&path);

        assert_eq!(
            degrees_of_freedom,
            vec![
                GeneNumericValue {
                    gene: "gene_1".to_string(),
                    value: 4.0
                },
                GeneNumericValue {
                    gene: "gene_2".to_string(),
                    value: 2.5
                },
            ]
        );
    }

    #[test]
    fn read_wald_t_degrees_of_freedom_tsv_validates_shape_and_values() {
        let bad_value = unique_test_path("bad_wald_t_df_value.tsv");
        fs::write(&bad_value, concat!("gene\tdf\n", "gene_1\tNA\n")).unwrap();
        assert!(read_wald_t_degrees_of_freedom_tsv(&bad_value).is_err());
        let _ = fs::remove_file(&bad_value);

        let bad_shape = unique_test_path("bad_wald_t_df_shape.tsv");
        fs::write(&bad_shape, concat!("gene\tdf\n", "gene_1\t4\textra\n")).unwrap();
        assert!(read_wald_t_degrees_of_freedom_tsv(&bad_shape).is_err());
        let _ = fs::remove_file(&bad_shape);
    }

    #[test]
    fn read_sample_levels_tsv_reads_labeled_string_levels() {
        let path = unique_test_path("read_sample_levels.tsv");
        fs::write(
            &path,
            concat!(
                "sample\tcondition\n",
                "sample_1\tA\n",
                "sample_2\tB\n",
                "sample_3\tA\n",
            ),
        )
        .unwrap();

        let levels = read_sample_levels_tsv(&path).unwrap();
        let _ = fs::remove_file(&path);

        assert_eq!(
            levels,
            vec![
                SampleLevel {
                    sample: "sample_1".to_string(),
                    level: "A".to_string()
                },
                SampleLevel {
                    sample: "sample_2".to_string(),
                    level: "B".to_string()
                },
                SampleLevel {
                    sample: "sample_3".to_string(),
                    level: "A".to_string()
                },
            ]
        );
    }

    #[test]
    fn read_sample_levels_tsv_validates_shape_and_values() {
        let bad_value = unique_test_path("bad_sample_level_value.tsv");
        fs::write(&bad_value, concat!("sample\tcondition\n", "sample_1\t\n")).unwrap();
        assert!(read_sample_levels_tsv(&bad_value).is_err());
        let _ = fs::remove_file(&bad_value);

        let bad_shape = unique_test_path("bad_sample_level_shape.tsv");
        fs::write(
            &bad_shape,
            concat!("sample\tcondition\n", "sample_1\tA\textra\n"),
        )
        .unwrap();
        assert!(read_sample_levels_tsv(&bad_shape).is_err());
        let _ = fs::remove_file(&bad_shape);
    }

    #[test]
    fn align_sample_levels_to_samples_uses_count_sample_order() {
        let levels = vec![
            SampleLevel {
                sample: "sample_3".to_string(),
                level: "C".to_string(),
            },
            SampleLevel {
                sample: "sample_1".to_string(),
                level: "A".to_string(),
            },
            SampleLevel {
                sample: "sample_2".to_string(),
                level: "B".to_string(),
            },
        ];
        let samples = vec![
            "sample_1".to_string(),
            "sample_2".to_string(),
            "sample_3".to_string(),
        ];

        assert_eq!(
            align_sample_levels_to_samples(&levels, &samples).unwrap(),
            vec!["A".to_string(), "B".to_string(), "C".to_string()]
        );
    }

    #[test]
    fn align_sample_levels_to_samples_rejects_missing_and_duplicate_samples() {
        let samples = vec!["sample_1".to_string(), "sample_2".to_string()];
        let missing = vec![
            SampleLevel {
                sample: "sample_1".to_string(),
                level: "A".to_string(),
            },
            SampleLevel {
                sample: "sample_3".to_string(),
                level: "C".to_string(),
            },
        ];
        assert!(align_sample_levels_to_samples(&missing, &samples).is_err());

        let duplicated = vec![
            SampleLevel {
                sample: "sample_1".to_string(),
                level: "A".to_string(),
            },
            SampleLevel {
                sample: "sample_1".to_string(),
                level: "B".to_string(),
            },
        ];
        assert!(align_sample_levels_to_samples(&duplicated, &samples).is_err());
    }

    #[test]
    fn write_count_matrices_tsv_preserves_names_and_fallbacks() {
        let raw_path = unique_test_path("counts.tsv");
        let normalized_path = unique_test_path("normalized_counts.tsv");
        let unnamed_path = unique_test_path("unnamed_normalized_counts.tsv");
        let counts = CountMatrix::from_row_major_u32_with_names(
            2,
            3,
            vec![10, 20, 30, 5, 10, 20],
            Some(vec!["gene_a".to_string(), "gene_b".to_string()]),
            Some(vec![
                "sample_1".to_string(),
                "sample_2".to_string(),
                "sample_3".to_string(),
            ]),
        )
        .unwrap();
        let normalized =
            RowMajorMatrix::from_row_major(2, 3, vec![10.0, 10.0, 6.0, 5.0, 5.0, 4.0]).unwrap();

        write_count_matrix_tsv(&raw_path, &counts).unwrap();
        write_normalized_counts_tsv(
            &normalized_path,
            counts.gene_names(),
            counts.sample_names(),
            &normalized,
        )
        .unwrap();
        write_normalized_counts_tsv(&unnamed_path, None, None, &normalized).unwrap();

        let raw = fs::read_to_string(&raw_path).unwrap();
        let normalized_text = fs::read_to_string(&normalized_path).unwrap();
        let unnamed = fs::read_to_string(&unnamed_path).unwrap();
        let _ = fs::remove_file(&raw_path);
        let _ = fs::remove_file(&normalized_path);
        let _ = fs::remove_file(&unnamed_path);

        assert_eq!(
            raw,
            concat!(
                "gene\tsample_1\tsample_2\tsample_3\n",
                "gene_a\t10\t20\t30\n",
                "gene_b\t5\t10\t20\n",
            )
        );
        assert_eq!(
            normalized_text,
            concat!(
                "gene\tsample_1\tsample_2\tsample_3\n",
                "gene_a\t10\t10\t6\n",
                "gene_b\t5\t5\t4\n",
            )
        );
        assert_eq!(
            unnamed,
            concat!(
                "gene\tsample1\tsample2\tsample3\n",
                "gene1\t10\t10\t6\n",
                "gene2\t5\t5\t4\n",
            )
        );
    }

    #[test]
    fn write_normalized_counts_tsv_validates_names_and_finite_values() {
        let path = unique_test_path("bad_normalized_counts.tsv");
        let normalized = RowMajorMatrix::from_row_major(1, 2, vec![1.0, f64::INFINITY]).unwrap();

        assert!(matches!(
            write_normalized_counts_tsv(&path, Some(&["gene_a".to_string()]), None, &normalized),
            Err(DeseqError::NonFiniteValue { .. })
        ));
        assert!(matches!(
            write_normalized_counts_tsv(
                &path,
                Some(&["gene_a".to_string(), "gene_b".to_string()]),
                None,
                &RowMajorMatrix::from_row_major(1, 2, vec![1.0, 2.0]).unwrap(),
            ),
            Err(DeseqError::InvalidDimensions { .. })
        ));
        let _ = fs::remove_file(&path);
    }

    #[test]
    fn write_normalization_factors_tsv_writes_deseq2_shape_and_validates_values() {
        let path = unique_test_path("normalization_factors.tsv");
        let bad_path = unique_test_path("bad_normalization_factors.tsv");
        let genes = vec!["gene_a".to_string(), "gene_b".to_string()];
        let samples = vec!["sample_1".to_string(), "sample_2".to_string()];
        let factors = RowMajorMatrix::from_row_major(2, 2, vec![1.0, 2.0, 0.5, 4.0]).unwrap();

        write_normalization_factors_tsv(&path, Some(&genes), Some(&samples), &factors).unwrap();

        let text = fs::read_to_string(&path).unwrap();
        let _ = fs::remove_file(&path);
        assert_eq!(
            text,
            concat!(
                "gene\tsample_1\tsample_2\n",
                "gene_a\t1\t2\n",
                "gene_b\t0.5\t4\n",
            )
        );
        assert!(matches!(
            write_normalization_factors_tsv(
                &bad_path,
                None,
                None,
                &RowMajorMatrix::from_row_major(1, 2, vec![1.0, 0.0]).unwrap(),
            ),
            Err(DeseqError::InvalidSizeFactors { .. })
        ));
        assert!(matches!(
            write_normalization_factors_tsv(
                &bad_path,
                Some(&genes),
                None,
                &RowMajorMatrix::from_row_major(1, 2, vec![1.0, 2.0]).unwrap(),
            ),
            Err(DeseqError::InvalidDimensions { .. })
        ));
        let _ = fs::remove_file(&bad_path);
    }

    #[test]
    fn write_size_factors_tsv_validates_names_and_values() {
        let path = unique_test_path("bad_size_factors.tsv");
        let sample_names = vec!["sample_1".to_string(), "sample_2".to_string()];

        assert!(matches!(
            write_size_factors_tsv(&path, Some(&sample_names), &[1.0]),
            Err(DeseqError::InvalidDimensions { .. })
        ));
        assert!(matches!(
            write_size_factors_tsv(&path, None, &[1.0, 0.0]),
            Err(DeseqError::InvalidSizeFactors { .. })
        ));
        assert!(matches!(
            write_size_factors_tsv(&path, None, &[1.0, f64::NAN]),
            Err(DeseqError::InvalidSizeFactors { .. })
        ));
        let _ = fs::remove_file(&path);
    }

    #[test]
    fn write_base_mean_tsv_validates_names_and_values() {
        let path = unique_test_path("bad_base_mean.tsv");
        let gene_names = vec!["gene_a".to_string(), "gene_b".to_string()];

        assert!(matches!(
            write_base_mean_tsv(&path, Some(&gene_names), &[1.0]),
            Err(DeseqError::InvalidDimensions { .. })
        ));
        assert!(matches!(
            write_base_mean_tsv(&path, None, &[1.0, -1.0]),
            Err(DeseqError::NonFiniteValue { .. })
        ));
        assert!(matches!(
            write_base_mean_tsv(&path, None, &[1.0, f64::INFINITY]),
            Err(DeseqError::NonFiniteValue { .. })
        ));
        let _ = fs::remove_file(&path);
    }

    #[test]
    fn write_base_metadata_tsv_writes_deseq2_shape_and_validates_lengths() {
        let path = unique_test_path("base_metadata.tsv");
        let bad_path = unique_test_path("bad_base_metadata.tsv");
        let genes = vec!["gene_a".to_string(), "gene_b".to_string()];

        write_base_metadata_tsv(
            &path,
            Some(&genes),
            &[10.0, f64::NAN],
            &[2.5, f64::INFINITY],
            &[false, true],
        )
        .unwrap();

        let text = fs::read_to_string(&path).unwrap();
        let _ = fs::remove_file(&path);
        assert_eq!(
            text,
            concat!(
                "gene\tbaseMean\tbaseVar\tallZero\n",
                "gene_a\t10\t2.5\tFALSE\n",
                "gene_b\tNA\tNA\tTRUE\n",
            )
        );
        assert!(matches!(
            write_base_metadata_tsv(&bad_path, Some(&genes), &[1.0], &[2.0], &[false]),
            Err(DeseqError::InvalidDimensions { .. })
        ));
        assert!(matches!(
            write_base_metadata_tsv(&bad_path, None, &[1.0, 2.0], &[3.0], &[false, true]),
            Err(DeseqError::InvalidDimensions { .. })
        ));
        let _ = fs::remove_file(&bad_path);
    }

    #[test]
    fn write_deseq_results_tsv_writes_columns_missing_values_and_logicals() {
        let path = unique_test_path("deseq_results.tsv");
        let tidy_path = unique_test_path("deseq_results_tidy.tsv");
        let metadata_path = unique_test_path("deseq_result_column_metadata.tsv");
        let table_metadata_path = unique_test_path("deseq_result_table_metadata.tsv");
        let results = DeseqResults {
            rows: vec![
                DeseqResultRow {
                    gene: Some("gene_a".to_string()),
                    base_mean: 10.0,
                    log2_fold_change: Some(1.25),
                    lfc_se: None,
                    stat: Some(-2.0),
                    pvalue: None,
                    padj: Some(0.05),
                    dispersion: Some(0.1),
                    converged: Some(true),
                    max_cooks: Some(4.0),
                    cooks_outlier: Some(false),
                    filtered: None,
                },
                DeseqResultRow {
                    gene: None,
                    base_mean: 20.0,
                    log2_fold_change: None,
                    lfc_se: Some(0.5),
                    stat: None,
                    pvalue: Some(0.8),
                    padj: None,
                    dispersion: None,
                    converged: Some(false),
                    max_cooks: None,
                    cooks_outlier: None,
                    filtered: None,
                },
            ],
            metadata: crate::results::DeseqResultsTableMetadata {
                test_type: Some(crate::options::TestType::Wald),
                result_name: Some("condition_B_vs_A".to_string()),
                comparison: Some("coefficient condition_B_vs_A".to_string()),
                lfc_threshold: 1.5,
                alt_hypothesis: Some("greater".to_string()),
                ..crate::results::DeseqResultsTableMetadata::default()
            },
            ..DeseqResults::default()
        };

        write_deseq_results_tsv(&path, &results).unwrap();
        write_deseq_results_tidy_tsv(&tidy_path, &results).unwrap();
        write_deseq_result_column_metadata_tsv(&metadata_path, &results).unwrap();
        write_deseq_result_table_metadata_tsv(&table_metadata_path, &results).unwrap();

        let text = fs::read_to_string(&path).unwrap();
        let tidy_text = fs::read_to_string(&tidy_path).unwrap();
        let metadata = fs::read_to_string(&metadata_path).unwrap();
        let table_metadata = fs::read_to_string(&table_metadata_path).unwrap();
        let _ = fs::remove_file(&path);
        let _ = fs::remove_file(&tidy_path);
        let _ = fs::remove_file(&metadata_path);
        let _ = fs::remove_file(&table_metadata_path);
        assert_eq!(
            text,
            concat!(
                "gene\tbaseMean\tlog2FoldChange\tlfcSE\tstat\tpvalue\tpadj\t",
                "dispersion\tconverged\tmaxCooks\tcooksOutlier\n",
                "gene_a\t10\t1.25\tNA\t-2\tNA\t0.05\t0.1\tTRUE\t4\tFALSE\n",
                "gene2\t20\tNA\t0.5\tNA\t0.8\tNA\tNA\tFALSE\tNA\tNA\n"
            )
        );
        assert_eq!(
            tidy_text,
            concat!(
                "row\tbaseMean\tlog2FoldChange\tlfcSE\tstat\tpvalue\tpadj\t",
                "dispersion\tconverged\tmaxCooks\tcooksOutlier\n",
                "gene_a\t10\t1.25\tNA\t-2\tNA\t0.05\t0.1\tTRUE\t4\tFALSE\n",
                "row2\t20\tNA\t0.5\tNA\t0.8\tNA\tNA\tFALSE\tNA\tNA\n"
            )
        );
        assert!(metadata.starts_with("name\ttype\tdescription\n"));
        assert!(metadata.contains("baseMean\tresults\tmean of normalized counts for all samples\n"));
        assert!(metadata.contains("converged\tdiagnostic\twhether beta fitting converged\n"));
        assert_eq!(
            table_metadata,
            concat!(
                "name\tvalue\n",
                "testType\tWald\n",
                "resultName\tcondition_B_vs_A\n",
                "comparison\tcoefficient condition_B_vs_A\n",
                "lfcThreshold\t1.5\n",
                "altHypothesis\tgreater\n",
                "pAdjustMethod\tBH\n",
            )
        );
    }

    #[test]
    fn write_deseq_results_tsv_rejects_invalid_numeric_rows() {
        let path = unique_test_path("bad_deseq_results.tsv");
        let mut results = DeseqResults {
            rows: vec![DeseqResultRow {
                gene: Some("gene_a".to_string()),
                base_mean: 10.0,
                log2_fold_change: Some(1.25),
                lfc_se: Some(0.5),
                stat: Some(-2.0),
                pvalue: Some(0.04),
                padj: Some(0.05),
                dispersion: Some(0.1),
                converged: Some(true),
                max_cooks: Some(4.0),
                cooks_outlier: Some(false),
                filtered: None,
            }],
            ..DeseqResults::default()
        };

        results.rows[0].pvalue = Some(1.2);
        assert!(write_deseq_results_tsv(&path, &results).is_err());
        results.rows[0].pvalue = Some(0.04);
        results.rows[0].max_cooks = Some(f64::NAN);
        assert!(write_deseq_results_tsv(&path, &results).is_err());
        results.rows[0].max_cooks = Some(4.0);
        results.rows[0].base_mean = -1.0;
        assert!(write_deseq_results_tsv(&path, &results).is_err());
        let _ = fs::remove_file(&path);
    }

    #[test]
    fn write_deseq_mcols_diagnostics_tsv_writes_present_columns_and_na_values() {
        let path = unique_test_path("mcols_diagnostics.tsv");
        let gene_names = vec!["gene_a".to_string(), "gene_b".to_string()];
        let diagnostics = Deseq2McolsDiagnostics {
            disp_gene_est: Some(vec![f64::NAN, 0.2]),
            disp_gene_iter: Some(vec![0, 3]),
            disp_map: Some(vec![f64::NAN, 0.25]),
            disp_outlier: Some(vec![false, true]),
            beta_optim_iter: Some(vec![f64::NAN, 8.0]),
            beta_optim_start_objective: Some(vec![f64::NAN, 12.5]),
            beta_optim_objective: Some(vec![f64::NAN, 4.25]),
            beta_optim_gradient_norm: Some(vec![f64::NAN, 1e-9]),
            max_cooks: Some(vec![None, Some(4.5)]),
            ..Deseq2McolsDiagnostics::default()
        };

        write_deseq_mcols_diagnostics_tsv(&path, Some(&gene_names), &diagnostics).unwrap();

        let text = fs::read_to_string(&path).unwrap();
        let _ = fs::remove_file(&path);
        assert_eq!(
            text,
            concat!(
                "gene\tdispGeneEst\tdispGeneIter\tdispMAP\tdispOutlier\t",
                "rustBetaOptimIter\trustBetaOptimStartObjective\t",
                "rustBetaOptimObjective\trustBetaOptimGradientNorm\tmaxCooks\n",
                "gene_a\tNA\t0\tNA\tFALSE\tNA\tNA\tNA\tNA\tNA\n",
                "gene_b\t0.2\t3\t0.25\tTRUE\t8\t12.5\t4.25\t0.000000001\t4.5\n"
            )
        );
    }

    #[test]
    fn write_cooks_diagnostics_tsv_writes_matrix_and_metadata() {
        let cooks_path = unique_test_path("cooks_distance.tsv");
        let row_path = unique_test_path("cooks_row_metadata.tsv");
        let sample_path = unique_test_path("cooks_sample_metadata.tsv");
        let gene_names = vec!["gene_a".to_string(), "gene_b".to_string()];
        let sample_names = vec!["sample_1".to_string(), "sample_2".to_string()];
        let cooks = CooksOutput {
            cooks: RowMajorMatrix::from_row_major(2, 2, vec![0.5, f64::NAN, 2.0, 4.0]).unwrap(),
            max_cooks: vec![None, Some(4.0)],
            robust_dispersion: vec![0.04, 0.25],
            samples_for_cooks: vec![true, false],
        };

        write_cooks_distance_tsv(&cooks_path, Some(&gene_names), Some(&sample_names), &cooks)
            .unwrap();
        write_cooks_row_metadata_tsv(&row_path, Some(&gene_names), &cooks).unwrap();
        write_cooks_sample_metadata_tsv(&sample_path, Some(&sample_names), &cooks).unwrap();

        let cooks_text = fs::read_to_string(&cooks_path).unwrap();
        let row_text = fs::read_to_string(&row_path).unwrap();
        let sample_text = fs::read_to_string(&sample_path).unwrap();
        let _ = fs::remove_file(&cooks_path);
        let _ = fs::remove_file(&row_path);
        let _ = fs::remove_file(&sample_path);

        assert_eq!(
            cooks_text,
            concat!(
                "gene\tsample_1\tsample_2\n",
                "gene_a\t0.5\tNA\n",
                "gene_b\t2\t4\n",
            )
        );
        assert_eq!(
            row_text,
            concat!(
                "gene\tmaxCooks\tcooksRobustDispersion\n",
                "gene_a\tNA\t0.04\n",
                "gene_b\t4\t0.25\n",
            )
        );
        assert_eq!(
            sample_text,
            concat!(
                "sample\tsamplesForCooks\n",
                "sample_1\tTRUE\n",
                "sample_2\tFALSE\n",
            )
        );
    }

    #[test]
    fn write_cooks_diagnostics_tsv_rejects_misaligned_metadata() {
        let path = unique_test_path("bad_cooks_metadata.tsv");
        let mut cooks = CooksOutput {
            cooks: RowMajorMatrix::from_row_major(2, 2, vec![0.5, 1.0, 2.0, 4.0]).unwrap(),
            max_cooks: vec![Some(1.0)],
            robust_dispersion: vec![0.04, 0.25],
            samples_for_cooks: vec![true, false],
        };

        assert!(write_cooks_row_metadata_tsv(&path, None, &cooks).is_err());
        cooks.max_cooks = vec![Some(1.0), Some(4.0)];
        cooks.robust_dispersion = vec![0.04];
        assert!(write_cooks_row_metadata_tsv(&path, None, &cooks).is_err());
        cooks.robust_dispersion = vec![0.04, 0.25];
        cooks.samples_for_cooks = vec![true];
        assert!(write_cooks_sample_metadata_tsv(&path, None, &cooks).is_err());
        let _ = fs::remove_file(&path);
    }

    #[test]
    fn write_cooks_replacement_metadata_tsv_writes_scalar_summary() {
        let path = unique_test_path("cooks_replacement_metadata.tsv");
        let replaced_path = unique_test_path("cooks_replaced_counts.tsv");
        let candidate_path = unique_test_path("cooks_candidate_counts.tsv");
        let outlier_path = unique_test_path("cooks_outlier_cells.tsv");
        let row_metadata_path = unique_test_path("cooks_replacement_row_metadata.tsv");
        let counts = CountMatrix::from_row_major_u32(
            2,
            4,
            vec![
                10, 30, 20, 50, //
                0, 0, 0, 0,
            ],
        )
        .unwrap();
        let size_factors = vec![1.0, 1.0, 1.0, 1.0];
        let normalized = crate::normalization::normalized_counts(&counts, &size_factors).unwrap();
        let cooks = RowMajorMatrix::from_row_major(
            2,
            4,
            vec![
                0.0, 9.0, 0.0, 0.0, //
                8.0, 0.0, 0.0, 0.0,
            ],
        )
        .unwrap();
        let design = DesignMatrix::from_row_major(4, 1, vec![1.0, 1.0, 1.0, 1.0], None).unwrap();
        let plan = crate::cooks::prepare_cooks_replacement_refit(
            &counts,
            &normalized,
            &size_factors,
            None,
            &cooks,
            &design,
            &crate::cooks::CooksReplacementOptions {
                trim: 0.25,
                cooks_cutoff: 5.0,
                min_replicates: 3,
                which_samples: None,
            },
        )
        .unwrap();

        write_cooks_replacement_metadata_tsv(&path, &plan).unwrap();
        write_cooks_replaced_counts_tsv(&replaced_path, &plan).unwrap();
        write_cooks_candidate_replacement_counts_tsv(&candidate_path, &plan).unwrap();
        write_cooks_outlier_cells_tsv(&outlier_path, &plan).unwrap();
        write_cooks_replacement_row_metadata_tsv(&row_metadata_path, &plan).unwrap();

        let text = fs::read_to_string(&path).unwrap();
        let replaced = fs::read_to_string(&replaced_path).unwrap();
        let candidate = fs::read_to_string(&candidate_path).unwrap();
        let outlier = fs::read_to_string(&outlier_path).unwrap();
        let row_metadata = fs::read_to_string(&row_metadata_path).unwrap();
        let _ = fs::remove_file(&path);
        let _ = fs::remove_file(&replaced_path);
        let _ = fs::remove_file(&candidate_path);
        let _ = fs::remove_file(&outlier_path);
        let _ = fs::remove_file(&row_metadata_path);
        assert_eq!(
            text,
            concat!(
                "name\tvalue\n",
                "nRefit\t2\n",
                "nRefitRows\t1\n",
                "nNewAllZero\t1\n",
                "nOutlierCells\t2\n",
                "nReplacedCells\t2\n",
                "nReplaceableSamples\t4\n",
                "shouldRefit\ttrue\n",
            )
        );
        assert_eq!(
            replaced,
            concat!(
                "gene\tsample1\tsample2\tsample3\tsample4\n",
                "gene1\t10\t25\t20\t50\n",
                "gene2\t0\t0\t0\t0\n",
            )
        );
        assert_eq!(
            candidate,
            concat!(
                "gene\tsample1\tsample2\tsample3\tsample4\n",
                "gene1\t25\t25\t25\t25\n",
                "gene2\t0\t0\t0\t0\n",
            )
        );
        assert_eq!(
            outlier,
            concat!(
                "gene\tsample1\tsample2\tsample3\tsample4\n",
                "gene1\tFALSE\tTRUE\tFALSE\tFALSE\n",
                "gene2\tTRUE\tFALSE\tFALSE\tFALSE\n",
            )
        );
        assert_eq!(
            row_metadata,
            concat!(
                "gene\treplace\trefitReplace\tnewAllZero\treplacedAllZero\t",
                "replacedBaseMean\treplacedBaseVar\tpostRefitMaxCooks\n",
                "gene1\tTRUE\tTRUE\tFALSE\tFALSE\t26.25\t289.58333333333337\tNA\n",
                "gene2\tTRUE\tFALSE\tTRUE\tTRUE\t0\t0\tNA\n",
            )
        );
    }

    #[test]
    fn write_deseq_mcols_diagnostics_tsv_rejects_misaligned_columns() {
        let path = unique_test_path("mcols_bad_diagnostics.tsv");
        let diagnostics = Deseq2McolsDiagnostics {
            disp_gene_est: Some(vec![0.1, 0.2]),
            disp_gene_iter: Some(vec![1]),
            ..Deseq2McolsDiagnostics::default()
        };

        assert!(write_deseq_mcols_diagnostics_tsv(&path, None, &diagnostics).is_err());
        let _ = fs::remove_file(&path);
    }

    #[test]
    fn write_deseq_mcols_diagnostics_tsv_rejects_misaligned_gene_names() {
        let path = unique_test_path("mcols_bad_names.tsv");
        let gene_names = vec!["gene_a".to_string()];
        let diagnostics = Deseq2McolsDiagnostics {
            disp_gene_est: Some(vec![0.1, 0.2]),
            ..Deseq2McolsDiagnostics::default()
        };

        assert!(write_deseq_mcols_diagnostics_tsv(&path, Some(&gene_names), &diagnostics).is_err());
        let _ = fs::remove_file(&path);
    }

    #[test]
    fn write_independent_filter_metadata_tsv_writes_deseq2_shapes() {
        let num_rej_path = unique_test_path("filter_num_rej.tsv");
        let lowess_path = unique_test_path("filter_lowess.tsv");
        let scalar_path = unique_test_path("filter_metadata.tsv");
        let filtering = IndependentFilteringOutput {
            enabled: true,
            theta: vec![0.0, 0.5],
            num_rejections: vec![1, 3],
            selected_index: Some(1),
            filter_theta: Some(0.5),
            filter_threshold: Some(10.0),
            lowess_fit: Some(vec![1.25, 2.75]),
            alpha: 0.1,
        };

        write_independent_filter_num_rej_tsv(&num_rej_path, &filtering).unwrap();
        write_independent_filter_lowess_tsv(&lowess_path, &filtering).unwrap();
        write_independent_filter_metadata_tsv(&scalar_path, &filtering).unwrap();

        let num_rej = fs::read_to_string(&num_rej_path).unwrap();
        let lowess = fs::read_to_string(&lowess_path).unwrap();
        let scalar = fs::read_to_string(&scalar_path).unwrap();
        let _ = fs::remove_file(&num_rej_path);
        let _ = fs::remove_file(&lowess_path);
        let _ = fs::remove_file(&scalar_path);
        assert_eq!(num_rej, "theta\tnumRej\n0\t1\n0.5\t3\n");
        assert_eq!(lowess, "x\ty\n0\t1.25\n0.5\t2.75\n");
        assert_eq!(
            scalar,
            "name\tvalue\nfilterThreshold\t10\nfilterTheta\t0.5\nalpha\t0.1\n"
        );
    }

    #[test]
    fn write_independent_filter_metadata_tsv_rejects_invalid_shapes() {
        let path = unique_test_path("bad_filter_metadata.tsv");
        let mut filtering = IndependentFilteringOutput {
            enabled: true,
            theta: vec![0.0, 0.5],
            num_rejections: vec![1],
            selected_index: Some(1),
            filter_theta: Some(0.5),
            filter_threshold: Some(10.0),
            lowess_fit: Some(vec![1.25, 2.75]),
            alpha: 0.1,
        };

        assert!(write_independent_filter_num_rej_tsv(&path, &filtering).is_err());
        filtering.num_rejections = vec![1, 3];
        filtering.lowess_fit = Some(vec![1.25]);
        assert!(write_independent_filter_lowess_tsv(&path, &filtering).is_err());
        filtering.lowess_fit = Some(vec![1.25, f64::NAN]);
        assert!(write_independent_filter_lowess_tsv(&path, &filtering).is_err());
        filtering.lowess_fit = Some(vec![1.25, 2.75]);
        filtering.theta[1] = 1.2;
        assert!(write_independent_filter_metadata_tsv(&path, &filtering).is_err());
        filtering.theta[1] = 0.5;
        filtering.alpha = f64::NAN;
        assert!(write_independent_filter_metadata_tsv(&path, &filtering).is_err());
        let _ = fs::remove_file(&path);
    }

    fn unique_test_path(name: &str) -> PathBuf {
        std::env::temp_dir().join(format!(
            "rsdeseq2_{}_{}_{}",
            std::process::id(),
            std::thread::current().name().unwrap_or("test"),
            name
        ))
    }
}
