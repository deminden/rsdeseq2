use super::*;

#[test]
fn fixed_dispersion_wald_results_contrast_canonicalizes_cleaned_sample_level_aliases() {
    let counts = CountMatrix::from_row_major_u32(
        2,
        8,
        vec![
            20, 22, 18, 24, 80, 84, 78, 88, //
            100, 98, 102, 96, 110, 112, 108, 114,
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
        Some(vec![
            "Intercept".into(),
            "cell.type_B.cell_vs_A.cell".into(),
        ]),
    )
    .unwrap();
    let levels = [
        "A cell", "A cell", "A cell", "A cell", "B cell", "B cell", "B cell", "B cell",
    ]
    .into_iter()
    .map(String::from)
    .collect::<Vec<_>>();
    let dispersions = [0.05, 0.08];
    let builder = DeseqBuilder::new()
        .size_factors(vec![1.0; 8])
        .disable_cooks_cutoff()
        .disable_independent_filtering();

    let (_expected_fit, expected_results) = builder
        .clone()
        .fit_fixed_dispersion_wald_factor_level_contrast(
            &counts,
            &design,
            &dispersions,
            FactorLevelContrast::new("cell.type", "B cell", "A cell", &levels),
        )
        .unwrap();
    let (_actual_fit, actual_results) = builder
        .clone()
        .fit_fixed_dispersion_wald_results_contrast(
            &counts,
            &design,
            &dispersions,
            &ResultsContrast::character_with_reference("cell.type", "B.cell", "A.cell", "A.cell"),
            Some(&levels),
        )
        .unwrap();

    assert_eq!(actual_results.rows, expected_results.rows);
    assert_eq!(
        actual_results.metadata.result_name.as_deref(),
        Some("cell.type_B.cell_vs_A.cell")
    );
    assert_eq!(
        actual_results.metadata.comparison.as_deref(),
        Some("factor-level contrast: cell.type B cell vs A cell")
    );

    let ambiguous_levels = [
        "A cell", "A-cell", "A cell", "A-cell", "B cell", "B cell", "B cell", "B cell",
    ];
    let err = builder
        .fit_fixed_dispersion_wald_results_contrast(
            &counts,
            &design,
            &dispersions,
            &ResultsContrast::character_with_reference("cell.type", "B.cell", "A.cell", "A.cell"),
            Some(&ambiguous_levels),
        )
        .unwrap_err()
        .to_string();
    assert!(err.contains("denominator level 'A.cell' resolves ambiguously"));
}

#[test]
fn fixed_dispersion_wald_list_contrast_accepts_cleaned_positive_and_negative_aliases() {
    let counts = CountMatrix::from_row_major_u32_with_names(
        2,
        6,
        vec![
            8, 9, 14, 16, 30, 33, //
            30, 36, 50, 58, 80, 91,
        ],
        Some(vec!["g1".into(), "g2".into()]),
        None,
    )
    .unwrap();
    let design = DesignMatrix::from_row_major(
        6,
        3,
        vec![
            1.0, 0.0, 0.0, //
            1.0, 0.0, 0.0, //
            1.0, 1.0, 0.0, //
            1.0, 1.0, 0.0, //
            1.0, 0.0, 1.0, //
            1.0, 0.0, 1.0,
        ],
        Some(vec![
            "Intercept".into(),
            "cell type_B cell_vs_T cell".into(),
            "cell type_NK cell_vs_T cell".into(),
        ]),
    )
    .unwrap();
    let builder = DeseqBuilder::new()
        .size_factors(vec![1.0; 6])
        .disable_cooks_cutoff()
        .disable_independent_filtering();

    let (_cleaned_fit, cleaned_results) = builder
        .clone()
        .fit_fixed_dispersion_wald_results_contrast(
            &counts,
            &design,
            &[0.1, 0.1],
            &ResultsContrast::list(
                vec!["cell.type_NK.cell_vs_T.cell".into()],
                vec!["cell.type_B.cell_vs_T.cell".into()],
            ),
            None::<&[String]>,
        )
        .unwrap();
    let (_raw_fit, raw_results) = builder
        .fit_fixed_dispersion_wald_results_contrast(
            &counts,
            &design,
            &[0.1, 0.1],
            &ResultsContrast::list(
                vec!["cell type_NK cell_vs_T cell".into()],
                vec!["cell type_B cell_vs_T cell".into()],
            ),
            None::<&[String]>,
        )
        .unwrap();

    assert_eq!(cleaned_results.rows, raw_results.rows);
    assert_eq!(
        cleaned_results.metadata.result_name.as_deref(),
        Some("contrast")
    );
    assert_eq!(
        cleaned_results.metadata.comparison.as_deref(),
        Some(
            "coefficient list contrast: cell.type_NK.cell_vs_T.cell vs cell.type_B.cell_vs_T.cell"
        )
    );
    assert_eq!(
        raw_results.metadata.comparison.as_deref(),
        Some(
            "coefficient list contrast: cell type_NK cell_vs_T cell vs cell type_B cell_vs_T cell"
        )
    );
}
