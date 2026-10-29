use super::*;

#[test]
fn fitted_lrt_object_uses_stored_model_frame_for_character_contrast_results() {
    let counts = counts_with_zero_row();
    let levels = [
        "T cell", "T cell", "T cell", "T cell", "B cell", "B cell", "B cell", "B cell",
    ]
    .into_iter()
    .map(String::from)
    .collect::<Vec<_>>();
    let model_frame = FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "cell type".to_string(),
            sample_levels: levels,
            levels: Some(vec!["T cell".to_string(), "B cell".to_string()]),
            reference: Some("T cell".to_string()),
        }],
        numeric_covariates: Vec::new(),
    };
    let builder = native_lrt_builder().model_frame(model_frame);

    let (fit, mut expected_results) = builder
        .fit_lrt_formula_with_results_contrast_request(
            &counts,
            "~ `cell type`",
            "~ 1",
            &ResultsContrast::character("cell.type", "B cell", "T cell"),
        )
        .unwrap();
    expected_results.independent_filtering = None;
    let object_results = fit
        .lrt_results_contrast_request(
            &counts,
            &ResultsContrast::character("cell.type", "B cell", "T cell"),
            counts.gene_names(),
        )
        .unwrap();
    let dispatched_results = fit
        .results_contrast_request(
            &counts,
            &ResultsContrast::character("cell.type", "B cell", "T cell"),
            counts.gene_names(),
        )
        .unwrap();
    let list_results = fit
        .results_contrast_request(
            &counts,
            &ResultsContrast::list(vec!["cell.type_B.cell_vs_T.cell".into()], Vec::new()),
            counts.gene_names(),
        )
        .unwrap();
    let numeric_results = fit
        .results_contrast_request(
            &counts,
            &ResultsContrast::numeric(vec![0.0, 1.0]),
            counts.gene_names(),
        )
        .unwrap();

    assert_eq!(object_results, expected_results);
    assert_eq!(dispatched_results, expected_results);
    assert_eq!(list_results.rows, expected_results.rows);
    assert_eq!(
        list_results.metadata.result_name.as_deref(),
        Some("contrast")
    );
    assert_eq!(
        list_results.metadata.comparison.as_deref(),
        Some("coefficient list contrast: cell.type_B.cell_vs_T.cell effect")
    );
    assert_eq!(list_results.metadata.contrast, Some(vec![0.0, 1.0]));
    assert_eq!(numeric_results.rows, expected_results.rows);
    assert_eq!(
        numeric_results.metadata.result_name.as_deref(),
        Some("contrast")
    );
    assert_eq!(
        numeric_results.metadata.comparison.as_deref(),
        Some("primitive numeric contrast")
    );
    assert_eq!(numeric_results.metadata.contrast, Some(vec![0.0, 1.0]));
}

#[test]
fn fitted_lrt_object_applies_numeric_contrast_all_zero_from_stored_design() {
    let counts = CountMatrix::from_row_major_u32_with_names(
        3,
        6,
        vec![
            0, 0, 0, 0, 45, 55, //
            12, 11, 24, 22, 18, 19, //
            0, 0, 0, 0, 0, 0,
        ],
        Some(vec![
            "outside_contrast".into(),
            "tested".into(),
            "all_zero".into(),
        ]),
        None,
    )
    .unwrap();
    let levels = ["A", "A", "B", "B", "C", "C"]
        .into_iter()
        .map(String::from)
        .collect::<Vec<_>>();
    let model_frame = FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "condition".to_string(),
            sample_levels: levels,
            levels: Some(vec!["A".into(), "B".into(), "C".into()]),
            reference: Some("A".into()),
        }],
        numeric_covariates: Vec::new(),
    };
    let builder = DeseqBuilder::new()
        .size_factors(vec![1.0; 6])
        .disable_cooks_cutoff()
        .disable_independent_filtering()
        .model_frame(model_frame);
    let full_design = builder
        .expanded_formula_design("~ 0 + condition")
        .unwrap()
        .standard_design;
    let reduced_design = DesignMatrix::from_row_major(
        6,
        1,
        vec![1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        Some(vec!["Intercept".into()]),
    )
    .unwrap();
    let dispersions = vec![0.2; counts.n_genes()];
    let (fit, character_results) = builder
        .fit_fixed_dispersion_lrt_results_contrast_from_model_frame(
            &counts,
            &full_design,
            &reduced_design,
            &dispersions,
            &ResultsContrast::character("condition", "B", "A"),
            builder.current_model_frame().unwrap(),
        )
        .unwrap();
    let list_results = fit
        .results_contrast_request(
            &counts,
            &ResultsContrast::list(vec!["conditionB".into()], vec!["conditionA".into()]),
            counts.gene_names(),
        )
        .unwrap();
    let numeric_results = fit
        .results_contrast_request(
            &counts,
            &ResultsContrast::numeric(vec![-1.0, 1.0, 0.0]),
            counts.gene_names(),
        )
        .unwrap();

    assert_eq!(list_results.rows, numeric_results.rows);
    assert_eq!(character_results.rows[0].log2_fold_change, Some(0.0));
    assert_eq!(list_results.rows[0].log2_fold_change, Some(0.0));
    assert_ne!(list_results.rows[0].stat, Some(0.0));
    assert_ne!(list_results.rows[0].pvalue, Some(1.0));
    assert_eq!(numeric_results.rows[0].log2_fold_change, Some(0.0));
    assert_ne!(numeric_results.rows[0].stat, Some(0.0));
    assert_ne!(numeric_results.rows[0].pvalue, Some(1.0));
    assert_eq!(numeric_results.rows[2].log2_fold_change, None);
}

#[test]
fn lrt_effect_projection_does_not_require_a_valid_wald_test() {
    let counts = counts_with_zero_row();
    let (mut fit, _) = native_lrt_builder()
        .fit_fixed_dispersion_lrt(
            &counts,
            &full_design(),
            &reduced_design(),
            &vec![0.2; counts.n_genes()],
            1,
        )
        .unwrap();
    let beta = fit.beta.as_ref().unwrap().row(1).unwrap().to_vec();
    let lrt = fit.lrt.as_ref().unwrap().clone();
    let se = fit.beta_se.as_ref().unwrap();
    let mut values = se.as_slice().to_vec();
    values[3] = 0.0;
    fit.beta_se = Some(RowMajorMatrix::from_row_major(counts.n_genes(), 2, values).unwrap());
    let covariance = fit.beta_covariance.as_ref().unwrap();
    let mut values = covariance.as_slice().to_vec();
    values[4..8].fill(f64::NAN);
    fit.beta_covariance =
        Some(RowMajorMatrix::from_row_major(counts.n_genes(), 4, values).unwrap());
    for (contrast, effect, se) in [
        (vec![0., 1.], beta[1], Some(0.)),
        (vec![1., 1.], beta[0] + beta[1], None),
    ] {
        let result = fit
            .lrt_results_contrast_request(
                &counts,
                &ResultsContrast::numeric(contrast),
                counts.gene_names(),
            )
            .unwrap();
        assert_relative_eq!(
            result.rows[1].log2_fold_change.unwrap(),
            effect,
            epsilon = 1e-14
        );
        assert_eq!(result.rows[1].lfc_se, se);
        for (gene, row) in result.rows.iter().enumerate() {
            assert_eq!(row.stat, lrt.deviance[gene]);
            assert_eq!(row.pvalue, lrt.pvalue[gene]);
        }
    }
}
