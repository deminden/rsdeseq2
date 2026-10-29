use super::*;

#[test]
fn fixed_dispersion_wald_pipeline_can_use_t_pvalues() {
    let counts = CountMatrix::from_row_major_u32_with_names(
        1,
        4,
        vec![10, 10, 20, 20],
        Some(vec!["gene_a".into()]),
        None,
    )
    .unwrap();
    let design = DesignMatrix::from_row_major(
        4,
        2,
        vec![1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0],
        Some(vec!["Intercept".into(), "condition_B_vs_A".into()]),
    )
    .unwrap();

    let (fit, results) = DeseqBuilder::new()
        .size_factors(vec![1.0, 1.0, 1.0, 1.0])
        .irls_options(IrlsOptions {
            ridge_lambda: 0.0,
            ..IrlsOptions::default()
        })
        .disable_cooks_cutoff()
        .disable_independent_filtering()
        .wald_t_degrees_of_freedom(4.0)
        .fit_fixed_dispersion_wald(&counts, &design, &[0.05], 1)
        .unwrap();

    let wald = fit.wald.as_ref().unwrap();
    assert_eq!(wald.degrees_of_freedom.as_ref().unwrap(), &vec![Some(4.0)]);
    assert_relative_eq!(
        wald.pvalue[0].unwrap(),
        two_sided_t_pvalue(wald.stat[0].unwrap(), 4.0).unwrap(),
        epsilon = 1e-15
    );
    assert_eq!(results.rows[0].pvalue, wald.pvalue[0]);
    assert_eq!(results.rows[0].padj, wald.pvalue[0]);
}

#[test]
fn fixed_dispersion_wald_pipeline_can_use_lfc_threshold() {
    let counts = CountMatrix::from_row_major_u32_with_names(
        1,
        4,
        vec![10, 10, 20, 20],
        Some(vec!["gene_a".into()]),
        None,
    )
    .unwrap();
    let design = DesignMatrix::from_row_major(
        4,
        2,
        vec![1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0],
        Some(vec!["Intercept".into(), "condition_B_vs_A".into()]),
    )
    .unwrap();

    let (fit, results) = DeseqBuilder::new()
        .size_factors(vec![1.0, 1.0, 1.0, 1.0])
        .irls_options(IrlsOptions {
            ridge_lambda: 0.0,
            ..IrlsOptions::default()
        })
        .disable_cooks_cutoff()
        .disable_independent_filtering()
        .wald_lfc_threshold(0.5, WaldAlternative::Greater)
        .fit_fixed_dispersion_wald(&counts, &design, &[0.05], 1)
        .unwrap();

    let wald = fit.wald.as_ref().unwrap();
    assert!(wald.stat[0].unwrap() >= 0.0);
    assert_eq!(results.rows[0].pvalue, wald.pvalue[0]);
    assert_eq!(results.rows[0].padj, wald.pvalue[0]);
    assert!(wald.pvalue[0].unwrap() > 0.0);
    assert!(wald.pvalue[0].unwrap() < 0.5);
    assert_eq!(results.metadata.lfc_threshold, 0.5);
    assert_eq!(results.metadata.alt_hypothesis.as_deref(), Some("greater"));
}

#[test]
fn fixed_dispersion_wald_pipeline_validates_per_gene_t_df_length() {
    let counts = CountMatrix::from_row_major_u32(2, 3, vec![2, 4, 6, 10, 10, 10]).unwrap();
    let design = DesignMatrix::from_row_major(3, 1, vec![1.0, 1.0, 1.0], None).unwrap();

    let err = DeseqBuilder::new()
        .size_factors(vec![1.0, 1.0, 1.0])
        .wald_t_per_gene_degrees_of_freedom(vec![4.0])
        .fit_fixed_dispersion_wald(&counts, &design, &[0.1, 0.2], 0)
        .unwrap_err();

    assert!(err.to_string().contains("degrees of freedom"));
}
