use approx::assert_relative_eq;
use rsdeseq2::prelude::*;

fn assert_optional_float_vec_eq_with_nan(actual: Option<&Vec<f64>>, expected: Option<&Vec<f64>>) {
    match (actual, expected) {
        (Some(actual), Some(expected)) => {
            assert_eq!(actual.len(), expected.len());
            for (left, right) in actual.iter().zip(expected) {
                if left.is_nan() || right.is_nan() {
                    assert!(left.is_nan() && right.is_nan());
                } else {
                    assert_relative_eq!(*left, *right, epsilon = 0.0);
                }
            }
        }
        (None, None) => {}
        _ => panic!("optional float vector presence differs"),
    }
}

fn assert_diagnostics_frame_excludes_matrix_state(frame: &Deseq2McolsDiagnosticsDataFrame) {
    let names = frame
        .columns
        .iter()
        .map(|column| column.name)
        .collect::<Vec<_>>();
    for matrix_name in [
        "mu",
        "hatDiagonal",
        "reducedMu",
        "reducedHatDiagonal",
        "beta",
        "betaSE",
        "betaCovariance",
    ] {
        assert!(
            !names.contains(&matrix_name),
            "mcols diagnostics unexpectedly exposed matrix column {matrix_name}"
        );
    }
}

#[test]
fn deseq2_mcols_diagnostics_are_empty_before_glm_stages() {
    let counts = CountMatrix::from_row_major_u32(2, 3, vec![10, 12, 14, 20, 22, 24]).unwrap();

    let fit = DeseqBuilder::new()
        .fit_size_factors_and_base_means(&counts)
        .unwrap();
    let diagnostics = fit.deseq2_mcols_diagnostics();

    assert_eq!(diagnostics, Deseq2McolsDiagnostics::default());
    assert!(diagnostics.present_column_names().is_empty());
}

#[test]
fn deseq2_mcols_diagnostics_include_gene_wise_dispersion_iterations() {
    let counts =
        CountMatrix::from_row_major_u32(2, 4, vec![10, 10, 20, 20, 10, 30, 10, 30]).unwrap();
    let design =
        DesignMatrix::from_row_major(4, 2, vec![1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0], None)
            .unwrap();

    let fit = DeseqBuilder::new()
        .size_factors(vec![1.0, 1.0, 1.0, 1.0])
        .gene_wise_dispersion_options(GeneWiseDispersionOptions {
            fit_method: GeneWiseDispersionFitMethod::Grid,
            use_cox_reid: false,
            ..GeneWiseDispersionOptions::default()
        })
        .fit_gene_wise_dispersions_linear_mu(&counts, &design)
        .unwrap();

    let diagnostics = fit.deseq2_mcols_diagnostics();
    assert_eq!(
        diagnostics.disp_gene_est.as_ref(),
        fit.disp_gene_est.as_ref()
    );
    assert_eq!(
        diagnostics.disp_gene_iter.as_ref(),
        fit.disp_gene_iter.as_ref()
    );
    assert!(
        diagnostics
            .disp_gene_iter
            .as_ref()
            .unwrap()
            .iter()
            .all(|iterations| *iterations > 0)
    );
    assert_eq!(diagnostics.beta_conv, None);
    assert_eq!(diagnostics.deviance, None);
    assert_eq!(
        diagnostics.present_column_names(),
        vec!["dispGeneEst", "dispGeneIter"]
    );
    let frame = diagnostics.data_frame();
    assert_diagnostics_frame_excludes_matrix_state(&frame);
    assert_eq!(frame.columns.len(), 2);
    assert_eq!(frame.columns[0].name, "dispGeneEst");
    assert!(matches!(
        &frame.columns[0].values,
        Deseq2McolsDiagnosticValues::Numeric(values)
            if values.len() == counts.n_genes()
    ));
    assert_eq!(frame.columns[1].name, "dispGeneIter");
    assert!(matches!(
        &frame.columns[1].values,
        Deseq2McolsDiagnosticValues::Integer(values)
            if values == fit.disp_gene_iter.as_ref().unwrap()
    ));
}

#[test]
fn deseq2_mcols_diagnostics_include_dispersion_fit_type() {
    let counts = CountMatrix::from_row_major_u32(
        8,
        8,
        vec![
            0, 0, 0, 0, 0, 0, 0, 0, //
            0, 20, 1, 19, 2, 18, 3, 17, //
            12, 28, 10, 30, 14, 26, 11, 29, //
            30, 50, 25, 55, 35, 45, 28, 52, //
            55, 105, 60, 100, 50, 110, 65, 95, //
            120, 200, 130, 190, 115, 205, 125, 195, //
            240, 400, 250, 390, 230, 410, 260, 380, //
            15, 18, 12, 17, 45, 50, 40, 55,
        ],
    )
    .unwrap();
    let design = DesignMatrix::from_row_major(
        8,
        2,
        vec![
            1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        ],
        None,
    )
    .unwrap();

    let fit = DeseqBuilder::new()
        .size_factors(vec![1.0; 8])
        .fit_type(FitType::Mean)
        .gene_wise_dispersion_options(GeneWiseDispersionOptions {
            use_cox_reid: false,
            fit_method: GeneWiseDispersionFitMethod::Grid,
            niter: 2,
            ..GeneWiseDispersionOptions::default()
        })
        .fit_dispersion_trend_glm_mu(&counts, &design)
        .unwrap();

    let diagnostics = fit.deseq2_mcols_diagnostics();
    assert_eq!(diagnostics.dispersion_fit_type, Some("mean"));
    assert_optional_float_vec_eq_with_nan(
        diagnostics.disp_gene_est.as_ref(),
        fit.disp_gene_est.as_ref(),
    );
    assert_optional_float_vec_eq_with_nan(diagnostics.disp_fit.as_ref(), fit.disp_fit.as_ref());
    assert_eq!(diagnostics.dispersion, None);
    assert_eq!(diagnostics.disp_iter, None);
    assert_eq!(diagnostics.disp_outlier, None);
    assert_eq!(
        diagnostics.dispersion_converged.as_ref(),
        fit.dispersion_converged.as_ref()
    );
    assert_eq!(
        diagnostics.present_column_names(),
        vec!["dispGeneEst", "dispGeneIter", "dispFit"]
    );
}

#[test]
fn deseq2_mcols_diagnostics_include_map_dispersion_columns() {
    let counts = CountMatrix::from_row_major_u32(
        8,
        8,
        vec![
            0, 0, 0, 0, 0, 0, 0, 0, //
            0, 20, 1, 19, 2, 18, 3, 17, //
            12, 28, 10, 30, 14, 26, 11, 29, //
            30, 50, 25, 55, 35, 45, 28, 52, //
            55, 105, 60, 100, 50, 110, 65, 95, //
            120, 200, 130, 190, 115, 205, 125, 195, //
            240, 400, 250, 390, 230, 410, 260, 380, //
            15, 18, 12, 17, 45, 50, 40, 55,
        ],
    )
    .unwrap();
    let design = DesignMatrix::from_row_major(
        8,
        2,
        vec![
            1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        ],
        None,
    )
    .unwrap();

    let fit = DeseqBuilder::new()
        .size_factors(vec![1.0; 8])
        .fit_type(FitType::Mean)
        .gene_wise_dispersion_options(GeneWiseDispersionOptions {
            use_cox_reid: false,
            fit_method: GeneWiseDispersionFitMethod::Grid,
            niter: 2,
            ..GeneWiseDispersionOptions::default()
        })
        .fit_map_dispersions_glm_mu(&counts, &design)
        .unwrap();

    let diagnostics = fit.deseq2_mcols_diagnostics();
    assert_optional_float_vec_eq_with_nan(
        diagnostics.disp_gene_est.as_ref(),
        fit.disp_gene_est.as_ref(),
    );
    assert_optional_float_vec_eq_with_nan(diagnostics.disp_fit.as_ref(), fit.disp_fit.as_ref());
    assert_optional_float_vec_eq_with_nan(diagnostics.disp_map.as_ref(), fit.disp_map.as_ref());
    assert_optional_float_vec_eq_with_nan(diagnostics.dispersion.as_ref(), fit.dispersion.as_ref());
    assert_eq!(diagnostics.disp_iter.as_ref(), fit.disp_iter.as_ref());
    assert_eq!(diagnostics.disp_outlier.as_ref(), fit.disp_outlier.as_ref());
    assert_eq!(
        diagnostics.dispersion_converged.as_ref(),
        fit.dispersion_converged.as_ref()
    );
    assert_eq!(diagnostics.dispersion_fit_type, Some("mean"));
    assert_eq!(
        diagnostics.present_column_names(),
        vec![
            "dispGeneEst",
            "dispGeneIter",
            "dispFit",
            "dispMAP",
            "dispersion",
            "dispIter",
            "dispOutlier",
        ]
    );
    let frame = diagnostics.data_frame();
    assert_diagnostics_frame_excludes_matrix_state(&frame);
    assert_eq!(
        frame
            .columns
            .iter()
            .map(|column| column.name)
            .collect::<Vec<_>>(),
        diagnostics.present_column_names()
    );
    assert!(matches!(
        &frame.columns[3].values,
        Deseq2McolsDiagnosticValues::Numeric(values)
            if values.len() == counts.n_genes()
    ));
    assert!(matches!(
        &frame.columns[6].values,
        Deseq2McolsDiagnosticValues::Logical(values)
            if values == fit.disp_outlier.as_ref().unwrap()
    ));
}

#[test]
fn deseq2_mcols_diagnostics_use_wald_beta_conv_shape() {
    let counts = CountMatrix::from_row_major_u32(1, 4, vec![10, 10, 20, 20]).unwrap();
    let design =
        DesignMatrix::from_row_major(4, 2, vec![1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0], None)
            .unwrap();

    let (fit, _results) = DeseqBuilder::new()
        .size_factors(vec![1.0, 1.0, 1.0, 1.0])
        .disable_cooks_cutoff()
        .disable_independent_filtering()
        .irls_options(IrlsOptions {
            ridge_lambda: 0.0,
            ..IrlsOptions::default()
        })
        .fit_fixed_dispersion_wald(&counts, &design, &[0.05], 1)
        .unwrap();

    let diagnostics = fit.deseq2_mcols_diagnostics();
    assert_eq!(diagnostics.beta_conv.as_ref(), fit.beta_converged.as_ref());
    assert_eq!(diagnostics.full_beta_conv, None);
    assert_eq!(diagnostics.reduced_beta_conv, None);
    assert_eq!(diagnostics.beta_iter.as_ref(), fit.beta_iter.as_ref());
    assert_optional_float_vec_eq_with_nan(
        diagnostics.beta_optim_iter.as_ref(),
        fit.beta_optim_iter.as_ref(),
    );
    assert_optional_float_vec_eq_with_nan(
        diagnostics.beta_optim_start_objective.as_ref(),
        fit.beta_optim_start_objective.as_ref(),
    );
    assert_optional_float_vec_eq_with_nan(
        diagnostics.beta_optim_objective.as_ref(),
        fit.beta_optim_objective.as_ref(),
    );
    assert_optional_float_vec_eq_with_nan(
        diagnostics.beta_optim_gradient_norm.as_ref(),
        fit.beta_optim_gradient_norm.as_ref(),
    );
    assert_eq!(diagnostics.reduced_beta_iter, None);
    assert_eq!(diagnostics.deviance.as_ref(), fit.full_deviance.as_ref());
    assert_eq!(diagnostics.max_cooks.as_ref(), fit.max_cooks.as_ref());
    assert_eq!(fit.mu.as_ref().unwrap().n_rows(), counts.n_genes());
    assert_eq!(fit.mu.as_ref().unwrap().n_cols(), counts.n_samples());
    assert_eq!(
        fit.hat_diagonal.as_ref().unwrap().n_rows(),
        counts.n_genes()
    );
    assert_eq!(
        fit.hat_diagonal.as_ref().unwrap().n_cols(),
        counts.n_samples()
    );
    assert_eq!(
        diagnostics.present_column_names(),
        vec![
            "dispersion",
            "betaConv",
            "betaIter",
            "rustBetaOptimIter",
            "rustBetaOptimStartObjective",
            "rustBetaOptimObjective",
            "rustBetaOptimGradientNorm",
            "deviance",
            "maxCooks",
        ]
    );
    let frame = diagnostics.data_frame();
    assert_diagnostics_frame_excludes_matrix_state(&frame);
    assert_eq!(
        frame
            .columns
            .iter()
            .map(|column| column.name)
            .collect::<Vec<_>>(),
        diagnostics.present_column_names()
    );
    assert!(matches!(
        &frame.columns[8].values,
        Deseq2McolsDiagnosticValues::OptionalNumeric(values)
            if values == fit.max_cooks.as_ref().unwrap()
    ));
    assert_relative_eq!(
        diagnostics.deviance.as_ref().unwrap()[0],
        -2.0 * fit.log_like.as_ref().unwrap()[0],
        epsilon = 1e-12
    );
}

#[test]
fn deseq2_mcols_diagnostics_use_lrt_full_and_reduced_shapes() {
    let counts = CountMatrix::from_row_major_u32(2, 4, vec![0, 0, 0, 0, 10, 10, 20, 20]).unwrap();
    let full =
        DesignMatrix::from_row_major(4, 2, vec![1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0], None)
            .unwrap();
    let reduced = DesignMatrix::from_row_major(4, 1, vec![1.0, 1.0, 1.0, 1.0], None).unwrap();

    let (fit, _results) = DeseqBuilder::new()
        .size_factors(vec![1.0, 1.0, 1.0, 1.0])
        .disable_cooks_cutoff()
        .disable_independent_filtering()
        .irls_options(IrlsOptions {
            ridge_lambda: 0.0,
            ..IrlsOptions::default()
        })
        .fit_fixed_dispersion_lrt(&counts, &full, &reduced, &[0.1, 0.05], 1)
        .unwrap();

    let diagnostics = fit.deseq2_mcols_diagnostics();
    assert_eq!(diagnostics.beta_conv, None);
    assert_eq!(
        diagnostics.full_beta_conv.as_ref(),
        fit.beta_converged.as_ref()
    );
    assert_eq!(
        diagnostics.reduced_beta_conv.as_ref(),
        fit.reduced_beta_converged.as_ref()
    );
    assert_eq!(diagnostics.beta_iter.as_ref(), fit.beta_iter.as_ref());
    assert_optional_float_vec_eq_with_nan(
        diagnostics.beta_optim_iter.as_ref(),
        fit.beta_optim_iter.as_ref(),
    );
    assert_optional_float_vec_eq_with_nan(
        diagnostics.beta_optim_start_objective.as_ref(),
        fit.beta_optim_start_objective.as_ref(),
    );
    assert_optional_float_vec_eq_with_nan(
        diagnostics.beta_optim_objective.as_ref(),
        fit.beta_optim_objective.as_ref(),
    );
    assert_optional_float_vec_eq_with_nan(
        diagnostics.beta_optim_gradient_norm.as_ref(),
        fit.beta_optim_gradient_norm.as_ref(),
    );
    assert_eq!(
        diagnostics.reduced_beta_iter.as_ref(),
        fit.reduced_beta_iter.as_ref()
    );
    let diagnostics_deviance = diagnostics.deviance.as_ref().unwrap();
    let fit_deviance = fit.full_deviance.as_ref().unwrap();
    assert_eq!(diagnostics_deviance.len(), fit_deviance.len());
    assert!(diagnostics_deviance[0].is_nan());
    assert!(fit_deviance[0].is_nan());
    assert_relative_eq!(diagnostics_deviance[1], fit_deviance[1], epsilon = 1e-12);
    assert_eq!(diagnostics.max_cooks.as_ref(), fit.max_cooks.as_ref());

    assert_eq!(diagnostics.full_beta_conv.as_ref().unwrap(), &[false, true]);
    assert_eq!(
        diagnostics.reduced_beta_conv.as_ref().unwrap(),
        &[false, true]
    );
    assert_eq!(diagnostics.beta_iter.as_ref().unwrap()[0], 0);
    assert_eq!(diagnostics.reduced_beta_iter.as_ref().unwrap()[0], 0);
    assert!(diagnostics_deviance[0].is_nan());
    assert_eq!(fit.mu.as_ref().unwrap().n_rows(), counts.n_genes());
    assert_eq!(fit.mu.as_ref().unwrap().n_cols(), counts.n_samples());
    assert_eq!(
        fit.hat_diagonal.as_ref().unwrap().n_rows(),
        counts.n_genes()
    );
    assert_eq!(
        fit.hat_diagonal.as_ref().unwrap().n_cols(),
        counts.n_samples()
    );
    assert_eq!(fit.reduced_mu.as_ref().unwrap().n_rows(), counts.n_genes());
    assert_eq!(
        fit.reduced_mu.as_ref().unwrap().n_cols(),
        counts.n_samples()
    );
    assert_eq!(
        fit.reduced_hat_diagonal.as_ref().unwrap().n_rows(),
        counts.n_genes()
    );
    assert_eq!(
        fit.reduced_hat_diagonal.as_ref().unwrap().n_cols(),
        counts.n_samples()
    );
    assert_eq!(
        diagnostics.present_column_names(),
        vec![
            "dispersion",
            "fullBetaConv",
            "reducedBetaConv",
            "betaIter",
            "rustBetaOptimIter",
            "rustBetaOptimStartObjective",
            "rustBetaOptimObjective",
            "rustBetaOptimGradientNorm",
            "reducedBetaIter",
            "deviance",
            "maxCooks",
        ]
    );
    assert_relative_eq!(
        diagnostics_deviance[1],
        -2.0 * fit.log_like.as_ref().unwrap()[1],
        epsilon = 1e-12
    );
    let frame = diagnostics.data_frame();
    assert_diagnostics_frame_excludes_matrix_state(&frame);
    assert_eq!(
        frame
            .columns
            .iter()
            .map(|column| column.name)
            .collect::<Vec<_>>(),
        diagnostics.present_column_names()
    );
}
