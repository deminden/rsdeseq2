use approx::assert_relative_eq;
use rsdeseq2::prelude::*;

#[test]
fn deseq2_mcols_diagnostics_are_empty_before_glm_stages() {
    let counts = CountMatrix::from_row_major_u32(2, 3, vec![10, 12, 14, 20, 22, 24]).unwrap();

    let fit = DeseqBuilder::new()
        .fit_size_factors_and_base_means(&counts)
        .unwrap();
    let diagnostics = fit.deseq2_mcols_diagnostics();

    assert_eq!(diagnostics, Deseq2McolsDiagnostics::default());
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
        diagnostics.disp_gene_iter.as_ref(),
        fit.disp_gene_iter.as_ref()
    );
    assert!(diagnostics
        .disp_gene_iter
        .as_ref()
        .unwrap()
        .iter()
        .all(|iterations| *iterations > 0));
    assert_eq!(diagnostics.beta_conv, None);
    assert_eq!(diagnostics.deviance, None);
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
    assert_eq!(diagnostics.reduced_beta_iter, None);
    assert_eq!(diagnostics.deviance.as_ref(), fit.full_deviance.as_ref());
    assert_eq!(diagnostics.max_cooks.as_ref(), fit.max_cooks.as_ref());
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
    assert_relative_eq!(
        diagnostics_deviance[1],
        -2.0 * fit.log_like.as_ref().unwrap()[1],
        epsilon = 1e-12
    );
}
