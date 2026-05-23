use approx::assert_relative_eq;
use rsdeseq2::prelude::*;

fn full_design() -> DesignMatrix {
    DesignMatrix::from_row_major(
        8,
        2,
        vec![
            1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        ],
        Some(vec!["Intercept".into(), "condition_B_vs_A".into()]),
    )
    .unwrap()
}

fn reduced_design() -> DesignMatrix {
    DesignMatrix::from_row_major(
        8,
        1,
        vec![1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        Some(vec!["Intercept".into()]),
    )
    .unwrap()
}

fn counts_with_zero_row() -> CountMatrix {
    CountMatrix::from_row_major_u32_with_names(
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
        Some(
            [
                "zero", "up", "flat", "variable", "high_up", "stable", "low_up", "broad",
            ]
            .into_iter()
            .map(String::from)
            .collect(),
        ),
        None,
    )
    .unwrap()
}

fn native_lrt_builder() -> DeseqBuilder {
    DeseqBuilder::new()
        .fit_type(FitType::Mean)
        .size_factors(vec![1.0; 8])
        .gene_wise_dispersion_options(GeneWiseDispersionOptions {
            use_cox_reid: false,
            fit_method: GeneWiseDispersionFitMethod::Grid,
            niter: 2,
            ..GeneWiseDispersionOptions::default()
        })
        .disable_cooks_cutoff()
        .disable_independent_filtering()
}

#[test]
fn native_glm_mu_lrt_preserves_diagnostics_and_all_zero_rows() {
    let counts = counts_with_zero_row();
    let full = full_design();
    let reduced = reduced_design();

    let (fit, results) = native_lrt_builder()
        .fit_lrt_glm_mu(&counts, &full, &reduced, 1)
        .unwrap();

    assert_eq!(fit.design.as_ref().unwrap(), &full);
    assert_eq!(fit.reduced_design.as_ref().unwrap(), &reduced);
    assert_eq!(
        fit.all_zero,
        vec![true, false, false, false, false, false, false, false]
    );
    assert!(fit.disp_prior_var.unwrap().is_finite());
    assert_eq!(fit.disp_gene_est.as_ref().unwrap().len(), counts.n_genes());
    assert_eq!(fit.disp_gene_iter.as_ref().unwrap().len(), counts.n_genes());
    assert_eq!(fit.disp_fit.as_ref().unwrap().len(), counts.n_genes());
    assert_eq!(fit.disp_map.as_ref().unwrap().len(), counts.n_genes());
    assert_eq!(fit.dispersion.as_ref().unwrap().len(), counts.n_genes());
    assert_eq!(fit.disp_gene_iter.as_ref().unwrap()[0], 0);
    assert!(fit.disp_gene_iter.as_ref().unwrap()[1..]
        .iter()
        .all(|iterations| *iterations > 0));

    assert_eq!(fit.beta.as_ref().unwrap().n_cols(), full.n_coefficients());
    assert_eq!(
        fit.reduced_log_like.as_ref().unwrap().len(),
        counts.n_genes()
    );
    assert_eq!(
        fit.reduced_beta_converged.as_deref(),
        Some(&[false, true, true, true, true, true, true, true][..])
    );
    assert_eq!(fit.lrt.as_ref().unwrap().degrees_of_freedom, 1);
    assert_eq!(fit.lrt.as_ref().unwrap().deviance[0], None);
    assert_eq!(fit.lrt.as_ref().unwrap().pvalue[0], None);
    assert_eq!(fit.cooks.as_ref().unwrap().n_rows(), counts.n_genes());

    assert_eq!(results.rows.len(), counts.n_genes());
    assert_eq!(results.rows[0].gene.as_deref(), Some("zero"));
    assert_eq!(results.rows[0].pvalue, None);
    assert_eq!(results.rows[0].padj, None);
    assert_eq!(results.rows[0].converged, None);
    assert!(results.rows[1].stat.unwrap().is_finite());
    assert!(results.rows[1].pvalue.unwrap().is_finite());
    assert_eq!(results.rows[1].pvalue, fit.lrt.as_ref().unwrap().pvalue[1]);
}

#[test]
fn native_glm_mu_lrt_matches_fixed_pipeline_when_reusing_final_dispersions() {
    let counts = counts_with_zero_row();
    let full = full_design();
    let reduced = reduced_design();
    let builder = native_lrt_builder();

    let (native_fit, native_results) = builder.fit_lrt_glm_mu(&counts, &full, &reduced, 1).unwrap();
    let final_dispersions = native_fit.dispersion.as_ref().unwrap().clone();
    let (fixed_fit, fixed_results) = builder
        .fit_fixed_dispersion_lrt(&counts, &full, &reduced, &final_dispersions, 1)
        .unwrap();

    for gene in 0..counts.n_genes() {
        let native_disp = native_fit.dispersion.as_ref().unwrap()[gene];
        let fixed_disp = fixed_fit.dispersion.as_ref().unwrap()[gene];
        if native_disp.is_nan() {
            assert!(fixed_disp.is_nan());
        } else {
            assert_relative_eq!(native_disp, fixed_disp, epsilon = 1e-12);
        }
        assert_eq!(
            native_fit.lrt.as_ref().unwrap().deviance[gene],
            fixed_fit.lrt.as_ref().unwrap().deviance[gene]
        );
        assert_eq!(
            native_results.rows[gene].pvalue,
            fixed_results.rows[gene].pvalue
        );
        assert_eq!(
            native_results.rows[gene].padj,
            fixed_results.rows[gene].padj
        );
        assert_eq!(
            native_results.rows[gene].log2_fold_change,
            fixed_results.rows[gene].log2_fold_change
        );
    }
}

#[test]
fn native_glm_mu_lrt_cooks_replacement_refit_merges_refit_rows() {
    let counts = counts_with_zero_row();
    let full = full_design();
    let reduced = reduced_design();
    let builder = native_lrt_builder();

    let output = builder
        .fit_lrt_glm_mu_with_cooks_replacement(
            &counts,
            &full,
            &reduced,
            1,
            &CooksReplacementOptions {
                trim: 0.2,
                cooks_cutoff: 0.0,
                min_replicates: 3,
                which_samples: Some(vec![true, false, false, false, false, false, false, false]),
            },
        )
        .unwrap();

    assert!(output.original_fit.lrt.is_some());
    assert!(output.refit_plan.n_refit > 0);
    assert!(output.refit_plan.should_refit);
    assert!(!output.refit_plan.refit_rows.is_empty());
    assert!(output.refit_fit.as_ref().unwrap().lrt.is_some());
    assert!(output.refit_results.is_some());
    assert_ne!(
        output.refit_plan.replacement.replaced_counts.as_slice(),
        counts.as_slice()
    );
    assert_eq!(
        output.refit_plan.replacement.replaceable_samples,
        vec![true, false, false, false, false, false, false, false]
    );

    let refit_results = output.refit_results.as_ref().unwrap();
    for gene in output.refit_plan.refit_rows.iter().copied() {
        assert_eq!(
            output.results.rows[gene].log2_fold_change,
            refit_results.rows[gene].log2_fold_change
        );
        assert_eq!(
            output.results.rows[gene].stat,
            refit_results.rows[gene].stat
        );
        assert_eq!(
            output.results.rows[gene].pvalue,
            refit_results.rows[gene].pvalue
        );
        assert_eq!(
            output.results.rows[gene].max_cooks,
            output.refit_plan.post_refit_max_cooks[gene]
        );
        assert_relative_eq!(
            output.results.rows[gene].base_mean,
            output.refit_plan.replaced_base_mean[gene],
            epsilon = 1e-12
        );
    }
    for gene in output.refit_plan.new_all_zero_rows.iter().copied() {
        assert_eq!(output.results.rows[gene].pvalue, None);
        assert_eq!(output.results.rows[gene].stat, None);
        assert_eq!(output.results.rows[gene].log2_fold_change, None);
        assert_eq!(output.results.rows[gene].dispersion, None);
    }
}

#[test]
fn native_glm_mu_lrt_cooks_replacement_refit_skips_when_no_rows_are_marked() {
    let counts = counts_with_zero_row();
    let full = full_design();
    let reduced = reduced_design();
    let builder = native_lrt_builder();

    let output = builder
        .fit_lrt_glm_mu_with_cooks_replacement(
            &counts,
            &full,
            &reduced,
            1,
            &CooksReplacementOptions {
                trim: 0.2,
                cooks_cutoff: f64::MAX,
                min_replicates: 3,
                which_samples: None,
            },
        )
        .unwrap();

    assert_eq!(output.refit_plan.n_refit, 0);
    assert!(!output.refit_plan.should_refit);
    assert!(output.refit_fit.is_none());
    assert!(output.refit_results.is_none());
    assert_eq!(output.results, output.original_results);
    assert_eq!(
        output.refit_plan.replacement.replaced_counts.as_slice(),
        counts.as_slice()
    );
}

#[test]
fn native_linear_mu_lrt_matches_fixed_pipeline_when_reusing_final_dispersions() {
    let counts = counts_with_zero_row();
    let full = full_design();
    let reduced = reduced_design();
    let builder = native_lrt_builder();

    let (native_fit, native_results) = builder
        .fit_lrt_linear_mu(&counts, &full, &reduced, 1)
        .unwrap();
    let final_dispersions = native_fit.dispersion.as_ref().unwrap().clone();
    let (fixed_fit, fixed_results) = builder
        .fit_fixed_dispersion_lrt(&counts, &full, &reduced, &final_dispersions, 1)
        .unwrap();

    assert!(native_fit.disp_prior_var.unwrap().is_finite());
    assert_eq!(native_fit.lrt.as_ref().unwrap().degrees_of_freedom, 1);
    assert_eq!(native_results.rows[0].pvalue, None);

    for gene in 0..counts.n_genes() {
        assert_eq!(
            native_fit.lrt.as_ref().unwrap().pvalue[gene],
            fixed_fit.lrt.as_ref().unwrap().pvalue[gene]
        );
        assert_eq!(
            native_results.rows[gene].pvalue,
            fixed_results.rows[gene].pvalue
        );
        assert_eq!(
            native_results.rows[gene].stat,
            fixed_results.rows[gene].stat
        );
    }
}
