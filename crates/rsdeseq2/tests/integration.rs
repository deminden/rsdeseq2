//! All integration tests, compiled as one binary.
//!
//! Each file in this directory is a module of this crate (loaded with
//! `#[path]`, so its own `mod` declarations resolve as before). Run one file's
//! tests with `cargo test --test integration <file>::`.
//!
//! Shared helpers (`common`, `support`, the R-interface evaluators) are still
//! declared by each file that uses them, exactly as when every file was its
//! own crate, so the same source is compiled once per user and a file module
//! may share its helper's name.
#![allow(clippy::duplicate_mod, clippy::module_inception)]

#[path = "adaptive_result.rs"]
mod adaptive_result;
#[path = "advanced_ashr_optimizers.rs"]
mod advanced_ashr_optimizers;
#[path = "ashr_mode.rs"]
mod ashr_mode;
#[path = "ashr_optimizers.rs"]
mod ashr_optimizers;
#[path = "ashr_output.rs"]
mod ashr_output;
#[path = "auto_dispersion.rs"]
mod auto_dispersion;
#[path = "base_mean.rs"]
mod base_mean;
#[path = "bounded_posterior.rs"]
mod bounded_posterior;
#[path = "builder_cooks_policy.rs"]
mod builder_cooks_policy;
#[path = "cli.rs"]
mod cli;
#[path = "cli_help.rs"]
mod cli_help;
#[path = "compiled_method_policy.rs"]
mod compiled_method_policy;
#[path = "continuous_likelihood.rs"]
mod continuous_likelihood;
#[path = "contrast_projection.rs"]
mod contrast_projection;
#[path = "contrast_zero_view.rs"]
mod contrast_zero_view;
#[path = "contrasts.rs"]
mod contrasts;
#[path = "cooks.rs"]
mod cooks;
#[path = "cooks_policy.rs"]
mod cooks_policy;
#[path = "cooks_policy_reference.rs"]
mod cooks_policy_reference;
#[path = "cooks_reference.rs"]
mod cooks_reference;
#[path = "count_likelihood.rs"]
mod count_likelihood;
#[path = "design.rs"]
mod design;
#[path = "diagnostics.rs"]
mod diagnostics;
#[path = "dispersion_gene_est.rs"]
mod dispersion_gene_est;
#[path = "dispersion_map.rs"]
mod dispersion_map;
#[path = "dispersion_prior.rs"]
mod dispersion_prior;
#[path = "dispersion_reference.rs"]
mod dispersion_reference;
#[path = "dispersion_trend.rs"]
mod dispersion_trend;
#[path = "edge_cases.rs"]
mod edge_cases;
#[path = "formula_expression.rs"]
mod formula_expression;
#[path = "formula_fitting.rs"]
mod formula_fitting;
#[path = "formula_term_expansion.rs"]
mod formula_term_expansion;
#[path = "frozen_local_vst.rs"]
mod frozen_local_vst;
#[path = "fuzz_smoke.rs"]
mod fuzz_smoke;
#[path = "glm_intercept.rs"]
mod glm_intercept;
#[path = "glm_irls.rs"]
mod glm_irls;
#[path = "gtex_real_wald_parity.rs"]
mod gtex_real_wald_parity;
#[path = "heavy_tail.rs"]
mod heavy_tail;
#[path = "independent_filtering.rs"]
mod independent_filtering;
#[path = "inference_primitives.rs"]
mod inference_primitives;
#[path = "large_shape_likelihood.rs"]
mod large_shape_likelihood;
#[path = "location_shrinkage.rs"]
mod location_shrinkage;
#[path = "lrt.rs"]
mod lrt;
#[path = "lrt_pipeline.rs"]
mod lrt_pipeline;
#[path = "lrt_reference.rs"]
mod lrt_reference;
#[path = "mixture_fit.rs"]
mod mixture_fit;
#[path = "mixture_grid.rs"]
mod mixture_grid;
#[path = "model_matrix.rs"]
mod model_matrix;
#[path = "model_validity.rs"]
mod model_validity;
#[path = "multiple_testing.rs"]
mod multiple_testing;
#[path = "native_factor.rs"]
mod native_factor;
#[path = "native_formula.rs"]
mod native_formula;
#[path = "native_lrt_pipeline.rs"]
mod native_lrt_pipeline;
#[path = "native_wald_pipeline.rs"]
mod native_wald_pipeline;
#[path = "nb_likelihood.rs"]
mod nb_likelihood;
#[path = "normal_log.rs"]
mod normal_log;
#[path = "normalmix_captures.rs"]
mod normalmix_captures;
#[path = "numeric_transform.rs"]
mod numeric_transform;
#[path = "observation_weights.rs"]
mod observation_weights;
#[path = "posterior_analysis.rs"]
mod posterior_analysis;
#[path = "posterior_primitives.rs"]
mod posterior_primitives;
#[path = "posterior_rng.rs"]
mod posterior_rng;
#[path = "prepared_contrasts.rs"]
mod prepared_contrasts;
#[path = "properties.rs"]
mod properties;
#[path = "reference_integrity.rs"]
mod reference_integrity;
#[path = "regression_method_policy.rs"]
mod regression_method_policy;
#[path = "regression_posterior.rs"]
mod regression_posterior;
#[path = "replacement_execution.rs"]
mod replacement_execution;
#[path = "results.rs"]
mod results;
#[path = "results_reference.rs"]
mod results_reference;
#[path = "row_preparation.rs"]
mod row_preparation;
#[path = "se_covariance_reference.rs"]
mod se_covariance_reference;
#[path = "size_factors.rs"]
mod size_factors;
#[path = "source_convergence.rs"]
mod source_convergence;
#[path = "sqp_controls.rs"]
mod sqp_controls;
#[path = "thread_invariance.rs"]
mod thread_invariance;
#[path = "transform.rs"]
mod transform;
#[path = "tsv_readers.rs"]
mod tsv_readers;
#[path = "wald.rs"]
mod wald;
#[path = "wald_pipeline.rs"]
mod wald_pipeline;
#[path = "wald_reference.rs"]
mod wald_reference;
#[path = "weighted_refit_metadata.rs"]
mod weighted_refit_metadata;

/// With `autotests = false`, a test file missing from this list would never
/// run; fail instead.
#[test]
fn every_test_file_is_listed() {
    let listed = include_str!("integration.rs");
    let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests");
    let mut missing = Vec::new();
    for entry in std::fs::read_dir(dir).unwrap() {
        let name = entry.unwrap().file_name().into_string().unwrap();
        if let Some(stem) = name.strip_suffix(".rs")
            && stem != "integration"
            && !listed.contains(&format!("\nmod {stem};\n"))
        {
            missing.push(name);
        }
    }
    assert!(
        missing.is_empty(),
        "add to tests/integration.rs: {missing:?}"
    );
}
