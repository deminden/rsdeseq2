//! Golden output harness for structural refactoring.
//!
//! Runs a fixed panel of inputs and workflows, then records one SHA-256 per
//! output field over the raw bit patterns of every value. A refactoring labelled
//! IDENTICAL must leave every hash unchanged at both 1 and 4 worker threads.
//!
//! ```sh
//! cargo run --release --locked -p rsdeseq2 --example golden -- --check
//! cargo run --release --locked -p rsdeseq2 --example golden -- --write
//! cargo run --release --locked -p rsdeseq2 --example golden -- --check --case pasilla
//! ```
//!
//! Hashing rules:
//!
//! - `f64` values contribute `to_bits()` (little endian), so signed zeros and
//!   NaN payloads such as R's `NA_real_` (payload 1954) are distinguished.
//! - Integers contribute a little-endian `u64`; booleans one byte.
//! - `Option`, vectors and matrices contribute a tag, length or shape first, so
//!   `None`, empty and differently shaped outputs never collide.
//! - Strings contribute their length and UTF-8 bytes.
//! - A few nested metadata types without public numeric fields
//!   (`FormulaModelFrame`, `DispersionTrendFit`, `DesignOrigin`, result-table
//!   metadata, mixture grids/fits) contribute their `Debug` text. Rust prints
//!   `f64` in `Debug` with the shortest round-trip representation, so this is
//!   exact for every non-NaN value. Those fields carry a `~debug` suffix.
//! - A workflow that returns an error records the hash of its `Display` text
//!   under the field `error`, so changed failure behaviour is also detected.
//!
//! Struct fields are destructured exhaustively where the type allows it, so
//! adding a public field to `DeseqFit` or the result types fails to compile here
//! until the harness hashes it.
//!
//! The manifest is committed per numeric profile and target, because platform
//! arithmetic is selected at compile time today (see docs/linux_numerical_parity.md).

use rsdeseq2::cooks::{CooksRefitPlan, CooksReplacementOptions, CooksReplacementOutput};
use rsdeseq2::core::{
    Analysis, CooksReplacementTestOutput, CooksReplacementWaldOutput, CountMatrix, CountsSummary,
    DeseqBuilder, DeseqFit, RunOutput,
};
use rsdeseq2::design::DesignMatrix;
use rsdeseq2::glm::{
    BetaPriorGlmFit, BetaPriorRefitOptions, GlmInput, LrtOutput, NbinomGlmFit, WaldOutput,
    fit_glms_with_estimated_beta_prior,
};
use rsdeseq2::independent_filtering::IndependentFilteringOutput;
use rsdeseq2::io::{read_count_matrix_tsv, read_design_matrix_tsv, read_size_factors_tsv};
use rsdeseq2::matrix::{MatrixOrder, MatrixView, RowMajorMatrix};
use rsdeseq2::normalization::Offsets;
use rsdeseq2::options::{FitType, SizeFactorMethod};
use rsdeseq2::results::{DeseqResultRow, DeseqResults};
use rsdeseq2::shrinkage::{
    AdaptiveShrinkageQuery, AdaptiveShrinkageResult, LocationShrinkage, LocationShrinkageOptions,
    NegativeBinomialPosterior, NegativeBinomialPosteriorInput, NegativeBinomialPosteriorOptions,
    PosteriorMapDiagnostic, PosteriorOffsets, adaptive_shrinkage, estimate_cauchy_prior_variance,
    fit_negative_binomial_posterior,
};
use rsdeseq2::transform::RlogOutput;
use rsdeseq2::{DeseqError, shrinkage::PosteriorAnalysis};
use sha2::{Digest, Sha256};
use std::collections::BTreeMap;
use std::error::Error;
use std::fmt::Debug;
use std::path::{Path, PathBuf};

const THREADS: [usize; 2] = [1, 4];
const MANIFEST_HEADER: &str = "# rsdeseq2 golden manifest v1";

// ---------------------------------------------------------------------------
// Canonical byte encoding

trait Golden {
    fn feed(&self, hasher: &mut Sha256);
}

impl Golden for f64 {
    fn feed(&self, hasher: &mut Sha256) {
        hasher.update(self.to_bits().to_le_bytes());
    }
}

impl Golden for bool {
    fn feed(&self, hasher: &mut Sha256) {
        hasher.update([u8::from(*self)]);
    }
}

impl Golden for usize {
    fn feed(&self, hasher: &mut Sha256) {
        hasher.update((*self as u64).to_le_bytes());
    }
}

impl Golden for u32 {
    fn feed(&self, hasher: &mut Sha256) {
        hasher.update(u64::from(*self).to_le_bytes());
    }
}

impl Golden for i32 {
    fn feed(&self, hasher: &mut Sha256) {
        hasher.update(i64::from(*self).to_le_bytes());
    }
}

impl Golden for str {
    fn feed(&self, hasher: &mut Sha256) {
        self.len().feed(hasher);
        hasher.update(self.as_bytes());
    }
}

impl Golden for String {
    fn feed(&self, hasher: &mut Sha256) {
        self.as_str().feed(hasher);
    }
}

impl<T: Golden> Golden for Option<T> {
    fn feed(&self, hasher: &mut Sha256) {
        match self {
            None => hasher.update([0u8]),
            Some(value) => {
                hasher.update([1u8]);
                value.feed(hasher);
            }
        }
    }
}

impl<T: Golden> Golden for [T] {
    fn feed(&self, hasher: &mut Sha256) {
        self.len().feed(hasher);
        for value in self {
            value.feed(hasher);
        }
    }
}

impl<T: Golden> Golden for Vec<T> {
    fn feed(&self, hasher: &mut Sha256) {
        self.as_slice().feed(hasher);
    }
}

impl<T: Golden> Golden for RowMajorMatrix<T> {
    fn feed(&self, hasher: &mut Sha256) {
        self.n_rows().feed(hasher);
        self.n_cols().feed(hasher);
        self.as_slice().feed(hasher);
    }
}

impl Golden for DesignMatrix {
    fn feed(&self, hasher: &mut Sha256) {
        self.matrix().feed(hasher);
        self.coefficient_names()
            .map(<[String]>::to_vec)
            .feed(hasher);
        format!("{:?}", self.origin()).feed(hasher);
    }
}

impl Golden for CountMatrix {
    fn feed(&self, hasher: &mut Sha256) {
        self.n_genes().feed(hasher);
        self.n_samples().feed(hasher);
        self.as_slice().feed(hasher);
        self.gene_names().map(<[String]>::to_vec).feed(hasher);
        self.sample_names().map(<[String]>::to_vec).feed(hasher);
    }
}

// ---------------------------------------------------------------------------
// Field records

#[derive(Default)]
struct Fields {
    entries: Vec<(String, String)>,
}

impl Fields {
    fn put<T: Golden + ?Sized>(&mut self, name: &str, value: &T) {
        let mut hasher = Sha256::new();
        value.feed(&mut hasher);
        self.entries
            .push((name.to_string(), hex(&hasher.finalize())));
    }

    fn put_debug<T: Debug + ?Sized>(&mut self, name: &str, value: &T) {
        let mut hasher = Sha256::new();
        hasher.update(format!("{value:?}").as_bytes());
        self.entries
            .push((format!("{name}~debug"), hex(&hasher.finalize())));
    }

    fn error(&mut self, prefix: &str, error: &dyn std::fmt::Display) {
        if std::env::var_os("GOLDEN_VERBOSE").is_some() {
            eprintln!("golden: workflow error: {error}");
        }
        self.put(&join(prefix, "error"), &error.to_string());
    }
}

fn hex(bytes: &[u8]) -> String {
    bytes.iter().map(|byte| format!("{byte:02x}")).collect()
}

fn join(prefix: &str, name: &str) -> String {
    if prefix.is_empty() {
        name.to_string()
    } else {
        format!("{prefix}.{name}")
    }
}

fn put_counts_summary(out: &mut Fields, p: &str, summary: &CountsSummary) {
    let CountsSummary {
        n_genes,
        n_samples,
        all_zero_genes,
    } = summary;
    out.put(&join(p, "counts_summary.n_genes"), n_genes);
    out.put(&join(p, "counts_summary.n_samples"), n_samples);
    out.put(&join(p, "counts_summary.all_zero_genes"), all_zero_genes);
}

fn put_wald(out: &mut Fields, p: &str, wald: &Option<WaldOutput>) {
    match wald {
        None => out.put(&join(p, "wald"), &None::<bool>),
        Some(WaldOutput {
            stat,
            pvalue,
            degrees_of_freedom,
        }) => {
            out.put(&join(p, "wald.stat"), stat);
            out.put(&join(p, "wald.pvalue"), pvalue);
            out.put(&join(p, "wald.degrees_of_freedom"), degrees_of_freedom);
        }
    }
}

fn put_lrt(out: &mut Fields, p: &str, lrt: &Option<LrtOutput>) {
    match lrt {
        None => out.put(&join(p, "lrt"), &None::<bool>),
        Some(LrtOutput {
            deviance,
            pvalue,
            degrees_of_freedom,
            reduced_converged,
        }) => {
            out.put(&join(p, "lrt.deviance"), deviance);
            out.put(&join(p, "lrt.pvalue"), pvalue);
            out.put(&join(p, "lrt.degrees_of_freedom"), degrees_of_freedom);
            out.put(&join(p, "lrt.reduced_converged"), reduced_converged);
        }
    }
}

fn put_fit(out: &mut Fields, p: &str, fit: &DeseqFit) {
    let DeseqFit {
        counts_summary,
        design,
        reduced_design,
        model_frame,
        size_factors,
        normalization_factors,
        observation_weights,
        weights_fail,
        weights_design_rank,
        base_mean,
        base_var,
        all_zero,
        disp_gene_est,
        disp_gene_iter,
        disp_fit,
        dispersion_trend,
        disp_map,
        dispersion,
        disp_iter,
        disp_outlier,
        disp_prior_var,
        var_log_disp_estimates,
        dispersion_converged,
        beta,
        beta_se,
        beta_optim_start,
        beta_covariance,
        beta_converged,
        beta_iter,
        beta_optim_iter,
        beta_optim_start_objective,
        beta_optim_objective,
        beta_optim_gradient_norm,
        log_like,
        full_deviance,
        reduced_log_like,
        reduced_beta_converged,
        reduced_beta_iter,
        reduced_mu,
        reduced_hat_diagonal,
        mu,
        cooks,
        max_cooks,
        hat_diagonal,
        wald,
        lrt,
    } = fit;
    put_counts_summary(out, p, counts_summary);
    out.put(&join(p, "design"), design);
    out.put(&join(p, "reduced_design"), reduced_design);
    out.put_debug(&join(p, "model_frame"), model_frame);
    out.put(&join(p, "size_factors"), size_factors);
    out.put(&join(p, "normalization_factors"), normalization_factors);
    out.put(&join(p, "observation_weights"), observation_weights);
    out.put(&join(p, "weights_fail"), weights_fail);
    out.put(&join(p, "weights_design_rank"), weights_design_rank);
    out.put(&join(p, "base_mean"), base_mean);
    out.put(&join(p, "base_var"), base_var);
    out.put(&join(p, "all_zero"), all_zero);
    out.put(&join(p, "disp_gene_est"), disp_gene_est);
    out.put(&join(p, "disp_gene_iter"), disp_gene_iter);
    out.put(&join(p, "disp_fit"), disp_fit);
    out.put_debug(&join(p, "dispersion_trend"), dispersion_trend);
    out.put(&join(p, "disp_map"), disp_map);
    out.put(&join(p, "dispersion"), dispersion);
    out.put(&join(p, "disp_iter"), disp_iter);
    out.put(&join(p, "disp_outlier"), disp_outlier);
    out.put(&join(p, "disp_prior_var"), disp_prior_var);
    out.put(&join(p, "var_log_disp_estimates"), var_log_disp_estimates);
    out.put(&join(p, "dispersion_converged"), dispersion_converged);
    out.put(&join(p, "beta"), beta);
    out.put(&join(p, "beta_se"), beta_se);
    out.put(&join(p, "beta_optim_start"), beta_optim_start);
    out.put(&join(p, "beta_covariance"), beta_covariance);
    out.put(&join(p, "beta_converged"), beta_converged);
    out.put(&join(p, "beta_iter"), beta_iter);
    out.put(&join(p, "beta_optim_iter"), beta_optim_iter);
    out.put(
        &join(p, "beta_optim_start_objective"),
        beta_optim_start_objective,
    );
    out.put(&join(p, "beta_optim_objective"), beta_optim_objective);
    out.put(
        &join(p, "beta_optim_gradient_norm"),
        beta_optim_gradient_norm,
    );
    out.put(&join(p, "log_like"), log_like);
    out.put(&join(p, "full_deviance"), full_deviance);
    out.put(&join(p, "reduced_log_like"), reduced_log_like);
    out.put(&join(p, "reduced_beta_converged"), reduced_beta_converged);
    out.put(&join(p, "reduced_beta_iter"), reduced_beta_iter);
    out.put(&join(p, "reduced_mu"), reduced_mu);
    out.put(&join(p, "reduced_hat_diagonal"), reduced_hat_diagonal);
    out.put(&join(p, "mu"), mu);
    out.put(&join(p, "cooks"), cooks);
    out.put(&join(p, "max_cooks"), max_cooks);
    out.put(&join(p, "hat_diagonal"), hat_diagonal);
    put_wald(out, p, wald);
    put_lrt(out, p, lrt);
}

fn put_filtering(out: &mut Fields, p: &str, filtering: &Option<IndependentFilteringOutput>) {
    let name = join(p, "independent_filtering");
    let Some(IndependentFilteringOutput {
        enabled,
        theta,
        num_rejections,
        selected_index,
        filter_theta,
        filter_threshold,
        lowess_fit,
        alpha,
    }) = filtering
    else {
        out.put(&name, &None::<bool>);
        return;
    };
    out.put(&format!("{name}.enabled"), enabled);
    out.put(&format!("{name}.theta"), theta);
    out.put(&format!("{name}.num_rejections"), num_rejections);
    out.put(&format!("{name}.selected_index"), selected_index);
    out.put(&format!("{name}.filter_theta"), filter_theta);
    out.put(&format!("{name}.filter_threshold"), filter_threshold);
    out.put(&format!("{name}.lowess_fit"), lowess_fit);
    out.put(&format!("{name}.alpha"), alpha);
}

fn put_results(out: &mut Fields, p: &str, results: &DeseqResults) {
    let DeseqResults {
        rows,
        metadata,
        independent_filtering,
    } = results;
    let column = |f: fn(&DeseqResultRow) -> Option<f64>| -> Vec<Option<f64>> {
        rows.iter().map(f).collect()
    };
    let flags = |f: fn(&DeseqResultRow) -> Option<bool>| -> Vec<Option<bool>> {
        rows.iter().map(f).collect()
    };
    // Exhaustive pattern: a new result column must be hashed here.
    if let Some(DeseqResultRow {
        gene: _,
        base_mean: _,
        log2_fold_change: _,
        lfc_se: _,
        stat: _,
        pvalue: _,
        padj: _,
        dispersion: _,
        converged: _,
        max_cooks: _,
        cooks_outlier: _,
        filtered: _,
    }) = rows.first()
    {}
    let genes: Vec<Option<String>> = rows.iter().map(|row| row.gene.clone()).collect();
    out.put(&join(p, "gene"), &genes);
    out.put(&join(p, "base_mean"), &column(|r| Some(r.base_mean)));
    out.put(
        &join(p, "log2_fold_change"),
        &column(|r| r.log2_fold_change),
    );
    out.put(&join(p, "lfc_se"), &column(|r| r.lfc_se));
    out.put(&join(p, "stat"), &column(|r| r.stat));
    out.put(&join(p, "pvalue"), &column(|r| r.pvalue));
    out.put(&join(p, "padj"), &column(|r| r.padj));
    out.put(&join(p, "dispersion"), &column(|r| r.dispersion));
    out.put(&join(p, "converged"), &flags(|r| r.converged));
    out.put(&join(p, "max_cooks"), &column(|r| r.max_cooks));
    out.put(&join(p, "cooks_outlier"), &flags(|r| r.cooks_outlier));
    out.put(&join(p, "filtered"), &flags(|r| r.filtered));
    out.put_debug(&join(p, "metadata"), metadata);
    put_filtering(out, p, independent_filtering);
}

fn put_replacement_output(out: &mut Fields, p: &str, replacement: &CooksReplacementOutput) {
    let CooksReplacementOutput {
        replaced_counts,
        candidate_replacement_counts,
        outlier_cells,
        replace,
        replaceable_samples,
    } = replacement;
    out.put(&join(p, "replaced_counts"), replaced_counts);
    out.put(
        &join(p, "candidate_replacement_counts"),
        candidate_replacement_counts,
    );
    out.put(&join(p, "outlier_cells"), outlier_cells);
    out.put(&join(p, "replace"), replace);
    out.put(&join(p, "replaceable_samples"), replaceable_samples);
}

fn put_refit_plan(out: &mut Fields, p: &str, plan: &CooksRefitPlan) {
    let CooksRefitPlan {
        replacement,
        replaced_normalized_counts,
        replaced_base_mean,
        replaced_base_var,
        replaced_all_zero,
        n_refit,
        refit_rows,
        new_all_zero_rows,
        should_refit,
        post_refit_max_cooks,
    } = plan;
    put_replacement_output(out, &join(p, "replacement"), replacement);
    out.put(
        &join(p, "replaced_normalized_counts"),
        replaced_normalized_counts,
    );
    out.put(&join(p, "replaced_base_mean"), replaced_base_mean);
    out.put(&join(p, "replaced_base_var"), replaced_base_var);
    out.put(&join(p, "replaced_all_zero"), replaced_all_zero);
    out.put(&join(p, "n_refit"), n_refit);
    out.put(&join(p, "refit_rows"), refit_rows);
    out.put(&join(p, "new_all_zero_rows"), new_all_zero_rows);
    out.put(&join(p, "should_refit"), should_refit);
    out.put(&join(p, "post_refit_max_cooks"), post_refit_max_cooks);
}

fn put_cooks_replacement(out: &mut Fields, output: &CooksReplacementWaldOutput) {
    let CooksReplacementWaldOutput {
        original_fit,
        original_results,
        refit_plan,
        refit_fit,
        refit_results,
        results,
    } = output;
    put_fit(out, "original_fit", original_fit);
    put_results(out, "original_results", original_results);
    put_refit_plan(out, "refit_plan", refit_plan);
    match refit_fit {
        Some(fit) => put_fit(out, "refit_fit", fit),
        None => out.put("refit_fit", &None::<bool>),
    }
    match refit_results {
        Some(refit) => put_results(out, "refit_results", refit),
        None => out.put("refit_results", &None::<bool>),
    }
    put_results(out, "results", results);
}

fn put_nbinom(out: &mut Fields, p: &str, fit: &NbinomGlmFit) {
    let NbinomGlmFit {
        log_like,
        beta_converged,
        beta,
        beta_se,
        beta_optim_start,
        beta_covariance,
        mu,
        beta_iter,
        beta_optim_iter,
        beta_optim_start_objective,
        beta_optim_objective,
        beta_optim_gradient_norm,
        model_matrix,
        n_terms,
        hat_diagonal,
    } = fit;
    out.put(&join(p, "log_like"), log_like);
    out.put(&join(p, "beta_converged"), beta_converged);
    out.put(&join(p, "beta"), beta);
    out.put(&join(p, "beta_se"), beta_se);
    out.put(&join(p, "beta_optim_start"), beta_optim_start);
    out.put(&join(p, "beta_covariance"), beta_covariance);
    out.put(&join(p, "mu"), mu);
    out.put(&join(p, "beta_iter"), beta_iter);
    out.put(&join(p, "beta_optim_iter"), beta_optim_iter);
    out.put(
        &join(p, "beta_optim_start_objective"),
        beta_optim_start_objective,
    );
    out.put(&join(p, "beta_optim_objective"), beta_optim_objective);
    out.put(
        &join(p, "beta_optim_gradient_norm"),
        beta_optim_gradient_norm,
    );
    out.put(&join(p, "model_matrix"), model_matrix);
    out.put(&join(p, "n_terms"), n_terms);
    out.put(&join(p, "hat_diagonal"), hat_diagonal);
}

fn put_beta_prior(out: &mut Fields, fit: &BetaPriorGlmFit) {
    let BetaPriorGlmFit {
        mle_fit,
        prior_fit,
        beta_prior_variance,
    } = fit;
    put_nbinom(out, "mle_fit", mle_fit);
    put_nbinom(out, "prior_fit", prior_fit);
    out.put("beta_prior_variance", beta_prior_variance);
}

fn put_rlog(out: &mut Fields, rlog: &RlogOutput) {
    let RlogOutput {
        transformed,
        intercept,
        sample_prior_variance,
        offset_mode,
    } = rlog;
    out.put("transformed", transformed);
    out.put("intercept", intercept);
    out.put("sample_prior_variance", sample_prior_variance);
    out.put_debug("offset_mode", offset_mode);
}

fn put_ashr(out: &mut Fields, result: &AdaptiveShrinkageResult) {
    let AdaptiveShrinkageResult {
        fit,
        posterior_mean,
        posterior_sd,
        svalues,
    } = result;
    let LocationShrinkage {
        fitted_grid,
        fitted_weights,
        fit: mixture_fit,
        posterior_mean: fit_posterior_mean,
        posterior_sd: fit_posterior_sd,
        negative_probability,
        positive_probability,
        local_false_discovery_rate,
        local_false_sign_rate,
        svalues: fit_svalues,
        qvalues,
        cdf,
        log_likelihood,
        warnings,
    } = fit;
    out.put_debug("fit.fitted_grid", fitted_grid);
    out.put("fit.fitted_weights", fitted_weights);
    out.put_debug("fit.fit", mixture_fit);
    out.put("fit.posterior_mean", fit_posterior_mean);
    out.put("fit.posterior_sd", fit_posterior_sd);
    out.put("fit.negative_probability", negative_probability);
    out.put("fit.positive_probability", positive_probability);
    out.put("fit.local_false_discovery_rate", local_false_discovery_rate);
    out.put("fit.local_false_sign_rate", local_false_sign_rate);
    out.put("fit.svalues", fit_svalues);
    out.put("fit.qvalues", qvalues);
    out.put("fit.cdf", cdf);
    out.put("fit.log_likelihood", log_likelihood);
    out.put("fit.warnings", warnings);
    out.put("posterior_mean", posterior_mean);
    out.put("posterior_sd", posterior_sd);
    out.put("svalues", svalues);
}

fn put_apeglm(out: &mut Fields, posterior: &NegativeBinomialPosterior) {
    // `summary_coefficient` is private; it is hashed through its accessor.
    let NegativeBinomialPosterior {
        mode,
        standard_deviation,
        false_sign_rate,
        sign_svalues,
        false_sign_or_small_rate,
        threshold_svalues,
        prior,
        diagnostics,
        warnings,
        random_state,
        analysis,
        ..
    } = posterior;
    out.put("summary_coefficient", &posterior.summary_coefficient());
    out.put("mode", mode);
    out.put("standard_deviation", standard_deviation);
    out.put("false_sign_rate", false_sign_rate);
    out.put("sign_svalues", sign_svalues);
    out.put("false_sign_or_small_rate", false_sign_or_small_rate);
    out.put("threshold_svalues", threshold_svalues);
    out.put("prior.nuisance", &prior.nuisance);
    out.put("prior.scale", &prior.scale);
    out.put("prior.nuisance_scale", &prior.nuisance_scale);
    out.put("prior.adapted_variance", &prior.adapted_variance);
    let convergence: Vec<Option<i32>> = diagnostics.iter().map(|d| d.convergence).collect();
    let evaluations: Vec<usize> = diagnostics.iter().map(|d| d.function_evaluations).collect();
    let opposite: Vec<bool> = diagnostics
        .iter()
        .map(|d| d.opposite_start_disagrees)
        .collect();
    let laplace: Vec<bool> = diagnostics.iter().map(|d| d.laplace_valid).collect();
    if let Some(PosteriorMapDiagnostic {
        convergence: _,
        function_evaluations: _,
        opposite_start_disagrees: _,
        laplace_valid: _,
    }) = diagnostics.first()
    {}
    out.put("diagnostics.convergence", &convergence);
    out.put("diagnostics.function_evaluations", &evaluations);
    out.put("diagnostics.opposite_start_disagrees", &opposite);
    out.put("diagnostics.laplace_valid", &laplace);
    out.put("warnings", warnings);
    out.put("random_state", random_state);
    let PosteriorAnalysis {
        intervals,
        contrast_mode,
        contrast_standard_deviation,
        grid_diagnostics,
    } = analysis;
    out.put("analysis.intervals", intervals);
    out.put("analysis.contrast_mode", contrast_mode);
    out.put(
        "analysis.contrast_standard_deviation",
        contrast_standard_deviation,
    );
    out.put_debug("analysis.grid_diagnostics", grid_diagnostics);
}

// ---------------------------------------------------------------------------
// Inputs

struct Input {
    name: &'static str,
    counts: CountMatrix,
    design: DesignMatrix,
    reduced: DesignMatrix,
    coefficient: usize,
    /// Size factors supplied by the fixture, used by every workflow except
    /// size-factor estimation itself.
    size_factors: Option<Vec<f64>>,
    /// Trend used by the weighted, normalization-factor, replacement,
    /// transform and shrinkage workflows. Inputs with too few genes for the
    /// parametric trend use the local trend; their parametric errors are still
    /// recorded by the per-trend workflows.
    primary_fit_type: FitType,
}

fn data_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data")
}

fn drop_last_column(design: &DesignMatrix) -> Result<DesignMatrix, DeseqError> {
    let n = design.n_samples();
    let k = design.n_coefficients();
    let mut values = Vec::with_capacity(n * (k - 1));
    for sample in 0..n {
        values.extend_from_slice(&design.matrix().row(sample)?[..k - 1]);
    }
    let names = design
        .coefficient_names()
        .map(|names| names[..k - 1].to_vec());
    DesignMatrix::from_row_major(n, k - 1, values, names)
}

fn load(
    name: &'static str,
    counts: &str,
    design: &str,
    reduced: Option<&str>,
    size_factors: Option<&str>,
    primary_fit_type: FitType,
) -> Result<Input, Box<dyn Error>> {
    let dir = data_dir();
    let counts = read_count_matrix_tsv(dir.join(counts))?;
    let design = read_design_matrix_tsv(dir.join(design))?;
    let reduced = match reduced {
        Some(path) => read_design_matrix_tsv(dir.join(path))?,
        None => drop_last_column(&design)?,
    };
    let coefficient = design.n_coefficients() - 1;
    let size_factors = size_factors
        .map(|path| read_size_factors_tsv(dir.join(path)))
        .transpose()?;
    Ok(Input {
        name,
        counts,
        design,
        reduced,
        coefficient,
        size_factors,
        primary_fit_type,
    })
}

fn inputs() -> Result<Vec<Input>, Box<dyn Error>> {
    Ok(vec![
        load(
            "toy4x4",
            "deseq2_reference/counts.tsv",
            "deseq2_reference/design_full.tsv",
            Some("deseq2_reference/design_reduced.tsv"),
            None,
            FitType::Local,
        )?,
        load(
            "gtex_auto_dispersion",
            "gtex_auto_dispersion/counts.tsv",
            "gtex_auto_dispersion/design_matrix.tsv",
            None,
            // Nine rows cannot support size-factor estimation; the fixture
            // supplies full-cohort factors (see its README).
            Some("gtex_auto_dispersion/reference_size_factors.tsv"),
            FitType::Local,
        )?,
        load(
            "pasilla_holdout",
            "pasilla_holdout/counts.tsv",
            "pasilla_holdout/design.tsv",
            None,
            None,
            FitType::Parametric,
        )?,
        load(
            "gtex_compact",
            "gtex_v8_real_wald/compact/counts.tsv",
            "gtex_v8_real_wald/design_matrix.tsv",
            None,
            None,
            FitType::Parametric,
        )?,
    ])
}

/// Deterministic observation weights in [0.5, 1.0].
fn observation_weights(counts: &CountMatrix) -> Result<RowMajorMatrix<f64>, DeseqError> {
    let (genes, samples) = (counts.n_genes(), counts.n_samples());
    let values = (0..genes)
        .flat_map(|g| (0..samples).map(move |s| 0.5 + ((g * 7 + s * 3) % 5) as f64 * 0.125))
        .collect();
    RowMajorMatrix::from_row_major(genes, samples, values)
}

/// Deterministic positive normalization factors with sample and gene structure.
fn normalization_factors(counts: &CountMatrix) -> Result<RowMajorMatrix<f64>, DeseqError> {
    let (genes, samples) = (counts.n_genes(), counts.n_samples());
    let values = (0..genes)
        .flat_map(|g| {
            (0..samples).map(move |s| {
                let gene_sample = 0.05 * (((g * 5 + s * 11) % 7) as f64 - 3.0);
                gene_sample.exp() * (1.0 + 0.1 * (s % 3) as f64)
            })
        })
        .collect();
    RowMajorMatrix::from_row_major(genes, samples, values)
}

// ---------------------------------------------------------------------------
// Workflows

type CaseResult = BTreeMap<String, String>;

fn record<T>(
    result: Result<T, DeseqError>,
    emit: impl FnOnce(&mut Fields, &T),
) -> Vec<(String, String)> {
    let mut out = Fields::default();
    match result {
        Ok(value) => emit(&mut out, &value),
        Err(error) => out.error("", &error),
    }
    out.entries
}

fn base_builder(threads: usize) -> DeseqBuilder {
    DeseqBuilder::new().threads(threads)
}

fn run_input(input: &Input, threads: usize, out: &mut CaseResult) -> Result<(), Box<dyn Error>> {
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()?;
    pool.install(|| run_input_in_pool(input, threads, out))
        .map_err(|error| error.to_string().into())
}

fn run_input_in_pool(
    input: &Input,
    threads: usize,
    out: &mut CaseResult,
) -> Result<(), Box<dyn Error + Send + Sync>> {
    let Input {
        name,
        counts,
        design,
        reduced,
        coefficient: _,
        size_factors,
        primary_fit_type,
    } = input;
    let builder = |threads: usize| match size_factors {
        Some(factors) => base_builder(threads).size_factors(factors.clone()),
        None => base_builder(threads),
    };
    let primary = |threads: usize| builder(threads).fit_type(*primary_fit_type);
    let mut add = |workflow: &str, entries: Vec<(String, String)>| {
        for (field, hash) in entries {
            out.insert(format!("{name}\t{workflow}\t{threads}\t{field}"), hash);
        }
    };

    for (label, method) in [
        ("size_factors_ratio", SizeFactorMethod::Ratio),
        ("size_factors_poscounts", SizeFactorMethod::PosCounts),
    ] {
        let fit = base_builder(threads)
            .size_factor_method(method)
            .run(&Analysis::without_design(counts))
            .and_then(stage_fit);
        add(label, record(fit, |o, fit| put_fit(o, "fit", fit)));
    }

    let wald_lrt = |o: &mut Fields, value: &(DeseqFit, DeseqResults)| {
        put_fit(o, "fit", &value.0);
        put_results(o, "results", &value.1);
    };
    let mut primary_wald = None;
    for (label, fit_type) in [
        ("parametric", FitType::Parametric),
        ("local", FitType::Local),
        ("mean", FitType::Mean),
    ] {
        let wald = builder(threads)
            .fit_type(fit_type)
            .run(&Analysis::new(counts, design))
            .and_then(results_pair);
        if fit_type == *primary_fit_type {
            primary_wald = wald.as_ref().ok().cloned();
        }
        add(&format!("wald_{label}"), record(wald, wald_lrt));
        let lrt = builder(threads)
            .fit_type(fit_type)
            .run(&Analysis::new(counts, design).lrt(reduced))
            .and_then(results_pair);
        add(&format!("lrt_{label}"), record(lrt, wald_lrt));
    }

    let weighted = primary(threads)
        .observation_weights(observation_weights(counts)?)
        .run(&Analysis::new(counts, design))
        .and_then(results_pair);
    add("wald_weights", record(weighted, wald_lrt));
    let weighted_lrt = primary(threads)
        .observation_weights(observation_weights(counts)?)
        .run(&Analysis::new(counts, design).lrt(reduced))
        .and_then(results_pair);
    add("lrt_weights", record(weighted_lrt, wald_lrt));

    let factors = primary(threads)
        .normalization_factors(normalization_factors(counts)?)
        .run(&Analysis::new(counts, design))
        .and_then(results_pair);
    add("wald_normalization_factors", record(factors, wald_lrt));
    let factors_lrt = primary(threads)
        .normalization_factors(normalization_factors(counts)?)
        .run(&Analysis::new(counts, design).lrt(reduced))
        .and_then(results_pair);
    add("lrt_normalization_factors", record(factors_lrt, wald_lrt));

    let options = CooksReplacementOptions::new(0.99);
    let cooks = primary(threads)
        .run(&Analysis::new(counts, design).cooks_replacement(&options))
        .and_then(|output| match output {
            RunOutput::Replacement(CooksReplacementTestOutput::Wald(output)) => Ok(output),
            _ => Err(shape_error()),
        });
    add(
        "wald_cooks_replacement",
        record(cooks, put_cooks_replacement),
    );

    let Some((fit, results)) = primary_wald else {
        return Ok(());
    };
    add(
        "vst",
        record(fit.vst(counts), |o, vst| o.put("transformed", vst)),
    );
    add("rlog", record(fit.rlog(counts), put_rlog));

    if let (Some(dispersion), Some(disp_fit)) = (&fit.dispersion, &fit.disp_fit) {
        // DESeq2 refits only rows that are not all zero.
        let rows: Vec<usize> = (0..counts.n_genes())
            .filter(|&g| !fit.all_zero[g])
            .collect();
        let pick = |values: &[f64]| -> Vec<f64> { rows.iter().map(|&g| values[g]).collect() };
        let normal = counts.select_rows(&rows).and_then(|nonzero| {
            let dispersions = pick(dispersion);
            fit_glms_with_estimated_beta_prior(
                GlmInput::new(
                    &nonzero,
                    Offsets::PerSample(&fit.size_factors),
                    &dispersions,
                ),
                design,
                &pick(&fit.base_mean),
                &pick(disp_fit),
                BetaPriorRefitOptions::default(),
            )
        });
        add("shrink_normal", record(normal, put_beta_prior));
    }

    let (estimates, standard_errors): (Vec<f64>, Vec<f64>) = results
        .rows
        .iter()
        .filter_map(|row| row.log2_fold_change.zip(row.lfc_se))
        .unzip();
    let ashr = adaptive_shrinkage(
        &estimates,
        &standard_errors,
        &LocationShrinkageOptions::default(),
        AdaptiveShrinkageQuery {
            estimates: &estimates,
            standard_errors: &standard_errors,
            lfc_threshold: 0.0,
            svalues: true,
        },
    );
    add("shrink_ashr", record(ashr, put_ashr));

    if let Some(dispersion) = &fit.dispersion {
        add(
            "shrink_apeglm",
            record(
                apeglm(input, &fit, dispersion, &estimates, &standard_errors),
                put_apeglm,
            ),
        );
    }
    Ok(())
}

fn apeglm(
    input: &Input,
    fit: &DeseqFit,
    dispersion: &[f64],
    estimates: &[f64],
    standard_errors: &[f64],
) -> Result<NegativeBinomialPosterior, DeseqError> {
    let counts: Vec<f64> = input
        .counts
        .as_slice()
        .iter()
        .map(|&c| f64::from(c))
        .collect();
    let genes = input.counts.n_genes();
    let samples = input.counts.n_samples();
    // DESeq2 passes natural-log MLEs to apeglm's adaptive prior.
    let ln2 = std::f64::consts::LN_2;
    let ln_estimates: Vec<f64> = estimates.iter().map(|x| x * ln2).collect();
    let ln_errors: Vec<f64> = standard_errors.iter().map(|x| x * ln2).collect();
    let variance = estimate_cauchy_prior_variance(&ln_estimates, &ln_errors)?;
    let mut options = NegativeBinomialPosteriorOptions {
        coefficient: Some(input.coefficient),
        ..Default::default()
    };
    options.prior.scale = variance.sqrt();
    let input = NegativeBinomialPosteriorInput {
        counts: MatrixView::new(genes, samples, &counts, MatrixOrder::RowMajor)?,
        design: MatrixView::new(
            samples,
            input.design.n_coefficients(),
            input.design.matrix().as_slice(),
            MatrixOrder::RowMajor,
        )?,
        dispersions: dispersion,
        offsets: PosteriorOffsets::SizeFactors(&fit.size_factors),
        weights: None,
        mle: None,
    };
    fit_negative_binomial_posterior(input, &options)
}

// ---------------------------------------------------------------------------
// Manifest

fn profile() -> &'static str {
    rsdeseq2::numeric_profile().name()
}

fn manifest_path() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tools/validation/golden")
        .join(format!(
            "{}.{}-{}.tsv",
            profile(),
            std::env::consts::ARCH,
            std::env::consts::OS
        ))
}

fn render(entries: &CaseResult) -> String {
    let mut text = format!(
        "{MANIFEST_HEADER}\n# profile: {}\n# target: {}-{}\n# columns: input, workflow, threads, field, sha256\n",
        profile(),
        std::env::consts::ARCH,
        std::env::consts::OS
    );
    for (key, hash) in entries {
        text.push_str(key);
        text.push('\t');
        text.push_str(hash);
        text.push('\n');
    }
    text
}

fn parse(text: &str) -> Result<CaseResult, Box<dyn Error>> {
    let mut entries = CaseResult::new();
    for line in text.lines().filter(|line| !line.starts_with('#')) {
        let (key, hash) = line
            .rsplit_once('\t')
            .ok_or_else(|| format!("malformed manifest line: {line}"))?;
        entries.insert(key.to_string(), hash.to_string());
    }
    Ok(entries)
}

struct Args {
    write: bool,
    only: Option<String>,
    manifest: PathBuf,
}

fn args() -> Result<Args, Box<dyn Error>> {
    let mut parsed = Args {
        write: false,
        only: None,
        manifest: manifest_path(),
    };
    let mut mode = None;
    let mut iter = std::env::args().skip(1);
    while let Some(arg) = iter.next() {
        match arg.as_str() {
            "--check" | "--write" => {
                if mode.replace(arg.clone()).is_some() {
                    return Err("choose one of --check or --write".into());
                }
            }
            "--case" => parsed.only = Some(iter.next().ok_or("--case needs an input name")?),
            "--manifest" => {
                parsed.manifest = PathBuf::from(iter.next().ok_or("--manifest needs a path")?)
            }
            other => return Err(format!("unknown argument {other}").into()),
        }
    }
    parsed.write = mode.as_deref() == Some("--write");
    if mode.is_none() {
        return Err("usage: golden (--check | --write) [--case INPUT] [--manifest PATH]".into());
    }
    Ok(parsed)
}

fn main() -> Result<(), Box<dyn Error>> {
    let args = args()?;
    let mut actual = CaseResult::new();
    for input in inputs()? {
        if args.only.as_deref().is_some_and(|only| only != input.name) {
            continue;
        }
        for threads in THREADS {
            let start = std::time::Instant::now();
            run_input(&input, threads, &mut actual)?;
            eprintln!(
                "golden: {} threads={threads} {:.2}s",
                input.name,
                start.elapsed().as_secs_f64()
            );
        }
    }

    // Thread invariance is part of the contract: report fields that differ.
    let mut thread_differences = Vec::new();
    for (key, hash) in &actual {
        let parts: Vec<&str> = key.split('\t').collect();
        if parts[2] == "1" {
            let other = format!("{}\t{}\t4\t{}", parts[0], parts[1], parts[3]);
            if actual.get(&other) != Some(hash) {
                thread_differences.push(format!("{}\t{}\t{}", parts[0], parts[1], parts[3]));
            }
        }
    }

    if args.write {
        if args.only.is_some() {
            return Err("--write requires the full panel".into());
        }
        if let Some(parent) = args.manifest.parent() {
            std::fs::create_dir_all(parent)?;
        }
        std::fs::write(&args.manifest, render(&actual))?;
        eprintln!(
            "golden: wrote {} fields to {}",
            actual.len(),
            args.manifest.display()
        );
        for difference in &thread_differences {
            eprintln!("golden: thread-dependent field (1 vs 4): {difference}");
        }
        return Ok(());
    }

    let text = std::fs::read_to_string(&args.manifest).map_err(|error| {
        format!(
            "cannot read golden manifest {}: {error}; generate it on this target with --write",
            args.manifest.display()
        )
    })?;
    let expected = parse(&text)?;
    let mut failures = Vec::new();
    for (key, hash) in &expected {
        if args
            .only
            .as_deref()
            .is_some_and(|only| !key.starts_with(&format!("{only}\t")))
        {
            continue;
        }
        match actual.get(key) {
            None => failures.push(format!("missing\t{key}")),
            Some(found) if found != hash => failures.push(format!("changed\t{key}")),
            Some(_) => {}
        }
    }
    for key in actual.keys() {
        if !expected.contains_key(key) {
            failures.push(format!("new\t{key}"));
        }
    }
    if failures.is_empty() {
        eprintln!(
            "golden: {} fields identical ({})",
            actual.len(),
            args.manifest.display()
        );
        Ok(())
    } else {
        for failure in &failures {
            println!("{failure}");
        }
        Err(format!("golden: {} fields differ from the manifest", failures.len()).into())
    }
}

fn shape_error() -> DeseqError {
    DeseqError::InvalidOptions {
        reason: "golden: unexpected run output shape".to_string(),
    }
}

fn stage_fit(output: RunOutput) -> Result<DeseqFit, DeseqError> {
    match output {
        RunOutput::Fit(fit) => Ok(fit),
        _ => Err(shape_error()),
    }
}

fn results_pair(output: RunOutput) -> Result<(DeseqFit, DeseqResults), DeseqError> {
    match output {
        RunOutput::Results { fit, results } => Ok((fit, results)),
        _ => Err(shape_error()),
    }
}
