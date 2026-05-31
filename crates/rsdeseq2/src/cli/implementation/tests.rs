#[cfg(test)]
mod tests {
    use std::fs;
    use std::path::PathBuf;
    use std::time::{SystemTime, UNIX_EPOCH};

    use super::*;
    use crate::core::CountsSummary;

    #[test]
    fn write_cli_result_sidecars_exports_fit_diagnostics() {
        let path = temp_tsv_path("fit_diagnostics");
        let paths = CliResultSidecarPaths {
            column_metadata: None,
            table_metadata: None,
            independent_filter_metadata: None,
            independent_filter_num_rej: None,
            independent_filter_lowess: None,
            fit_diagnostics: Some(path.clone()),
            refit_diagnostics: None,
            fit_beta: None,
            fit_beta_se: None,
            fit_beta_optim_start: None,
            refit_beta: None,
            refit_beta_se: None,
            refit_beta_optim_start: None,
        };
        let analysis = CliAnalysisOutput {
            results: DeseqResults::default(),
            fit: Some(one_gene_fit(0.2)),
            refit: None,
            cooks: None,
            refit_plan: None,
        };
        let gene_names = vec!["geneA".to_string()];

        write_cli_result_sidecars(&paths, Some(&gene_names), &analysis).unwrap();

        let got = fs::read_to_string(&path).unwrap();
        fs::remove_file(&path).unwrap();
        assert!(got.starts_with("gene\tdispGeneEst\tdispFit\tdispersion\n"));
        assert!(got.contains("geneA\t0.2\t0.1\t0.2\n"));
    }

    #[test]
    fn write_cli_result_sidecars_rejects_missing_refit_diagnostics() {
        let path = temp_tsv_path("missing_refit_diagnostics");
        let paths = CliResultSidecarPaths {
            column_metadata: None,
            table_metadata: None,
            independent_filter_metadata: None,
            independent_filter_num_rej: None,
            independent_filter_lowess: None,
            fit_diagnostics: None,
            refit_diagnostics: Some(path),
            fit_beta: None,
            fit_beta_se: None,
            fit_beta_optim_start: None,
            refit_beta: None,
            refit_beta_se: None,
            refit_beta_optim_start: None,
        };
        let analysis = CliAnalysisOutput {
            results: DeseqResults::default(),
            fit: Some(one_gene_fit(0.2)),
            refit: None,
            cooks: None,
            refit_plan: None,
        };

        let err = write_cli_result_sidecars(&paths, None, &analysis).unwrap_err();

        assert!(matches!(
            err,
            DeseqError::InvalidOptions { reason }
                if reason.contains("replacement refit diagnostics sidecar output requires rows to be refit")
        ));
    }

    fn one_gene_fit(dispersion: f64) -> DeseqFit {
        DeseqFit {
            counts_summary: CountsSummary {
                n_genes: 1,
                n_samples: 2,
                all_zero_genes: 0,
            },
            design: None,
            reduced_design: None,
            size_factors: vec![1.0, 1.0],
            normalization_factors: None,
            observation_weights: None,
            weights_fail: None,
            weights_design_rank: None,
            base_mean: vec![1.0],
            base_var: vec![0.0],
            all_zero: vec![false],
            disp_gene_est: Some(vec![dispersion]),
            disp_gene_iter: None,
            disp_fit: Some(vec![0.1]),
            dispersion_trend: None,
            disp_map: None,
            dispersion: Some(vec![dispersion]),
            disp_iter: None,
            disp_outlier: None,
            disp_prior_var: None,
            var_log_disp_estimates: None,
            dispersion_converged: None,
            beta: None,
            beta_se: None,
            beta_optim_start: None,
            beta_covariance: None,
            beta_converged: None,
            beta_iter: None,
            beta_optim_iter: None,
            beta_optim_start_objective: None,
            beta_optim_objective: None,
            beta_optim_gradient_norm: None,
            log_like: None,
            full_deviance: None,
            reduced_log_like: None,
            reduced_beta_converged: None,
            reduced_beta_iter: None,
            reduced_mu: None,
            reduced_hat_diagonal: None,
            mu: None,
            cooks: None,
            max_cooks: None,
            hat_diagonal: None,
            wald: None,
            lrt: None,
        }
    }

    fn temp_tsv_path(name: &str) -> PathBuf {
        let nonce = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        std::env::temp_dir().join(format!(
            "rsdeseq2_cli_{name}_{}_{}.tsv",
            std::process::id(),
            nonce
        ))
    }
}
