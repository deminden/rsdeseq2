use super::*;

#[test]
fn native_weighted_glm_mu_map_matches_deseq2_reference() {
    let reference_rows = read_reference_tsv("native_weighted_glm_mu_reference.tsv");
    let mu_rows = read_reference_tsv("native_weighted_glm_mu_dispersion_mu.tsv");
    let size_factors = read_size_factors("size_factors_ratio.tsv");
    let observation_weights = read_reference_matrix("observation_weights.tsv");

    let counts = reference_counts();
    let design = reference_full_design();
    let fit = DeseqBuilder::new()
        .size_factors(size_factors)
        .observation_weights(observation_weights)
        .fit_type(FitType::Mean)
        .gene_wise_dispersion_options(GeneWiseDispersionOptions {
            use_cox_reid: false,
            niter: 2,
            ..GeneWiseDispersionOptions::default()
        })
        .fit_map_dispersions_glm_mu(&counts, &design)
        .unwrap();

    let disp_gene_est = fit.disp_gene_est.as_ref().unwrap();
    let disp_gene_iter = fit.disp_gene_iter.as_ref().unwrap();
    let disp_fit = fit.disp_fit.as_ref().unwrap();
    let disp_map = fit.disp_map.as_ref().unwrap();
    let dispersion = fit.dispersion.as_ref().unwrap();
    let disp_iter = fit.disp_iter.as_ref().unwrap();
    let disp_outlier = fit.disp_outlier.as_ref().unwrap();
    let mu = fit.mu.as_ref().unwrap();
    let samples = reference_sample_names();

    assert_eq!(reference_rows.len(), counts.n_genes());
    assert_eq!(mu_rows.len(), counts.n_genes());
    for (gene, row) in reference_rows.iter().enumerate() {
        assert_eq!(
            row.get("gene").map(String::as_str),
            Some(reference_gene_names()[gene].as_str())
        );
        assert_eq!(fit.all_zero[gene], parse_required_bool(row, "allZero"));
        assert_eq!(
            fit.weights_fail.as_ref().unwrap()[gene],
            parse_required_bool(row, "weightsFail")
        );
        assert_float_close(
            fit.base_mean[gene],
            parse_required_f64(row, "baseMean"),
            1e-10,
            1e-10,
            &format!("weighted GLM-mu baseMean gene {gene}"),
        );
        assert_float_close(
            fit.base_var[gene],
            parse_required_f64(row, "baseVar"),
            1e-10,
            1e-10,
            &format!("weighted GLM-mu baseVar gene {gene}"),
        );
        assert_f64_or_missing(
            disp_gene_est[gene],
            parse_optional_f64(row, "dispGeneEst"),
            1e-4,
            1e-4,
            &format!("weighted GLM-mu dispGeneEst gene {gene}"),
        );
        assert_f64_or_missing(
            disp_fit[gene],
            parse_optional_f64(row, "dispFit"),
            1e-4,
            1e-4,
            &format!("weighted GLM-mu dispFit gene {gene}"),
        );
        assert_f64_or_missing(
            disp_map[gene],
            parse_optional_f64(row, "dispMAP"),
            1e-4,
            1e-4,
            &format!("weighted GLM-mu dispMAP gene {gene}"),
        );
        assert_f64_or_missing(
            dispersion[gene],
            parse_optional_f64(row, "dispersion"),
            1e-4,
            1e-4,
            &format!("weighted GLM-mu final dispersion gene {gene}"),
        );
        if !fit.all_zero[gene] {
            // The Armijo path can take a different number of equivalent
            // line-search steps while landing on the same dispersion and mu.
            assert!(
                disp_gene_iter[gene] > 0,
                "weighted GLM-mu gene-wise iterations gene {gene}"
            );
            assert!(
                parse_required_f64(row, "dispGeneIter") > 0.0,
                "DESeq2 weighted GLM-mu gene-wise iterations gene {gene}"
            );
            assert_eq!(
                disp_iter[gene],
                parse_required_f64(row, "dispIter") as usize,
                "weighted GLM-mu MAP iterations gene {gene}"
            );
            assert_eq!(
                disp_outlier[gene],
                parse_required_bool(row, "dispOutlier"),
                "weighted GLM-mu MAP outlier gene {gene}"
            );
        }
    }

    assert_float_close(
        fit.disp_prior_var.unwrap(),
        parse_required_f64(&reference_rows[0], "dispPriorVar"),
        1e-4,
        1e-4,
        "weighted GLM-mu dispersion prior variance",
    );

    for (gene, row) in mu_rows.iter().enumerate() {
        for (sample, sample_name) in samples.iter().enumerate() {
            assert_f64_or_missing(
                mu.row(gene).unwrap()[sample],
                parse_optional_f64(row, sample_name),
                1e-4,
                1e-4,
                &format!("weighted GLM-mu dispersion mu gene {gene} sample {sample}"),
            );
        }
    }
}

#[test]
fn native_weighted_glm_mu_local_map_matches_deseq2_reference() {
    let reference_rows = read_reference_tsv("native_weighted_glm_mu_local_reference.tsv");
    let mu_rows = read_reference_tsv("native_weighted_glm_mu_local_dispersion_mu.tsv");
    let size_factors = read_size_factors("size_factors_ratio.tsv");
    let observation_weights = read_reference_matrix("observation_weights.tsv");

    let counts = reference_counts();
    let design = reference_full_design();
    let fit = DeseqBuilder::new()
        .size_factors(size_factors)
        .observation_weights(observation_weights)
        .fit_type(FitType::Local)
        .gene_wise_dispersion_options(GeneWiseDispersionOptions {
            use_cox_reid: false,
            niter: 2,
            ..GeneWiseDispersionOptions::default()
        })
        .fit_map_dispersions_glm_mu(&counts, &design)
        .unwrap();

    let disp_gene_est = fit.disp_gene_est.as_ref().unwrap();
    let disp_gene_iter = fit.disp_gene_iter.as_ref().unwrap();
    let disp_fit = fit.disp_fit.as_ref().unwrap();
    let disp_map = fit.disp_map.as_ref().unwrap();
    let dispersion = fit.dispersion.as_ref().unwrap();
    let disp_iter = fit.disp_iter.as_ref().unwrap();
    let disp_outlier = fit.disp_outlier.as_ref().unwrap();
    let mu = fit.mu.as_ref().unwrap();
    let samples = reference_sample_names();

    assert_eq!(reference_rows.len(), counts.n_genes());
    assert_eq!(mu_rows.len(), counts.n_genes());
    for (gene, row) in reference_rows.iter().enumerate() {
        assert_eq!(
            row.get("gene").map(String::as_str),
            Some(reference_gene_names()[gene].as_str())
        );
        let all_zero = parse_required_bool(row, "allZero");
        let weights_fail = parse_required_bool(row, "weightsFail");
        let skipped = all_zero || weights_fail;
        assert_eq!(fit.all_zero[gene], skipped);
        assert_eq!(fit.weights_fail.as_ref().unwrap()[gene], weights_fail);
        assert_float_close(
            fit.base_mean[gene],
            parse_required_f64(row, "baseMean"),
            1e-10,
            1e-10,
            &format!("weighted GLM-mu local MAP baseMean gene {gene}"),
        );
        assert_float_close(
            fit.base_var[gene],
            parse_required_f64(row, "baseVar"),
            1e-10,
            1e-10,
            &format!("weighted GLM-mu local MAP baseVar gene {gene}"),
        );
        assert_f64_or_missing(
            disp_gene_est[gene],
            parse_optional_f64(row, "dispGeneEst"),
            1e-4,
            1e-4,
            &format!("weighted GLM-mu local MAP dispGeneEst gene {gene}"),
        );
        assert_f64_or_missing(
            disp_fit[gene],
            parse_optional_f64(row, "dispFit"),
            1e-4,
            1e-4,
            &format!("weighted GLM-mu local MAP dispFit gene {gene}"),
        );
        assert_f64_or_missing(
            disp_map[gene],
            parse_optional_f64(row, "dispMAP"),
            1e-4,
            1e-4,
            &format!("weighted GLM-mu local MAP dispMAP gene {gene}"),
        );
        assert_f64_or_missing(
            dispersion[gene],
            parse_optional_f64(row, "dispersion"),
            1e-4,
            1e-4,
            &format!("weighted GLM-mu local MAP final dispersion gene {gene}"),
        );
        if !skipped {
            assert!(
                disp_gene_iter[gene] > 0,
                "weighted GLM-mu local MAP gene-wise iterations gene {gene}"
            );
            assert!(
                parse_required_f64(row, "dispGeneIter") > 0.0,
                "DESeq2 weighted GLM-mu local MAP gene-wise iterations gene {gene}"
            );
            assert_eq!(
                disp_iter[gene],
                parse_required_f64(row, "dispIter") as usize,
                "weighted GLM-mu local MAP iterations gene {gene}"
            );
            assert_eq!(
                disp_outlier[gene],
                parse_required_bool(row, "dispOutlier"),
                "weighted GLM-mu local MAP outlier gene {gene}"
            );
        }
    }

    assert_float_close(
        fit.disp_prior_var.unwrap(),
        parse_required_f64(&reference_rows[0], "dispPriorVar"),
        1e-4,
        1e-4,
        "weighted GLM-mu local MAP dispersion prior variance",
    );

    for (gene, row) in mu_rows.iter().enumerate() {
        for (sample, sample_name) in samples.iter().enumerate() {
            assert_f64_or_missing(
                mu.row(gene).unwrap()[sample],
                parse_optional_f64(row, sample_name),
                1e-4,
                1e-4,
                &format!("weighted GLM-mu local MAP dispersion mu gene {gene} sample {sample}"),
            );
        }
    }
}

#[test]
fn native_weighted_glm_mu_local_cox_reid_map_matches_deseq2_reference() {
    let reference_rows = read_reference_tsv("native_weighted_glm_mu_local_cr_map_reference.tsv");
    let mu_rows = read_reference_tsv("native_weighted_glm_mu_local_cr_map_dispersion_mu.tsv");
    let size_factors = read_size_factors("size_factors_ratio.tsv");
    let observation_weights = read_reference_matrix("observation_weights.tsv");

    let counts = reference_counts();
    let design = reference_full_design();
    let fit = DeseqBuilder::new()
        .size_factors(size_factors)
        .observation_weights(observation_weights)
        .fit_type(FitType::Local)
        .gene_wise_dispersion_options(GeneWiseDispersionOptions {
            niter: 2,
            ..GeneWiseDispersionOptions::default()
        })
        .fit_map_dispersions_glm_mu(&counts, &design)
        .unwrap();

    let disp_gene_est = fit.disp_gene_est.as_ref().unwrap();
    let disp_gene_iter = fit.disp_gene_iter.as_ref().unwrap();
    let disp_fit = fit.disp_fit.as_ref().unwrap();
    let disp_map = fit.disp_map.as_ref().unwrap();
    let dispersion = fit.dispersion.as_ref().unwrap();
    let disp_iter = fit.disp_iter.as_ref().unwrap();
    let disp_outlier = fit.disp_outlier.as_ref().unwrap();
    let mu = fit.mu.as_ref().unwrap();
    let samples = reference_sample_names();

    assert_eq!(reference_rows.len(), counts.n_genes());
    assert_eq!(mu_rows.len(), counts.n_genes());
    for (gene, row) in reference_rows.iter().enumerate() {
        assert_eq!(
            row.get("gene").map(String::as_str),
            Some(reference_gene_names()[gene].as_str())
        );
        let all_zero = parse_required_bool(row, "allZero");
        let weights_fail = parse_required_bool(row, "weightsFail");
        let skipped = all_zero || weights_fail;
        assert_eq!(fit.all_zero[gene], skipped);
        assert_eq!(fit.weights_fail.as_ref().unwrap()[gene], weights_fail);
        assert_float_close(
            fit.base_mean[gene],
            parse_required_f64(row, "baseMean"),
            1e-10,
            1e-10,
            &format!("weighted GLM-mu Cox-Reid local MAP baseMean gene {gene}"),
        );
        assert_float_close(
            fit.base_var[gene],
            parse_required_f64(row, "baseVar"),
            1e-10,
            1e-10,
            &format!("weighted GLM-mu Cox-Reid local MAP baseVar gene {gene}"),
        );
        assert_f64_or_missing(
            disp_gene_est[gene],
            parse_optional_f64(row, "dispGeneEst"),
            1e-4,
            1e-4,
            &format!("weighted GLM-mu Cox-Reid local MAP dispGeneEst gene {gene}"),
        );
        assert_f64_or_missing(
            disp_fit[gene],
            parse_optional_f64(row, "dispFit"),
            1e-4,
            1e-4,
            &format!("weighted GLM-mu Cox-Reid local MAP dispFit gene {gene}"),
        );
        assert_f64_or_missing(
            disp_map[gene],
            parse_optional_f64(row, "dispMAP"),
            1e-4,
            1e-4,
            &format!("weighted GLM-mu Cox-Reid local MAP dispMAP gene {gene}"),
        );
        assert_f64_or_missing(
            dispersion[gene],
            parse_optional_f64(row, "dispersion"),
            1e-4,
            1e-4,
            &format!("weighted GLM-mu Cox-Reid local MAP final dispersion gene {gene}"),
        );
        if !skipped {
            assert!(
                disp_gene_iter[gene] > 0,
                "weighted GLM-mu Cox-Reid local MAP gene-wise iterations gene {gene}"
            );
            assert!(
                parse_required_f64(row, "dispGeneIter") > 0.0,
                "DESeq2 weighted GLM-mu Cox-Reid local MAP gene-wise iterations gene {gene}"
            );
            assert_eq!(
                disp_iter[gene],
                parse_required_f64(row, "dispIter") as usize,
                "weighted GLM-mu Cox-Reid local MAP iterations gene {gene}"
            );
            assert_eq!(
                disp_outlier[gene],
                parse_required_bool(row, "dispOutlier"),
                "weighted GLM-mu Cox-Reid local MAP outlier gene {gene}"
            );
        }
    }

    assert_float_close(
        fit.disp_prior_var.unwrap(),
        parse_required_f64(&reference_rows[0], "dispPriorVar"),
        1e-4,
        1e-4,
        "weighted GLM-mu Cox-Reid local MAP dispersion prior variance",
    );

    for (gene, row) in mu_rows.iter().enumerate() {
        for (sample, sample_name) in samples.iter().enumerate() {
            assert_f64_or_missing(
                mu.row(gene).unwrap()[sample],
                parse_optional_f64(row, sample_name),
                1e-4,
                1e-4,
                &format!(
                    "weighted GLM-mu Cox-Reid local MAP dispersion mu gene {gene} sample {sample}"
                ),
            );
        }
    }
}

#[test]
fn native_weighted_glm_mu_mean_cox_reid_map_matches_deseq2_reference() {
    let reference_rows = read_reference_tsv("native_weighted_glm_mu_mean_cr_map_reference.tsv");
    let mu_rows = read_reference_tsv("native_weighted_glm_mu_mean_cr_map_dispersion_mu.tsv");
    let size_factors = read_size_factors("size_factors_ratio.tsv");
    let observation_weights = read_reference_matrix("observation_weights.tsv");

    let counts = reference_counts();
    let design = reference_full_design();
    let fit = DeseqBuilder::new()
        .size_factors(size_factors)
        .observation_weights(observation_weights)
        .fit_type(FitType::Mean)
        .gene_wise_dispersion_options(GeneWiseDispersionOptions {
            niter: 2,
            ..GeneWiseDispersionOptions::default()
        })
        .fit_map_dispersions_glm_mu(&counts, &design)
        .unwrap();

    let disp_gene_est = fit.disp_gene_est.as_ref().unwrap();
    let disp_gene_iter = fit.disp_gene_iter.as_ref().unwrap();
    let disp_fit = fit.disp_fit.as_ref().unwrap();
    let disp_map = fit.disp_map.as_ref().unwrap();
    let dispersion = fit.dispersion.as_ref().unwrap();
    let disp_iter = fit.disp_iter.as_ref().unwrap();
    let disp_outlier = fit.disp_outlier.as_ref().unwrap();
    let mu = fit.mu.as_ref().unwrap();
    let samples = reference_sample_names();

    assert_eq!(reference_rows.len(), counts.n_genes());
    assert_eq!(mu_rows.len(), counts.n_genes());
    for (gene, row) in reference_rows.iter().enumerate() {
        assert_eq!(
            row.get("gene").map(String::as_str),
            Some(reference_gene_names()[gene].as_str())
        );
        assert_eq!(fit.all_zero[gene], parse_required_bool(row, "allZero"));
        assert_eq!(
            fit.weights_fail.as_ref().unwrap()[gene],
            parse_required_bool(row, "weightsFail")
        );
        assert_float_close(
            fit.base_mean[gene],
            parse_required_f64(row, "baseMean"),
            1e-10,
            1e-10,
            &format!("weighted GLM-mu Cox-Reid mean MAP baseMean gene {gene}"),
        );
        assert_float_close(
            fit.base_var[gene],
            parse_required_f64(row, "baseVar"),
            1e-10,
            1e-10,
            &format!("weighted GLM-mu Cox-Reid mean MAP baseVar gene {gene}"),
        );
        assert_f64_or_missing(
            disp_gene_est[gene],
            parse_optional_f64(row, "dispGeneEst"),
            1e-4,
            1e-4,
            &format!("weighted GLM-mu Cox-Reid mean MAP dispGeneEst gene {gene}"),
        );
        assert_f64_or_missing(
            disp_fit[gene],
            parse_optional_f64(row, "dispFit"),
            1e-4,
            1e-4,
            &format!("weighted GLM-mu Cox-Reid mean MAP dispFit gene {gene}"),
        );
        assert_f64_or_missing(
            disp_map[gene],
            parse_optional_f64(row, "dispMAP"),
            1e-4,
            1e-4,
            &format!("weighted GLM-mu Cox-Reid mean MAP dispMAP gene {gene}"),
        );
        assert_f64_or_missing(
            dispersion[gene],
            parse_optional_f64(row, "dispersion"),
            1e-4,
            1e-4,
            &format!("weighted GLM-mu Cox-Reid mean MAP final dispersion gene {gene}"),
        );
        if !fit.all_zero[gene] {
            assert!(
                disp_gene_iter[gene] > 0,
                "weighted GLM-mu Cox-Reid mean MAP gene-wise iterations gene {gene}"
            );
            assert!(
                parse_required_f64(row, "dispGeneIter") > 0.0,
                "DESeq2 weighted GLM-mu Cox-Reid mean MAP gene-wise iterations gene {gene}"
            );
            assert_eq!(
                disp_iter[gene],
                parse_required_f64(row, "dispIter") as usize,
                "weighted GLM-mu Cox-Reid mean MAP iterations gene {gene}"
            );
            assert_eq!(
                disp_outlier[gene],
                parse_required_bool(row, "dispOutlier"),
                "weighted GLM-mu Cox-Reid mean MAP outlier gene {gene}"
            );
        }
    }

    assert_float_close(
        fit.disp_prior_var.unwrap(),
        parse_required_f64(&reference_rows[0], "dispPriorVar"),
        1e-4,
        1e-4,
        "weighted GLM-mu Cox-Reid mean MAP dispersion prior variance",
    );

    for (gene, row) in mu_rows.iter().enumerate() {
        for (sample, sample_name) in samples.iter().enumerate() {
            assert_f64_or_missing(
                mu.row(gene).unwrap()[sample],
                parse_optional_f64(row, sample_name),
                1e-4,
                1e-4,
                &format!(
                    "weighted GLM-mu Cox-Reid mean MAP dispersion mu gene {gene} sample {sample}"
                ),
            );
        }
    }
}
