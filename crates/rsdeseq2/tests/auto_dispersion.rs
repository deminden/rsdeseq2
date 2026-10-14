use approx::assert_relative_eq;
use rsdeseq2::prelude::*;
use std::path::PathBuf;

fn fixture() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/gtex_auto_dispersion")
}

fn table(name: &str) -> Vec<csv::StringRecord> {
    csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(fixture().join(name))
        .unwrap()
        .records()
        .collect::<Result<Vec<_>, _>>()
        .unwrap()
}

fn data() -> (CountMatrix, DesignMatrix, Vec<f64>) {
    let counts = rsdeseq2::io::read_count_matrix_tsv(fixture().join("counts.tsv")).unwrap();
    let design = rsdeseq2::io::read_design_matrix_tsv(fixture().join("design_matrix.tsv")).unwrap();
    let factors = table("reference_size_factors.tsv")
        .iter()
        .map(|row| row[1].parse().unwrap())
        .collect();
    (counts, design, factors)
}

fn check_reference(fit: &DeseqFit, mode: &str, check_mu: bool) {
    let reference = table(&format!("reference_{mode}.tsv"));
    let estimates = fit.disp_gene_est.as_ref().unwrap();
    for (gene, row) in reference.iter().enumerate() {
        if row[1] == *"NA" {
            assert!(estimates[gene].is_nan());
        } else {
            assert_relative_eq!(
                estimates[gene],
                row[1].parse::<f64>().unwrap(),
                epsilon = 1e-9,
                max_relative = 1e-8
            );
        }
    }
    if check_mu {
        let expected = table(&format!("reference_{mode}_mu.tsv"));
        let mu = fit.mu.as_ref().unwrap();
        for (gene, row) in expected.iter().enumerate() {
            for (sample, value) in row.iter().skip(1).enumerate() {
                let actual = mu.row(gene).unwrap()[sample];
                if value == "NA" {
                    assert!(actual.is_nan());
                } else {
                    assert_relative_eq!(
                        actual,
                        value.parse::<f64>().unwrap(),
                        epsilon = 1e-9,
                        max_relative = 1e-9
                    );
                }
            }
        }
    }
}

#[test]
fn automatic_gene_estimates_match_deseq2_and_preserve_explicit_glm_mode() {
    let (counts, design, factors) = data();
    let builder = DeseqBuilder::new().size_factors(factors);
    let automatic = builder
        .clone()
        .automatic_dispersion_means()
        .fit_gene_wise_dispersions_glm_mu(&counts, &design)
        .unwrap();
    let explicit = builder
        .fit_gene_wise_dispersions_glm_mu(&counts, &design)
        .unwrap();
    check_reference(&automatic, "auto", true);
    check_reference(&explicit, "glm", true);
    let difference = automatic
        .disp_gene_est
        .as_ref()
        .unwrap()
        .iter()
        .zip(explicit.disp_gene_est.as_ref().unwrap())
        .filter_map(|(a, b)| (a.is_finite() && b.is_finite()).then_some((a - b).abs()))
        .fold(0.0_f64, f64::max);
    assert!(
        difference > 0.1,
        "fixture must distinguish DESeq2's two mean-fitting branches"
    );
}

#[test]
fn supplied_weights_disable_automatic_linear_means() {
    let (counts, design, factors) = data();
    let weights = RowMajorMatrix::from_row_major(
        counts.n_genes(),
        counts.n_samples(),
        vec![1.0; counts.n_genes() * counts.n_samples()],
    )
    .unwrap();
    let fit = DeseqBuilder::new()
        .automatic_dispersion_means()
        .size_factors(factors)
        .observation_weights(weights)
        .fit_gene_wise_dispersions_glm_mu(&counts, &design)
        .unwrap();
    check_reference(&fit, "weighted", true);
}

#[test]
fn automatic_linear_means_use_gene_sample_normalization_factors() {
    let (counts, design, factors) = data();
    let normalization = RowMajorMatrix::from_row_major(
        counts.n_genes(),
        counts.n_samples(),
        factors.repeat(counts.n_genes()),
    )
    .unwrap();
    let fit = DeseqBuilder::new()
        .automatic_dispersion_means()
        .normalization_factors(normalization)
        .fit_gene_wise_dispersions_glm_mu(&counts, &design)
        .unwrap();
    check_reference(&fit, "auto", true);
}

#[test]
fn top_level_workflows_select_deseq2_automatic_means() {
    let (counts, design, factors) = data();
    let builder = DeseqBuilder::new()
        .size_factors(factors)
        .fit_type(FitType::Mean)
        .disable_cooks_cutoff()
        .disable_independent_filtering();
    let (fit, _) = builder.fit_with_results(&counts, &design).unwrap();
    check_reference(&fit, "auto", false);
    let (named, _) = builder
        .fit_with_results_name(&counts, &design, "sexmale")
        .unwrap();
    check_reference(&named, "auto", false);
    let (contrast, _) = builder
        .fit_with_results_contrast(&counts, &design, &[0.0, 1.0])
        .unwrap();
    check_reference(&contrast, "auto", false);
}

#[test]
fn automatic_selection_recognizes_saturated_groups_in_other_parameterizations() {
    let intercept = DesignMatrix::from_row_major(4, 1, vec![1.0; 4], None).unwrap();
    let no_intercept =
        DesignMatrix::from_row_major(4, 2, vec![1., 0., 1., 0., 0., 1., 0., 1.], None).unwrap();
    let additive = DesignMatrix::from_row_major(
        4,
        3,
        vec![1., 0., 0., 1., 0., 1., 1., 1., 0., 1., 1., 1.],
        None,
    )
    .unwrap();
    let interaction = DesignMatrix::from_row_major(
        4,
        4,
        vec![
            1., 0., 0., 0., 1., 0., 1., 0., 1., 1., 0., 0., 1., 1., 1., 1.,
        ],
        None,
    )
    .unwrap();
    let numeric =
        DesignMatrix::from_row_major(4, 2, vec![1., 0., 1., 1., 1., 2., 1., 3.], None).unwrap();
    assert!(design_uses_linear_dispersion_means(&intercept).unwrap());
    assert!(design_uses_linear_dispersion_means(&no_intercept).unwrap());
    assert!(design_uses_linear_dispersion_means(&interaction).unwrap());
    assert!(!design_uses_linear_dispersion_means(&additive).unwrap());
    assert!(!design_uses_linear_dispersion_means(&numeric).unwrap());
}
