use super::*;

#[test]
fn dispersion_prior_variance_matches_deseq2_reference() {
    let rows = read_reference_tsv("dispersion_prior_variance_reference.tsv");
    let disp_gene_est = rows
        .iter()
        .map(|row| parse_required_f64(row, "dispGeneEst"))
        .collect::<Vec<_>>();
    let disp_fit = rows
        .iter()
        .map(|row| parse_required_f64(row, "dispFit"))
        .collect::<Vec<_>>();
    let n_samples = parse_required_f64(&rows[0], "nSamples") as usize;
    let n_coefficients = parse_required_f64(&rows[0], "nCoefficients") as usize;

    let output = estimate_dispersion_prior_variance(
        &disp_gene_est,
        &disp_fit,
        1e-8,
        n_samples,
        n_coefficients,
    )
    .unwrap();

    assert_eq!(
        output.residual_degrees_of_freedom,
        parse_required_f64(&rows[0], "residualDf") as usize
    );
    assert_float_close(
        output.var_log_disp_estimates,
        parse_required_f64(&rows[0], "varLogDispEsts"),
        1e-10,
        1e-10,
        "dispersion prior varLogDispEsts",
    );
    assert_float_close(
        output.expected_log_dispersion_variance,
        parse_required_f64(&rows[0], "expectedLogDispVariance"),
        1e-10,
        1e-10,
        "dispersion prior expected sampling variance",
    );
    assert_float_close(
        output.disp_prior_var,
        parse_required_f64(&rows[0], "dispPriorVar"),
        1e-10,
        1e-10,
        "dispersion prior variance",
    );
    for (idx, row) in rows.iter().enumerate() {
        assert_eq!(
            output.above_min_disp[idx],
            parse_required_bool(row, "aboveMinDisp")
        );
    }
}

#[test]
fn low_df_dispersion_prior_variance_matches_deseq2_reference_inputs() {
    let rows = read_reference_tsv("dispersion_prior_variance_low_df_reference.tsv");
    let disp_gene_est = rows
        .iter()
        .map(|row| parse_required_f64(row, "dispGeneEst"))
        .collect::<Vec<_>>();
    let disp_fit = rows
        .iter()
        .map(|row| parse_required_f64(row, "dispFit"))
        .collect::<Vec<_>>();
    let n_samples = parse_required_f64(&rows[0], "nSamples") as usize;
    let n_coefficients = parse_required_f64(&rows[0], "nCoefficients") as usize;

    let output = estimate_dispersion_prior_variance(
        &disp_gene_est,
        &disp_fit,
        1e-8,
        n_samples,
        n_coefficients,
    )
    .unwrap();

    assert_eq!(output.residual_degrees_of_freedom, 2);
    for (idx, row) in rows.iter().enumerate() {
        assert_eq!(
            output.above_min_disp[idx],
            parse_required_bool(row, "aboveMinDisp")
        );
    }
    assert_float_close(
        output.var_log_disp_estimates,
        parse_required_f64(&rows[0], "varLogDispEsts"),
        1e-12,
        1e-12,
        "low-df dispersion prior varLogDispEsts",
    );
    assert_float_close(
        output.expected_log_dispersion_variance,
        parse_required_f64(&rows[0], "expectedLogDispVariance"),
        1e-12,
        1e-12,
        "low-df dispersion prior expected sampling variance",
    );
    let deseq2_prior_var = parse_required_f64(&rows[0], "dispPriorVar");
    assert!(deseq2_prior_var >= 0.25);
    assert!(deseq2_prior_var <= 8.0);
    assert!(output.disp_prior_var >= 0.25);
    assert!(output.disp_prior_var <= 8.0);
}
