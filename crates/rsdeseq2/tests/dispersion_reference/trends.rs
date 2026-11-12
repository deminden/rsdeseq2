use super::*;

#[test]
fn parametric_dispersion_trend_matches_deseq2_reference() {
    let rows = read_reference_tsv("parametric_trend_reference.tsv");
    let means = rows
        .iter()
        .map(|row| parse_required_f64(row, "baseMean"))
        .collect::<Vec<_>>();
    let disps = rows
        .iter()
        .map(|row| parse_required_f64(row, "dispGeneEst"))
        .collect::<Vec<_>>();

    let fit = fit_parametric_dispersion_trend(
        &means,
        &disps,
        ParametricDispersionTrendOptions::default(),
    )
    .unwrap();
    let expected_asympt = parse_required_f64(&rows[0], "asymptDisp");
    let expected_extra = parse_required_f64(&rows[0], "extraPois");

    assert_float_close(
        fit.trend.asympt_disp,
        expected_asympt,
        1e-8,
        1e-8,
        "parametric trend asymptDisp",
    );
    assert_float_close(
        fit.trend.extra_pois,
        expected_extra,
        1e-8,
        1e-8,
        "parametric trend extraPois",
    );
    for (idx, row) in rows.iter().enumerate() {
        assert_eq!(fit.use_for_fit[idx], parse_required_bool(row, "useForFit"));
        assert_float_close(
            fit.disp_fit[idx],
            parse_required_f64(row, "dispFit"),
            1e-8,
            1e-8,
            &format!("parametric trend dispFit row {idx}"),
        );
    }
}

#[test]
fn parametric_dispersion_trend_predicts_deseq2_reference_means() {
    let fit_rows = read_reference_tsv("parametric_trend_reference.tsv");
    let prediction_rows = read_reference_tsv("parametric_trend_prediction_reference.tsv");
    let means = fit_rows
        .iter()
        .map(|row| parse_required_f64(row, "baseMean"))
        .collect::<Vec<_>>();
    let disps = fit_rows
        .iter()
        .map(|row| parse_required_f64(row, "dispGeneEst"))
        .collect::<Vec<_>>();
    let prediction_means = prediction_rows
        .iter()
        .map(|row| parse_required_f64(row, "mean"))
        .collect::<Vec<_>>();

    let fit = fit_parametric_dispersion_trend(
        &means,
        &disps,
        ParametricDispersionTrendOptions::default(),
    )
    .unwrap();
    let predicted = fit
        .trend
        .evaluate_many_allow_missing(&prediction_means)
        .unwrap();

    for (idx, row) in prediction_rows.iter().enumerate() {
        assert_float_close(
            predicted[idx],
            parse_required_f64(row, "dispFit"),
            1e-8,
            1e-8,
            &format!("parametric trend prediction row {idx}"),
        );
    }
}

#[test]
fn mean_dispersion_trend_matches_deseq2_reference() {
    let rows = read_reference_tsv("mean_trend_reference.tsv");
    let means = rows
        .iter()
        .map(|row| parse_required_f64(row, "baseMean"))
        .collect::<Vec<_>>();
    let disps = rows
        .iter()
        .map(|row| parse_required_f64(row, "dispGeneEst"))
        .collect::<Vec<_>>();

    let fit =
        fit_mean_dispersion_trend(&means, &disps, MeanDispersionTrendOptions::default()).unwrap();
    let expected_mean = parse_required_f64(&rows[0], "meanDisp");

    assert_float_close(
        fit.trend.mean_disp,
        expected_mean,
        1e-12,
        1e-12,
        "mean trend meanDisp",
    );
    for (idx, row) in rows.iter().enumerate() {
        assert_eq!(fit.use_for_fit[idx], parse_required_bool(row, "useForFit"));
        assert_eq!(
            fit.use_for_mean[idx],
            parse_required_bool(row, "useForMean")
        );
        assert_float_close(
            fit.disp_fit[idx],
            parse_required_f64(row, "dispFit"),
            1e-12,
            1e-12,
            &format!("mean trend dispFit row {idx}"),
        );
    }
}

#[test]
fn local_dispersion_trend_matches_deseq2_reference_shape() {
    let rows = read_reference_tsv("local_trend_reference.tsv");
    let means = rows
        .iter()
        .map(|row| parse_required_f64(row, "baseMean"))
        .collect::<Vec<_>>();
    let disps = rows
        .iter()
        .map(|row| parse_required_f64(row, "dispGeneEst"))
        .collect::<Vec<_>>();

    let fit =
        fit_local_dispersion_trend(&means, &disps, LocalDispersionTrendOptions::default()).unwrap();

    assert_eq!(
        fit.used_min_disp_floor,
        parse_required_bool(&rows[0], "usedMinDispFloor")
    );
    for (idx, row) in rows.iter().enumerate() {
        assert_eq!(fit.use_for_fit[idx], parse_required_bool(row, "useForFit"));
        assert_float_close(
            fit.disp_fit[idx],
            parse_required_f64(row, "dispFit"),
            1e-3,
            1e-3,
            &format!("local trend dispFit row {idx}"),
        );
    }
}

#[test]
fn local_dispersion_trend_predicts_deseq2_reference_means() {
    let fit_rows = read_reference_tsv("local_trend_reference.tsv");
    let prediction_rows = read_reference_tsv("local_trend_prediction_reference.tsv");
    let means = fit_rows
        .iter()
        .map(|row| parse_required_f64(row, "baseMean"))
        .collect::<Vec<_>>();
    let disps = fit_rows
        .iter()
        .map(|row| parse_required_f64(row, "dispGeneEst"))
        .collect::<Vec<_>>();
    let prediction_means = prediction_rows
        .iter()
        .map(|row| parse_required_f64(row, "mean"))
        .collect::<Vec<_>>();

    let fit =
        fit_local_dispersion_trend(&means, &disps, LocalDispersionTrendOptions::default()).unwrap();
    let predicted = fit
        .trend
        .evaluate_many_allow_missing(&prediction_means)
        .unwrap();

    for (idx, row) in prediction_rows.iter().enumerate() {
        assert_float_close(
            predicted[idx],
            parse_required_f64(row, "dispFit"),
            2e-3,
            2e-3,
            &format!("local trend prediction row {idx}"),
        );
    }
}

#[test]
fn local_dispersion_trend_floor_matches_deseq2_reference() {
    let rows = read_reference_tsv("local_trend_floor_reference.tsv");
    let means = rows
        .iter()
        .map(|row| parse_required_f64(row, "baseMean"))
        .collect::<Vec<_>>();
    let disps = rows
        .iter()
        .map(|row| parse_required_f64(row, "dispGeneEst"))
        .collect::<Vec<_>>();

    let fit =
        fit_local_dispersion_trend(&means, &disps, LocalDispersionTrendOptions::default()).unwrap();

    assert_eq!(
        fit.used_min_disp_floor,
        parse_required_bool(&rows[0], "usedMinDispFloor")
    );
    assert_eq!(fit.genes_used, 0);
    for (idx, row) in rows.iter().enumerate() {
        assert_eq!(fit.use_for_fit[idx], parse_required_bool(row, "useForFit"));
        assert_float_close(
            fit.disp_fit[idx],
            parse_required_f64(row, "dispFit"),
            1e-15,
            1e-15,
            &format!("local trend floor dispFit row {idx}"),
        );
    }
}

#[test]
fn local_dispersion_trend_mixed_threshold_matches_deseq2_reference() {
    let rows = read_reference_tsv("local_trend_mixed_threshold_reference.tsv");
    let means = rows
        .iter()
        .map(|row| parse_required_f64(row, "baseMean"))
        .collect::<Vec<_>>();
    let disps = rows
        .iter()
        .map(|row| parse_required_f64(row, "dispGeneEst"))
        .collect::<Vec<_>>();

    let fit =
        fit_local_dispersion_trend(&means, &disps, LocalDispersionTrendOptions::default()).unwrap();

    assert_eq!(
        fit.used_min_disp_floor,
        parse_required_bool(&rows[0], "usedMinDispFloor")
    );
    assert!(fit.genes_used > 0);
    assert!(fit.genes_used < rows.len());
    for (idx, row) in rows.iter().enumerate() {
        assert_eq!(fit.use_for_fit[idx], parse_required_bool(row, "useForFit"));
        assert_float_close(
            fit.disp_fit[idx],
            parse_required_f64(row, "dispFit"),
            2e-2,
            2e-2,
            &format!("local trend mixed threshold dispFit row {idx}"),
        );
    }
}
