use rsdeseq2::prelude::*;

use approx::assert_abs_diff_eq;

#[test]
fn default_theta_starts_at_fraction_of_zero_filter_values() {
    let theta = default_theta(&[0.0, 0.0, 10.0, 20.0]).unwrap();

    assert_eq!(theta.len(), 50);
    assert_eq!(theta[0], 0.5);
    assert_eq!(theta[49], 0.95);
}

#[test]
fn filtered_p_adjustments_match_filtered_bh_shape() {
    let filter = vec![0.0, 1.0, 2.0];
    let pvalues = vec![Some(0.01), Some(0.02), Some(0.03)];
    let (columns, num_rej) = filtered_p_adjustments(&filter, &pvalues, &[1.0], 0.05).unwrap();

    assert_eq!(columns.len(), 1);
    assert_eq!(columns[0][0], None);
    assert_eq!(columns[0][1], Some(0.03));
    assert_eq!(columns[0][2], Some(0.03));
    assert_eq!(num_rej, vec![2]);
}

#[test]
fn independent_filtering_applies_selected_threshold_and_metadata() {
    let fit = toy_fit(vec![1.0, 1.0, 1.0], vec![1.0, 1.0, 1.0], vec![true; 3]);
    let mut results = build_wald_results(&[0.0, 1.0, 2.0], &fit, 0, None, None).unwrap();
    let output = apply_independent_filtering(
        &mut results,
        &IndependentFilteringOptions {
            enabled: true,
            alpha: 0.5,
            theta: Some(vec![1.0, 1.0]),
        },
    )
    .unwrap();

    assert!(output.enabled);
    assert_eq!(output.selected_index, Some(0));
    assert_eq!(output.filter_theta, Some(1.0));
    assert_eq!(output.filter_threshold, Some(2.0));
    assert_eq!(results.independent_filtering, Some(output));
    assert_eq!(results.rows[0].filtered, Some(true));
    assert_eq!(results.rows[1].filtered, Some(true));
    assert_eq!(results.rows[2].filtered, Some(false));
    assert_eq!(results.rows[0].padj, None);
    assert_eq!(results.rows[1].padj, None);
    assert!(results.rows[2].padj.is_some());
    assert_eq!(
        results
            .independent_filtering
            .as_ref()
            .unwrap()
            .lowess_fit
            .as_ref()
            .unwrap()
            .len(),
        2
    );
}

#[test]
fn independent_filtering_disabled_recomputes_regular_bh() {
    let fit = toy_fit(vec![2.0, 1.0], vec![0.5, 1.0], vec![true; 2]);
    let mut results = build_wald_results(&[0.0, 10.0], &fit, 0, None, None).unwrap();
    results.rows[0].padj = None;
    results.rows[0].filtered = Some(true);

    let output = apply_independent_filtering(
        &mut results,
        &IndependentFilteringOptions {
            enabled: false,
            ..IndependentFilteringOptions::default()
        },
    )
    .unwrap();

    assert!(!output.enabled);
    assert_eq!(results.independent_filtering, Some(output));
    assert_eq!(results.rows[0].filtered, None);
    assert_eq!(results.rows[1].filtered, None);
    assert!(results.rows[0].padj.is_some());
    assert!(results.rows[1].padj.is_some());
    assert_eq!(
        results.independent_filtering.as_ref().unwrap().lowess_fit,
        None
    );
}

#[test]
#[allow(clippy::excessive_precision)]
fn lowess_fit_matches_r_stats_lowess_for_default_filter_grid_fixture() {
    let theta = (0..50)
        .map(|idx| 0.95 * idx as f64 / 49.0)
        .collect::<Vec<_>>();
    let mut num_rejections = vec![0; 10];
    num_rejections.extend(1..=20);
    num_rejections.extend(vec![30; 20]);
    let expected = vec![
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.020076251918791414,
        0.12340003422977887,
        0.3683849514569546,
        0.78200893423769502,
        1.3683849514569535,
        2.1234000342297783,
        3.0200762519187925,
        3.9999999999999987,
        4.9999999999999982,
        6.0,
        6.9999999999999982,
        7.9999999999999982,
        9.0000000000000018,
        9.9999999999999982,
        11.000000000000002,
        11.999999999999998,
        12.999999999999998,
        13.999999999999996,
        15.000000000000005,
        15.999999999999998,
        17.180686267269127,
        18.909837788880086,
        21.081464220814784,
        23.354230893569678,
        25.495375220735649,
        27.426750793498471,
        28.946685924971323,
        29.799237480812085,
        30.000000000000004,
        30.0,
        30.0,
        29.999999999999993,
        30.000000000000007,
        30.000000000000004,
        30.000000000000004,
        30.0,
        30.0,
        29.999999999999996,
        30.000000000000004,
        29.999999999999996,
        30.000000000000011,
        29.999999999999993,
        29.999999999999961,
        30.000000000000263,
    ];

    let fitted = lowess_fitted_values(&theta, &num_rejections, 0.2).unwrap();
    assert_eq!(fitted.len(), expected.len());
    for (actual, expected) in fitted.into_iter().zip(expected) {
        assert_abs_diff_eq!(actual, expected, epsilon = 1e-9);
    }
}

#[test]
fn selected_filter_index_matches_deseq2_lowess_rule_fixture() {
    let theta = (0..50)
        .map(|idx| 0.95 * idx as f64 / 49.0)
        .collect::<Vec<_>>();
    let mut num_rejections = vec![0; 10];
    num_rejections.extend(1..=20);
    num_rejections.extend(vec![30; 20]);

    let (selected_index, lowess_fit) =
        select_filter_index_with_lowess(&theta, &num_rejections).unwrap();

    assert_eq!(selected_index, 30);
    assert_eq!(select_filter_index(&theta, &num_rejections), 30);
    assert_eq!(lowess_fit.len(), theta.len());
}

#[test]
#[allow(clippy::excessive_precision)]
fn lowess_fit_matches_r_stats_lowess_for_dense_custom_theta_fixture() {
    let theta = (0..200)
        .map(|idx| 0.95 * idx as f64 / 199.0)
        .collect::<Vec<_>>();
    let mut num_rejections = (0..200)
        .map(|idx| ((70.0 * idx as f64 / 199.0).floor() as usize).min(40))
        .collect::<Vec<_>>();
    for value in num_rejections.iter_mut().take(30) {
        *value = 0;
    }
    for value in num_rejections.iter_mut().skip(120) {
        *value = 40;
    }
    let expected = [
        (0, 0.0),
        (1, 0.0),
        (2, 0.0),
        (3, 0.0),
        (4, 0.0),
        (5, 0.0),
        (6, 0.0),
        (7, 0.0),
        (24, 3.7888632462679142),
        (25, 4.7417355283352274),
        (26, 5.4626941789122077),
        (27, 6.0444099585334268),
        (28, 6.5544558252245286),
        (29, 7.0310065934337649),
        (30, 7.4969326772329419),
        (31, 7.9687390603411865),
        (79, 27.290618749886633),
        (80, 27.642230661313068),
        (81, 27.993415579683759),
        (82, 28.344177994321473),
        (83, 28.69458092245068),
        (84, 29.04448976205753),
        (85, 29.393902011800801),
        (86, 29.742812475621754),
        (87, 30.090999210496296),
        (119, 38.643807673194146),
        (120, 38.840932813678222),
        (121, 39.054516100171881),
        (122, 39.287711287885969),
        (123, 39.529510153052243),
        (124, 39.741374126571941),
        (125, 39.890326961689254),
        (126, 39.970392181969849),
        (127, 39.996757741446025),
        (189, 40.00000000000005),
        (190, 40.000000000000085),
        (191, 39.999999999999957),
        (192, 39.999999999999936),
        (193, 39.999999999999943),
        (194, 40.000000000000206),
        (195, 40.000000000000171),
        (196, 40.000000000000249),
        (197, 39.999999999999893),
        (198, 40.000000000000036),
        (199, 39.999999999999872),
    ];

    let fitted = lowess_fitted_values(&theta, &num_rejections, 0.2).unwrap();

    assert_eq!(fitted.len(), theta.len());
    for (idx, expected) in expected {
        assert_abs_diff_eq!(fitted[idx], expected, epsilon = 1e-9);
    }
}

#[test]
fn selected_filter_index_matches_deseq2_dense_custom_theta_fixture() {
    let theta = (0..200)
        .map(|idx| 0.95 * idx as f64 / 199.0)
        .collect::<Vec<_>>();
    let mut num_rejections = (0..200)
        .map(|idx| ((70.0 * idx as f64 / 199.0).floor() as usize).min(40))
        .collect::<Vec<_>>();
    for value in num_rejections.iter_mut().take(30) {
        *value = 0;
    }
    for value in num_rejections.iter_mut().skip(120) {
        *value = 40;
    }

    let (selected_index, lowess_fit) =
        select_filter_index_with_lowess(&theta, &num_rejections).unwrap();

    assert_eq!(selected_index, 114);
    assert_eq!(select_filter_index(&theta, &num_rejections), 114);
    assert_eq!(lowess_fit.len(), theta.len());
}

#[test]
fn independent_filtering_rejects_invalid_options() {
    let fit = toy_fit(vec![1.0], vec![1.0], vec![true]);
    let mut results = build_wald_results(&[1.0], &fit, 0, None, None).unwrap();

    assert!(apply_independent_filtering(
        &mut results,
        &IndependentFilteringOptions {
            alpha: 1.0,
            ..IndependentFilteringOptions::default()
        },
    )
    .is_err());
    assert!(apply_independent_filtering(
        &mut results,
        &IndependentFilteringOptions {
            theta: Some(vec![0.5]),
            ..IndependentFilteringOptions::default()
        },
    )
    .is_err());
}

fn toy_fit(beta: Vec<f64>, beta_se: Vec<f64>, beta_converged: Vec<bool>) -> NbinomGlmFit {
    let n_genes = beta_converged.len();
    let n_samples = 2;
    NbinomGlmFit {
        log_like: vec![0.0; n_genes],
        beta_converged,
        beta: RowMajorMatrix::from_row_major(n_genes, 1, beta).unwrap(),
        beta_se: RowMajorMatrix::from_row_major(n_genes, 1, beta_se).unwrap(),
        beta_covariance: None,
        mu: RowMajorMatrix::from_row_major(n_genes, n_samples, vec![1.0; n_genes * n_samples])
            .unwrap(),
        beta_iter: vec![1; n_genes],
        model_matrix: DesignMatrix::from_row_major(n_samples, 1, vec![1.0, 1.0], None).unwrap(),
        n_terms: 1,
        hat_diagonal: RowMajorMatrix::from_row_major(
            n_genes,
            n_samples,
            vec![0.5; n_genes * n_samples],
        )
        .unwrap(),
    }
}
