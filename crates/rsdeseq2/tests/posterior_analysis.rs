//! Optional posterior queries, numerical references and validation behavior.
#[path = "common/posterior.rs"]
mod support;
use support::*;

#[test]
fn posterior_analysis_matches_independent_intervals_and_correlated_contrasts() {
    let reference: serde_json::Value =
        serde_json::from_str(include_str!("data/posterior_analysis_r462.json")).unwrap();
    let mut failures = vec![];
    for c in reference["cases"].as_array().unwrap() {
        let input: regression::Case = serde_json::from_value(c.clone()).unwrap();
        let r = match regression::run(&input) {
            Ok(r) => r,
            Err(e) => {
                failures.push(format!("{}: {e}", input.name));
                continue;
            }
        };
        let reorder = |m: &rsdeseq2::matrix::RowMajorMatrix<f64>| {
            (0..m.n_cols())
                .flat_map(|j| (0..m.n_rows()).map(move |i| *m.get(i, j).unwrap()))
                .collect::<Vec<_>>()
        };
        let gate = 0.01 * std::f64::consts::LN_2;
        for (field, m) in [
            ("mode", &r.mode),
            ("sd", r.standard_deviation.as_ref().unwrap()),
            ("contrast_mode", r.analysis.contrast_mode.as_ref().unwrap()),
            (
                "contrast_sd",
                r.analysis.contrast_standard_deviation.as_ref().unwrap(),
            ),
        ] {
            compare(
                &reorder(m),
                &c[field],
                gate,
                0.0,
                &format!("{} {field}", input.name),
                &mut failures,
            );
        }
        for (field, v) in [
            ("fsr", &r.false_sign_rate),
            ("sign_svalues", &r.sign_svalues),
        ] {
            compare(
                &v.as_ref()
                    .unwrap()
                    .iter()
                    .map(|v| v.unwrap_or(f64::from_bits(0x7ff00000000007a2)))
                    .collect::<Vec<_>>(),
                &c[field],
                0.02,
                0.0,
                &format!("{} {field}", input.name),
                &mut failures,
            );
        }
        let analysis = input.analysis.as_ref().unwrap();
        let interval = r.analysis.intervals.as_ref().unwrap();
        for gene in 0..input.observations.rows {
            // New mesh intervals have a predeclared one-cell approximation
            // gate. Existing MAP/SD/sign gates above are unchanged. HPD cutoff
            // membership can jump a cell for nearly tied density heights.
            let width = if analysis.kind == "laplace" {
                0.0
            } else {
                2.0 * parse_numeric_input(&analysis.grid_standard_deviations).unwrap()
                    * r.standard_deviation
                        .as_ref()
                        .unwrap()
                        .get(gene, input.coefficient)
                        .unwrap()
                    / analysis.grid_intervals as f64
            };
            for column in 0..2 {
                let repair = c["interval_repairs"].as_array().and_then(|rows| {
                    rows.iter().find(|v| {
                        v["gene"].as_u64() == Some(gene as u64)
                            && v["column"].as_u64() == Some(column as u64)
                    })
                });
                let expected = serde_json::json!([repair
                    .map(|v| &v["expected"])
                    .unwrap_or(&c["intervals"][column * input.observations.rows + gene])]);
                if repair.is_some() {
                    assert!(
                        *interval.get(gene, 0).unwrap()
                            < *r.mode.get(gene, input.coefficient).unwrap()
                    );
                }
                compare(
                    &[*interval.get(gene, column).unwrap()],
                    &expected,
                    gate + width,
                    0.0,
                    &format!("{} interval[{gene},{column}]", input.name),
                    &mut failures,
                );
            }
        }
        // Full covariance matters: the zero contrast must have zero variance,
        // and its estimate is zero even when the nonzero contrasts are correlated.
        for gene in 0..input.observations.rows {
            if r.diagnostics[gene].laplace_valid {
                assert_eq!(
                    *r.analysis
                        .contrast_mode
                        .as_ref()
                        .unwrap()
                        .get(gene, 2)
                        .unwrap(),
                    0.0
                );
                assert_eq!(
                    *r.analysis
                        .contrast_standard_deviation
                        .as_ref()
                        .unwrap()
                        .get(gene, 2)
                        .unwrap(),
                    0.0
                );
            }
        }
    }
    assert!(
        failures.is_empty(),
        "{} mismatches:\n{}",
        failures.len(),
        failures.join("\n")
    );
}

#[test]
fn posterior_analysis_controls_fail_early_and_map_only_queries_keep_missing_summaries() {
    let reference: serde_json::Value =
        serde_json::from_str(include_str!("data/posterior_analysis_r462.json")).unwrap();
    let base = reference["cases"]
        .as_array()
        .unwrap()
        .iter()
        .find(|v| v["name"] == "normal_credible")
        .unwrap();
    let mut variants = vec![];
    let mut add = |path: &str, value: serde_json::Value, message: &str| {
        let mut c = base.clone();
        *c.pointer_mut(path).unwrap() = value;
        variants.push((c, message.to_owned()));
    };
    add(
        "/analysis/grid_intervals",
        serde_json::json!(51),
        "integration grid",
    );
    add(
        "/analysis/level",
        serde_json::json!("-0.1"),
        "interval level",
    );
    add(
        "/analysis/contrasts/rows",
        serde_json::json!(3),
        "dimensions",
    );
    add(
        "/analysis/contrasts/values/0",
        serde_json::json!("NaN"),
        "posterior contrasts",
    );
    add(
        "/analysis/nuisance_points",
        serde_json::json!(0),
        "integration grid",
    );
    add(
        "/analysis/grid_intervals",
        serde_json::json!(usize::MAX - 1),
        "overflows",
    );
    add(
        "/analysis/parameter_standard_errors",
        serde_json::json!(["-1", "0", "0"]),
        "uncertainties",
    );
    for (c, message) in variants {
        let input: regression::Case = serde_json::from_value(c).unwrap();
        let error = regression::run(&input).unwrap_err().to_string();
        assert!(error.contains(&message), "expected {message}, got {error}");
    }
    let mut c = base.clone();
    c["analysis"]["parameter_standard_errors"] = serde_json::json!(["0.1", "0.1", "0.1"]);
    c["threshold"] = serde_json::json!("0.2");
    let error = regression::run(&serde_json::from_value(c.clone()).unwrap())
        .unwrap_err()
        .to_string();
    assert!(error.contains("threshold"));
    c["threshold"] = serde_json::Value::Null;
    c["analysis"]["kind"] = serde_json::json!("laplace");
    assert!(
        regression::run(&serde_json::from_value(c).unwrap())
            .unwrap_err()
            .to_string()
            .contains("grid interval")
    );
    let mut c = reference["cases"]
        .as_array()
        .unwrap()
        .iter()
        .find(|v| v["name"] == "negative_binomial_credible")
        .unwrap()
        .clone();
    c["method"] = serde_json::json!("nbinomC");
    c["analysis"]["grid_intervals"] = serde_json::json!(1);
    c["analysis"]["parameter_standard_errors"] = serde_json::json!(["0.1", "0.1", "0.1"]);
    let r = regression::run(&serde_json::from_value(c.clone()).unwrap()).unwrap();
    c["analysis"] = serde_json::Value::Null;
    let baseline = regression::run(&serde_json::from_value(c).unwrap()).unwrap();
    assert!(
        r.mode
            .as_slice()
            .iter()
            .zip(baseline.mode.as_slice())
            .all(|(a, b)| a.to_bits() == b.to_bits())
    );
    assert!(r.standard_deviation.is_none());
    for m in [
        r.analysis.intervals.unwrap(),
        r.analysis.contrast_mode.unwrap(),
        r.analysis.contrast_standard_deviation.unwrap(),
    ] {
        assert!(
            m.as_slice()
                .iter()
                .all(|v| v.to_bits() == 0x7ff00000000007a2)
        );
    }
    assert!(r.analysis.grid_diagnostics.is_empty());
}

#[test]
fn one_coefficient_grid_has_a_valid_zero_dimensional_nuisance_integral() {
    let reference: serde_json::Value =
        serde_json::from_str(include_str!("data/posterior_analysis_r462.json")).unwrap();
    let mut c = reference["cases"]
        .as_array()
        .unwrap()
        .iter()
        .find(|v| v["name"] == "normal_credible")
        .unwrap()
        .clone();
    c["design"]["columns"] = serde_json::json!(1);
    c["design"]["values"].as_array_mut().unwrap().truncate(8);
    c["coefficient"] = serde_json::json!(0);
    c["analysis"]["contrasts"] = serde_json::Value::Null;
    let input: regression::Case = serde_json::from_value(c).unwrap();
    let r = regression::run(&input).unwrap();
    let y = parse_numeric_inputs(&input.observations.values).unwrap();
    let variance = 1.0 / (8.0 / (0.7 * 0.7) + 1.0 / 225.0);
    for gene in 0..3 {
        let mean = (0..8).map(|i| y[i * 3 + gene]).sum::<f64>() / (0.7 * 0.7) * variance;
        assert!((r.mode.get(gene, 0).unwrap() - mean).abs() < 1e-6);
        assert!(
            (r.standard_deviation.as_ref().unwrap().get(gene, 0).unwrap() - variance.sqrt()).abs()
                < 1e-6
        );
        let interval = r.analysis.intervals.as_ref().unwrap();
        let width = 10.0 * variance.sqrt() / 50.0;
        assert!(
            (interval.get(gene, 0).unwrap() - (mean - 1.6448536269514722 * variance.sqrt())).abs()
                < width
        );
        assert!(
            (interval.get(gene, 1).unwrap() - (mean + 1.6448536269514722 * variance.sqrt())).abs()
                < width
        );
        assert_eq!(
            r.analysis.grid_diagnostics[gene]
                .as_ref()
                .unwrap()
                .evaluations,
            51
        );
    }
}

#[test]
fn posterior_work_preflight_is_batch_wide_and_cancellation_is_cooperative() {
    use rsdeseq2::shrinkage::{
        PosteriorAnalysisOptions, PosteriorIntervalKind, PosteriorIntervalOptions,
    };
    let mut analysis = PosteriorAnalysisOptions {
        interval: Some(PosteriorIntervalOptions {
            kind: PosteriorIntervalKind::EqualTail,
            ..Default::default()
        }),
        ..Default::default()
    };
    assert_eq!(
        analysis
            .estimated_grid_evaluations(2, 10, Some(1), false)
            .unwrap(),
        199_218_750
    );
    assert_eq!(
        analysis
            .estimated_grid_evaluations(2, 10, Some(1), true)
            .unwrap(),
        0
    );
    assert_eq!(
        analysis
            .estimated_grid_evaluations(2, 10, None, false)
            .unwrap(),
        0
    );
    analysis.parameter_standard_errors = Some(vec![0.1; 2]);
    assert_eq!(
        analysis
            .estimated_grid_evaluations(2, 10, Some(1), false)
            .unwrap(),
        996_093_750
    );
    // Exercise policy through the actual fit, rather than only a count helper.
    use rsdeseq2::matrix::{MatrixOrder, MatrixView};
    use rsdeseq2::shrinkage::{
        NegativeBinomialPosteriorInput, NegativeBinomialPosteriorOptions, PosteriorOffsets,
        fit_negative_binomial_posterior,
    };
    let counts = [10., 12., 20., 22.];
    let design = [1., 0., 1., 0., 1., 1., 1., 1.];
    let input = NegativeBinomialPosteriorInput {
        counts: MatrixView::new(1, 4, &counts, MatrixOrder::RowMajor).unwrap(),
        design: MatrixView::new(4, 2, &design, MatrixOrder::RowMajor).unwrap(),
        dispersions: &[0.1],
        offsets: PosteriorOffsets::None,
        weights: None,
        mle: None,
    };
    analysis.parameter_standard_errors = None;
    analysis.work.max_grid_evaluations = Some(254); // 51 * 5 = 255
    let mut options = NegativeBinomialPosteriorOptions {
        coefficient: Some(1),
        analysis,
        ..Default::default()
    };
    assert!(
        fit_negative_binomial_posterior(input, &options)
            .unwrap_err()
            .to_string()
            .contains("255 evaluations")
    );
    options.analysis.work.max_grid_evaluations = None;
    let accepted = fit_negative_binomial_posterior(input, &options).unwrap();
    assert!(accepted.diagnostics[0].laplace_valid);
    options.analysis.work.cancellation.cancel();
    assert!(
        fit_negative_binomial_posterior(input, &options)
            .unwrap_err()
            .to_string()
            .contains("cancelled")
    );
    // A frontend callback (such as R's interrupt poll) must remain confined to
    // its originating thread, while explicit cancellation reaches workers.
    static AFFINITY_POLLS: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
    fn observe_owner() -> bool {
        AFFINITY_POLLS.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        false
    }
    options.analysis.work =
        rsdeseq2::shrinkage::PosteriorWorkControl::default().with_interrupt_check(observe_owner);
    std::thread::scope(|scope| {
        let worker = scope.spawn(|| fit_negative_binomial_posterior(input, &options));
        assert!(worker.join().unwrap().is_ok());
    });
    assert_eq!(AFFINITY_POLLS.load(std::sync::atomic::Ordering::Relaxed), 0);
    assert!(fit_negative_binomial_posterior(input, &options).is_ok());
    let owner_polls = AFFINITY_POLLS.load(std::sync::atomic::Ordering::Relaxed);
    assert!(owner_polls > 0);
    let shared = options.analysis.work.cancellation.clone();
    std::thread::scope(|scope| {
        scope.spawn(|| shared.cancel()).join().unwrap();
        let worker = scope.spawn(|| fit_negative_binomial_posterior(input, &options));
        assert!(
            worker
                .join()
                .unwrap()
                .unwrap_err()
                .to_string()
                .contains("cancelled")
        );
    });
    assert_eq!(
        AFFINITY_POLLS.load(std::sync::atomic::Ordering::Relaxed),
        owner_polls
    );

    // Poll cancellation during an active grid, beyond entry and row preflight.
    static POLLS: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
    fn interrupt() -> bool {
        POLLS.fetch_add(1, std::sync::atomic::Ordering::Relaxed) >= 2
    }
    options.analysis.work =
        rsdeseq2::shrinkage::PosteriorWorkControl::default().with_interrupt_check(interrupt);
    assert!(
        fit_negative_binomial_posterior(input, &options)
            .unwrap_err()
            .to_string()
            .contains("cancelled")
    );
}
