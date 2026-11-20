use crate::support as common;
use common::hex_bits as bits;
use rsdeseq2::DeseqError;
use rsdeseq2::errors::{PosteriorArray, PosteriorArrayAxis};
use std::error::Error;
#[allow(dead_code)]
#[path = "../tools/validation/evaluate_r_interface/heavy_tail.rs"]
mod heavy_tail;
#[allow(dead_code)]
#[path = "../tools/validation/evaluate_r_interface/random.rs"]
mod random;
use rsdeseq2::shrinkage::estimate_cauchy_prior_variance;
fn parse_numeric_input(s: &str) -> Result<f64, Box<dyn Error>> {
    Ok(if let Some(s) = s.strip_prefix("bits:") {
        f64::from_bits(u64::from_str_radix(s, 16)?)
    } else if s == "NA" {
        f64::from_bits(0x7ff0_0000_0000_07a2)
    } else {
        s.parse()?
    })
}
fn parse_numeric_inputs(v: &[String]) -> Result<Vec<f64>, Box<dyn Error>> {
    v.iter().map(|s| parse_numeric_input(s)).collect()
}
fn compare(
    actual: &[f64],
    expected: &serde_json::Value,
    bound: f64,
    context: &str,
    failures: &mut Vec<String>,
) {
    let expected = expected.as_array().unwrap();
    assert_eq!(actual.len(), expected.len(), "{context}");
    for (i, (&a, b)) in actual.iter().zip(expected).enumerate() {
        let b = parse_numeric_input(b.as_str().unwrap()).unwrap();
        if a == b
            || (a.is_nan()
                && b.is_nan()
                && (a.to_bits() as u32 == 1954) == (b.to_bits() as u32 == 1954))
        {
            continue;
        }
        if a.is_finite() && b.is_finite() && (a - b).abs() <= bound {
            continue;
        }
        failures.push(format!(
            "{context}[{i}]: {a:.17e} != {b:.17e}; error {:.3e}",
            (a - b).abs()
        ));
    }
}
#[test]
fn negative_binomial_cauchy_posteriors_match_original_apeglm() {
    let reference: serde_json::Value =
        serde_json::from_str(include_str!("data/cauchy_posterior_r462.json")).unwrap();
    let mut failures = vec![];
    let cases = reference["cases"].as_array().unwrap();
    assert_eq!(cases.len(), 13);
    for c in cases {
        let input: heavy_tail::Case = serde_json::from_value(c.clone()).unwrap();
        let r = heavy_tail::run(&input).unwrap();
        let g = input.counts.rows;
        let p = input.design.columns;
        // Original matrices are stored column-major; native outputs row-major.
        let reorder = |m: &rsdeseq2::matrix::RowMajorMatrix<f64>| {
            (0..p)
                .flat_map(|j| (0..g).map(move |i| *m.get(i, j).unwrap()))
                .collect::<Vec<_>>()
        };
        // Existing public posterior bounds are 0.01 log2 units. Inputs and
        // numerical outputs here use natural logs; the same gate is rescaled.
        compare(
            &reorder(&r.mode),
            &c["mode"],
            0.01 * std::f64::consts::LN_2,
            &format!("{} mode", input.name),
            &mut failures,
        );
        let sd = r
            .standard_deviation
            .as_ref()
            .map(reorder)
            .unwrap_or_else(|| vec![f64::from_bits(0x7ff0_0000_0000_07a2); g * p]);
        compare(
            &sd,
            &c["sd"],
            0.01 * std::f64::consts::LN_2,
            &format!("{} sd", input.name),
            &mut failures,
        );
        for (field, values, bound) in [
            ("fsr", r.false_sign_rate, 0.02),
            ("sign_svalues", r.sign_svalues, 0.02),
            ("fsos", r.false_sign_or_small_rate, 0.05),
            ("threshold_svalues", r.threshold_svalues, 0.05),
        ] {
            let values = values
                .unwrap()
                .into_iter()
                .map(|v| v.unwrap_or(f64::from_bits(0x7ff0_0000_0000_07a2)))
                .collect::<Vec<_>>();
            compare(
                &values,
                &c[field],
                bound,
                &format!("{} {field}", input.name),
                &mut failures,
            );
        }
        let expected = parse_numeric_input(c["fitted_scale"].as_str().unwrap()).unwrap();
        assert!(
            (r.prior.scale - expected).abs() <= 1e-12,
            "{} prior scale {} != {}",
            input.name,
            r.prior.scale,
            expected
        );
    }
    assert!(
        failures.is_empty(),
        "{} numerical mismatches:\n{}",
        failures.len(),
        failures.join("\n")
    );
}
#[test]
fn adaptive_prior_variance_matches_original_bracket_and_missing_policy() {
    let reference: serde_json::Value =
        serde_json::from_str(include_str!("data/cauchy_posterior_r462.json")).unwrap();
    let cases = reference["priors"].as_array().unwrap();
    assert_eq!(cases.len(), 6);
    let mut failures = vec![];
    for c in cases {
        let input: heavy_tail::PriorCase = serde_json::from_value(c.clone()).unwrap();
        let actual = estimate_cauchy_prior_variance(
            &parse_numeric_inputs(&input.estimates).unwrap(),
            &parse_numeric_inputs(&input.standard_errors).unwrap(),
        )
        .unwrap();
        let expected = parse_numeric_input(c["variance"].as_str().unwrap()).unwrap();
        compare(
            &[actual],
            &serde_json::json!([c["variance"]]),
            1e-12 * (1.0 + expected.abs()),
            &input.name,
            &mut failures,
        );
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
    assert!(estimate_cauchy_prior_variance(&[f64::NAN], &[1.0]).is_err());
    assert!(estimate_cauchy_prior_variance(&[1.0], &[f64::NAN]).is_err());
}

#[test]
fn invalid_rows_and_degenerate_scales_error_while_zero_rows_remain_missing() {
    use rsdeseq2::{
        matrix::{MatrixOrder, MatrixView},
        shrinkage::{
            NegativeBinomialPosteriorInput, NegativeBinomialPosteriorOptions, PosteriorOffsets,
            fit_negative_binomial_posterior,
        },
    };
    let counts = [1., 2.];
    let dispersions = [0.1];
    let design = [1., 0., 1., 1.];
    let input = NegativeBinomialPosteriorInput {
        counts: MatrixView::new(1, 2, &counts, MatrixOrder::RowMajor).unwrap(),
        design: MatrixView::new(2, 2, &design, MatrixOrder::RowMajor).unwrap(),
        dispersions: &dispersions,
        offsets: PosteriorOffsets::None,
        weights: None,
        mle: None,
    };
    let options = NegativeBinomialPosteriorOptions {
        coefficient: Some(1),
        ..Default::default()
    };
    for scale in [0., f64::NAN, f64::INFINITY, 1e-300, 1e300] {
        let mut bad = options.clone();
        bad.prior.scale = scale;
        assert!(fit_negative_binomial_posterior(input, &bad).is_err());
    }
    for dispersion in [f64::NAN, 0., -0.1] {
        assert!(
            fit_negative_binomial_posterior(
                NegativeBinomialPosteriorInput {
                    dispersions: &[dispersion],
                    ..input
                },
                &options
            )
            .is_err()
        );
    }
    for y in [[-1., 2.], [f64::NAN, 2.], [f64::MAX, f64::MAX]] {
        assert!(
            fit_negative_binomial_posterior(
                NegativeBinomialPosteriorInput {
                    counts: MatrixView::new(1, 2, &y, MatrixOrder::RowMajor).unwrap(),
                    ..input
                },
                &options
            )
            .is_err()
        );
    }
    let zero = fit_negative_binomial_posterior(
        NegativeBinomialPosteriorInput {
            counts: MatrixView::new(1, 2, &[0., 0.], MatrixOrder::RowMajor).unwrap(),
            dispersions: &[f64::NAN],
            ..input
        },
        &options,
    )
    .unwrap();
    assert!(
        zero.mode
            .as_slice()
            .iter()
            .all(|x| x.is_nan() && x.to_bits() as u32 == 1954)
    );
    assert_eq!(zero.false_sign_rate, Some(vec![None]));
    assert!(zero.diagnostics[0].convergence.is_none());
}

#[test]
fn posterior_summaries_keep_coefficient_identity_and_reject_inconsistent_state() {
    use rsdeseq2::matrix::{MatrixOrder, MatrixView, RowMajorMatrix};
    use rsdeseq2::shrinkage::{
        CoefficientScale, NegativeBinomialPosteriorInput, NegativeBinomialPosteriorOptions,
        PosteriorOffsets, fit_negative_binomial_posterior,
    };
    let counts = [
        100., 105., 98., 110., 108., 102., 10., 12., 9., 80., 85., 78.,
    ];
    let design = [1., 0., 1., 0., 1., 0., 1., 1., 1., 1., 1., 1.];
    let input = NegativeBinomialPosteriorInput {
        counts: MatrixView::new(2, 6, &counts, MatrixOrder::RowMajor).unwrap(),
        design: MatrixView::new(6, 2, &design, MatrixOrder::RowMajor).unwrap(),
        dispersions: &[0.1, 0.1],
        offsets: PosteriorOffsets::None,
        weights: None,
        mle: None,
    };
    let options = NegativeBinomialPosteriorOptions {
        coefficient: Some(1),
        ..Default::default()
    };
    let fit = fit_negative_binomial_posterior(input, &options).unwrap();
    assert_eq!(fit.summary_coefficient(), Some(1));
    let summary = fit
        .coefficient_summary(1, CoefficientScale::Natural)
        .unwrap();
    assert_eq!(summary.sign_svalues, fit.sign_svalues);
    std::assert_matches!(
        fit.coefficient_summary(0, CoefficientScale::Natural)
            .unwrap_err(),
        DeseqError::PosteriorSummaryCoefficientMismatch {
            requested: 0,
            fitted: Some(1)
        }
    );
    // The data distinguish the two scientific summaries; a wrong-column test
    // would be ineffective if both coefficients happened to have equal s-values.
    let other = fit_negative_binomial_posterior(
        input,
        &NegativeBinomialPosteriorOptions {
            coefficient: Some(0),
            ..options
        },
    )
    .unwrap();
    assert_eq!(fit.mode, other.mode);
    assert_ne!(fit.sign_svalues, other.sign_svalues);
    use PosteriorArrayAxis::{Columns, Rows};
    let failures = [
        (PosteriorArray::StandardDeviation, Rows, 1),
        (PosteriorArray::SignSvalues, Rows, 0),
        (PosteriorArray::ThresholdSvalues, Rows, 0),
        (PosteriorArray::FalseSignRate, Rows, 0),
        (PosteriorArray::FalseSignOrSmallRate, Rows, 0),
        (PosteriorArray::Diagnostics, Rows, 0),
        (PosteriorArray::Intervals, Columns, 1),
        (PosteriorArray::StandardDeviation, Columns, 1),
        (PosteriorArray::ContrastStandardDeviation, Columns, 1),
        (PosteriorArray::GridDiagnostics, Rows, 1),
    ];
    for (damage, (expected_array, expected_axis, expected_actual)) in
        failures.into_iter().enumerate()
    {
        let mut broken = fit.clone();
        match damage {
            0 => {
                broken.standard_deviation =
                    Some(RowMajorMatrix::from_row_major(1, 1, vec![1.]).unwrap())
            }
            1 => broken.sign_svalues = Some(vec![]),
            2 => broken.threshold_svalues = Some(vec![]),
            3 => broken.false_sign_rate = Some(vec![]),
            4 => broken.false_sign_or_small_rate = Some(vec![]),
            5 => broken.diagnostics.clear(),
            6 => {
                broken.analysis.intervals =
                    Some(RowMajorMatrix::from_row_major(2, 1, vec![1.; 2]).unwrap())
            }
            7 => broken.standard_deviation = Some(RowMajorMatrix::from_elem(2, 1, 1.).unwrap()),
            8 => {
                broken.analysis.contrast_mode = Some(RowMajorMatrix::from_elem(2, 2, 1.).unwrap());
                broken.analysis.contrast_standard_deviation =
                    Some(RowMajorMatrix::from_elem(2, 1, 1.).unwrap());
            }
            _ => broken.analysis.grid_diagnostics = vec![None],
        }
        let error = broken
            .coefficient_summary(1, CoefficientScale::Natural)
            .unwrap_err();
        if damage == 8 {
            assert_eq!(
                error.to_string(),
                "invalid options: heavy-tailed shrinkage: posterior contrast arrays have inconsistent shapes",
            );
        }
        std::assert_matches!(
            error,
            DeseqError::PosteriorShapeMismatch { array, axis, expected: 2, actual, .. }
                if (array, axis, actual) == (expected_array, expected_axis, expected_actual)
        );
    }
}
