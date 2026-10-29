//! Independent probability references and fitted/primitive inference conformance.
use rsdeseq2::prelude::*;
use serde::Deserialize;

#[derive(Deserialize)]
struct Reference {
    family: String,
    input_bits: String,
    df_bits: String,
    samples: usize,
    coefficients: usize,
    expected_bits: String,
}
fn value(bits: &str) -> f64 {
    f64::from_bits(u64::from_str_radix(bits, 16).unwrap())
}

#[test]
fn standalone_probability_primitives_match_pinned_r462_references() {
    let mut table = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(include_bytes!("data/inference_r462.tsv").as_slice());
    let mut families = std::collections::BTreeMap::new();
    for row in table.deserialize::<Reference>() {
        let r = row.unwrap();
        let q = value(&r.input_bits);
        let df = value(&r.df_bits);
        let expected = value(&r.expected_bits);
        let actual = match r.family.as_str() {
            "normal_upper" | "normal_two_sided" | "t_upper" => {
                let options = if r.family == "t_upper" {
                    WaldTestOptions::t_degrees_of_freedom(df)
                } else {
                    WaldTestOptions::normal()
                };
                let options = options.with_lfc_threshold(
                    0.0,
                    if r.family == "normal_two_sided" {
                        WaldAlternative::GreaterAbs
                    } else {
                        WaldAlternative::Greater
                    },
                );
                wald_test_from_estimates(&[q], &[1.0], &options)
                    .unwrap()
                    .pvalue[0]
                    .unwrap()
            }
            "chi_square" => lrt_pvalues_from_statistics(&[q], &[df]).unwrap()[0].unwrap(),
            "cooks" => rsdeseq2::results::default_cooks_cutoff(r.samples, r.coefficients)
                .unwrap()
                .unwrap(),
            other => panic!("unrecognized reference family {other}"),
        };
        // Relative tail accuracy matters even far below ordinary p-value
        // cutoffs. The absolute floor only permits subnormal rounding.
        let relative = if r.family.starts_with("normal") {
            8.0 * f64::EPSILON
        } else {
            2e-12
        };
        assert!(
            (actual - expected).abs() <= (relative * expected.abs()).max(4.0 * f64::from_bits(1)),
            "{} q={q} df={df} m={} p={}: actual={actual:.17e}, R={expected:.17e}",
            r.family,
            r.samples,
            r.coefficients
        );
        *families.entry(r.family).or_insert(0usize) += 1;
    }
    assert_eq!(families.len(), 5);
}

#[test]
fn primitive_estimates_and_fitted_coefficients_have_identical_bits() {
    let counts =
        CountMatrix::from_row_major_u32(3, 4, vec![10, 12, 20, 22, 40, 38, 20, 18, 4, 8, 7, 9])
            .unwrap();
    let design =
        DesignMatrix::from_row_major(4, 2, vec![1., 0., 1., 0., 1., 1., 1., 1.], None).unwrap();
    let fit = fit_irls(
        &counts,
        &design,
        &[1.; 4],
        &[0.1; 3],
        IrlsOptions::default(),
    )
    .unwrap();
    let beta = (0..3)
        .map(|g| fit.beta.row(g).unwrap()[1])
        .collect::<Vec<_>>();
    let se = (0..3)
        .map(|g| fit.beta_se.row(g).unwrap()[1])
        .collect::<Vec<_>>();
    let bits = |v: Vec<Option<f64>>| {
        v.into_iter()
            .map(|x| x.map(f64::to_bits))
            .collect::<Vec<_>>()
    };
    for alternative in [
        WaldAlternative::GreaterAbs,
        WaldAlternative::GreaterAbsUpshot,
        WaldAlternative::GreaterAbs2014,
        WaldAlternative::LessAbs,
        WaldAlternative::Greater,
        WaldAlternative::Less,
    ] {
        for t in [false, true] {
            if t && alternative == WaldAlternative::GreaterAbsUpshot {
                continue;
            }
            let options = if t {
                WaldTestOptions::t_per_gene_degrees_of_freedom(vec![10., f64::NAN, 2.])
            } else {
                WaldTestOptions::normal()
            }
            .with_lfc_threshold(0.75, alternative);
            let fitted = wald_test_coefficient_with_options(&fit, 1, &options).unwrap();
            let primitive = wald_test_from_estimates(&beta, &se, &options).unwrap();
            assert_eq!(bits(fitted.stat), bits(primitive.stat));
            assert_eq!(bits(fitted.pvalue), bits(primitive.pvalue));
            assert_eq!(fitted.degrees_of_freedom, primitive.degrees_of_freedom);
        }
    }
}

#[test]
fn missing_t_df_retains_the_thresholded_statistic_and_invalid_estimates_stay_missing() {
    let alternatives = [
        (WaldAlternative::GreaterAbs, -3.0),
        (WaldAlternative::GreaterAbs2014, -1.0),
        (WaldAlternative::LessAbs, 0.0),
        (WaldAlternative::Greater, 0.0),
        (WaldAlternative::Less, -1.0),
    ];
    for (alternative, expected) in alternatives {
        let options = WaldTestOptions::t_per_gene_degrees_of_freedom(vec![f64::NAN, 0., 10., 10.])
            .with_lfc_threshold(1.0, alternative);
        let result = wald_test_from_estimates(
            &[-1.5, -1.5, f64::INFINITY, 1.],
            &[0.5, 0.5, 1., 0.],
            &options,
        )
        .unwrap();
        assert_eq!(
            result.stat,
            vec![Some(expected), Some(expected), None, None]
        );
        assert_eq!(result.pvalue, vec![None; 4]);
    }
    assert_eq!(
        lrt_pvalues_from_statistics(&[f64::NAN, 6., 6., -1.], &[2., f64::NAN, 0., 2.]).unwrap(),
        vec![None, None, None, Some(1.)]
    );
}

#[test]
fn primitive_shapes_and_model_dependent_df_are_validated_before_evaluation() {
    assert!(wald_test_from_estimates(&[1.], &[], &WaldTestOptions::normal()).is_err());
    assert!(
        wald_test_from_estimates(
            &[1., 2.],
            &[1., 1.],
            &WaldTestOptions::t_per_gene_degrees_of_freedom(vec![4.])
        )
        .is_err()
    );
    assert!(
        wald_test_from_estimates(
            &[1.],
            &[1.],
            &WaldTestOptions::t_residual_degrees_of_freedom()
        )
        .is_err()
    );
    assert!(lrt_pvalues_from_statistics(&[1., 2.], &[]).is_err());
    assert!(lrt_pvalues_from_statistics(&[1., 2.], &[1., 2., 3.]).is_err());
    assert!(
        wald_test_from_estimates(&[], &[], &WaldTestOptions::normal())
            .unwrap()
            .pvalue
            .is_empty()
    );
    assert!(lrt_pvalues_from_statistics(&[], &[]).unwrap().is_empty());
}
