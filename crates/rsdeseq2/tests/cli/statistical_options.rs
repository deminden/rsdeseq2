use super::numerical_assertions::*;
use super::*;

fn distinct_matrix(a: &RowMajorMatrix<f64>, b: &RowMajorMatrix<f64>) -> bool {
    numerically_distinct(
        &a.as_slice().iter().copied().map(Some).collect::<Vec<_>>(),
        &b.as_slice().iter().copied().map(Some).collect::<Vec<_>>(),
    )
}

#[test]
fn frozen_rlog_options_match_library_and_change_transformed_values() {
    let dir = temp_dir("frozen-rlog-numerical");
    let intercept_path = dir.join("intercept.tsv");
    write_frozen_intercept_fixture(&intercept_path);
    let counts = common::reference_counts();
    let intercept = [3.5, -0.5, 2.75, 4.25];
    let output = dir.join("rlog.tsv");
    let builder = test_builder();
    let mut blind_expected = None;
    for blind in [true, false] {
        let design = if blind {
            DesignMatrix::intercept_only(4).unwrap()
        } else {
            common::reference_full_design()
        };
        let fit = builder
            .fit_map_dispersions_glm_mu(&counts, &design)
            .unwrap();
        let expected = fit
            .frozen_rlog(&counts, &intercept, 2.5)
            .unwrap()
            .transformed;
        let ordinary = fit.rlog(&counts).unwrap().transformed;
        assert!(
            distinct_matrix(&expected, &ordinary),
            "ignored frozen intercept must change answers"
        );
        let wrong_prior = fit
            .frozen_rlog(&counts, &intercept, 0.1)
            .unwrap()
            .transformed;
        assert!(
            distinct_matrix(&expected, &wrong_prior),
            "ignored prior variance must change answers"
        );
        if let Some(blind_expected) = &blind_expected {
            assert!(
                distinct_matrix(&expected, blind_expected),
                "ignored blind=false must change answers"
            );
        } else {
            blind_expected = Some(expected.clone());
        }
        let label = format!("frozen rlog blind={blind}");
        let mut args = basic_args("rlog", &reference_data_path("counts.tsv"), &output);
        args.extend([
            "--fit-type".into(),
            "mean".into(),
            "--frozen-intercept".into(),
            intercept_path.to_str().unwrap().into(),
            "--rlog-prior-variance".into(),
            "2.5".into(),
        ]);
        if !blind {
            args.extend([
                "--blind=false".into(),
                "--design".into(),
                reference_data_path("design_full.tsv")
                    .to_str()
                    .unwrap()
                    .into(),
            ]);
        }
        invoke(&label, &args);
        assert_matrix_export(&output, &counts, &expected, &label);
    }
}

#[test]
fn wald_threshold_and_t_tail_options_match_library_numerically() {
    let dir = temp_dir("wald-statistical-numerical");
    let output = dir.join("wald.tsv");
    let df_path = dir.join("df.tsv");
    // Different degrees of freedom in shuffled gene order detect row misalignment.
    fs::write(
        &df_path,
        "gene\tdf\ngene3\t7\ngene1\t3\ngene4\t9\ngene2\t5\n",
    )
    .unwrap();
    let counts = common::reference_counts();
    let design = common::reference_full_design();
    let result = |builder: &DeseqBuilder| builder.fit_with_results(&counts, &design).unwrap().1;
    let normal = result(&test_builder());
    let scalar = result(&test_builder().wald_t_degrees_of_freedom(4.0));
    let cases = [
        (
            "threshold and greater alternative",
            vec!["--lfc-threshold", "0.5", "--alternative", "greater"],
            test_builder().wald_lfc_threshold(0.5, WaldAlternative::Greater),
        ),
        (
            "residual t degrees of freedom",
            vec!["--use-t"],
            test_builder().wald_t_residual_degrees_of_freedom(),
        ),
        (
            "scalar t degrees of freedom",
            vec!["--t-degrees-of-freedom", "4"],
            test_builder().wald_t_degrees_of_freedom(4.0),
        ),
        (
            "shuffled per-gene t degrees of freedom",
            vec!["--t-degrees-of-freedom-file", df_path.to_str().unwrap()],
            test_builder().wald_t_per_gene_degrees_of_freedom(vec![3.0, 5.0, 7.0, 9.0]),
        ),
    ];
    for (label, flags, builder) in cases {
        let expected = result(&builder);
        assert!(
            numerically_distinct(&result_values(&expected), &result_values(&normal)),
            "{label}: ignored flag undetectable"
        );
        if label.starts_with("threshold") {
            for ignored in [
                test_builder().wald_lfc_threshold(0.0, WaldAlternative::Greater),
                test_builder().wald_lfc_threshold(0.5, WaldAlternative::GreaterAbs),
            ] {
                assert!(
                    numerically_distinct(
                        &result_values(&expected),
                        &result_values(&result(&ignored))
                    ),
                    "{label}: each option must matter"
                );
            }
        } else if label.starts_with("residual") || label.starts_with("shuffled") {
            assert!(
                numerically_distinct(&result_values(&expected), &result_values(&scalar)),
                "{label}: must detect a wrong scalar df"
            );
        }
        let mut args = basic_args("wald", &reference_data_path("counts.tsv"), &output);
        model_args(
            &mut args,
            "wald",
            &reference_data_path("design_full.tsv"),
            &reference_data_path("design_reduced.tsv"),
        );
        args.extend(["--coefficient".into(), "1".into()]);
        args.extend(flags.into_iter().map(str::to_owned));
        invoke(label, &args);
        assert_result_export(&output, &expected, label);
    }
}
