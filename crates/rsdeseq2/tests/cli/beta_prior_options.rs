use super::numerical_assertions::*;
use super::*;

#[derive(Clone, Copy, Debug)]
enum Family {
    Expanded,
    Factor,
    Additive,
}

struct PriorFixture {
    counts: CountMatrix,
    standard: DesignMatrix,
    expanded: DesignMatrix,
    groups: Vec<Vec<usize>>,
    factors: RowMajorMatrix<f64>,
    weights: RowMajorMatrix<f64>,
    size_factors: Vec<f64>,
}

impl PriorFixture {
    fn new(additive: bool) -> Self {
        let counts = common::reference_counts();
        let names = if additive {
            vec!["Intercept", "condition_B_vs_A", "batch_Y_vs_X"]
        } else {
            vec!["Intercept", "condition_B_vs_A"]
        };
        let expanded_names = if additive {
            vec!["Intercept", "conditionA", "conditionB", "batchX", "batchY"]
        } else {
            vec!["Intercept", "conditionA", "conditionB"]
        };
        // Explicit matrices keep the oracle independent of CLI factor parsing
        // and the generated-design entry points being exercised.
        let standard = DesignMatrix::from_row_major(
            4,
            names.len(),
            (0..4)
                .flat_map(|i| {
                    let mut row = vec![1.0, f64::from(i >= 2)];
                    if additive {
                        row.push(f64::from(i == 1 || i == 2));
                    }
                    row
                })
                .collect(),
            Some(names.into_iter().map(str::to_owned).collect()),
        )
        .unwrap();
        let expanded = DesignMatrix::from_row_major(
            4,
            expanded_names.len(),
            (0..4)
                .flat_map(|i| {
                    let mut row = vec![1.0, f64::from(i < 2), f64::from(i >= 2)];
                    if additive {
                        row.extend([f64::from(i == 0 || i == 3), f64::from(i == 1 || i == 2)]);
                    }
                    row
                })
                .collect(),
            Some(expanded_names.into_iter().map(str::to_owned).collect()),
        )
        .unwrap();
        Self {
            size_factors: estimate_size_factors_ratio(&counts).unwrap(),
            counts,
            standard,
            expanded,
            groups: if additive {
                vec![vec![0], vec![2], vec![4]]
            } else {
                vec![vec![0], vec![2]]
            },
            factors: rsdeseq2::io::read_normalization_factors_tsv(reference_data_path(
                "normalization_factors.tsv",
            ))
            .unwrap(),
            weights: rsdeseq2::io::read_observation_weights_tsv(reference_data_path(
                "observation_weights.tsv",
            ))
            .unwrap(),
        }
    }

    fn expected(&self, case: &PriorCase, normalized: bool, weighted: bool) -> DeseqResults {
        let design = ExpandedModelBetaPriorDesignInput {
            standard_design: &self.standard,
            expanded_design: &self.expanded,
            coefficient_groups: &self.groups,
        };
        let dispersions = [0.08, 0.12, 0.05, 0.2];
        let weights = weighted.then_some(&self.weights);
        if normalized {
            let input = ExpandedBetaPriorWaldNormalizedResultsInput {
                counts: &self.counts,
                design,
                normalization_factors: &self.factors,
                weights,
                dispersions: &dispersions,
                base_mean: &case.base_mean,
                disp_fit: &dispersions,
                gene_names: self.counts.gene_names(),
                options: BetaPriorRefitOptions::default(),
            };
            match &case.selection {
                Selection::Coefficient(index) => fit_expanded_beta_prior_wald_results_with_normalization_factors_and_weights(input, *index),
                Selection::Contrast(values) => fit_expanded_beta_prior_wald_contrast_results_with_normalization_factors_and_weights(input, values),
            }.unwrap().results
        } else {
            let input = ExpandedBetaPriorWaldResultsInput {
                counts: &self.counts,
                design,
                size_factors: &self.size_factors,
                weights,
                dispersions: &dispersions,
                base_mean: &case.base_mean,
                disp_fit: &dispersions,
                gene_names: self.counts.gene_names(),
                options: BetaPriorRefitOptions::default(),
            };
            match &case.selection {
                Selection::Coefficient(index) => {
                    fit_expanded_beta_prior_wald_results(input, *index)
                }
                Selection::Contrast(values) => {
                    fit_expanded_beta_prior_wald_contrast_results(input, values)
                }
            }
            .unwrap()
            .results
        }
    }
}

enum Selection {
    Coefficient(usize),
    Contrast(Vec<f64>),
}

struct PriorCase {
    name: &'static str,
    args: Vec<String>,
    selection: Selection,
    normalized: bool,
    weighted: bool,
    base_mean: [f64; 4],
}

fn check_prior_options(family: Family) {
    let dir = temp_dir("beta-prior-numerical");
    let additive = matches!(family, Family::Additive);
    let fixture = PriorFixture::new(additive);
    let standard = dir.join("standard.tsv");
    let expanded = dir.join("expanded.tsv");
    let levels = dir.join("condition.tsv");
    let batch = dir.join("batch.tsv");
    if additive {
        write_additive_contrast_design_fixture(&standard);
    } else {
        write_standard_contrast_design_fixture(&standard);
    }
    write_expanded_factor_design_fixture(&expanded);
    write_sample_level_fixture(&levels);
    write_batch_sample_level_fixture(&batch);
    let dispersions = dir.join("dispersions.tsv");
    let base_mean = dir.join("base_mean.tsv");
    let disp_fit = dir.join("disp_fit.tsv");
    write_gene_numeric_fixture(&dispersions, "dispersion", &[0.08, 0.12, 0.05, 0.2]);
    write_gene_numeric_fixture(&disp_fit, "dispFit", &[0.08, 0.12, 0.05, 0.2]);
    let output = dir.join("results.tsv");
    let mut args = basic_args("wald", &reference_data_path("counts.tsv"), &output);
    args.extend([
        "--design".into(),
        standard.to_str().unwrap().into(),
        "--fit-type".into(),
        "mean".into(),
        "--disable-cooks-cutoff".into(),
        "--beta-prior-dispersions".into(),
        dispersions.to_str().unwrap().into(),
        "--beta-prior-base-mean".into(),
        base_mean.to_str().unwrap().into(),
        "--beta-prior-disp-fit".into(),
        disp_fit.to_str().unwrap().into(),
    ]);
    args.extend(match family {
        Family::Expanded => vec![
            "--beta-prior-expanded-design".into(),
            expanded.to_str().unwrap().into(),
            "--beta-prior-coefficient-groups".into(),
            "0|2".into(),
        ],
        Family::Factor => vec![
            "--beta-prior-factor".into(),
            "condition".into(),
            "--beta-prior-reference".into(),
            "A".into(),
            "--beta-prior-sample-levels".into(),
            levels.to_str().unwrap().into(),
        ],
        Family::Additive => vec![
            "--beta-prior-additive-factors".into(),
            "condition,batch".into(),
            "--beta-prior-additive-references".into(),
            "A,X".into(),
            "--beta-prior-additive-sample-levels".into(),
            format!("{},{}", levels.display(), batch.display()),
        ],
    });
    let target = if additive {
        "batch_Y_vs_X"
    } else {
        "condition_B_vs_A"
    };
    let last = fixture.standard.n_coefficients() - 1;
    let make = |name, selector: Vec<String>, selection, normalized, weighted| PriorCase {
        name,
        args: selector,
        selection,
        normalized,
        weighted,
        base_mean: if normalized {
            [12.375, 2.5, 60.9375, 7.03125]
        } else {
            [16.5, 3.0, 97.5, 7.5]
        },
    };
    let mut contrast = vec![0.0; last + 1];
    contrast[last] = 1.0;
    let mut cases = vec![
        make(
            "coefficient",
            vec!["--coefficient-name".into(), target.into()],
            Selection::Coefficient(last),
            false,
            false,
        ),
        make(
            "named contrast",
            vec!["--contrast-name".into(), target.into()],
            Selection::Contrast(contrast.clone()),
            false,
            false,
        ),
        make(
            "normalization factors",
            vec![
                if additive {
                    "--contrast-name"
                } else {
                    "--coefficient-name"
                }
                .into(),
                target.into(),
                "--normalization-factors".into(),
                reference_data_path("normalization_factors.tsv")
                    .to_str()
                    .unwrap()
                    .into(),
            ],
            if additive {
                Selection::Contrast(contrast.clone())
            } else {
                Selection::Coefficient(last)
            },
            true,
            false,
        ),
        make(
            "nondefault coefficient",
            vec![
                "--coefficient-name".into(),
                if additive {
                    "condition_B_vs_A"
                } else {
                    "Intercept"
                }
                .into(),
            ],
            Selection::Coefficient(if additive { 1 } else { 0 }),
            false,
            false,
        ),
    ];
    if additive {
        cases.extend([
            make(
                "observation weights",
                vec![
                    "--coefficient-name".into(),
                    target.into(),
                    "--observation-weights".into(),
                    reference_data_path("observation_weights.tsv")
                        .to_str()
                        .unwrap()
                        .into(),
                ],
                Selection::Coefficient(last),
                false,
                true,
            ),
            make(
                "positive list",
                vec![
                    "--contrast-positive".into(),
                    "condition_B_vs_A,batch_Y_vs_X".into(),
                    "--contrast-negative-weight".into(),
                    "-0.5".into(),
                ],
                Selection::Contrast(vec![0.0, 1.0, 1.0]),
                false,
                false,
            ),
            make(
                "weighted positive and negative list",
                vec![
                    "--contrast-positive".into(),
                    "condition_B_vs_A".into(),
                    "--contrast-negative".into(),
                    "batch_Y_vs_X".into(),
                    "--contrast-positive-weight".into(),
                    "0.75".into(),
                    "--contrast-negative-weight".into(),
                    "-0.5".into(),
                ],
                Selection::Contrast(vec![0.0, 0.75, -0.5]),
                false,
                false,
            ),
        ]);
    }
    let baseline = fixture.expected(&cases[0], false, false);
    for case in cases {
        let label = format!("{family:?}: {}", case.name);
        write_gene_numeric_fixture(&base_mean, "baseMean", &case.base_mean);
        let expected = fixture.expected(&case, case.normalized, case.weighted);
        assert!(
            expected.rows.iter().any(|r| r.stat.is_some()),
            "{label}: no testable genes"
        );
        if case.normalized || case.weighted {
            let ignored = fixture.expected(&case, false, false);
            assert!(
                numerically_distinct(&result_values(&expected), &result_values(&ignored)),
                "{label}: ignored option undetectable"
            );
        }
        if case.name == "nondefault coefficient" || case.name.contains("list") {
            assert!(
                numerically_distinct(&result_values(&expected), &result_values(&baseline)),
                "{label}: selector undetectable"
            );
        }
        let mut selected = args.clone();
        selected.extend(case.args);
        invoke(&label, &selected);
        assert_result_export(&output, &expected, &label);
    }
}

#[test]
fn expanded_beta_prior_options_match_explicit_library_inputs() {
    check_prior_options(Family::Expanded);
}
#[test]
fn factor_beta_prior_options_match_explicit_library_inputs() {
    check_prior_options(Family::Factor);
}
#[test]
fn additive_beta_prior_options_match_explicit_library_inputs() {
    check_prior_options(Family::Additive);
}
