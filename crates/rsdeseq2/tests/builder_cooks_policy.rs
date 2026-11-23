//! Compare builder presentation with shared stored-query policy on the same fit.
//! The independent original-DESeq2 origin/row oracle lives in cooks_policy_reference.
use crate::support::run::run_results;
use rsdeseq2::errors::NameScope;
use rsdeseq2::prelude::*;

fn frame() -> FormulaModelFrame {
    FormulaModelFrame {
        factors: vec![FormulaFactorColumn {
            name: "condition".into(),
            sample_levels: ["A", "A", "A", "B", "B", "B"].map(String::from).to_vec(),
            levels: Some(vec!["A".into(), "B".into()]),
            reference: Some("A".into()),
        }],
        ..Default::default()
    }
}
fn counts() -> CountMatrix {
    CountMatrix::from_row_major_u32(
        3,
        6,
        vec![
            1, 20, 21, 20, 20, 20, 10, 11, 12, 15, 16, 17, 30, 32, 31, 50, 52, 51,
        ],
    )
    .unwrap()
}
fn builder(frame: FormulaModelFrame) -> DeseqBuilder {
    DeseqBuilder::new()
        .threads(1)
        .size_factors(vec![1.; 6])
        .fit_type(FitType::Mean)
        .gene_wise_dispersion_options(GeneWiseDispersionOptions {
            use_cox_reid: false,
            fit_method: GeneWiseDispersionFitMethod::Grid,
            niter: 2,
            ..Default::default()
        })
        .model_frame(frame)
        .cooks_cutoff_threshold(0.)
        .disable_independent_filtering()
}
fn request(contrast: Vec<f64>) -> ResultRequest {
    ResultRequest::new(
        ResultsContrast::Numeric(contrast),
        WaldTestOptions::normal(),
        CooksCutoff::Threshold(0.),
        IndependentFilteringOptions {
            enabled: false,
            ..Default::default()
        },
    )
    .unwrap()
}
#[derive(Clone)]
enum Query {
    Coefficient(usize),
    Numeric(Vec<f64>),
}
impl Query {
    fn contrast(&self, design: &DesignMatrix) -> Vec<f64> {
        match self {
            Self::Coefficient(i) => {
                let mut v = vec![0.; design.n_coefficients()];
                v[*i] = 1.;
                v
            }
            Self::Numeric(v) => v.clone(),
        }
    }
}
struct Case {
    name: &'static str,
    design: DesignMatrix,
    frame: FormulaModelFrame,
    query: Query,
    eligible: bool,
}
fn cases() -> Vec<Case> {
    let original = frame();
    let standard =
        expanded_formula_design_with_offsets_from_model_frame("~condition", &original).unwrap();
    let no_intercept =
        expanded_formula_design_with_offsets_from_model_frame("~0+condition", &original).unwrap();
    let transformed = expanded_formula_design_with_offsets_from_model_frame(
        "~relevel(condition,ref='B')",
        &original,
    )
    .unwrap();
    let mut unused = original.clone();
    let mut extra = unused.factors[0].clone();
    extra.name = "unused".into();
    unused.factors.push(extra);
    let mut numeric = original.clone();
    numeric.numeric_covariates.push(FormulaNumericColumn {
        name: "numeric".into(),
        values: vec![0., 0., 0., 1., 1., 1.],
    });
    let mut removed = original.clone();
    removed.numeric_covariates.push(FormulaNumericColumn {
        name: "nuisance".into(),
        values: vec![0., 1., 2., 3., 4., 5.],
    });
    vec![
        Case {
            name: "intercept",
            design: standard.design.standard_design.clone(),
            frame: original.clone(),
            query: Query::Coefficient(0),
            eligible: true,
        },
        Case {
            name: "scaled",
            design: standard.design.standard_design.clone(),
            frame: original.clone(),
            query: Query::Numeric(vec![0., 2.]),
            eligible: true,
        },
        Case {
            name: "mixed",
            design: standard.design.standard_design.clone(),
            frame: original.clone(),
            query: Query::Numeric(vec![1., 1.]),
            eligible: true,
        },
        Case {
            name: "unused_factor",
            design: standard.design.standard_design.clone(),
            frame: unused,
            query: Query::Coefficient(1),
            eligible: true,
        },
        Case {
            name: "transformed_factor",
            design: transformed.design.standard_design,
            frame: transformed.model_frame,
            query: Query::Numeric(vec![0., 1.]),
            eligible: true,
        },
        Case {
            name: "no_intercept",
            design: no_intercept.design.standard_design,
            frame: no_intercept.model_frame,
            query: Query::Coefficient(0),
            eligible: true,
        },
        Case {
            name: "matrix",
            design: standard
                .design
                .standard_design
                .clone()
                .with_origin(DesignOrigin::matrix()),
            frame: original,
            query: Query::Coefficient(1),
            eligible: false,
        },
        Case {
            name: "numeric_source",
            design: standard
                .design
                .standard_design
                .clone()
                .with_origin(DesignOrigin::formula(vec!["numeric".into()]).unwrap()),
            frame: numeric,
            query: Query::Coefficient(1),
            eligible: false,
        },
        Case {
            name: "removed_variable",
            design: standard.design.standard_design.with_origin(
                DesignOrigin::formula(vec!["condition".into(), "nuisance".into()]).unwrap(),
            ),
            frame: removed,
            query: Query::Coefficient(1),
            eligible: false,
        },
    ]
}
#[derive(Clone, Copy, Debug)]
enum Engine {
    FixedWald,
    FixedLrt,
    NativeWald,
    NativeLrt,
}
fn evaluate(
    engine: Engine,
    replacement: bool,
    b: &DeseqBuilder,
    counts: &CountMatrix,
    design: &DesignMatrix,
    query: &Query,
) -> (DeseqFit, DeseqResults) {
    let reduced = DesignMatrix::intercept_only(counts.n_samples()).unwrap();
    let d = [0.2; 3];
    let options = CooksReplacementOptions::new(f64::MAX);
    let column = match query {
        Query::Coefficient(index) => ResultColumn::Coefficient(*index),
        Query::Numeric(contrast) => ResultColumn::Numeric(contrast),
    };
    let request = Analysis::new(counts, design).result(column);
    let request = match engine {
        Engine::FixedWald => request.fixed_dispersion(&d),
        Engine::FixedLrt => request.lrt(&reduced).fixed_dispersion(&d),
        Engine::NativeWald => request.means(MeanModel::GlmMu),
        Engine::NativeLrt => request.lrt(&reduced).means(MeanModel::GlmMu),
    };
    let request = if replacement {
        request.cooks_replacement(&options)
    } else {
        request
    };
    match b.run(&request).unwrap() {
        RunOutput::Results { fit, results } => (fit, results),
        RunOutput::Replacement(CooksReplacementTestOutput::Wald(out)) => {
            assert_eq!(out.refit_plan.n_refit, 0);
            (out.original_fit, out.results)
        }
        RunOutput::Replacement(CooksReplacementTestOutput::Lrt(out)) => {
            assert_eq!(out.refit_plan.n_refit, 0);
            (out.original_fit, out.results)
        }
        _ => panic!("expected fitted results"),
    }
}
// Derive-based PartialEq deliberately follows IEEE NaN != NaN. Compare the
// optimizer's missing diagnostic payload bits before comparing remaining state.
fn assert_same_fit_bits(mut a: DeseqFit, mut b: DeseqFit, context: &str) {
    let bits = |v: &[f64]| v.iter().map(|v| v.to_bits()).collect::<Vec<_>>();
    assert_eq!(
        a.beta_optim_start
            .as_ref()
            .map(|m| (m.n_rows(), m.n_cols(), bits(m.as_slice()))),
        b.beta_optim_start
            .as_ref()
            .map(|m| (m.n_rows(), m.n_cols(), bits(m.as_slice()))),
        "{context}"
    );
    a.beta_optim_start = None;
    b.beta_optim_start = None;
    macro_rules! diagnostics {
        ($($field:ident),*) => {$(
            assert_eq!(a.$field.as_deref().map(bits),b.$field.as_deref().map(bits),"{} {context}",stringify!($field));
            a.$field=None;b.$field=None;
        )*}
    }
    diagnostics!(
        beta_optim_iter,
        beta_optim_start_objective,
        beta_optim_objective,
        beta_optim_gradient_norm
    );
    assert_eq!(a, b, "fitted numerical state must be unchanged: {context}");
}
fn check(case: &Case, engine: Engine, replacement: bool) {
    let counts = counts();
    let b = builder(case.frame.clone());
    let (fit, actual) = evaluate(engine, replacement, &b, &counts, &case.design, &case.query);
    let (raw_fit, _) = evaluate(
        engine,
        replacement,
        &b.disable_cooks_cutoff(),
        &counts,
        &case.design,
        &case.query,
    );
    let expected = fit
        .results_for_request(&counts, &request(case.query.contrast(&case.design)), None)
        .unwrap();
    let context = format!("{} {engine:?} replacement={replacement}", case.name);
    assert_eq!(
        actual.rows[0].cooks_outlier,
        Some(!case.eligible),
        "{context}"
    );
    assert_eq!(actual.rows, expected.rows, "{context}");
    assert_same_fit_bits(fit, raw_fit, &context);
}
#[test]
fn fixed_builders_and_no_refit_replacement_use_original_scope_for_every_display() {
    for case in cases() {
        for engine in [Engine::FixedWald, Engine::FixedLrt] {
            for replacement in [false, true] {
                check(&case, engine, replacement);
            }
        }
    }
}
#[test]
fn native_builders_and_no_refit_replacement_share_the_same_policy() {
    for case in cases().into_iter().filter(|c| {
        matches!(
            c.name,
            "intercept" | "scaled" | "unused_factor" | "transformed_factor" | "no_intercept"
        )
    }) {
        for engine in [Engine::NativeWald, Engine::NativeLrt] {
            for replacement in [false, true] {
                check(&case, engine, replacement);
            }
        }
    }
}

#[test]
fn invalid_frame_is_rejected_before_normalization_even_when_filtering_is_disabled() {
    let counts = counts();
    let original = frame();
    let design = expanded_formula_design_from_model_frame("~condition", &original)
        .unwrap()
        .standard_design;
    let mut ambiguous = original.clone();
    ambiguous.factors[0].name = "cell type".into();
    ambiguous.numeric_covariates.push(FormulaNumericColumn {
        name: "cell-type".into(),
        values: vec![0.; 6],
    });
    let mut wrong_samples = original;
    wrong_samples.factors[0].sample_levels.pop();
    // Invalid size factors would fail first if frame preflight were lost.
    let reject = |bad| {
        run_results(
            builder(bad)
                .size_factors(vec![0.; 6])
                .disable_cooks_cutoff()
                .run(
                    &Analysis::new(&counts, &design)
                        .fixed_dispersion(&[0.2; 3])
                        .result(ResultColumn::Coefficient(0)),
                ),
        )
        .unwrap_err()
    };
    crate::support::assert_ambiguous_name(
        reject(ambiguous),
        ErrorKind::Formula,
        NameScope::FormulaVariable,
        "cell-type",
        &["cell type", "cell-type"],
    );
    std::assert_matches!(
        reject(wrong_samples),
        DeseqError::InvalidDimensions { context, expected: 6, actual: 5 }
            if context == "formula model frame samples"
    );
    std::assert_matches!(
        reject(FormulaModelFrame::default()),
        DeseqError::MissingPrerequisite { domain: ErrorKind::Formula, operation, prerequisite }
            if operation == "formula model frame" && prerequisite == "at least one column"
    );
}

#[test]
fn explicit_factor_alias_without_a_frame_preserves_original_formula_scope() {
    let counts = counts();
    let design = DesignMatrix::from_row_major(
        6,
        2,
        vec![1., 0., 1., 0., 1., 0., 1., 1., 1., 1., 1., 1.],
        Some(vec!["Intercept".into(), "cell.type_B_vs_A".into()]),
    )
    .unwrap()
    .with_origin(DesignOrigin::formula(vec!["cell type".into()]).unwrap());
    let levels = frame().factors.remove(0).sample_levels;
    for factor in ["cell type", "cell.type"] {
        let (_, result) = run_results(
            DeseqBuilder::new()
                .size_factors(vec![1.; 6])
                .cooks_cutoff_threshold(0.)
                .disable_independent_filtering()
                .run(
                    &Analysis::new(&counts, &design)
                        .fixed_dispersion(&[0.2; 3])
                        .result(ResultColumn::FactorLevel(FactorLevelContrast::new(
                            factor, "B", "A", &levels,
                        ))),
                ),
        )
        .unwrap();
        assert_eq!(result.rows[0].cooks_outlier, Some(false), "{factor}");
        assert!(result.rows[0].pvalue.is_some());
    }
}
