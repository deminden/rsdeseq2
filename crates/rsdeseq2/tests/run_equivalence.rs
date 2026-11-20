//! Equivalence of every public `DeseqBuilder::fit*` method with its
//! `DeseqBuilder::run` request (design note 4.2, step 1).
//!
//! Each method and its `Analysis` run on the same builder and inputs; the
//! outputs (or error messages) must be identical to the bit, compared through
//! their `Debug` form so NaN fields compare equal.

use rsdeseq2::prelude::*;

struct Fixture {
    counts: CountMatrix,
    design: DesignMatrix,
    reduced: DesignMatrix,
    frame: FormulaModelFrame,
    levels: Vec<String>,
    dispersions: Vec<f64>,
    numeric: Vec<f64>,
    spec: ContrastSpec,
    cooks: CooksReplacementOptions,
    full_formula: &'static str,
    reduced_formula: &'static str,
}

impl Fixture {
    fn new() -> Self {
        let n = 16;
        let levels = (0..n)
            .map(|j| if j < 8 { "A" } else { "B" }.to_string())
            .collect::<Vec<_>>();
        let batch = (0..n)
            .map(|j| if j % 2 == 0 { "x" } else { "y" }.to_string())
            .collect::<Vec<_>>();
        let frame = FormulaModelFrame {
            factors: vec![
                FormulaFactorColumn {
                    name: "batch".to_string(),
                    sample_levels: batch,
                    levels: None,
                    reference: None,
                },
                FormulaFactorColumn {
                    name: "condition".to_string(),
                    sample_levels: levels.clone(),
                    levels: None,
                    reference: None,
                },
            ],
            numeric_covariates: Vec::new(),
        };
        let builder = DeseqBuilder::new().model_frame(frame.clone());
        let design = builder
            .expanded_formula_design("~ batch + condition")
            .unwrap()
            .standard_design;
        let reduced = builder
            .expanded_formula_design("~ batch")
            .unwrap()
            .standard_design;
        // 14 genes: deterministic spread of means, a condition effect on some
        // rows, one all-zero row and one Cook's outlier.
        let genes = 14;
        let mut values = Vec::with_capacity(genes * n);
        for g in 0..genes {
            for j in 0..n {
                let base = 5.0 * (1.6f64).powi(g as i32);
                let effect = if j >= 8 && g % 3 == 0 { 2.5 } else { 1.0 };
                let wobble = 1.0 + 0.15 * (((g * 7 + j * 3) % 5) as f64 - 2.0);
                let mut value = (base * effect * wobble).round() as u32;
                if g == 0 {
                    value = 0;
                }
                if g == 5 && j == 3 {
                    value *= 40;
                }
                values.push(value);
            }
        }
        let counts = CountMatrix::from_row_major_u32(genes, n, values).unwrap();
        let cutoff = default_cooks_cutoff(n, design.n_coefficients())
            .unwrap()
            .unwrap();
        Self {
            dispersions: (0..genes).map(|g| 0.05 + 0.01 * g as f64).collect(),
            counts,
            design,
            reduced,
            frame,
            levels,
            numeric: vec![0.0, 0.0, 1.0],
            spec: ContrastSpec::factor_level("condition", "B", "A"),
            cooks: CooksReplacementOptions::new(cutoff),
            full_formula: "~ batch + condition",
            reduced_formula: "~ batch",
        }
    }

    fn builder(&self, test: &str) -> DeseqBuilder {
        if test == "bare" {
            // No stored model frame: methods that take one must record it.
            return DeseqBuilder::new();
        }
        let builder = DeseqBuilder::new().model_frame(self.frame.clone());
        match test {
            "lrt" => builder
                .test(TestType::Lrt)
                .reduced_design(self.reduced.clone()),
            _ => builder,
        }
    }

    fn factor_level(&self) -> FactorLevelContrast<'_> {
        FactorLevelContrast::new("condition", "B", "A", &self.levels)
    }

    fn character(&self) -> ResultsContrast {
        ResultsContrast::character("condition", "B", "A")
    }

    fn numeric_request(&self) -> ResultsContrast {
        ResultsContrast::numeric(vec![0.0, 0.0, 1.0])
    }
}

thread_local! {
    static ERRORS: std::cell::RefCell<Vec<String>> = const { std::cell::RefCell::new(Vec::new()) };
}

fn shown<T: std::fmt::Debug>(value: Result<T, DeseqError>) -> String {
    match value {
        Ok(value) => format!("{value:?}"),
        Err(error) => format!("error: {error}"),
    }
}

fn fr(value: Result<(DeseqFit, DeseqResults), DeseqError>) -> String {
    shown(value)
}

fn fo(value: Result<DeseqFit, DeseqError>) -> String {
    shown(value)
}

fn rw(value: Result<CooksReplacementWaldOutput, DeseqError>) -> String {
    shown(value)
}

fn rl(value: Result<CooksReplacementLrtOutput, DeseqError>) -> String {
    shown(value)
}

fn rt(value: Result<CooksReplacementTestOutput, DeseqError>) -> String {
    shown(value)
}

fn run_fr(value: Result<RunOutput, DeseqError>) -> String {
    shown(value.and_then(|output| match output {
        RunOutput::Results { fit, results } => Ok((fit, results)),
        other => Err(DeseqError::InvalidOptions {
            reason: format!("unexpected run output shape: {other:?}"),
        }),
    }))
}

fn run_fo(value: Result<RunOutput, DeseqError>) -> String {
    shown(value.map(|output| output.fit().clone()))
}

fn run_rw(value: Result<RunOutput, DeseqError>) -> String {
    shown(value.and_then(|output| match output {
        RunOutput::Replacement(output) => match output {
            CooksReplacementTestOutput::Wald(output) => Ok(output),
            other => Err(DeseqError::InvalidOptions {
                reason: format!("unexpected test: {:?}", other.test_type()),
            }),
        },
        other => Err(DeseqError::InvalidOptions {
            reason: format!("unexpected run output shape: {}", short(&other)),
        }),
    }))
}

fn run_rl(value: Result<RunOutput, DeseqError>) -> String {
    shown(value.and_then(|output| match output {
        RunOutput::Replacement(output) => match output {
            CooksReplacementTestOutput::Lrt(output) => Ok(output),
            other => Err(DeseqError::InvalidOptions {
                reason: format!("unexpected test: {:?}", other.test_type()),
            }),
        },
        other => Err(DeseqError::InvalidOptions {
            reason: format!("unexpected run output shape: {}", short(&other)),
        }),
    }))
}

fn run_rt(value: Result<RunOutput, DeseqError>) -> String {
    shown(value.and_then(|output| match output {
        RunOutput::Replacement(output) => Ok(output),
        other => Err(DeseqError::InvalidOptions {
            reason: format!("unexpected run output shape: {}", short(&other)),
        }),
    }))
}

fn short(output: &RunOutput) -> &'static str {
    match output {
        RunOutput::Fit(_) => "Fit",
        RunOutput::Results { .. } => "Results",
        RunOutput::Replacement(_) => "Replacement",
        _ => "other",
    }
}

fn check(report: &mut Vec<String>, tag: &str, method: String, run: String) {
    // RSDESEQ2_RUN_EQUIVALENCE_DUMP=<file> records every method output, so
    // two commits can be compared method by method.
    if let Some(path) = std::env::var_os("RSDESEQ2_RUN_EQUIVALENCE_DUMP") {
        use std::io::Write;
        let mut file = std::fs::OpenOptions::new()
            .create(true)
            .append(true)
            .open(path)
            .unwrap();
        writeln!(file, "{tag}\t{}", method.replace(['\n', '\t'], " ")).unwrap();
    }
    if method.starts_with("error: ") {
        ERRORS.with(|errors| errors.borrow_mut().push(format!("{tag}: {method}")));
    }
    if method != run {
        let first = |s: &str| s.chars().take(160).collect::<String>();
        report.push(format!(
            "{tag}\n    method: {}\n    run:    {}",
            first(&method),
            first(&run)
        ));
    }
}

#[test]
fn every_fit_method_matches_its_run_request() {
    // One large generated function: give it a large stack in debug builds.
    let report = std::thread::Builder::new()
        .stack_size(256 << 20)
        .spawn(|| {
            let report = compare_all(&Fixture::new());
            let errors = ERRORS.with(|errors| errors.take());
            // Every pair must exercise a successful run, not a shared error.
            for error in &errors {
                eprintln!("both error: {error}");
            }
            assert!(
                errors.is_empty(),
                "{} pairs failed on both sides",
                errors.len()
            );
            report
        })
        .unwrap()
        .join()
        .unwrap();
    assert!(
        report.is_empty(),
        "{} of the method/run pairs differ:\n{}",
        report.len(),
        report.join("\n")
    );
}

#[test]
fn run_rejects_invalid_combinations_before_work() {
    use rsdeseq2::errors::RunConflict;
    let f = Fixture::new();
    let b = f.builder("wald");
    let reject = |analysis: Analysis<'_>| b.run(&analysis).unwrap_err();
    let base = || Analysis::new(&f.counts, &f.design);
    let checks = [
        (
            base()
                .fixed_dispersion(&f.dispersions)
                .stop_after(Stage::Map),
            RunConflict::FixedDispersionEarlyStop,
        ),
        (
            base()
                .fixed_dispersion(&f.dispersions)
                .means(MeanModel::LinearMu),
            RunConflict::FixedDispersionMeanModel,
        ),
        (
            base().stop_after(Stage::Trend).cooks_replacement(&f.cooks),
            RunConflict::ReplacementEarlyStop,
        ),
        (
            base()
                .stop_after(Stage::Trend)
                .result(ResultColumn::Coefficient(1)),
            RunConflict::ResultColumnEarlyStop,
        ),
        (
            base().lrt(&f.reduced).stop_after(Stage::Map),
            RunConflict::LrtEarlyStop,
        ),
        (
            base().lrt_formula(f.reduced_formula),
            RunConflict::ReducedFormulaNeedsFullFormula,
        ),
        (
            Analysis::without_design(&f.counts).stop_after(Stage::Results),
            RunConflict::DesignRequired,
        ),
    ];
    for (analysis, expected) in checks {
        std::assert_matches!(
            reject(analysis), DeseqError::InvalidRun { conflict } if conflict == expected
        );
    }
    std::assert_matches!(
        reject(base().fixed_dispersion(&f.dispersions[1..])),
        DeseqError::InvalidDimensions { context, expected, actual }
            if context == "run fixed dispersions" && expected == f.counts.n_genes()
                && actual == f.dispersions.len() - 1
    );
    for analysis in [base(), Analysis::formula(&f.counts, "~ invalid_metadata")] {
        std::assert_matches!(
            reject(analysis.means(MeanModel::LinearMu).cooks_replacement(&f.cooks)),
            DeseqError::UnsupportedFeature { feature }
                if feature == "run: Cook's replacement with linear-mu dispersions"
        );
    }
    let character = f.character();
    let both = ResultColumn::Results {
        contrast: &character,
        sample_levels: Some(&f.levels),
        model_frame: Some(&f.frame),
    };
    for analysis in [
        base().result(both),
        Analysis::formula(&f.counts, "~ invalid_metadata").result(both),
    ] {
        std::assert_matches!(
            reject(analysis),
            DeseqError::InvalidRun {
                conflict: RunConflict::DuplicateContrastMetadata
            }
        );
    }
}

#[test]
fn run_trend_override_matches_a_builder_with_that_trend() {
    // The dedicated leaves and the builder-copy path must agree with a
    // builder whose own fit type is the override.
    let f = Fixture::new();
    let b = f.builder("wald");
    for fit_type in [FitType::Mean, FitType::Local] {
        let typed = b.clone().fit_type(fit_type);
        for stage in [Stage::Trend, Stage::Map, Stage::Results] {
            let analysis = Analysis::new(&f.counts, &f.design).stop_after(stage);
            assert_eq!(
                shown(b.run(&analysis.fit_type(fit_type))),
                shown(typed.run(&analysis)),
                "{fit_type:?} {stage:?}"
            );
        }
    }
}
// @generated below by tools/validation/gen_run_equivalence.py from the builder
// signatures: one block per method and builder/request variant.
// 224 comparisons.
fn compare_all(f: &Fixture) -> Vec<String> {
    let mut report = Vec::new();
    {
        let b = f.builder("wald");
        let rc = f.character();
        let method = fr(b.fit_fixed_dispersion_lrt_results_contrast(
            &f.counts,
            &f.design,
            &f.reduced,
            &f.dispersions,
            &rc,
            Some(f.levels.as_slice()),
        ));
        let run = run_fr(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .fixed_dispersion(&f.dispersions)
                    .result(ResultColumn::Results {
                        contrast: &rc,
                        sample_levels: Some(f.levels.as_slice()),
                        model_frame: None,
                    }),
            ),
        );
        check(
            &mut report,
            "fit_fixed_dispersion_lrt_results_contrast/wald/character/levels",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");
        let rc = f.character();
        let method = fr(b.fit_fixed_dispersion_lrt_results_contrast(
            &f.counts,
            &f.design,
            &f.reduced,
            &f.dispersions,
            &rc,
            None::<&[String]>,
        ));
        let run = run_fr(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .fixed_dispersion(&f.dispersions)
                    .result(ResultColumn::Results {
                        contrast: &rc,
                        sample_levels: None,
                        model_frame: None,
                    }),
            ),
        );
        check(
            &mut report,
            "fit_fixed_dispersion_lrt_results_contrast/wald/character/nolevels",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");
        let rc = f.numeric_request();
        let method = fr(b.fit_fixed_dispersion_lrt_results_contrast(
            &f.counts,
            &f.design,
            &f.reduced,
            &f.dispersions,
            &rc,
            Some(f.levels.as_slice()),
        ));
        let run = run_fr(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .fixed_dispersion(&f.dispersions)
                    .result(ResultColumn::Results {
                        contrast: &rc,
                        sample_levels: Some(f.levels.as_slice()),
                        model_frame: None,
                    }),
            ),
        );
        check(
            &mut report,
            "fit_fixed_dispersion_lrt_results_contrast/wald/numeric/levels",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");
        let rc = f.numeric_request();
        let method = fr(b.fit_fixed_dispersion_lrt_results_contrast(
            &f.counts,
            &f.design,
            &f.reduced,
            &f.dispersions,
            &rc,
            None::<&[String]>,
        ));
        let run = run_fr(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .fixed_dispersion(&f.dispersions)
                    .result(ResultColumn::Results {
                        contrast: &rc,
                        sample_levels: None,
                        model_frame: None,
                    }),
            ),
        );
        check(
            &mut report,
            "fit_fixed_dispersion_lrt_results_contrast/wald/numeric/nolevels",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");
        let rc = f.character();
        let method = fr(
            b.fit_fixed_dispersion_lrt_results_contrast_from_model_frame(
                &f.counts,
                &f.design,
                &f.reduced,
                &f.dispersions,
                &rc,
                &f.frame,
            ),
        );
        let run = run_fr(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .fixed_dispersion(&f.dispersions)
                    .result(ResultColumn::Results {
                        contrast: &rc,
                        sample_levels: None,
                        model_frame: Some(&f.frame),
                    }),
            ),
        );
        check(
            &mut report,
            "fit_fixed_dispersion_lrt_results_contrast_from_model_frame/wald/character",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");
        let rc = f.numeric_request();
        let method = fr(
            b.fit_fixed_dispersion_lrt_results_contrast_from_model_frame(
                &f.counts,
                &f.design,
                &f.reduced,
                &f.dispersions,
                &rc,
                &f.frame,
            ),
        );
        let run = run_fr(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .fixed_dispersion(&f.dispersions)
                    .result(ResultColumn::Results {
                        contrast: &rc,
                        sample_levels: None,
                        model_frame: Some(&f.frame),
                    }),
            ),
        );
        check(
            &mut report,
            "fit_fixed_dispersion_lrt_results_contrast_from_model_frame/wald/numeric",
            method,
            run,
        );
    }
    {
        let b = f.builder("bare");
        let rc = f.character();
        let method = fr(
            b.fit_fixed_dispersion_lrt_results_contrast_from_model_frame(
                &f.counts,
                &f.design,
                &f.reduced,
                &f.dispersions,
                &rc,
                &f.frame,
            ),
        );
        let run = run_fr(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .fixed_dispersion(&f.dispersions)
                    .result(ResultColumn::Results {
                        contrast: &rc,
                        sample_levels: None,
                        model_frame: Some(&f.frame),
                    }),
            ),
        );
        check(
            &mut report,
            "fit_fixed_dispersion_lrt_results_contrast_from_model_frame/bare/character",
            method,
            run,
        );
    }
    {
        let b = f.builder("bare");
        let rc = f.numeric_request();
        let method = fr(
            b.fit_fixed_dispersion_lrt_results_contrast_from_model_frame(
                &f.counts,
                &f.design,
                &f.reduced,
                &f.dispersions,
                &rc,
                &f.frame,
            ),
        );
        let run = run_fr(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .fixed_dispersion(&f.dispersions)
                    .result(ResultColumn::Results {
                        contrast: &rc,
                        sample_levels: None,
                        model_frame: Some(&f.frame),
                    }),
            ),
        );
        check(
            &mut report,
            "fit_fixed_dispersion_lrt_results_contrast_from_model_frame/bare/numeric",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");
        let rc = f.character();
        let method = rl(
            b.fit_fixed_dispersion_lrt_results_contrast_with_cooks_replacement(
                &f.counts,
                &f.design,
                &f.reduced,
                &f.dispersions,
                &rc,
                Some(f.levels.as_slice()),
                &f.cooks,
            ),
        );
        let run = run_rl(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .fixed_dispersion(&f.dispersions)
                    .result(ResultColumn::Results {
                        contrast: &rc,
                        sample_levels: Some(f.levels.as_slice()),
                        model_frame: None,
                    })
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_fixed_dispersion_lrt_results_contrast_with_cooks_replacement/wald/character/levels",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");
        let rc = f.character();
        let method = rl(
            b.fit_fixed_dispersion_lrt_results_contrast_with_cooks_replacement(
                &f.counts,
                &f.design,
                &f.reduced,
                &f.dispersions,
                &rc,
                None::<&[String]>,
                &f.cooks,
            ),
        );
        let run = run_rl(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .fixed_dispersion(&f.dispersions)
                    .result(ResultColumn::Results {
                        contrast: &rc,
                        sample_levels: None,
                        model_frame: None,
                    })
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_fixed_dispersion_lrt_results_contrast_with_cooks_replacement/wald/character/nolevels",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");
        let rc = f.numeric_request();
        let method = rl(
            b.fit_fixed_dispersion_lrt_results_contrast_with_cooks_replacement(
                &f.counts,
                &f.design,
                &f.reduced,
                &f.dispersions,
                &rc,
                Some(f.levels.as_slice()),
                &f.cooks,
            ),
        );
        let run = run_rl(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .fixed_dispersion(&f.dispersions)
                    .result(ResultColumn::Results {
                        contrast: &rc,
                        sample_levels: Some(f.levels.as_slice()),
                        model_frame: None,
                    })
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_fixed_dispersion_lrt_results_contrast_with_cooks_replacement/wald/numeric/levels",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");
        let rc = f.numeric_request();
        let method = rl(
            b.fit_fixed_dispersion_lrt_results_contrast_with_cooks_replacement(
                &f.counts,
                &f.design,
                &f.reduced,
                &f.dispersions,
                &rc,
                None::<&[String]>,
                &f.cooks,
            ),
        );
        let run = run_rl(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .fixed_dispersion(&f.dispersions)
                    .result(ResultColumn::Results {
                        contrast: &rc,
                        sample_levels: None,
                        model_frame: None,
                    })
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_fixed_dispersion_lrt_results_contrast_with_cooks_replacement/wald/numeric/nolevels",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");
        let rc = f.character();
        let method = rl(b
            .fit_fixed_dispersion_lrt_results_contrast_from_model_frame_with_cooks_replacement(
                &f.counts,
                &f.design,
                &f.reduced,
                &f.dispersions,
                &rc,
                &f.frame,
                &f.cooks,
            ));
        let run = run_rl(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .fixed_dispersion(&f.dispersions)
                    .result(ResultColumn::Results {
                        contrast: &rc,
                        sample_levels: None,
                        model_frame: Some(&f.frame),
                    })
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_fixed_dispersion_lrt_results_contrast_from_model_frame_with_cooks_replacement/wald/character",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");
        let rc = f.numeric_request();
        let method = rl(b
            .fit_fixed_dispersion_lrt_results_contrast_from_model_frame_with_cooks_replacement(
                &f.counts,
                &f.design,
                &f.reduced,
                &f.dispersions,
                &rc,
                &f.frame,
                &f.cooks,
            ));
        let run = run_rl(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .fixed_dispersion(&f.dispersions)
                    .result(ResultColumn::Results {
                        contrast: &rc,
                        sample_levels: None,
                        model_frame: Some(&f.frame),
                    })
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_fixed_dispersion_lrt_results_contrast_from_model_frame_with_cooks_replacement/wald/numeric",
            method,
            run,
        );
    }
    {
        let b = f.builder("bare");
        let rc = f.character();
        let method = rl(b
            .fit_fixed_dispersion_lrt_results_contrast_from_model_frame_with_cooks_replacement(
                &f.counts,
                &f.design,
                &f.reduced,
                &f.dispersions,
                &rc,
                &f.frame,
                &f.cooks,
            ));
        let run = run_rl(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .fixed_dispersion(&f.dispersions)
                    .result(ResultColumn::Results {
                        contrast: &rc,
                        sample_levels: None,
                        model_frame: Some(&f.frame),
                    })
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_fixed_dispersion_lrt_results_contrast_from_model_frame_with_cooks_replacement/bare/character",
            method,
            run,
        );
    }
    {
        let b = f.builder("bare");
        let rc = f.numeric_request();
        let method = rl(b
            .fit_fixed_dispersion_lrt_results_contrast_from_model_frame_with_cooks_replacement(
                &f.counts,
                &f.design,
                &f.reduced,
                &f.dispersions,
                &rc,
                &f.frame,
                &f.cooks,
            ));
        let run = run_rl(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .fixed_dispersion(&f.dispersions)
                    .result(ResultColumn::Results {
                        contrast: &rc,
                        sample_levels: None,
                        model_frame: Some(&f.frame),
                    })
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_fixed_dispersion_lrt_results_contrast_from_model_frame_with_cooks_replacement/bare/numeric",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");
        let rc = f.character();
        let method = fr(b.fit_fixed_dispersion_wald_results_contrast(
            &f.counts,
            &f.design,
            &f.dispersions,
            &rc,
            Some(f.levels.as_slice()),
        ));
        let run = run_fr(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .fixed_dispersion(&f.dispersions)
                    .result(ResultColumn::Results {
                        contrast: &rc,
                        sample_levels: Some(f.levels.as_slice()),
                        model_frame: None,
                    }),
            ),
        );
        check(
            &mut report,
            "fit_fixed_dispersion_wald_results_contrast/wald/character/levels",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");
        let rc = f.character();
        let method = fr(b.fit_fixed_dispersion_wald_results_contrast(
            &f.counts,
            &f.design,
            &f.dispersions,
            &rc,
            None::<&[String]>,
        ));
        let run = run_fr(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .fixed_dispersion(&f.dispersions)
                    .result(ResultColumn::Results {
                        contrast: &rc,
                        sample_levels: None,
                        model_frame: None,
                    }),
            ),
        );
        check(
            &mut report,
            "fit_fixed_dispersion_wald_results_contrast/wald/character/nolevels",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");
        let rc = f.numeric_request();
        let method = fr(b.fit_fixed_dispersion_wald_results_contrast(
            &f.counts,
            &f.design,
            &f.dispersions,
            &rc,
            Some(f.levels.as_slice()),
        ));
        let run = run_fr(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .fixed_dispersion(&f.dispersions)
                    .result(ResultColumn::Results {
                        contrast: &rc,
                        sample_levels: Some(f.levels.as_slice()),
                        model_frame: None,
                    }),
            ),
        );
        check(
            &mut report,
            "fit_fixed_dispersion_wald_results_contrast/wald/numeric/levels",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");
        let rc = f.numeric_request();
        let method = fr(b.fit_fixed_dispersion_wald_results_contrast(
            &f.counts,
            &f.design,
            &f.dispersions,
            &rc,
            None::<&[String]>,
        ));
        let run = run_fr(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .fixed_dispersion(&f.dispersions)
                    .result(ResultColumn::Results {
                        contrast: &rc,
                        sample_levels: None,
                        model_frame: None,
                    }),
            ),
        );
        check(
            &mut report,
            "fit_fixed_dispersion_wald_results_contrast/wald/numeric/nolevels",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");
        let rc = f.character();
        let method = fr(
            b.fit_fixed_dispersion_wald_results_contrast_from_model_frame(
                &f.counts,
                &f.design,
                &f.dispersions,
                &rc,
                &f.frame,
            ),
        );
        let run = run_fr(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .fixed_dispersion(&f.dispersions)
                    .result(ResultColumn::Results {
                        contrast: &rc,
                        sample_levels: None,
                        model_frame: Some(&f.frame),
                    }),
            ),
        );
        check(
            &mut report,
            "fit_fixed_dispersion_wald_results_contrast_from_model_frame/wald/character",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");
        let rc = f.numeric_request();
        let method = fr(
            b.fit_fixed_dispersion_wald_results_contrast_from_model_frame(
                &f.counts,
                &f.design,
                &f.dispersions,
                &rc,
                &f.frame,
            ),
        );
        let run = run_fr(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .fixed_dispersion(&f.dispersions)
                    .result(ResultColumn::Results {
                        contrast: &rc,
                        sample_levels: None,
                        model_frame: Some(&f.frame),
                    }),
            ),
        );
        check(
            &mut report,
            "fit_fixed_dispersion_wald_results_contrast_from_model_frame/wald/numeric",
            method,
            run,
        );
    }
    {
        let b = f.builder("bare");
        let rc = f.character();
        let method = fr(
            b.fit_fixed_dispersion_wald_results_contrast_from_model_frame(
                &f.counts,
                &f.design,
                &f.dispersions,
                &rc,
                &f.frame,
            ),
        );
        let run = run_fr(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .fixed_dispersion(&f.dispersions)
                    .result(ResultColumn::Results {
                        contrast: &rc,
                        sample_levels: None,
                        model_frame: Some(&f.frame),
                    }),
            ),
        );
        check(
            &mut report,
            "fit_fixed_dispersion_wald_results_contrast_from_model_frame/bare/character",
            method,
            run,
        );
    }
    {
        let b = f.builder("bare");
        let rc = f.numeric_request();
        let method = fr(
            b.fit_fixed_dispersion_wald_results_contrast_from_model_frame(
                &f.counts,
                &f.design,
                &f.dispersions,
                &rc,
                &f.frame,
            ),
        );
        let run = run_fr(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .fixed_dispersion(&f.dispersions)
                    .result(ResultColumn::Results {
                        contrast: &rc,
                        sample_levels: None,
                        model_frame: Some(&f.frame),
                    }),
            ),
        );
        check(
            &mut report,
            "fit_fixed_dispersion_wald_results_contrast_from_model_frame/bare/numeric",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");
        let rc = f.character();
        let method = rw(
            b.fit_fixed_dispersion_wald_results_contrast_with_cooks_replacement(
                &f.counts,
                &f.design,
                &f.dispersions,
                &rc,
                Some(f.levels.as_slice()),
                &f.cooks,
            ),
        );
        let run = run_rw(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .fixed_dispersion(&f.dispersions)
                    .result(ResultColumn::Results {
                        contrast: &rc,
                        sample_levels: Some(f.levels.as_slice()),
                        model_frame: None,
                    })
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_fixed_dispersion_wald_results_contrast_with_cooks_replacement/wald/character/levels",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");
        let rc = f.character();
        let method = rw(
            b.fit_fixed_dispersion_wald_results_contrast_with_cooks_replacement(
                &f.counts,
                &f.design,
                &f.dispersions,
                &rc,
                None::<&[String]>,
                &f.cooks,
            ),
        );
        let run = run_rw(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .fixed_dispersion(&f.dispersions)
                    .result(ResultColumn::Results {
                        contrast: &rc,
                        sample_levels: None,
                        model_frame: None,
                    })
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_fixed_dispersion_wald_results_contrast_with_cooks_replacement/wald/character/nolevels",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");
        let rc = f.numeric_request();
        let method = rw(
            b.fit_fixed_dispersion_wald_results_contrast_with_cooks_replacement(
                &f.counts,
                &f.design,
                &f.dispersions,
                &rc,
                Some(f.levels.as_slice()),
                &f.cooks,
            ),
        );
        let run = run_rw(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .fixed_dispersion(&f.dispersions)
                    .result(ResultColumn::Results {
                        contrast: &rc,
                        sample_levels: Some(f.levels.as_slice()),
                        model_frame: None,
                    })
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_fixed_dispersion_wald_results_contrast_with_cooks_replacement/wald/numeric/levels",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");
        let rc = f.numeric_request();
        let method = rw(
            b.fit_fixed_dispersion_wald_results_contrast_with_cooks_replacement(
                &f.counts,
                &f.design,
                &f.dispersions,
                &rc,
                None::<&[String]>,
                &f.cooks,
            ),
        );
        let run = run_rw(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .fixed_dispersion(&f.dispersions)
                    .result(ResultColumn::Results {
                        contrast: &rc,
                        sample_levels: None,
                        model_frame: None,
                    })
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_fixed_dispersion_wald_results_contrast_with_cooks_replacement/wald/numeric/nolevels",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");
        let rc = f.character();
        let method = rw(b
            .fit_fixed_dispersion_wald_results_contrast_from_model_frame_with_cooks_replacement(
                &f.counts,
                &f.design,
                &f.dispersions,
                &rc,
                &f.frame,
                &f.cooks,
            ));
        let run = run_rw(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .fixed_dispersion(&f.dispersions)
                    .result(ResultColumn::Results {
                        contrast: &rc,
                        sample_levels: None,
                        model_frame: Some(&f.frame),
                    })
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_fixed_dispersion_wald_results_contrast_from_model_frame_with_cooks_replacement/wald/character",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");
        let rc = f.numeric_request();
        let method = rw(b
            .fit_fixed_dispersion_wald_results_contrast_from_model_frame_with_cooks_replacement(
                &f.counts,
                &f.design,
                &f.dispersions,
                &rc,
                &f.frame,
                &f.cooks,
            ));
        let run = run_rw(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .fixed_dispersion(&f.dispersions)
                    .result(ResultColumn::Results {
                        contrast: &rc,
                        sample_levels: None,
                        model_frame: Some(&f.frame),
                    })
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_fixed_dispersion_wald_results_contrast_from_model_frame_with_cooks_replacement/wald/numeric",
            method,
            run,
        );
    }
    {
        let b = f.builder("bare");
        let rc = f.character();
        let method = rw(b
            .fit_fixed_dispersion_wald_results_contrast_from_model_frame_with_cooks_replacement(
                &f.counts,
                &f.design,
                &f.dispersions,
                &rc,
                &f.frame,
                &f.cooks,
            ));
        let run = run_rw(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .fixed_dispersion(&f.dispersions)
                    .result(ResultColumn::Results {
                        contrast: &rc,
                        sample_levels: None,
                        model_frame: Some(&f.frame),
                    })
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_fixed_dispersion_wald_results_contrast_from_model_frame_with_cooks_replacement/bare/character",
            method,
            run,
        );
    }
    {
        let b = f.builder("bare");
        let rc = f.numeric_request();
        let method = rw(b
            .fit_fixed_dispersion_wald_results_contrast_from_model_frame_with_cooks_replacement(
                &f.counts,
                &f.design,
                &f.dispersions,
                &rc,
                &f.frame,
                &f.cooks,
            ));
        let run = run_rw(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .fixed_dispersion(&f.dispersions)
                    .result(ResultColumn::Results {
                        contrast: &rc,
                        sample_levels: None,
                        model_frame: Some(&f.frame),
                    })
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_fixed_dispersion_wald_results_contrast_from_model_frame_with_cooks_replacement/bare/numeric",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method =
            rw(b.fit_wald_glm_mu_with_cooks_replacement(&f.counts, &f.design, 2, &f.cooks));
        let run = run_rw(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .means(MeanModel::GlmMu)
                    .result(ResultColumn::Coefficient(2))
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_wald_glm_mu_with_cooks_replacement/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method = rw(b.fit_wald_glm_mu_contrast_with_cooks_replacement(
            &f.counts, &f.design, &f.numeric, &f.cooks,
        ));
        let run = run_rw(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .means(MeanModel::GlmMu)
                    .result(ResultColumn::Numeric(&f.numeric))
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_wald_glm_mu_contrast_with_cooks_replacement/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method = rw(b.fit_wald_glm_mu_contrast_spec_with_cooks_replacement(
            &f.counts, &f.design, &f.spec, &f.cooks,
        ));
        let run = run_rw(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .means(MeanModel::GlmMu)
                    .result(ResultColumn::Spec(&f.spec))
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_wald_glm_mu_contrast_spec_with_cooks_replacement/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method = rw(
            b.fit_wald_glm_mu_factor_level_contrast_with_cooks_replacement(
                &f.counts,
                &f.design,
                f.factor_level(),
                &f.cooks,
            ),
        );
        let run = run_rw(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .means(MeanModel::GlmMu)
                    .result(ResultColumn::FactorLevel(f.factor_level()))
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_wald_glm_mu_factor_level_contrast_with_cooks_replacement/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method = rl(
            b.fit_lrt_glm_mu_with_cooks_replacement(&f.counts, &f.design, &f.reduced, 2, &f.cooks)
        );
        let run = run_rl(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .means(MeanModel::GlmMu)
                    .result(ResultColumn::Coefficient(2))
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_lrt_glm_mu_with_cooks_replacement/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method = rl(b.fit_lrt_glm_mu_contrast_with_cooks_replacement(
            &f.counts, &f.design, &f.reduced, &f.numeric, &f.cooks,
        ));
        let run = run_rl(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .means(MeanModel::GlmMu)
                    .result(ResultColumn::Numeric(&f.numeric))
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_lrt_glm_mu_contrast_with_cooks_replacement/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method = rl(
            b.fit_lrt_glm_mu_factor_level_contrast_with_cooks_replacement(
                &f.counts,
                &f.design,
                &f.reduced,
                f.factor_level(),
                &f.cooks,
            ),
        );
        let run = run_rl(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .means(MeanModel::GlmMu)
                    .result(ResultColumn::FactorLevel(f.factor_level()))
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_lrt_glm_mu_factor_level_contrast_with_cooks_replacement/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method = rl(b.fit_lrt_glm_mu_contrast_spec_with_cooks_replacement(
            &f.counts, &f.design, &f.reduced, &f.spec, &f.cooks,
        ));
        let run = run_rl(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .means(MeanModel::GlmMu)
                    .result(ResultColumn::Spec(&f.spec))
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_lrt_glm_mu_contrast_spec_with_cooks_replacement/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method =
            fr(b.fit_fixed_dispersion_lrt(&f.counts, &f.design, &f.reduced, &f.dispersions, 2));
        let run = run_fr(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .fixed_dispersion(&f.dispersions)
                    .result(ResultColumn::Coefficient(2)),
            ),
        );
        check(&mut report, "fit_fixed_dispersion_lrt/wald", method, run);
    }
    {
        let b = f.builder("wald");

        let method = rl(b.fit_fixed_dispersion_lrt_with_cooks_replacement(
            &f.counts,
            &f.design,
            &f.reduced,
            &f.dispersions,
            2,
            &f.cooks,
        ));
        let run = run_rl(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .fixed_dispersion(&f.dispersions)
                    .result(ResultColumn::Coefficient(2))
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_fixed_dispersion_lrt_with_cooks_replacement/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method = fr(b.fit_fixed_dispersion_lrt_contrast(
            &f.counts,
            &f.design,
            &f.reduced,
            &f.dispersions,
            &f.numeric,
        ));
        let run = run_fr(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .fixed_dispersion(&f.dispersions)
                    .result(ResultColumn::Numeric(&f.numeric)),
            ),
        );
        check(
            &mut report,
            "fit_fixed_dispersion_lrt_contrast/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method = rl(b.fit_fixed_dispersion_lrt_contrast_with_cooks_replacement(
            &f.counts,
            &f.design,
            &f.reduced,
            &f.dispersions,
            &f.numeric,
            &f.cooks,
        ));
        let run = run_rl(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .fixed_dispersion(&f.dispersions)
                    .result(ResultColumn::Numeric(&f.numeric))
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_fixed_dispersion_lrt_contrast_with_cooks_replacement/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method = fr(b.fit_fixed_dispersion_lrt_contrast_spec(
            &f.counts,
            &f.design,
            &f.reduced,
            &f.dispersions,
            &f.spec,
        ));
        let run = run_fr(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .fixed_dispersion(&f.dispersions)
                    .result(ResultColumn::Spec(&f.spec)),
            ),
        );
        check(
            &mut report,
            "fit_fixed_dispersion_lrt_contrast_spec/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method = rl(
            b.fit_fixed_dispersion_lrt_contrast_spec_with_cooks_replacement(
                &f.counts,
                &f.design,
                &f.reduced,
                &f.dispersions,
                &f.spec,
                &f.cooks,
            ),
        );
        let run = run_rl(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .fixed_dispersion(&f.dispersions)
                    .result(ResultColumn::Spec(&f.spec))
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_fixed_dispersion_lrt_contrast_spec_with_cooks_replacement/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method = fr(b.fit_fixed_dispersion_lrt_factor_level_contrast(
            &f.counts,
            &f.design,
            &f.reduced,
            &f.dispersions,
            f.factor_level(),
        ));
        let run = run_fr(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .fixed_dispersion(&f.dispersions)
                    .result(ResultColumn::FactorLevel(f.factor_level())),
            ),
        );
        check(
            &mut report,
            "fit_fixed_dispersion_lrt_factor_level_contrast/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method = rl(
            b.fit_fixed_dispersion_lrt_factor_level_contrast_with_cooks_replacement(
                &f.counts,
                &f.design,
                &f.reduced,
                &f.dispersions,
                f.factor_level(),
                &f.cooks,
            ),
        );
        let run = run_rl(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .fixed_dispersion(&f.dispersions)
                    .result(ResultColumn::FactorLevel(f.factor_level()))
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_fixed_dispersion_lrt_factor_level_contrast_with_cooks_replacement/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method = fr(b.fit_fixed_dispersion_wald(&f.counts, &f.design, &f.dispersions, 2));
        let run = run_fr(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .fixed_dispersion(&f.dispersions)
                    .result(ResultColumn::Coefficient(2)),
            ),
        );
        check(&mut report, "fit_fixed_dispersion_wald/wald", method, run);
    }
    {
        let b = f.builder("wald");

        let method = fr(b.fit_fixed_dispersion_wald_contrast(
            &f.counts,
            &f.design,
            &f.dispersions,
            &f.numeric,
        ));
        let run = run_fr(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .fixed_dispersion(&f.dispersions)
                    .result(ResultColumn::Numeric(&f.numeric)),
            ),
        );
        check(
            &mut report,
            "fit_fixed_dispersion_wald_contrast/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method = rw(b.fit_fixed_dispersion_wald_with_cooks_replacement(
            &f.counts,
            &f.design,
            &f.dispersions,
            2,
            &f.cooks,
        ));
        let run = run_rw(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .fixed_dispersion(&f.dispersions)
                    .result(ResultColumn::Coefficient(2))
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_fixed_dispersion_wald_with_cooks_replacement/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method = rw(b.fit_fixed_dispersion_wald_contrast_with_cooks_replacement(
            &f.counts,
            &f.design,
            &f.dispersions,
            &f.numeric,
            &f.cooks,
        ));
        let run = run_rw(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .fixed_dispersion(&f.dispersions)
                    .result(ResultColumn::Numeric(&f.numeric))
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_fixed_dispersion_wald_contrast_with_cooks_replacement/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method = rw(
            b.fit_fixed_dispersion_wald_contrast_spec_with_cooks_replacement(
                &f.counts,
                &f.design,
                &f.dispersions,
                &f.spec,
                &f.cooks,
            ),
        );
        let run = run_rw(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .fixed_dispersion(&f.dispersions)
                    .result(ResultColumn::Spec(&f.spec))
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_fixed_dispersion_wald_contrast_spec_with_cooks_replacement/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method = fr(b.fit_fixed_dispersion_wald_contrast_spec(
            &f.counts,
            &f.design,
            &f.dispersions,
            &f.spec,
        ));
        let run = run_fr(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .fixed_dispersion(&f.dispersions)
                    .result(ResultColumn::Spec(&f.spec)),
            ),
        );
        check(
            &mut report,
            "fit_fixed_dispersion_wald_contrast_spec/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method = fr(b.fit_fixed_dispersion_wald_factor_level_contrast(
            &f.counts,
            &f.design,
            &f.dispersions,
            f.factor_level(),
        ));
        let run = run_fr(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .fixed_dispersion(&f.dispersions)
                    .result(ResultColumn::FactorLevel(f.factor_level())),
            ),
        );
        check(
            &mut report,
            "fit_fixed_dispersion_wald_factor_level_contrast/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method = rw(
            b.fit_fixed_dispersion_wald_factor_level_contrast_with_cooks_replacement(
                &f.counts,
                &f.design,
                &f.dispersions,
                f.factor_level(),
                &f.cooks,
            ),
        );
        let run = run_rw(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .fixed_dispersion(&f.dispersions)
                    .result(ResultColumn::FactorLevel(f.factor_level()))
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_fixed_dispersion_wald_factor_level_contrast_with_cooks_replacement/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method = fr(b.fit_wald_glm_mu_parametric(&f.counts, &f.design, 2));
        let run = run_fr(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .means(MeanModel::GlmMu)
                    .fit_type(FitType::Parametric)
                    .result(ResultColumn::Coefficient(2)),
            ),
        );
        check(&mut report, "fit_wald_glm_mu_parametric/wald", method, run);
    }
    {
        let b = f.builder("wald");

        let method = fr(b.fit_wald_glm_mu_contrast_parametric(&f.counts, &f.design, &f.numeric));
        let run = run_fr(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .means(MeanModel::GlmMu)
                    .fit_type(FitType::Parametric)
                    .result(ResultColumn::Numeric(&f.numeric)),
            ),
        );
        check(
            &mut report,
            "fit_wald_glm_mu_contrast_parametric/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method = fr(b.fit_wald_glm_mu_contrast_spec_parametric(&f.counts, &f.design, &f.spec));
        let run = run_fr(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .means(MeanModel::GlmMu)
                    .fit_type(FitType::Parametric)
                    .result(ResultColumn::Spec(&f.spec)),
            ),
        );
        check(
            &mut report,
            "fit_wald_glm_mu_contrast_spec_parametric/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method = fr(b.fit_wald_glm_mu_factor_level_contrast_parametric(
            &f.counts,
            &f.design,
            f.factor_level(),
        ));
        let run = run_fr(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .means(MeanModel::GlmMu)
                    .fit_type(FitType::Parametric)
                    .result(ResultColumn::FactorLevel(f.factor_level())),
            ),
        );
        check(
            &mut report,
            "fit_wald_glm_mu_factor_level_contrast_parametric/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method = fo(b.fit_gene_wise_dispersions_linear_mu(&f.counts, &f.design));
        let run = run_fo(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .means(MeanModel::LinearMu)
                    .stop_after(Stage::GeneWise),
            ),
        );
        check(
            &mut report,
            "fit_gene_wise_dispersions_linear_mu/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method = fo(b.fit_gene_wise_dispersions_glm_mu(&f.counts, &f.design));
        let run = run_fo(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .means(MeanModel::GlmMu)
                    .stop_after(Stage::GeneWise),
            ),
        );
        check(
            &mut report,
            "fit_gene_wise_dispersions_glm_mu/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method = fo(b.fit_map_dispersions_linear_mu(&f.counts, &f.design));
        let run = run_fo(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .means(MeanModel::LinearMu)
                    .stop_after(Stage::Map),
            ),
        );
        check(
            &mut report,
            "fit_map_dispersions_linear_mu/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method = fo(b.fit_map_dispersions_linear_mu_parametric(&f.counts, &f.design));
        let run = run_fo(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .means(MeanModel::LinearMu)
                    .fit_type(FitType::Parametric)
                    .stop_after(Stage::Map),
            ),
        );
        check(
            &mut report,
            "fit_map_dispersions_linear_mu_parametric/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method = fo(b.fit_map_dispersions_glm_mu(&f.counts, &f.design));
        let run = run_fo(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .means(MeanModel::GlmMu)
                    .stop_after(Stage::Map),
            ),
        );
        check(&mut report, "fit_map_dispersions_glm_mu/wald", method, run);
    }
    {
        let b = f.builder("wald");

        let method = fo(b.fit_map_dispersions_glm_mu_parametric(&f.counts, &f.design));
        let run = run_fo(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .means(MeanModel::GlmMu)
                    .fit_type(FitType::Parametric)
                    .stop_after(Stage::Map),
            ),
        );
        check(
            &mut report,
            "fit_map_dispersions_glm_mu_parametric/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method = fo(b.fit_size_factors_and_base_means(&f.counts));
        let run = run_fo(b.run(&Analysis::without_design(&f.counts)));
        check(
            &mut report,
            "fit_size_factors_and_base_means/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method = fo(b.fit_size_factors_and_base_means_with_design(&f.counts, &f.design));
        let run =
            run_fo(b.run(&Analysis::new(&f.counts, &f.design).stop_after(Stage::SizeFactors)));
        check(
            &mut report,
            "fit_size_factors_and_base_means_with_design/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method = fo(b.fit_parametric_dispersion_trend_linear_mu(&f.counts, &f.design));
        let run = run_fo(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .means(MeanModel::LinearMu)
                    .fit_type(FitType::Parametric)
                    .stop_after(Stage::Trend),
            ),
        );
        check(
            &mut report,
            "fit_parametric_dispersion_trend_linear_mu/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method = fo(b.fit_mean_dispersion_trend_linear_mu(&f.counts, &f.design));
        let run = run_fo(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .means(MeanModel::LinearMu)
                    .fit_type(FitType::Mean)
                    .stop_after(Stage::Trend),
            ),
        );
        check(
            &mut report,
            "fit_mean_dispersion_trend_linear_mu/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method = fo(b.fit_parametric_dispersion_trend_glm_mu(&f.counts, &f.design));
        let run = run_fo(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .means(MeanModel::GlmMu)
                    .fit_type(FitType::Parametric)
                    .stop_after(Stage::Trend),
            ),
        );
        check(
            &mut report,
            "fit_parametric_dispersion_trend_glm_mu/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method = fo(b.fit_mean_dispersion_trend_glm_mu(&f.counts, &f.design));
        let run = run_fo(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .means(MeanModel::GlmMu)
                    .fit_type(FitType::Mean)
                    .stop_after(Stage::Trend),
            ),
        );
        check(
            &mut report,
            "fit_mean_dispersion_trend_glm_mu/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method = fo(b.fit_local_dispersion_trend_linear_mu(&f.counts, &f.design));
        let run = run_fo(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .means(MeanModel::LinearMu)
                    .fit_type(FitType::Local)
                    .stop_after(Stage::Trend),
            ),
        );
        check(
            &mut report,
            "fit_local_dispersion_trend_linear_mu/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method = fo(b.fit_local_dispersion_trend_glm_mu(&f.counts, &f.design));
        let run = run_fo(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .means(MeanModel::GlmMu)
                    .fit_type(FitType::Local)
                    .stop_after(Stage::Trend),
            ),
        );
        check(
            &mut report,
            "fit_local_dispersion_trend_glm_mu/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method = fo(b.fit_dispersion_trend_linear_mu(&f.counts, &f.design));
        let run = run_fo(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .means(MeanModel::LinearMu)
                    .stop_after(Stage::Trend),
            ),
        );
        check(
            &mut report,
            "fit_dispersion_trend_linear_mu/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method = fo(b.fit_dispersion_trend_glm_mu(&f.counts, &f.design));
        let run = run_fo(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .means(MeanModel::GlmMu)
                    .stop_after(Stage::Trend),
            ),
        );
        check(&mut report, "fit_dispersion_trend_glm_mu/wald", method, run);
    }
    {
        let b = f.builder("wald");

        let method = fr(b.fit_wald_linear_mu(&f.counts, &f.design, 2));
        let run = run_fr(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .means(MeanModel::LinearMu)
                    .result(ResultColumn::Coefficient(2)),
            ),
        );
        check(&mut report, "fit_wald_linear_mu/wald", method, run);
    }
    {
        let b = f.builder("wald");

        let method = fr(b.fit_wald_linear_mu_contrast(&f.counts, &f.design, &f.numeric));
        let run = run_fr(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .means(MeanModel::LinearMu)
                    .result(ResultColumn::Numeric(&f.numeric)),
            ),
        );
        check(&mut report, "fit_wald_linear_mu_contrast/wald", method, run);
    }
    {
        let b = f.builder("wald");

        let method = fr(b.fit_wald_linear_mu_contrast_spec(&f.counts, &f.design, &f.spec));
        let run = run_fr(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .means(MeanModel::LinearMu)
                    .result(ResultColumn::Spec(&f.spec)),
            ),
        );
        check(
            &mut report,
            "fit_wald_linear_mu_contrast_spec/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method =
            fr(b.fit_wald_linear_mu_factor_level_contrast(&f.counts, &f.design, f.factor_level()));
        let run = run_fr(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .means(MeanModel::LinearMu)
                    .result(ResultColumn::FactorLevel(f.factor_level())),
            ),
        );
        check(
            &mut report,
            "fit_wald_linear_mu_factor_level_contrast/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method = fr(b.fit_wald_linear_mu_parametric(&f.counts, &f.design, 2));
        let run = run_fr(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .means(MeanModel::LinearMu)
                    .fit_type(FitType::Parametric)
                    .result(ResultColumn::Coefficient(2)),
            ),
        );
        check(
            &mut report,
            "fit_wald_linear_mu_parametric/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method = fr(b.fit_wald_linear_mu_contrast_parametric(&f.counts, &f.design, &f.numeric));
        let run = run_fr(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .means(MeanModel::LinearMu)
                    .fit_type(FitType::Parametric)
                    .result(ResultColumn::Numeric(&f.numeric)),
            ),
        );
        check(
            &mut report,
            "fit_wald_linear_mu_contrast_parametric/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method =
            fr(b.fit_wald_linear_mu_contrast_spec_parametric(&f.counts, &f.design, &f.spec));
        let run = run_fr(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .means(MeanModel::LinearMu)
                    .fit_type(FitType::Parametric)
                    .result(ResultColumn::Spec(&f.spec)),
            ),
        );
        check(
            &mut report,
            "fit_wald_linear_mu_contrast_spec_parametric/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method = fr(b.fit_wald_linear_mu_factor_level_contrast_parametric(
            &f.counts,
            &f.design,
            f.factor_level(),
        ));
        let run = run_fr(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .means(MeanModel::LinearMu)
                    .fit_type(FitType::Parametric)
                    .result(ResultColumn::FactorLevel(f.factor_level())),
            ),
        );
        check(
            &mut report,
            "fit_wald_linear_mu_factor_level_contrast_parametric/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method = fr(b.fit_wald_glm_mu(&f.counts, &f.design, 2));
        let run = run_fr(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .means(MeanModel::GlmMu)
                    .result(ResultColumn::Coefficient(2)),
            ),
        );
        check(&mut report, "fit_wald_glm_mu/wald", method, run);
    }
    {
        let b = f.builder("wald");

        let method = fr(b.fit_wald_glm_mu_contrast(&f.counts, &f.design, &f.numeric));
        let run = run_fr(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .means(MeanModel::GlmMu)
                    .result(ResultColumn::Numeric(&f.numeric)),
            ),
        );
        check(&mut report, "fit_wald_glm_mu_contrast/wald", method, run);
    }
    {
        let b = f.builder("wald");

        let method = fr(b.fit_wald_glm_mu_contrast_spec(&f.counts, &f.design, &f.spec));
        let run = run_fr(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .means(MeanModel::GlmMu)
                    .result(ResultColumn::Spec(&f.spec)),
            ),
        );
        check(
            &mut report,
            "fit_wald_glm_mu_contrast_spec/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method =
            fr(b.fit_wald_glm_mu_factor_level_contrast(&f.counts, &f.design, f.factor_level()));
        let run = run_fr(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .means(MeanModel::GlmMu)
                    .result(ResultColumn::FactorLevel(f.factor_level())),
            ),
        );
        check(
            &mut report,
            "fit_wald_glm_mu_factor_level_contrast/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method = fr(b.fit_lrt_linear_mu(&f.counts, &f.design, &f.reduced, 2));
        let run = run_fr(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .means(MeanModel::LinearMu)
                    .result(ResultColumn::Coefficient(2)),
            ),
        );
        check(&mut report, "fit_lrt_linear_mu/wald", method, run);
    }
    {
        let b = f.builder("wald");

        let method = fr(b.fit_lrt_linear_mu_contrast(&f.counts, &f.design, &f.reduced, &f.numeric));
        let run = run_fr(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .means(MeanModel::LinearMu)
                    .result(ResultColumn::Numeric(&f.numeric)),
            ),
        );
        check(&mut report, "fit_lrt_linear_mu_contrast/wald", method, run);
    }
    {
        let b = f.builder("wald");

        let method =
            fr(b.fit_lrt_linear_mu_contrast_spec(&f.counts, &f.design, &f.reduced, &f.spec));
        let run = run_fr(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .means(MeanModel::LinearMu)
                    .result(ResultColumn::Spec(&f.spec)),
            ),
        );
        check(
            &mut report,
            "fit_lrt_linear_mu_contrast_spec/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method = fr(b.fit_lrt_linear_mu_factor_level_contrast(
            &f.counts,
            &f.design,
            &f.reduced,
            f.factor_level(),
        ));
        let run = run_fr(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .means(MeanModel::LinearMu)
                    .result(ResultColumn::FactorLevel(f.factor_level())),
            ),
        );
        check(
            &mut report,
            "fit_lrt_linear_mu_factor_level_contrast/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method = fr(b.fit_lrt_linear_mu_parametric(&f.counts, &f.design, &f.reduced, 2));
        let run = run_fr(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .means(MeanModel::LinearMu)
                    .fit_type(FitType::Parametric)
                    .result(ResultColumn::Coefficient(2)),
            ),
        );
        check(
            &mut report,
            "fit_lrt_linear_mu_parametric/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method = fr(
            b.fit_lrt_linear_mu_contrast_parametric(&f.counts, &f.design, &f.reduced, &f.numeric)
        );
        let run = run_fr(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .means(MeanModel::LinearMu)
                    .fit_type(FitType::Parametric)
                    .result(ResultColumn::Numeric(&f.numeric)),
            ),
        );
        check(
            &mut report,
            "fit_lrt_linear_mu_contrast_parametric/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method = fr(
            b.fit_lrt_linear_mu_contrast_spec_parametric(&f.counts, &f.design, &f.reduced, &f.spec)
        );
        let run = run_fr(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .means(MeanModel::LinearMu)
                    .fit_type(FitType::Parametric)
                    .result(ResultColumn::Spec(&f.spec)),
            ),
        );
        check(
            &mut report,
            "fit_lrt_linear_mu_contrast_spec_parametric/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method = fr(b.fit_lrt_linear_mu_factor_level_contrast_parametric(
            &f.counts,
            &f.design,
            &f.reduced,
            f.factor_level(),
        ));
        let run = run_fr(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .means(MeanModel::LinearMu)
                    .fit_type(FitType::Parametric)
                    .result(ResultColumn::FactorLevel(f.factor_level())),
            ),
        );
        check(
            &mut report,
            "fit_lrt_linear_mu_factor_level_contrast_parametric/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method = fr(b.fit_lrt_glm_mu(&f.counts, &f.design, &f.reduced, 2));
        let run = run_fr(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .means(MeanModel::GlmMu)
                    .result(ResultColumn::Coefficient(2)),
            ),
        );
        check(&mut report, "fit_lrt_glm_mu/wald", method, run);
    }
    {
        let b = f.builder("wald");

        let method = fr(b.fit_lrt_glm_mu_contrast(&f.counts, &f.design, &f.reduced, &f.numeric));
        let run = run_fr(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .means(MeanModel::GlmMu)
                    .result(ResultColumn::Numeric(&f.numeric)),
            ),
        );
        check(&mut report, "fit_lrt_glm_mu_contrast/wald", method, run);
    }
    {
        let b = f.builder("wald");

        let method = fr(b.fit_lrt_glm_mu_contrast_spec(&f.counts, &f.design, &f.reduced, &f.spec));
        let run = run_fr(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .means(MeanModel::GlmMu)
                    .result(ResultColumn::Spec(&f.spec)),
            ),
        );
        check(
            &mut report,
            "fit_lrt_glm_mu_contrast_spec/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method = fr(b.fit_lrt_glm_mu_factor_level_contrast(
            &f.counts,
            &f.design,
            &f.reduced,
            f.factor_level(),
        ));
        let run = run_fr(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .means(MeanModel::GlmMu)
                    .result(ResultColumn::FactorLevel(f.factor_level())),
            ),
        );
        check(
            &mut report,
            "fit_lrt_glm_mu_factor_level_contrast/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method = fr(b.fit_lrt_glm_mu_parametric(&f.counts, &f.design, &f.reduced, 2));
        let run = run_fr(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .means(MeanModel::GlmMu)
                    .fit_type(FitType::Parametric)
                    .result(ResultColumn::Coefficient(2)),
            ),
        );
        check(&mut report, "fit_lrt_glm_mu_parametric/wald", method, run);
    }
    {
        let b = f.builder("wald");

        let method =
            fr(b.fit_lrt_glm_mu_contrast_parametric(&f.counts, &f.design, &f.reduced, &f.numeric));
        let run = run_fr(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .means(MeanModel::GlmMu)
                    .fit_type(FitType::Parametric)
                    .result(ResultColumn::Numeric(&f.numeric)),
            ),
        );
        check(
            &mut report,
            "fit_lrt_glm_mu_contrast_parametric/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method = fr(
            b.fit_lrt_glm_mu_contrast_spec_parametric(&f.counts, &f.design, &f.reduced, &f.spec)
        );
        let run = run_fr(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .means(MeanModel::GlmMu)
                    .fit_type(FitType::Parametric)
                    .result(ResultColumn::Spec(&f.spec)),
            ),
        );
        check(
            &mut report,
            "fit_lrt_glm_mu_contrast_spec_parametric/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method = fr(b.fit_lrt_glm_mu_factor_level_contrast_parametric(
            &f.counts,
            &f.design,
            &f.reduced,
            f.factor_level(),
        ));
        let run = run_fr(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .means(MeanModel::GlmMu)
                    .fit_type(FitType::Parametric)
                    .result(ResultColumn::FactorLevel(f.factor_level())),
            ),
        );
        check(
            &mut report,
            "fit_lrt_glm_mu_factor_level_contrast_parametric/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method = fr(b.fit_with_results_contrast(&f.counts, &f.design, &f.numeric));
        let run = run_fr(
            b.run(&Analysis::new(&f.counts, &f.design).result(ResultColumn::Numeric(&f.numeric))),
        );
        check(&mut report, "fit_with_results_contrast/wald", method, run);
    }
    {
        let b = f.builder("lrt");

        let method = fr(b.fit_with_results_contrast(&f.counts, &f.design, &f.numeric));
        let run = run_fr(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .result(ResultColumn::Numeric(&f.numeric)),
            ),
        );
        check(&mut report, "fit_with_results_contrast/lrt", method, run);
    }
    {
        let b = f.builder("wald");

        let method = fo(b.fit_contrast(&f.counts, &f.design, &f.numeric));
        let run = run_fo(
            b.run(&Analysis::new(&f.counts, &f.design).result(ResultColumn::Numeric(&f.numeric))),
        );
        check(&mut report, "fit_contrast/wald", method, run);
    }
    {
        let b = f.builder("lrt");

        let method = fo(b.fit_contrast(&f.counts, &f.design, &f.numeric));
        let run = run_fo(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .result(ResultColumn::Numeric(&f.numeric)),
            ),
        );
        check(&mut report, "fit_contrast/lrt", method, run);
    }
    {
        let b = f.builder("wald");

        let method = fr(b.fit_with_results_contrast_spec(&f.counts, &f.design, &f.spec));
        let run =
            run_fr(b.run(&Analysis::new(&f.counts, &f.design).result(ResultColumn::Spec(&f.spec))));
        check(
            &mut report,
            "fit_with_results_contrast_spec/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("lrt");

        let method = fr(b.fit_with_results_contrast_spec(&f.counts, &f.design, &f.spec));
        let run = run_fr(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .result(ResultColumn::Spec(&f.spec)),
            ),
        );
        check(
            &mut report,
            "fit_with_results_contrast_spec/lrt",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");
        let rc = f.character();
        let method = fr(b.fit_with_results_contrast_request(
            &f.counts,
            &f.design,
            &rc,
            Some(f.levels.as_slice()),
        ));
        let run = run_fr(b.run(&Analysis::new(&f.counts, &f.design).result(
            ResultColumn::Results {
                contrast: &rc,
                sample_levels: Some(f.levels.as_slice()),
                model_frame: None,
            },
        )));
        check(
            &mut report,
            "fit_with_results_contrast_request/wald/character/levels",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");
        let rc = f.character();
        let method =
            fr(b.fit_with_results_contrast_request(&f.counts, &f.design, &rc, None::<&[String]>));
        let run = run_fr(b.run(&Analysis::new(&f.counts, &f.design).result(
            ResultColumn::Results {
                contrast: &rc,
                sample_levels: None,
                model_frame: None,
            },
        )));
        check(
            &mut report,
            "fit_with_results_contrast_request/wald/character/nolevels",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");
        let rc = f.numeric_request();
        let method = fr(b.fit_with_results_contrast_request(
            &f.counts,
            &f.design,
            &rc,
            Some(f.levels.as_slice()),
        ));
        let run = run_fr(b.run(&Analysis::new(&f.counts, &f.design).result(
            ResultColumn::Results {
                contrast: &rc,
                sample_levels: Some(f.levels.as_slice()),
                model_frame: None,
            },
        )));
        check(
            &mut report,
            "fit_with_results_contrast_request/wald/numeric/levels",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");
        let rc = f.numeric_request();
        let method =
            fr(b.fit_with_results_contrast_request(&f.counts, &f.design, &rc, None::<&[String]>));
        let run = run_fr(b.run(&Analysis::new(&f.counts, &f.design).result(
            ResultColumn::Results {
                contrast: &rc,
                sample_levels: None,
                model_frame: None,
            },
        )));
        check(
            &mut report,
            "fit_with_results_contrast_request/wald/numeric/nolevels",
            method,
            run,
        );
    }
    {
        let b = f.builder("lrt");
        let rc = f.character();
        let method = fr(b.fit_with_results_contrast_request(
            &f.counts,
            &f.design,
            &rc,
            Some(f.levels.as_slice()),
        ));
        let run = run_fr(
            b.run(&Analysis::new(&f.counts, &f.design).lrt(&f.reduced).result(
                ResultColumn::Results {
                    contrast: &rc,
                    sample_levels: Some(f.levels.as_slice()),
                    model_frame: None,
                },
            )),
        );
        check(
            &mut report,
            "fit_with_results_contrast_request/lrt/character/levels",
            method,
            run,
        );
    }
    {
        let b = f.builder("lrt");
        let rc = f.character();
        let method =
            fr(b.fit_with_results_contrast_request(&f.counts, &f.design, &rc, None::<&[String]>));
        let run = run_fr(
            b.run(&Analysis::new(&f.counts, &f.design).lrt(&f.reduced).result(
                ResultColumn::Results {
                    contrast: &rc,
                    sample_levels: None,
                    model_frame: None,
                },
            )),
        );
        check(
            &mut report,
            "fit_with_results_contrast_request/lrt/character/nolevels",
            method,
            run,
        );
    }
    {
        let b = f.builder("lrt");
        let rc = f.numeric_request();
        let method = fr(b.fit_with_results_contrast_request(
            &f.counts,
            &f.design,
            &rc,
            Some(f.levels.as_slice()),
        ));
        let run = run_fr(
            b.run(&Analysis::new(&f.counts, &f.design).lrt(&f.reduced).result(
                ResultColumn::Results {
                    contrast: &rc,
                    sample_levels: Some(f.levels.as_slice()),
                    model_frame: None,
                },
            )),
        );
        check(
            &mut report,
            "fit_with_results_contrast_request/lrt/numeric/levels",
            method,
            run,
        );
    }
    {
        let b = f.builder("lrt");
        let rc = f.numeric_request();
        let method =
            fr(b.fit_with_results_contrast_request(&f.counts, &f.design, &rc, None::<&[String]>));
        let run = run_fr(
            b.run(&Analysis::new(&f.counts, &f.design).lrt(&f.reduced).result(
                ResultColumn::Results {
                    contrast: &rc,
                    sample_levels: None,
                    model_frame: None,
                },
            )),
        );
        check(
            &mut report,
            "fit_with_results_contrast_request/lrt/numeric/nolevels",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");
        let rc = f.character();
        let method =
            fr(b.fit_formula_with_results_contrast_request(&f.counts, f.full_formula, &rc));
        let run = run_fr(b.run(&Analysis::formula(&f.counts, f.full_formula).result(
            ResultColumn::Results {
                contrast: &rc,
                sample_levels: None,
                model_frame: None,
            },
        )));
        check(
            &mut report,
            "fit_formula_with_results_contrast_request/wald/character",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");
        let rc = f.numeric_request();
        let method =
            fr(b.fit_formula_with_results_contrast_request(&f.counts, f.full_formula, &rc));
        let run = run_fr(b.run(&Analysis::formula(&f.counts, f.full_formula).result(
            ResultColumn::Results {
                contrast: &rc,
                sample_levels: None,
                model_frame: None,
            },
        )));
        check(
            &mut report,
            "fit_formula_with_results_contrast_request/wald/numeric",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");
        let rc = f.character();
        let method = fr(b.fit_with_results_contrast_request_from_model_frame(
            &f.counts, &f.design, &rc, &f.frame,
        ));
        let run = run_fr(b.run(&Analysis::new(&f.counts, &f.design).result(
            ResultColumn::Results {
                contrast: &rc,
                sample_levels: None,
                model_frame: Some(&f.frame),
            },
        )));
        check(
            &mut report,
            "fit_with_results_contrast_request_from_model_frame/wald/character",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");
        let rc = f.numeric_request();
        let method = fr(b.fit_with_results_contrast_request_from_model_frame(
            &f.counts, &f.design, &rc, &f.frame,
        ));
        let run = run_fr(b.run(&Analysis::new(&f.counts, &f.design).result(
            ResultColumn::Results {
                contrast: &rc,
                sample_levels: None,
                model_frame: Some(&f.frame),
            },
        )));
        check(
            &mut report,
            "fit_with_results_contrast_request_from_model_frame/wald/numeric",
            method,
            run,
        );
    }
    {
        let b = f.builder("lrt");
        let rc = f.character();
        let method = fr(b.fit_with_results_contrast_request_from_model_frame(
            &f.counts, &f.design, &rc, &f.frame,
        ));
        let run = run_fr(
            b.run(&Analysis::new(&f.counts, &f.design).lrt(&f.reduced).result(
                ResultColumn::Results {
                    contrast: &rc,
                    sample_levels: None,
                    model_frame: Some(&f.frame),
                },
            )),
        );
        check(
            &mut report,
            "fit_with_results_contrast_request_from_model_frame/lrt/character",
            method,
            run,
        );
    }
    {
        let b = f.builder("lrt");
        let rc = f.numeric_request();
        let method = fr(b.fit_with_results_contrast_request_from_model_frame(
            &f.counts, &f.design, &rc, &f.frame,
        ));
        let run = run_fr(
            b.run(&Analysis::new(&f.counts, &f.design).lrt(&f.reduced).result(
                ResultColumn::Results {
                    contrast: &rc,
                    sample_levels: None,
                    model_frame: Some(&f.frame),
                },
            )),
        );
        check(
            &mut report,
            "fit_with_results_contrast_request_from_model_frame/lrt/numeric",
            method,
            run,
        );
    }
    {
        let b = f.builder("bare");
        let rc = f.character();
        let method = fr(b.fit_with_results_contrast_request_from_model_frame(
            &f.counts, &f.design, &rc, &f.frame,
        ));
        let run = run_fr(b.run(&Analysis::new(&f.counts, &f.design).result(
            ResultColumn::Results {
                contrast: &rc,
                sample_levels: None,
                model_frame: Some(&f.frame),
            },
        )));
        check(
            &mut report,
            "fit_with_results_contrast_request_from_model_frame/bare/character",
            method,
            run,
        );
    }
    {
        let b = f.builder("bare");
        let rc = f.numeric_request();
        let method = fr(b.fit_with_results_contrast_request_from_model_frame(
            &f.counts, &f.design, &rc, &f.frame,
        ));
        let run = run_fr(b.run(&Analysis::new(&f.counts, &f.design).result(
            ResultColumn::Results {
                contrast: &rc,
                sample_levels: None,
                model_frame: Some(&f.frame),
            },
        )));
        check(
            &mut report,
            "fit_with_results_contrast_request_from_model_frame/bare/numeric",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method = fo(b.fit_contrast_spec(&f.counts, &f.design, &f.spec));
        let run =
            run_fo(b.run(&Analysis::new(&f.counts, &f.design).result(ResultColumn::Spec(&f.spec))));
        check(&mut report, "fit_contrast_spec/wald", method, run);
    }
    {
        let b = f.builder("lrt");

        let method = fo(b.fit_contrast_spec(&f.counts, &f.design, &f.spec));
        let run = run_fo(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .result(ResultColumn::Spec(&f.spec)),
            ),
        );
        check(&mut report, "fit_contrast_spec/lrt", method, run);
    }
    {
        let b = f.builder("wald");

        let method =
            fr(b.fit_with_results_factor_level_contrast(&f.counts, &f.design, f.factor_level()));
        let run = run_fr(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .result(ResultColumn::FactorLevel(f.factor_level())),
            ),
        );
        check(
            &mut report,
            "fit_with_results_factor_level_contrast/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("lrt");

        let method =
            fr(b.fit_with_results_factor_level_contrast(&f.counts, &f.design, f.factor_level()));
        let run = run_fr(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .result(ResultColumn::FactorLevel(f.factor_level())),
            ),
        );
        check(
            &mut report,
            "fit_with_results_factor_level_contrast/lrt",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method = fo(b.fit_factor_level_contrast(&f.counts, &f.design, f.factor_level()));
        let run = run_fo(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .result(ResultColumn::FactorLevel(f.factor_level())),
            ),
        );
        check(&mut report, "fit_factor_level_contrast/wald", method, run);
    }
    {
        let b = f.builder("lrt");

        let method = fo(b.fit_factor_level_contrast(&f.counts, &f.design, f.factor_level()));
        let run = run_fo(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .result(ResultColumn::FactorLevel(f.factor_level())),
            ),
        );
        check(&mut report, "fit_factor_level_contrast/lrt", method, run);
    }
    {
        let b = f.builder("wald");

        let method = fo(b.fit(&f.counts, &f.design));
        let run = run_fo(b.run(&Analysis::new(&f.counts, &f.design)));
        check(&mut report, "fit/wald", method, run);
    }
    {
        let b = f.builder("lrt");

        let method = fo(b.fit(&f.counts, &f.design));
        let run = run_fo(b.run(&Analysis::new(&f.counts, &f.design).lrt(&f.reduced)));
        check(&mut report, "fit/lrt", method, run);
    }
    {
        let b = f.builder("wald");

        let method = fo(b.fit_formula(&f.counts, f.full_formula));
        let run = run_fo(b.run(&Analysis::formula(&f.counts, f.full_formula)));
        check(&mut report, "fit_formula/wald", method, run);
    }
    {
        let b = f.builder("lrt");

        let method = fo(b.fit_formula(&f.counts, f.full_formula));
        let run = run_fo(b.run(&Analysis::formula(&f.counts, f.full_formula).lrt(&f.reduced)));
        check(&mut report, "fit_formula/lrt", method, run);
    }
    {
        let b = f.builder("wald");

        let method = fr(b.fit_with_results(&f.counts, &f.design));
        let run = run_fr(b.run(&Analysis::new(&f.counts, &f.design)));
        check(&mut report, "fit_with_results/wald", method, run);
    }
    {
        let b = f.builder("lrt");

        let method = fr(b.fit_with_results(&f.counts, &f.design));
        let run = run_fr(b.run(&Analysis::new(&f.counts, &f.design).lrt(&f.reduced)));
        check(&mut report, "fit_with_results/lrt", method, run);
    }
    {
        let b = f.builder("wald");

        let method = fr(b.fit_formula_with_results(&f.counts, f.full_formula));
        let run = run_fr(b.run(&Analysis::formula(&f.counts, f.full_formula)));
        check(&mut report, "fit_formula_with_results/wald", method, run);
    }
    {
        let b = f.builder("wald");

        let method = fr(b.fit_with_results_name(&f.counts, &f.design, "condition_B_vs_A"));
        let run = run_fr(b.run(
            &Analysis::new(&f.counts, &f.design).result(ResultColumn::Name("condition_B_vs_A")),
        ));
        check(&mut report, "fit_with_results_name/wald", method, run);
    }
    {
        let b = f.builder("lrt");

        let method = fr(b.fit_with_results_name(&f.counts, &f.design, "condition_B_vs_A"));
        let run = run_fr(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .result(ResultColumn::Name("condition_B_vs_A")),
            ),
        );
        check(&mut report, "fit_with_results_name/lrt", method, run);
    }
    {
        let b = f.builder("wald");

        let method =
            fr(b.fit_formula_with_results_name(&f.counts, f.full_formula, "condition_B_vs_A"));
        let run = run_fr(
            b.run(
                &Analysis::formula(&f.counts, f.full_formula)
                    .result(ResultColumn::Name("condition_B_vs_A")),
            ),
        );
        check(
            &mut report,
            "fit_formula_with_results_name/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method = fo(b.fit_formula_name(&f.counts, f.full_formula, "condition_B_vs_A"));
        let run = run_fo(
            b.run(
                &Analysis::formula(&f.counts, f.full_formula)
                    .result(ResultColumn::Name("condition_B_vs_A")),
            ),
        );
        check(&mut report, "fit_formula_name/wald", method, run);
    }
    {
        let b = f.builder("lrt");

        let method = fo(b.fit_formula_name(&f.counts, f.full_formula, "condition_B_vs_A"));
        let run = run_fo(
            b.run(
                &Analysis::formula(&f.counts, f.full_formula)
                    .lrt(&f.reduced)
                    .result(ResultColumn::Name("condition_B_vs_A")),
            ),
        );
        check(&mut report, "fit_formula_name/lrt", method, run);
    }
    {
        let b = f.builder("wald");

        let method = fo(b.fit_lrt(&f.counts, &f.design, &f.reduced));
        let run = run_fo(b.run(&Analysis::new(&f.counts, &f.design).lrt(&f.reduced)));
        check(&mut report, "fit_lrt/wald", method, run);
    }
    {
        let b = f.builder("wald");

        let method = fo(b.fit_lrt_formula(&f.counts, f.full_formula, f.reduced_formula));
        let run = run_fo(
            b.run(&Analysis::formula(&f.counts, f.full_formula).lrt_formula(f.reduced_formula)),
        );
        check(&mut report, "fit_lrt_formula/wald", method, run);
    }
    {
        let b = f.builder("wald");

        let method = fr(b.fit_lrt_with_results(&f.counts, &f.design, &f.reduced));
        let run = run_fr(b.run(&Analysis::new(&f.counts, &f.design).lrt(&f.reduced)));
        check(&mut report, "fit_lrt_with_results/wald", method, run);
    }
    {
        let b = f.builder("wald");

        let method =
            fr(b.fit_lrt_formula_with_results(&f.counts, f.full_formula, f.reduced_formula));
        let run = run_fr(
            b.run(&Analysis::formula(&f.counts, f.full_formula).lrt_formula(f.reduced_formula)),
        );
        check(
            &mut report,
            "fit_lrt_formula_with_results/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method =
            fr(b.fit_lrt_with_results_name(&f.counts, &f.design, &f.reduced, "condition_B_vs_A"));
        let run = run_fr(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .result(ResultColumn::Name("condition_B_vs_A")),
            ),
        );
        check(&mut report, "fit_lrt_with_results_name/wald", method, run);
    }
    {
        let b = f.builder("wald");

        let method = fr(b.fit_lrt_formula_with_results_name(
            &f.counts,
            f.full_formula,
            f.reduced_formula,
            "condition_B_vs_A",
        ));
        let run = run_fr(
            b.run(
                &Analysis::formula(&f.counts, f.full_formula)
                    .lrt_formula(f.reduced_formula)
                    .result(ResultColumn::Name("condition_B_vs_A")),
            ),
        );
        check(
            &mut report,
            "fit_lrt_formula_with_results_name/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method = fo(b.fit_lrt_name(&f.counts, &f.design, &f.reduced, "condition_B_vs_A"));
        let run = run_fo(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .result(ResultColumn::Name("condition_B_vs_A")),
            ),
        );
        check(&mut report, "fit_lrt_name/wald", method, run);
    }
    {
        let b = f.builder("wald");

        let method = fo(b.fit_lrt_formula_name(
            &f.counts,
            f.full_formula,
            f.reduced_formula,
            "condition_B_vs_A",
        ));
        let run = run_fo(
            b.run(
                &Analysis::formula(&f.counts, f.full_formula)
                    .lrt_formula(f.reduced_formula)
                    .result(ResultColumn::Name("condition_B_vs_A")),
            ),
        );
        check(&mut report, "fit_lrt_formula_name/wald", method, run);
    }
    {
        let b = f.builder("wald");

        let method =
            fr(b.fit_lrt_with_results_contrast(&f.counts, &f.design, &f.reduced, &f.numeric));
        let run = run_fr(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .result(ResultColumn::Numeric(&f.numeric)),
            ),
        );
        check(
            &mut report,
            "fit_lrt_with_results_contrast/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method = fo(b.fit_lrt_contrast(&f.counts, &f.design, &f.reduced, &f.numeric));
        let run = run_fo(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .result(ResultColumn::Numeric(&f.numeric)),
            ),
        );
        check(&mut report, "fit_lrt_contrast/wald", method, run);
    }
    {
        let b = f.builder("wald");

        let method =
            fr(b.fit_lrt_with_results_contrast_spec(&f.counts, &f.design, &f.reduced, &f.spec));
        let run = run_fr(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .result(ResultColumn::Spec(&f.spec)),
            ),
        );
        check(
            &mut report,
            "fit_lrt_with_results_contrast_spec/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method = fr(b.fit_lrt_with_results_factor_level_contrast(
            &f.counts,
            &f.design,
            &f.reduced,
            f.factor_level(),
        ));
        let run = run_fr(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .result(ResultColumn::FactorLevel(f.factor_level())),
            ),
        );
        check(
            &mut report,
            "fit_lrt_with_results_factor_level_contrast/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");
        let rc = f.character();
        let method = fr(b.fit_lrt_with_results_contrast_request(
            &f.counts,
            &f.design,
            &f.reduced,
            &rc,
            Some(f.levels.as_slice()),
        ));
        let run = run_fr(
            b.run(&Analysis::new(&f.counts, &f.design).lrt(&f.reduced).result(
                ResultColumn::Results {
                    contrast: &rc,
                    sample_levels: Some(f.levels.as_slice()),
                    model_frame: None,
                },
            )),
        );
        check(
            &mut report,
            "fit_lrt_with_results_contrast_request/wald/character/levels",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");
        let rc = f.character();
        let method = fr(b.fit_lrt_with_results_contrast_request(
            &f.counts,
            &f.design,
            &f.reduced,
            &rc,
            None::<&[String]>,
        ));
        let run = run_fr(
            b.run(&Analysis::new(&f.counts, &f.design).lrt(&f.reduced).result(
                ResultColumn::Results {
                    contrast: &rc,
                    sample_levels: None,
                    model_frame: None,
                },
            )),
        );
        check(
            &mut report,
            "fit_lrt_with_results_contrast_request/wald/character/nolevels",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");
        let rc = f.numeric_request();
        let method = fr(b.fit_lrt_with_results_contrast_request(
            &f.counts,
            &f.design,
            &f.reduced,
            &rc,
            Some(f.levels.as_slice()),
        ));
        let run = run_fr(
            b.run(&Analysis::new(&f.counts, &f.design).lrt(&f.reduced).result(
                ResultColumn::Results {
                    contrast: &rc,
                    sample_levels: Some(f.levels.as_slice()),
                    model_frame: None,
                },
            )),
        );
        check(
            &mut report,
            "fit_lrt_with_results_contrast_request/wald/numeric/levels",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");
        let rc = f.numeric_request();
        let method = fr(b.fit_lrt_with_results_contrast_request(
            &f.counts,
            &f.design,
            &f.reduced,
            &rc,
            None::<&[String]>,
        ));
        let run = run_fr(
            b.run(&Analysis::new(&f.counts, &f.design).lrt(&f.reduced).result(
                ResultColumn::Results {
                    contrast: &rc,
                    sample_levels: None,
                    model_frame: None,
                },
            )),
        );
        check(
            &mut report,
            "fit_lrt_with_results_contrast_request/wald/numeric/nolevels",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");
        let rc = f.character();
        let method = fr(b.fit_lrt_with_results_contrast_request_from_model_frame(
            &f.counts, &f.design, &f.reduced, &rc, &f.frame,
        ));
        let run = run_fr(
            b.run(&Analysis::new(&f.counts, &f.design).lrt(&f.reduced).result(
                ResultColumn::Results {
                    contrast: &rc,
                    sample_levels: None,
                    model_frame: Some(&f.frame),
                },
            )),
        );
        check(
            &mut report,
            "fit_lrt_with_results_contrast_request_from_model_frame/wald/character",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");
        let rc = f.numeric_request();
        let method = fr(b.fit_lrt_with_results_contrast_request_from_model_frame(
            &f.counts, &f.design, &f.reduced, &rc, &f.frame,
        ));
        let run = run_fr(
            b.run(&Analysis::new(&f.counts, &f.design).lrt(&f.reduced).result(
                ResultColumn::Results {
                    contrast: &rc,
                    sample_levels: None,
                    model_frame: Some(&f.frame),
                },
            )),
        );
        check(
            &mut report,
            "fit_lrt_with_results_contrast_request_from_model_frame/wald/numeric",
            method,
            run,
        );
    }
    {
        let b = f.builder("bare");
        let rc = f.character();
        let method = fr(b.fit_lrt_with_results_contrast_request_from_model_frame(
            &f.counts, &f.design, &f.reduced, &rc, &f.frame,
        ));
        let run = run_fr(
            b.run(&Analysis::new(&f.counts, &f.design).lrt(&f.reduced).result(
                ResultColumn::Results {
                    contrast: &rc,
                    sample_levels: None,
                    model_frame: Some(&f.frame),
                },
            )),
        );
        check(
            &mut report,
            "fit_lrt_with_results_contrast_request_from_model_frame/bare/character",
            method,
            run,
        );
    }
    {
        let b = f.builder("bare");
        let rc = f.numeric_request();
        let method = fr(b.fit_lrt_with_results_contrast_request_from_model_frame(
            &f.counts, &f.design, &f.reduced, &rc, &f.frame,
        ));
        let run = run_fr(
            b.run(&Analysis::new(&f.counts, &f.design).lrt(&f.reduced).result(
                ResultColumn::Results {
                    contrast: &rc,
                    sample_levels: None,
                    model_frame: Some(&f.frame),
                },
            )),
        );
        check(
            &mut report,
            "fit_lrt_with_results_contrast_request_from_model_frame/bare/numeric",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");
        let rc = f.character();
        let method = fr(b.fit_lrt_formula_with_results_contrast_request(
            &f.counts,
            f.full_formula,
            f.reduced_formula,
            &rc,
        ));
        let run = run_fr(
            b.run(
                &Analysis::formula(&f.counts, f.full_formula)
                    .lrt_formula(f.reduced_formula)
                    .result(ResultColumn::Results {
                        contrast: &rc,
                        sample_levels: None,
                        model_frame: None,
                    }),
            ),
        );
        check(
            &mut report,
            "fit_lrt_formula_with_results_contrast_request/wald/character",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");
        let rc = f.numeric_request();
        let method = fr(b.fit_lrt_formula_with_results_contrast_request(
            &f.counts,
            f.full_formula,
            f.reduced_formula,
            &rc,
        ));
        let run = run_fr(
            b.run(
                &Analysis::formula(&f.counts, f.full_formula)
                    .lrt_formula(f.reduced_formula)
                    .result(ResultColumn::Results {
                        contrast: &rc,
                        sample_levels: None,
                        model_frame: None,
                    }),
            ),
        );
        check(
            &mut report,
            "fit_lrt_formula_with_results_contrast_request/wald/numeric",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method = fo(b.fit_lrt_contrast_spec(&f.counts, &f.design, &f.reduced, &f.spec));
        let run = run_fo(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .result(ResultColumn::Spec(&f.spec)),
            ),
        );
        check(&mut report, "fit_lrt_contrast_spec/wald", method, run);
    }
    {
        let b = f.builder("wald");

        let method =
            fo(b.fit_lrt_factor_level_contrast(&f.counts, &f.design, &f.reduced, f.factor_level()));
        let run = run_fo(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .result(ResultColumn::FactorLevel(f.factor_level())),
            ),
        );
        check(
            &mut report,
            "fit_lrt_factor_level_contrast/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method = rl(b.fit_lrt_formula_with_results_with_cooks_replacement(
            &f.counts,
            f.full_formula,
            f.reduced_formula,
            &f.cooks,
        ));
        let run = run_rl(
            b.run(
                &Analysis::formula(&f.counts, f.full_formula)
                    .lrt_formula(f.reduced_formula)
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_lrt_formula_with_results_with_cooks_replacement/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method = rl(b.fit_lrt_formula_with_results_name_with_cooks_replacement(
            &f.counts,
            f.full_formula,
            f.reduced_formula,
            "condition_B_vs_A",
            &f.cooks,
        ));
        let run = run_rl(
            b.run(
                &Analysis::formula(&f.counts, f.full_formula)
                    .lrt_formula(f.reduced_formula)
                    .result(ResultColumn::Name("condition_B_vs_A"))
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_lrt_formula_with_results_name_with_cooks_replacement/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method = rl(b.fit_lrt_with_results_with_cooks_replacement(
            &f.counts, &f.design, &f.reduced, &f.cooks,
        ));
        let run = run_rl(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_lrt_with_results_with_cooks_replacement/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method = rl(b.fit_lrt_with_results_name_with_cooks_replacement(
            &f.counts,
            &f.design,
            &f.reduced,
            "condition_B_vs_A",
            &f.cooks,
        ));
        let run = run_rl(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .result(ResultColumn::Name("condition_B_vs_A"))
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_lrt_with_results_name_with_cooks_replacement/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method = rl(b.fit_lrt_with_results_contrast_with_cooks_replacement(
            &f.counts, &f.design, &f.reduced, &f.numeric, &f.cooks,
        ));
        let run = run_rl(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .result(ResultColumn::Numeric(&f.numeric))
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_lrt_with_results_contrast_with_cooks_replacement/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method = rl(b.fit_lrt_with_results_contrast_spec_with_cooks_replacement(
            &f.counts, &f.design, &f.reduced, &f.spec, &f.cooks,
        ));
        let run = run_rl(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .result(ResultColumn::Spec(&f.spec))
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_lrt_with_results_contrast_spec_with_cooks_replacement/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method = rl(
            b.fit_lrt_with_results_factor_level_contrast_with_cooks_replacement(
                &f.counts,
                &f.design,
                &f.reduced,
                f.factor_level(),
                &f.cooks,
            ),
        );
        let run = run_rl(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .result(ResultColumn::FactorLevel(f.factor_level()))
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_lrt_with_results_factor_level_contrast_with_cooks_replacement/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");
        let rc = f.character();
        let method = rl(
            b.fit_lrt_with_results_contrast_request_with_cooks_replacement(
                &f.counts,
                &f.design,
                &f.reduced,
                &rc,
                Some(f.levels.as_slice()),
                &f.cooks,
            ),
        );
        let run = run_rl(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .result(ResultColumn::Results {
                        contrast: &rc,
                        sample_levels: Some(f.levels.as_slice()),
                        model_frame: None,
                    })
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_lrt_with_results_contrast_request_with_cooks_replacement/wald/character/levels",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");
        let rc = f.character();
        let method = rl(
            b.fit_lrt_with_results_contrast_request_with_cooks_replacement(
                &f.counts,
                &f.design,
                &f.reduced,
                &rc,
                None::<&[String]>,
                &f.cooks,
            ),
        );
        let run = run_rl(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .result(ResultColumn::Results {
                        contrast: &rc,
                        sample_levels: None,
                        model_frame: None,
                    })
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_lrt_with_results_contrast_request_with_cooks_replacement/wald/character/nolevels",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");
        let rc = f.numeric_request();
        let method = rl(
            b.fit_lrt_with_results_contrast_request_with_cooks_replacement(
                &f.counts,
                &f.design,
                &f.reduced,
                &rc,
                Some(f.levels.as_slice()),
                &f.cooks,
            ),
        );
        let run = run_rl(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .result(ResultColumn::Results {
                        contrast: &rc,
                        sample_levels: Some(f.levels.as_slice()),
                        model_frame: None,
                    })
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_lrt_with_results_contrast_request_with_cooks_replacement/wald/numeric/levels",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");
        let rc = f.numeric_request();
        let method = rl(
            b.fit_lrt_with_results_contrast_request_with_cooks_replacement(
                &f.counts,
                &f.design,
                &f.reduced,
                &rc,
                None::<&[String]>,
                &f.cooks,
            ),
        );
        let run = run_rl(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .result(ResultColumn::Results {
                        contrast: &rc,
                        sample_levels: None,
                        model_frame: None,
                    })
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_lrt_with_results_contrast_request_with_cooks_replacement/wald/numeric/nolevels",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");
        let rc = f.character();
        let method = rl(b
            .fit_lrt_with_results_contrast_request_from_model_frame_with_cooks_replacement(
                &f.counts, &f.design, &f.reduced, &rc, &f.frame, &f.cooks,
            ));
        let run = run_rl(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .result(ResultColumn::Results {
                        contrast: &rc,
                        sample_levels: None,
                        model_frame: Some(&f.frame),
                    })
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_lrt_with_results_contrast_request_from_model_frame_with_cooks_replacement/wald/character",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");
        let rc = f.numeric_request();
        let method = rl(b
            .fit_lrt_with_results_contrast_request_from_model_frame_with_cooks_replacement(
                &f.counts, &f.design, &f.reduced, &rc, &f.frame, &f.cooks,
            ));
        let run = run_rl(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .result(ResultColumn::Results {
                        contrast: &rc,
                        sample_levels: None,
                        model_frame: Some(&f.frame),
                    })
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_lrt_with_results_contrast_request_from_model_frame_with_cooks_replacement/wald/numeric",
            method,
            run,
        );
    }
    {
        let b = f.builder("bare");
        let rc = f.character();
        let method = rl(b
            .fit_lrt_with_results_contrast_request_from_model_frame_with_cooks_replacement(
                &f.counts, &f.design, &f.reduced, &rc, &f.frame, &f.cooks,
            ));
        let run = run_rl(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .result(ResultColumn::Results {
                        contrast: &rc,
                        sample_levels: None,
                        model_frame: Some(&f.frame),
                    })
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_lrt_with_results_contrast_request_from_model_frame_with_cooks_replacement/bare/character",
            method,
            run,
        );
    }
    {
        let b = f.builder("bare");
        let rc = f.numeric_request();
        let method = rl(b
            .fit_lrt_with_results_contrast_request_from_model_frame_with_cooks_replacement(
                &f.counts, &f.design, &f.reduced, &rc, &f.frame, &f.cooks,
            ));
        let run = run_rl(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .result(ResultColumn::Results {
                        contrast: &rc,
                        sample_levels: None,
                        model_frame: Some(&f.frame),
                    })
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_lrt_with_results_contrast_request_from_model_frame_with_cooks_replacement/bare/numeric",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");
        let rc = f.character();
        let method = rl(
            b.fit_lrt_formula_with_results_contrast_request_with_cooks_replacement(
                &f.counts,
                f.full_formula,
                f.reduced_formula,
                &rc,
                &f.cooks,
            ),
        );
        let run = run_rl(
            b.run(
                &Analysis::formula(&f.counts, f.full_formula)
                    .lrt_formula(f.reduced_formula)
                    .result(ResultColumn::Results {
                        contrast: &rc,
                        sample_levels: None,
                        model_frame: None,
                    })
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_lrt_formula_with_results_contrast_request_with_cooks_replacement/wald/character",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");
        let rc = f.numeric_request();
        let method = rl(
            b.fit_lrt_formula_with_results_contrast_request_with_cooks_replacement(
                &f.counts,
                f.full_formula,
                f.reduced_formula,
                &rc,
                &f.cooks,
            ),
        );
        let run = run_rl(
            b.run(
                &Analysis::formula(&f.counts, f.full_formula)
                    .lrt_formula(f.reduced_formula)
                    .result(ResultColumn::Results {
                        contrast: &rc,
                        sample_levels: None,
                        model_frame: None,
                    })
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_lrt_formula_with_results_contrast_request_with_cooks_replacement/wald/numeric",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method = rw(b.fit_with_results_with_cooks_replacement(&f.counts, &f.design, &f.cooks));
        let run = run_rw(b.run(&Analysis::new(&f.counts, &f.design).cooks_replacement(&f.cooks)));
        check(
            &mut report,
            "fit_with_results_with_cooks_replacement/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method = rw(b.fit_formula_with_results_with_cooks_replacement(
            &f.counts,
            f.full_formula,
            &f.cooks,
        ));
        let run = run_rw(
            b.run(&Analysis::formula(&f.counts, f.full_formula).cooks_replacement(&f.cooks)),
        );
        check(
            &mut report,
            "fit_formula_with_results_with_cooks_replacement/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method = rw(b.fit_with_results_name_with_cooks_replacement(
            &f.counts,
            &f.design,
            "condition_B_vs_A",
            &f.cooks,
        ));
        let run = run_rw(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .result(ResultColumn::Name("condition_B_vs_A"))
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_with_results_name_with_cooks_replacement/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method = rw(b.fit_formula_with_results_name_with_cooks_replacement(
            &f.counts,
            f.full_formula,
            "condition_B_vs_A",
            &f.cooks,
        ));
        let run = run_rw(
            b.run(
                &Analysis::formula(&f.counts, f.full_formula)
                    .result(ResultColumn::Name("condition_B_vs_A"))
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_formula_with_results_name_with_cooks_replacement/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method =
            rt(b.fit_with_test_results_with_cooks_replacement(&f.counts, &f.design, &f.cooks));
        let run = run_rt(b.run(&Analysis::new(&f.counts, &f.design).cooks_replacement(&f.cooks)));
        check(
            &mut report,
            "fit_with_test_results_with_cooks_replacement/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("lrt");

        let method =
            rt(b.fit_with_test_results_with_cooks_replacement(&f.counts, &f.design, &f.cooks));
        let run = run_rt(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_with_test_results_with_cooks_replacement/lrt",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method = rt(b.fit_formula_with_test_results_with_cooks_replacement(
            &f.counts,
            f.full_formula,
            &f.cooks,
        ));
        let run = run_rt(
            b.run(&Analysis::formula(&f.counts, f.full_formula).cooks_replacement(&f.cooks)),
        );
        check(
            &mut report,
            "fit_formula_with_test_results_with_cooks_replacement/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("lrt");

        let method = rt(b.fit_formula_with_test_results_with_cooks_replacement(
            &f.counts,
            f.full_formula,
            &f.cooks,
        ));
        let run = run_rt(
            b.run(
                &Analysis::formula(&f.counts, f.full_formula)
                    .lrt(&f.reduced)
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_formula_with_test_results_with_cooks_replacement/lrt",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method = rt(b.fit_with_test_results_name_with_cooks_replacement(
            &f.counts,
            &f.design,
            "condition_B_vs_A",
            &f.cooks,
        ));
        let run = run_rt(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .result(ResultColumn::Name("condition_B_vs_A"))
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_with_test_results_name_with_cooks_replacement/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("lrt");

        let method = rt(b.fit_with_test_results_name_with_cooks_replacement(
            &f.counts,
            &f.design,
            "condition_B_vs_A",
            &f.cooks,
        ));
        let run = run_rt(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .result(ResultColumn::Name("condition_B_vs_A"))
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_with_test_results_name_with_cooks_replacement/lrt",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method = rt(b.fit_formula_with_test_results_name_with_cooks_replacement(
            &f.counts,
            f.full_formula,
            "condition_B_vs_A",
            &f.cooks,
        ));
        let run = run_rt(
            b.run(
                &Analysis::formula(&f.counts, f.full_formula)
                    .result(ResultColumn::Name("condition_B_vs_A"))
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_formula_with_test_results_name_with_cooks_replacement/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("lrt");

        let method = rt(b.fit_formula_with_test_results_name_with_cooks_replacement(
            &f.counts,
            f.full_formula,
            "condition_B_vs_A",
            &f.cooks,
        ));
        let run = run_rt(
            b.run(
                &Analysis::formula(&f.counts, f.full_formula)
                    .lrt(&f.reduced)
                    .result(ResultColumn::Name("condition_B_vs_A"))
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_formula_with_test_results_name_with_cooks_replacement/lrt",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method = rw(b.fit_with_results_contrast_with_cooks_replacement(
            &f.counts, &f.design, &f.numeric, &f.cooks,
        ));
        let run = run_rw(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .result(ResultColumn::Numeric(&f.numeric))
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_with_results_contrast_with_cooks_replacement/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method = rt(b.fit_with_test_results_contrast_with_cooks_replacement(
            &f.counts, &f.design, &f.numeric, &f.cooks,
        ));
        let run = run_rt(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .result(ResultColumn::Numeric(&f.numeric))
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_with_test_results_contrast_with_cooks_replacement/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("lrt");

        let method = rt(b.fit_with_test_results_contrast_with_cooks_replacement(
            &f.counts, &f.design, &f.numeric, &f.cooks,
        ));
        let run = run_rt(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .result(ResultColumn::Numeric(&f.numeric))
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_with_test_results_contrast_with_cooks_replacement/lrt",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method = rw(b.fit_with_results_contrast_spec_with_cooks_replacement(
            &f.counts, &f.design, &f.spec, &f.cooks,
        ));
        let run = run_rw(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .result(ResultColumn::Spec(&f.spec))
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_with_results_contrast_spec_with_cooks_replacement/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method = rt(
            b.fit_with_test_results_contrast_spec_with_cooks_replacement(
                &f.counts, &f.design, &f.spec, &f.cooks,
            ),
        );
        let run = run_rt(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .result(ResultColumn::Spec(&f.spec))
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_with_test_results_contrast_spec_with_cooks_replacement/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("lrt");

        let method = rt(
            b.fit_with_test_results_contrast_spec_with_cooks_replacement(
                &f.counts, &f.design, &f.spec, &f.cooks,
            ),
        );
        let run = run_rt(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .result(ResultColumn::Spec(&f.spec))
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_with_test_results_contrast_spec_with_cooks_replacement/lrt",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");
        let rc = f.character();
        let method = rt(
            b.fit_with_test_results_contrast_request_with_cooks_replacement(
                &f.counts,
                &f.design,
                &rc,
                Some(f.levels.as_slice()),
                &f.cooks,
            ),
        );
        let run = run_rt(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .result(ResultColumn::Results {
                        contrast: &rc,
                        sample_levels: Some(f.levels.as_slice()),
                        model_frame: None,
                    })
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_with_test_results_contrast_request_with_cooks_replacement/wald/character/levels",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");
        let rc = f.character();
        let method = rt(
            b.fit_with_test_results_contrast_request_with_cooks_replacement(
                &f.counts,
                &f.design,
                &rc,
                None::<&[String]>,
                &f.cooks,
            ),
        );
        let run = run_rt(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .result(ResultColumn::Results {
                        contrast: &rc,
                        sample_levels: None,
                        model_frame: None,
                    })
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_with_test_results_contrast_request_with_cooks_replacement/wald/character/nolevels",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");
        let rc = f.numeric_request();
        let method = rt(
            b.fit_with_test_results_contrast_request_with_cooks_replacement(
                &f.counts,
                &f.design,
                &rc,
                Some(f.levels.as_slice()),
                &f.cooks,
            ),
        );
        let run = run_rt(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .result(ResultColumn::Results {
                        contrast: &rc,
                        sample_levels: Some(f.levels.as_slice()),
                        model_frame: None,
                    })
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_with_test_results_contrast_request_with_cooks_replacement/wald/numeric/levels",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");
        let rc = f.numeric_request();
        let method = rt(
            b.fit_with_test_results_contrast_request_with_cooks_replacement(
                &f.counts,
                &f.design,
                &rc,
                None::<&[String]>,
                &f.cooks,
            ),
        );
        let run = run_rt(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .result(ResultColumn::Results {
                        contrast: &rc,
                        sample_levels: None,
                        model_frame: None,
                    })
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_with_test_results_contrast_request_with_cooks_replacement/wald/numeric/nolevels",
            method,
            run,
        );
    }
    {
        let b = f.builder("lrt");
        let rc = f.character();
        let method = rt(
            b.fit_with_test_results_contrast_request_with_cooks_replacement(
                &f.counts,
                &f.design,
                &rc,
                Some(f.levels.as_slice()),
                &f.cooks,
            ),
        );
        let run = run_rt(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .result(ResultColumn::Results {
                        contrast: &rc,
                        sample_levels: Some(f.levels.as_slice()),
                        model_frame: None,
                    })
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_with_test_results_contrast_request_with_cooks_replacement/lrt/character/levels",
            method,
            run,
        );
    }
    {
        let b = f.builder("lrt");
        let rc = f.character();
        let method = rt(
            b.fit_with_test_results_contrast_request_with_cooks_replacement(
                &f.counts,
                &f.design,
                &rc,
                None::<&[String]>,
                &f.cooks,
            ),
        );
        let run = run_rt(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .result(ResultColumn::Results {
                        contrast: &rc,
                        sample_levels: None,
                        model_frame: None,
                    })
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_with_test_results_contrast_request_with_cooks_replacement/lrt/character/nolevels",
            method,
            run,
        );
    }
    {
        let b = f.builder("lrt");
        let rc = f.numeric_request();
        let method = rt(
            b.fit_with_test_results_contrast_request_with_cooks_replacement(
                &f.counts,
                &f.design,
                &rc,
                Some(f.levels.as_slice()),
                &f.cooks,
            ),
        );
        let run = run_rt(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .result(ResultColumn::Results {
                        contrast: &rc,
                        sample_levels: Some(f.levels.as_slice()),
                        model_frame: None,
                    })
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_with_test_results_contrast_request_with_cooks_replacement/lrt/numeric/levels",
            method,
            run,
        );
    }
    {
        let b = f.builder("lrt");
        let rc = f.numeric_request();
        let method = rt(
            b.fit_with_test_results_contrast_request_with_cooks_replacement(
                &f.counts,
                &f.design,
                &rc,
                None::<&[String]>,
                &f.cooks,
            ),
        );
        let run = run_rt(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .result(ResultColumn::Results {
                        contrast: &rc,
                        sample_levels: None,
                        model_frame: None,
                    })
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_with_test_results_contrast_request_with_cooks_replacement/lrt/numeric/nolevels",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");
        let rc = f.character();
        let method = rt(
            b.fit_formula_with_test_results_contrast_request_with_cooks_replacement(
                &f.counts,
                f.full_formula,
                &rc,
                &f.cooks,
            ),
        );
        let run = run_rt(
            b.run(
                &Analysis::formula(&f.counts, f.full_formula)
                    .result(ResultColumn::Results {
                        contrast: &rc,
                        sample_levels: None,
                        model_frame: None,
                    })
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_formula_with_test_results_contrast_request_with_cooks_replacement/wald/character",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");
        let rc = f.numeric_request();
        let method = rt(
            b.fit_formula_with_test_results_contrast_request_with_cooks_replacement(
                &f.counts,
                f.full_formula,
                &rc,
                &f.cooks,
            ),
        );
        let run = run_rt(
            b.run(
                &Analysis::formula(&f.counts, f.full_formula)
                    .result(ResultColumn::Results {
                        contrast: &rc,
                        sample_levels: None,
                        model_frame: None,
                    })
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_formula_with_test_results_contrast_request_with_cooks_replacement/wald/numeric",
            method,
            run,
        );
    }
    {
        let b = f.builder("lrt");
        let rc = f.character();
        let method = rt(
            b.fit_formula_with_test_results_contrast_request_with_cooks_replacement(
                &f.counts,
                f.full_formula,
                &rc,
                &f.cooks,
            ),
        );
        let run = run_rt(
            b.run(
                &Analysis::formula(&f.counts, f.full_formula)
                    .lrt(&f.reduced)
                    .result(ResultColumn::Results {
                        contrast: &rc,
                        sample_levels: None,
                        model_frame: None,
                    })
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_formula_with_test_results_contrast_request_with_cooks_replacement/lrt/character",
            method,
            run,
        );
    }
    {
        let b = f.builder("lrt");
        let rc = f.numeric_request();
        let method = rt(
            b.fit_formula_with_test_results_contrast_request_with_cooks_replacement(
                &f.counts,
                f.full_formula,
                &rc,
                &f.cooks,
            ),
        );
        let run = run_rt(
            b.run(
                &Analysis::formula(&f.counts, f.full_formula)
                    .lrt(&f.reduced)
                    .result(ResultColumn::Results {
                        contrast: &rc,
                        sample_levels: None,
                        model_frame: None,
                    })
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_formula_with_test_results_contrast_request_with_cooks_replacement/lrt/numeric",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");
        let rc = f.character();
        let method = rt(b
            .fit_with_test_results_contrast_request_from_model_frame_with_cooks_replacement(
                &f.counts, &f.design, &rc, &f.frame, &f.cooks,
            ));
        let run = run_rt(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .result(ResultColumn::Results {
                        contrast: &rc,
                        sample_levels: None,
                        model_frame: Some(&f.frame),
                    })
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_with_test_results_contrast_request_from_model_frame_with_cooks_replacement/wald/character",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");
        let rc = f.numeric_request();
        let method = rt(b
            .fit_with_test_results_contrast_request_from_model_frame_with_cooks_replacement(
                &f.counts, &f.design, &rc, &f.frame, &f.cooks,
            ));
        let run = run_rt(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .result(ResultColumn::Results {
                        contrast: &rc,
                        sample_levels: None,
                        model_frame: Some(&f.frame),
                    })
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_with_test_results_contrast_request_from_model_frame_with_cooks_replacement/wald/numeric",
            method,
            run,
        );
    }
    {
        let b = f.builder("lrt");
        let rc = f.character();
        let method = rt(b
            .fit_with_test_results_contrast_request_from_model_frame_with_cooks_replacement(
                &f.counts, &f.design, &rc, &f.frame, &f.cooks,
            ));
        let run = run_rt(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .result(ResultColumn::Results {
                        contrast: &rc,
                        sample_levels: None,
                        model_frame: Some(&f.frame),
                    })
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_with_test_results_contrast_request_from_model_frame_with_cooks_replacement/lrt/character",
            method,
            run,
        );
    }
    {
        let b = f.builder("lrt");
        let rc = f.numeric_request();
        let method = rt(b
            .fit_with_test_results_contrast_request_from_model_frame_with_cooks_replacement(
                &f.counts, &f.design, &rc, &f.frame, &f.cooks,
            ));
        let run = run_rt(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .result(ResultColumn::Results {
                        contrast: &rc,
                        sample_levels: None,
                        model_frame: Some(&f.frame),
                    })
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_with_test_results_contrast_request_from_model_frame_with_cooks_replacement/lrt/numeric",
            method,
            run,
        );
    }
    {
        let b = f.builder("bare");
        let rc = f.character();
        let method = rt(b
            .fit_with_test_results_contrast_request_from_model_frame_with_cooks_replacement(
                &f.counts, &f.design, &rc, &f.frame, &f.cooks,
            ));
        let run = run_rt(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .result(ResultColumn::Results {
                        contrast: &rc,
                        sample_levels: None,
                        model_frame: Some(&f.frame),
                    })
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_with_test_results_contrast_request_from_model_frame_with_cooks_replacement/bare/character",
            method,
            run,
        );
    }
    {
        let b = f.builder("bare");
        let rc = f.numeric_request();
        let method = rt(b
            .fit_with_test_results_contrast_request_from_model_frame_with_cooks_replacement(
                &f.counts, &f.design, &rc, &f.frame, &f.cooks,
            ));
        let run = run_rt(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .result(ResultColumn::Results {
                        contrast: &rc,
                        sample_levels: None,
                        model_frame: Some(&f.frame),
                    })
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_with_test_results_contrast_request_from_model_frame_with_cooks_replacement/bare/numeric",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method = rw(
            b.fit_with_results_factor_level_contrast_with_cooks_replacement(
                &f.counts,
                &f.design,
                f.factor_level(),
                &f.cooks,
            ),
        );
        let run = run_rw(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .result(ResultColumn::FactorLevel(f.factor_level()))
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_with_results_factor_level_contrast_with_cooks_replacement/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("wald");

        let method = rt(
            b.fit_with_test_results_factor_level_contrast_with_cooks_replacement(
                &f.counts,
                &f.design,
                f.factor_level(),
                &f.cooks,
            ),
        );
        let run = run_rt(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .result(ResultColumn::FactorLevel(f.factor_level()))
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_with_test_results_factor_level_contrast_with_cooks_replacement/wald",
            method,
            run,
        );
    }
    {
        let b = f.builder("lrt");

        let method = rt(
            b.fit_with_test_results_factor_level_contrast_with_cooks_replacement(
                &f.counts,
                &f.design,
                f.factor_level(),
                &f.cooks,
            ),
        );
        let run = run_rt(
            b.run(
                &Analysis::new(&f.counts, &f.design)
                    .lrt(&f.reduced)
                    .result(ResultColumn::FactorLevel(f.factor_level()))
                    .cooks_replacement(&f.cooks),
            ),
        );
        check(
            &mut report,
            "fit_with_test_results_factor_level_contrast_with_cooks_replacement/lrt",
            method,
            run,
        );
    }

    report
}
