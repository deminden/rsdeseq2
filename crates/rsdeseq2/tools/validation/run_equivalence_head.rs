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
            other => Err(DeseqError::InvalidOptions { reason: format!("unexpected test: {:?}", other.test_type()) }),
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
            other => Err(DeseqError::InvalidOptions { reason: format!("unexpected test: {:?}", other.test_type()) }),
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
            assert!(errors.is_empty(), "{} pairs failed on both sides", errors.len());
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
