//! Validation and per-run configuration of the analysis entry point.
use rsdeseq2::prelude::*;

struct Fixture {
    counts: CountMatrix,
    design: DesignMatrix,
    reduced: DesignMatrix,
    frame: FormulaModelFrame,
    levels: Vec<String>,
    dispersions: Vec<f64>,
    cooks: CooksReplacementOptions,
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
            cooks: CooksReplacementOptions::new(cutoff),
        }
    }
}

#[test]
fn run_rejects_invalid_combinations_before_work() {
    use rsdeseq2::errors::RunConflict;
    let f = Fixture::new();
    let b = DeseqBuilder::new().model_frame(f.frame.clone());
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
            base().lrt_formula("~ batch"),
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
    let character = ResultsContrast::character("condition", "B", "A");
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
    let b = DeseqBuilder::new().model_frame(f.frame.clone());
    for fit_type in [FitType::Mean, FitType::Local] {
        let typed = b.clone().fit_type(fit_type);
        for stage in [Stage::Trend, Stage::Map, Stage::Results] {
            let analysis = Analysis::new(&f.counts, &f.design).stop_after(stage);
            assert_eq!(
                format!("{:?}", b.run(&analysis.fit_type(fit_type)).unwrap()),
                format!("{:?}", typed.run(&analysis).unwrap()),
                "{fit_type:?} {stage:?}"
            );
        }
    }
}
