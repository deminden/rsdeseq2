# Analysis entry point

Use `DeseqBuilder::run(&Analysis)` for an end-to-end fit or an early stage.
The Rust CLI and R bridge use the same workflow. The builder configures
threading, size-factor estimation and other fitting settings. Each `Analysis`
selects the design, test and outputs for one run.

## Request and output

| Setting | Current API |
| --- | --- |
| Design | `Analysis::new(counts, design)` for a numeric matrix, `Analysis::formula(counts, formula)` for the builder's model frame, or `Analysis::without_design(counts)` for size factors |
| Test and dispersion | Wald with estimated dispersion by default. Use `lrt` or `lrt_formula` to supply a reduced design, and `fixed_dispersion` to supply dispersions |
| Mean model and trend | `means(MeanModel)` selects automatic, GLM or linear means. `fit_type(FitType)` overrides the builder's trend for one run |
| Result column | `result(ResultColumn)` selects a coefficient or contrast |
| Other output controls | `cooks_replacement` requests replacement and refitting, `stop_after(Stage)` selects an early stage, and `beta_prior(true)` selects the supported formula beta-prior Wald fit |

`run` validates the request before fitting. Its `RunOutput` distinguishes
fit-only, result, replacement and beta-prior output. Use `fit()` to borrow
the fitted state or `into_fit_and_results()` to consume the output and
assemble the final fit and result, including merged replacement rows.

Use `Analysis::lrt(reduced)` or `Analysis::lrt_formula(reduced)`
to supply the reduced design as part of the request.

## Preferred public API

These are the preferred public paths for complete analyses, fitted-model
queries and individual DESeq2 stages in 0.6.1.

| Area | Public paths |
| --- | --- |
| Workflow | `core::DeseqBuilder::run`, `Analysis`, its matrix-design selectors, `RunOutput`, `DeseqFit` and `Stage` |
| Stored fits | `results::ResultRequest`, `DeseqFit::results_for_request` and the replacement-aware query with its refit plan |
| Data and policy | Count, design and row-major matrices, labeled readers in `io`, observation weights, `Offsets`, contrasts, result rows, Cook's replacement settings and fit diagnostics |
| Standalone stages | Size-factor options and estimation, gene-wise/trend/MAP dispersion estimation, `DispersionFitInput` for lower-level search, `GlmInput` and the canonical fixed-dispersion fits, `vst_with_dispersion_trend_and_offsets` and the offset-based rlog functions |
| Configuration | `DeseqError`, `ErrorKind`, runtime `NumericProfile` and `DeseqBuilder::threads` |

The prelude exports the preferred workflow types. Formula evaluation and
ashr and apeglm shrinkage are experimental. `PosteriorWorkControl` provides
cooperative cancellation for shrinkage.
