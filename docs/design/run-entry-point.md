# Design note 4.2 — one `DeseqBuilder::run` entry point

Status: **shape approved** (2026-11-17). One `run`, as the plan's target table
says: "1 `run` plus deprecated wrappers, then 1"; stopping early is part of
`run`. Implementation follows the migration steps below.

## Problem

`DeseqBuilder` exposes 143 `pub fn fit*` methods. Each method is one point
in a product of independent choices, spelled out in its name:

| Axis | Values in method names | Methods using it |
| --- | --- | ---: |
| Test | `wald`, `lrt` (reduced design as an extra argument) | 32 / 59 |
| Mean model for dispersion | automatic, `glm_mu`, `linear_mu` | 33 / 23 |
| Dispersion | estimated, `fixed_dispersion` (caller-supplied vector) | 24 fixed |
| Trend | builder `fit_type`, or `_parametric` suffix | 20 |
| Design input | matrix, `formula` (string plus model frame) | 18 formula |
| Result selection | default last coefficient, `_name`, `_contrast` (numeric), `_contrast_spec`, `_factor_level_contrast`, `_results_contrast`, `_contrast_request`, `_from_model_frame` | 83 |
| Cook's replacement | off, `_with_cooks_replacement` | 47 |
| Output | fit only, `_with_results`, `_with_test_results` | — |
| Stage | full pipeline, or stop after size factors, gene-wise, trend or MAP | ~10 |

Adding one axis value, for example a new contrast form, doubles the methods
on some axes. The R bridge and CLI call 25 distinct methods, and docs and
tests cannot cover every combination.

## Proposal

Add one entry point that takes a request value. Keep the existing methods as
thin wrappers around it, first deprecated and then removed.

```rust
pub struct Analysis<'a> {
    pub counts: &'a CountMatrix,
    pub design: DesignInput<'a>,          // Matrix(&DesignMatrix) | Formula { formula: &str, frame: &FormulaModelFrame }
    pub test: TestRequest<'a>,            // Wald | Lrt { reduced: DesignInput<'a> }
    pub dispersion: DispersionRequest<'a>,// Estimate | Fixed(&[f64])
    pub means: MeanModel,                 // Auto (default) | GlmMu | LinearMu
    pub result: ResultRequest<'a>,        // FitOnly | Default | Coefficient(usize) | Name(&str)
                                          // | Numeric(&[f64]) | Spec(&ContrastSpec)
                                          // | Results { contrast: &ResultsContrast, sample_levels: Option<&[&str]> }
    pub replacement: Replacement<'a>,     // Off | Cooks(&CooksReplacementOptions)
    pub stop_after: Stage,                // Results (default) | SizeFactors | GeneWise | Trend | Map
}

impl DeseqBuilder {
    pub fn run(&self, analysis: &Analysis<'_>) -> Result<RunOutput, DeseqError>;
}

#[non_exhaustive]
pub struct RunOutput {
    pub fit: DeseqFit,                    // merged fit when replacement refit rows
    pub results: Option<DeseqResults>,    // None for FitOnly or an early Stage
    pub replacement: Option<ReplacementOutput>, // original fit and results, refit plan and refit
}
```

- The `_parametric` suffix becomes `fit_type(FitType::Parametric)` on the
  builder, which already exists. Options stay on the builder: thread count,
  Cook's cutoff, filtering, IRLS and size-factor options.
- `Analysis::new(counts, design)` fills the defaults. These match today's
  `fit_with_results`: Wald, estimated dispersion, automatic means, the default
  coefficient, no replacement. Each axis then has a `with_*` setter.
- Every combination is validated once, in `run`. For example, `Fixed`
  dispersion with `Stage::Trend` becomes an `Options` error.

## Migration (one PR per step, IDENTICAL unless noted)

1. **Add `run`.** Implement it by dispatching to the existing private
   workflow functions, not the public wrappers. Golden harness, plus a new
   test that runs each of the 143 methods and its `run` equivalent on the toy
   inputs and compares the outputs bit for bit.
2. **Rewrite every public `fit*` as a one-line call to `run`** and mark them
   `#[deprecated(since = "0.6.0", note = "use DeseqBuilder::run with …")]`.
   Each note names the exact `Analysis` that replaces the method.
3. **Move the R bridge and CLI to `run`.** This removes 25 call sites and
   leaves no deprecated calls inside the workspace.
4. **Release 0.6.0** with the deprecations, and remove the wrappers in 0.7.0.

The metric `builder_fit_methods` goes from 143 to 1 plus the deprecated
wrappers, then to 1, the target in `hardening-progress.md`.

## Risks

- **Semver.** Step 2 adds deprecation warnings for downstream users. The R
  package is in-tree and moves in step 3. No other dependents are known.
- **Hidden differences between near-duplicate methods.** Some pairs may
  differ subtly, for example which contrast forms consult the model frame.
  The step-1 equivalence test exists to find these before anything is
  rewritten. Any real difference becomes an explicit `ResultRequest` variant,
  not a silent unification.
- **Diff size.** The wrappers are mechanical (about 1,500 lines across
  `builder/`). Steps 2 and 3 are large but contain no logic changes.

## Implementation notes (step 1, 2026-11-18)

Three small changes from the sketch above, found while writing step 1:

- `ResultRequest` is named `ResultColumn`, because `rsdeseq2::results`
  already exports a `ResultRequest`. The early-stop "fit only" form needs no
  variant: methods that return only a fit still choose a column, and the fit
  records that column's test, so their wrappers take `run(..)?.fit()`.
- `RunOutput` is an enum (`Fit`, `Results`, `Replacement`) rather than a
  struct. With replacement, the merged fit is built by consuming the original
  fit, so a struct holding both would need a full copy of the fit.
  `RunOutput::into_fit_and_results` merges on request.
- `DesignInput::Formula` takes the formula only and uses the builder's stored
  model frame, as the formula methods do today; a different frame for one run
  is `ResultColumn::Results { model_frame, .. }` for contrast resolution. The
  `_parametric` suffix is `Analysis::fit_type`, an optional per-run override
  of the builder's trend type, so wrappers do not clone the builder to change
  it. `DesignInput::None` covers the design-free size-factor stage.

Step 1 dispatches to the existing leaf methods. The generated test
`tests/run_equivalence.rs` (from `tools/validation/gen_run_equivalence.py`)
checks 208 method/request pairs, covering all 140 `fit*` methods except the
two fast-VST trend helpers, which stay outside `run` and will be renamed in
4.3. All 208 pairs are identical, and none fails on both sides.

## Decisions

1. One `run(&Analysis)` returning `RunOutput`, not separate `run_wald` and
   `run_lrt`.
2. Deprecated in 0.6.0 and all removed in 0.7.0, leaving one entry point.
3. Stopping early is `Analysis::stop_after: Stage`; the staged methods become
   deprecated wrappers too.
