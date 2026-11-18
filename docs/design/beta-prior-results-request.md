# Design note 4.3d — one beta-prior Wald results entry point

Status: **beta prior on `run` implemented** (2026-11-18). Denis chose the
`run` option ("yes, if original DESeq offers it": DESeq2 has
`DESeq(betaPrior = TRUE)`). The request-type collapse below is not done.

## Problem

`rsdeseq2::results` has 40 public `fit_expanded_*_beta_prior_wald_*`
functions, the `betaPrior = TRUE` Wald path. Each is one point in a product
of four independent choices:

| Axis | Values | Today |
| --- | --- | --- |
| Design source | expanded matrices, one factor, additive factors (+ numeric, interactions), formula + specs, formula + model frame | 5 name stems, 10 input structs |
| Offsets | size factors, normalization factors (+ weights either way) | `_with_normalization_factors_and_weights` suffix and a second input struct per stem |
| Reported column | coefficient, numeric contrast | `_contrast` infix |
| Cook's replacement | off, on | `_and_cooks_replacement` / `_with_cooks_replacement` suffix |

5 × 2 × 2 × 2 = 40 functions and 10 near-identical input structs (every
`*NormalizedResultsInput` differs from its twin only in
`normalization_factors` versus `size_factors`). The CLI calls them through
one dispatcher in `cli/implementation/analysis.rs`; the R bridge does not
use them.

## Proposal

```rust
pub enum BetaPriorDesign<'a> {
    Expanded(ExpandedModelBetaPriorDesignInput<'a>),
    Factor { factor: &'a str, sample_levels: &'a [String], reference: &'a str },
    Additive { factors, numeric_covariates, interactions,
               factor_numeric_interactions, numeric_interactions },
    Formula { formula: &'a str, factors, numeric_covariates },
    FormulaModelFrame { formula: &'a str, model_frame: &'a FormulaModelFrame },
}

pub enum BetaPriorColumn<'a> { Coefficient(usize), Contrast(&'a [f64]) }

#[non_exhaustive]
pub struct BetaPriorWaldRequest<'a> {
    pub glm: GlmInput<'a>,             // counts, Offsets, dispersions, weights (exists)
    pub design: BetaPriorDesign<'a>,
    pub base_mean: &'a [f64],
    pub disp_fit: &'a [f64],
    pub gene_names: Option<&'a [String]>,
    pub column: BetaPriorColumn<'a>,
    pub replacement: Option<&'a CooksReplacementOptions>,
    pub options: BetaPriorRefitOptions,
}

pub fn beta_prior_wald_results(
    request: &BetaPriorWaldRequest<'_>,
) -> Result<BetaPriorWaldOutput, DeseqError>;
```

`BetaPriorWaldOutput` is an enum `{ Results(ExpandedBetaPriorWaldResults),
Replacement(ExpandedBetaPriorWaldReplacementResults) }`, the same split as
`RunOutput`; the five per-stem result structs (`ExpandedFactor…Results`,
`ExpandedAdditive…Results`, …) differ only in which design they echo back, so
the output would carry the resolved `ExpandedModelBetaPriorDesignInput`
data once.

## Migration (as for `run`)

1. Add the request, output and function, dispatching to the existing bodies;
   a generated test runs all 40 functions and their requests and compares the
   outputs bit for bit.
2. Make the 40 functions one-line wrappers and deprecate them with the
   10 input structs (0.6.0), removal in 0.7.0.
3. Move the CLI dispatcher to the request.

## Questions for Denis

1. Is a single `beta_prior_wald_results(&BetaPriorWaldRequest)` the right
   shape, or should the beta prior become an `Analysis` option on
   `DeseqBuilder::run` instead (DESeq2's `DESeq(betaPrior = TRUE)`), with
   these functions kept only as low-level primitives?
2. Keep one output enum, or keep the five per-design result structs?

## Implemented: `Analysis::beta_prior`

`Analysis::formula(&counts, formula).beta_prior(true)` runs DESeq2's
`DESeq(betaPrior = TRUE)` and returns `RunOutput::BetaPrior`
(`BetaPriorRunOutput`: fit, results, expanded design, `betaPriorVar`,
expanded MAP fit and contrast). It follows `nbinomWaldTest(betaPrior = TRUE)`
with the expanded model matrix (`fitGLMsWithPrior`, `estimateBetaPriorVar`,
`addAllContrasts`, `averagePriorsOverLevels`); see
`core/implementation/builder/beta_prior.rs`.

- Prior variances come from the MLE fold changes of the *standard* design,
  plus every pairwise difference between a factor's non-reference levels,
  averaged over each factor and assigned to all its expanded level columns.
- Cook's distances use the standard-design MLE fit's means and hat
  diagonals, as in DESeq2.
- The reported column is the expanded contrast equal to the requested
  standard column (`condition_B_vs_A` is `conditionB - conditionA`).
- Not supported, with an error: matrix designs, LRT, fixed dispersions,
  early stages, interaction terms (DESeq2 rejects these too), Cook's
  replacement, and sample-level factor contrasts.

Against R 4.6.2 / DESeq2 1.52.0 `DESeq(betaPrior = TRUE,
minReplicatesForReplace = Inf)` on `three_vs_three` (`~ condition`) and
`many_samples` (`~ batch + condition`, 500 samples): `betaPriorVar` within
1e-12 relative, every result field within 7e-12 relative, identical
missingness, 0 flips at padj < 0.1 (`tests/beta_prior_run.rs`).

## Finding: the expanded primitives do not reproduce `DESeq(betaPrior = TRUE)`

The 40 `fit_expanded_*_beta_prior_wald_*` functions (and the CLI's
`--beta-prior-expanded-design` path) estimate the prior variances from an
MLE fit on the *expanded* design and report coefficients "collapsed" onto the
standard design by averaging expanded columns. DESeq2 estimates the
variances from the standard-design MLE (above) and reports expanded
contrasts. On `three_vs_three` the primitives give a prior variance of
2.97/3.13 where DESeq2 gives 0.53, and log2 fold changes 40 % smaller than
DESeq2's. They have no end-to-end DESeq2 reference. Options: deprecate them in
favour of `run`, or re-base them on the `run` implementation; this needs
Denis's decision.
