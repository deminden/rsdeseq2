# Beta-prior Wald analyses

## Formula workflow

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

Independent R 4.6.2 / DESeq2 1.52.0 fixtures cover `~ condition` on
`three_vs_three` and `~ batch + condition` on the 500-sample design.
[The beta-prior tests](../../crates/rsdeseq2/tests/beta_prior_run.rs)
check prior variances, full result fields, missingness and significance.

## Expanded-matrix fits

The `fit_expanded_*_beta_prior_wald_*` primitives and the CLI's
`--beta-prior-expanded-design` path take an explicit expanded design.
They estimate prior variances from that design's MLE fit and average
expanded coefficients when reporting a collapsed column. Use the formula
workflow above for DESeq2's standard-design prior estimation and expanded
contrasts.
