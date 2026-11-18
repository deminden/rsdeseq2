# Changelog

Public API and user-visible behaviour changes. Release notes for published
versions live in [docs/releases](docs/releases/0.5.0.md).

## Unreleased

### Added

- `local_vst_inverse_factor_mean(Offsets)` and
  `vst_with_dispersion_trend_and_offsets`.

- `rlog_transform`, `rlog_fit`, `rlog_frozen` and `rlog_with_estimated_prior`
  take `Offsets` (size factors or normalization factors).

- `fit_fixed_dispersion_glm(GlmInput, design, options)` and
  `fit_intercept_only_glm(GlmInput)`: the fixed-dispersion IRLS fit and the
  intercept-only closed form with size-factor or normalization-factor
  `Offsets` and optional weights.

- `Offsets` (per-sample size factors or a per-cell normalization-factor
  matrix) is now public, with `GlmInput` (counts, offsets, dispersions,
  optional weights) and three beta-prior fits that take it:
  `fit_glms_with_beta_prior`, `fit_glms_with_estimated_beta_prior` and
  `fit_expanded_glms_with_estimated_beta_prior`.

- `DispersionFitInput` with `line_search_dispersion` and
  `grid_search_dispersion`: one input struct (counts, means, optional design,
  prior and weights) for a single gene's dispersion fit, in place of twelve
  argument-combination functions.

- `DeseqBuilder::run(&Analysis) -> RunOutput`: one entry point for the
  whole pipeline. `Analysis` holds the design (matrix, formula, or none for
  size factors only), the test, fixed or estimated dispersions, the mean
  model, the reported column (`ResultColumn`), Cook's replacement, an early
  `Stage` and an optional trend type. Invalid combinations are rejected before
  any work. Every existing `fit*` method gives bit-identical output to its
  `run` request; those methods will be deprecated in 0.6.0.

- End-to-end parity test against `DESeq()` + `results()` defaults on a
  3,745-gene GTEx kidney-cortex panel (size factors, parametric dispersion
  trend, Cook's replacement, independent filtering). 3,744 genes match R to
  about 1e-13; one very high-count gene (GPX3) ends its MAP line search at a
  different point because the Armijo test there compares log posteriors that
  agree to the last bit. Its dispersion differs by 2.6e-5 on the log scale and
  its p-value by 3.7e-6; no significance call changes.

- `NumericProfile`, `numeric_profile()` and `set_numeric_profile()`, with
  three profiles chosen at run time:
  - `portable` (the library default): pure-Rust elementary and special
    functions and fixed-order reductions, so the same input gives the same
    bits on every OS and architecture. One golden manifest covers all targets.
  - `r-linux-x86_64`: R on Linux x86_64 (x87 extended reductions, the OpenBLAS
    Haswell dot schedule, glibc libm); available on any x86_64 host.
  - `r-macos`: R on macOS (host libm); available on macOS only.
  `RSDESEQ2_NUMERIC_PROFILE` selects one by name, and `r-host` picks the R
  profile for the current host. The R package selects `r-host` unless the
  variable is set, so its output is unchanged. Library users who relied on the
  old host-dependent default should set `RSDESEQ2_NUMERIC_PROFILE=r-host`;
  On full GTEx (56,200 genes, two contrasts, parametric and local trends)
  `portable` differs from `r-linux-x86_64` by at most 1.9e-4 relative in
  p-values and 1.2e-4 in dispersions, with no significance call changing at
  padj 0.01, 0.05 or 0.1; VST values differ by at most 1e-15 relative.

- `_from_reader` forms of the TSV readers that open their own file:
  `read_count_matrix_tsv`, `read_labeled_design_matrix_tsv`,
  `read_labeled_size_factors_tsv`, `read_labeled_geometric_means_tsv`,
  `read_labeled_sample_numeric_tsv`, `read_labeled_gene_numeric_tsv`,
  `read_labeled_wald_t_degrees_of_freedom_tsv` and `read_sample_levels_tsv`.
  They accept any `std::io::Read` (for example a byte slice) and share the
  parsing code of the path forms. `read_count_matrix_tsv_from_reader` is in
  the prelude.

- Cargo features `cli` (the binary and `rsdeseq2::cli`; pulls in clap, sha2 and
  same-file), `formula` (ICU locale collation for formula factor levels and
  string comparisons) and `shrinkage` (`rsdeseq2::shrinkage`). All three are
  on by default, so existing users see no change. Without `formula`, formula
  evaluation accepts only `"C"` collation and reports `UnsupportedFeature`
  otherwise. The R package builds with `formula` and `shrinkage` only.

- `ErrorKind` and `DeseqError::kind()`: a stable coarse classification
  (shape, counts, normalization, dispersion, options, numeric, unsupported,
  I/O, parse) so callers can branch without matching message text.

### Removed

- `rsdeseq2::bindings` and its `language_bindings()` placeholder, which only
  returned `UnsupportedFeature`.
- `FitType::GlmGamPoi`. It was never implemented (every path returned
  `UnsupportedFeature`) and neither the CLI nor the R package could select it.
- `ExecutionMode`, `DeseqBuilder::execution_mode` and
  `DeseqBuilder::current_execution_mode`, and the R `DESeq(parityMode = )` /
  `fitRust(parityMode = )` argument. Only the strict mode ever existed in the
  numerical code; "fast" and "diagnostic" were accepted and ignored. Passing
  `parityMode` to the R functions is now an unused-argument error.

### Changed

- Errors are more structured; every message is unchanged. Out-of-domain
  parameters ("*x* must be finite and positive", and the other phrases in
  `Requirement`) are `DeseqError::InvalidParameter { domain, parameter,
  requirement }` (92 sites), formula parse and evaluation errors are
  `InvalidFormula` (203 sites), missing coefficients are
  `CoefficientNotFound { candidates, .. }`, and internal invariant failures
  are `Internal`. They were `InvalidOptions`, `InvalidDispersion`,
  `InvalidSizeFactors` or `InvalidCounts` before, so code that matches those
  variants for these cases should match the new variants or `kind()` instead.
  `kind()` is unchanged except for formula errors (`ErrorKind::Formula`) and
  invariants (`ErrorKind::Internal`). Design note:
  `docs/design/structured-errors.md`.

- The numeric foundation (`errors`, `math`, `matrix`, `glm::nb`, numeric
  profiles, special functions and distributions) moved into a new
  `rsdeseq2-base` crate, formula evaluation and design matrices
  (`rsdeseq2::design`) into `rsdeseq2-formula`, and the empirical-Bayes
  shrinkage engines (`rsdeseq2::shrinkage`) into `rsdeseq2-shrink`.
  `rsdeseq2` depends on them and re-exports every public path unchanged, so no
  user code changes, and its `formula` and `shrinkage` features work as
  before. Releases publish the four crates together, and the R package
  vendors all four.

- rlog now fits its sample-effect GLM with the settings DESeq2's `rlogData`
  passes to `fitNbinomGLMs` (`betaTol = 1e-4`, no optim fallback), available as
  `IrlsOptions::deseq2_rlog()`, and rank-deficient designs (every rlog design)
  start IRLS as `fitNbinomGLMs` does (NUMERIC). On pasilla_holdout the largest
  difference from `DESeq2::rlog(blind = TRUE)` drops from 1.06e-2 to 2.5e-14.

- `DesignMatrix::rank`, `is_full_rank` and `validate_full_rank` (used by every
  fit's full-rank check) now follow R's default `qr(x)$rank`, the rule DESeq2's
  `checkFullRank` uses: LINPACK `dqrdc2` with tolerance 1e-7, exposed as
  `math::rank::matrix_rank_r_qr` and `R_QR_TOLERANCE` (NUMERIC). The former
  column-relative elimination at 1e-10 accepted near-collinear designs that
  DESeq2 rejects: 8 of 24 new R 4.6.2 reference cases, plus 2 cases of the
  existing `rank_r462.tsv`. `rank_with_tolerance` is unchanged.

- Dispersion trend fitting follows DESeq2 1.52 `estimateDispersionsFit`
  (NUMERIC): a failed parametric trend falls back to the local trend instead of
  returning an error, and the local trend is fitted only on rows with
  `dispGeneEst > 100 * minDisp` (it previously also used rows from 10 to 100
  times `minDisp`). New `estimate_dispersion_trend(base_mean, disp_gene_est,
  fit_type, min_disp)` exposes this dispatch; `fit_dispersion_trend` now uses
  it. `DeseqFit::dispersion_trend` records the trend actually fitted. The
  golden harness changed only for the toy 4x4 parametric cases, which
  previously errored and now match the local fit bit for bit.

- `NonFiniteValue` errors for a finite value that is out of range (for
  example a negative dispersion) now read "out-of-range value in ..." instead
  of "non-finite value in ...". The variant, its fields and its `ErrorKind`
  are unchanged.

- Error variants that named the wrong input: `bh_adjust_f64` now reports
  invalid p-values as `InvalidOptions` (was `InvalidSizeFactors`);
  `math::median` reports empty input as `InvalidDimensions` with context
  `median input` (was `InvalidSizeFactors`); an impossible LRT chi-square
  distribution is `InvalidOptions` (was `InvalidCounts`). Messages are otherwise
  unchanged.
- `NonFiniteValue` errors report the rejected value instead of `NaN` in 34
  places (normalization, base means and variances, IRLS working weights,
  responses and linear predictors, beta-prior helpers, contrast all-zero
  checks, independent filtering and local VST). Numerical results are
  unchanged.

- Release builds of the workspace use thin LTO and one codegen unit. Outputs
  are bit-identical; the CLI binary is about 16% smaller. The R package build
  keeps its own profile and unwinding panics.

- `DeseqBuilder`, `DeseqFit`, `DeseqResults`, `NbinomGlmFit`, `RlogOutput` and
  the Cook's replacement outputs are `#[must_use]`, so discarding a fit or
  result now warns. Impossible internal states return `InvalidOptions`
  ("internal invariant violated: ...") instead of panicking.

- `DeseqError` is `#[non_exhaustive]`; exhaustive matches need a wildcard arm.

### Fixed

- `poscounts` size factors in the pipeline (`DeseqBuilder`, `DeseqBuilder::run`,
  the CLI beta-prior paths) now match `DESeq(dds, sfType = "poscounts")`: the
  `DESeqDataSet` method passes its own non-zero geometric means and so
  rescales the size factors to a geometric mean of 1, which the builder did
  not do. `estimate_size_factors_for_dataset` exposes that rule; the
  matrix-level `estimate_size_factors_poscounts*` functions still match
  `estimateSizeFactorsForMatrix` and are unchanged.

- A formula containing a call to an unknown function, such as `~ q(age)` or
  `~ f() + condition`, recursed until the stack overflowed and aborted the
  process (CLI or R session). It now returns the usual "formula variable ...
  is not present" error. Found by the formula fuzz target.

### Deprecated

- `local_vst_inverse_size_factor_mean{,_from_normalization_factors}` and
  `vst_with_dispersion_trend_and_{size,normalization}_factors`: use the
  `Offsets` forms above. Output is unchanged; removal in 0.7.0.

- The eight `rlog_*_with_size_factors` / `rlog_*_with_normalization_factors`
  functions (including `rlog_with_estimated_prior_and_*`): use the four
  `Offsets` forms above. Output is unchanged; removal in 0.7.0.

- `fit_fixed_dispersion_irls{,_with_weights,_with_normalization_factors,_with_normalization_factors_and_weights}`
  and `fit_intercept_only_fixed_dispersion{,_with_weights,_with_normalization_factors}`:
  use `fit_fixed_dispersion_glm` or `fit_intercept_only_glm` with a
  `GlmInput`. Output is unchanged; removal in 0.7.0. `fit_irls` stays.

- The twelve `fit_glms_with_beta_prior_variance*`,
  `fit_glms_with_estimated_beta_prior_variance*` and
  `fit_expanded_glms_with_estimated_beta_prior_variance*` variants: use the
  three `GlmInput` functions above; the size-factor or normalization-factor
  and weights suffixes become `Offsets` and `GlmInput::weights`. Output is
  unchanged; removal in 0.7.0.

- The twelve `fit_dispersion_{line_search,grid}{,_with_prior,_with_prior_and_weights}`
  and `_no_cr` variants: use `line_search_dispersion` or
  `grid_search_dispersion` with a `DispersionFitInput`. `_no_cr` is a
  `DispersionFitInput` without a design. Output is unchanged; removal in 0.7.0.

- All 140 `DeseqBuilder::fit*` methods (for example `fit_with_results`,
  `fit_lrt_with_results`, `fit_wald_glm_mu_contrast_spec`,
  `fit_fixed_dispersion_wald_with_cooks_replacement`,
  `fit_map_dispersions_glm_mu`): use `DeseqBuilder::run` with an `Analysis`.
  Each deprecation note gives the exact replacement request. The methods are
  one-line wrappers with unchanged output and will be removed in 0.7.0. The
  fast-VST trend helpers are not affected.

- `rsdeseq2::math::qr`, renamed to `rsdeseq2::math::rank`: the module computes
  matrix ranks and contains no QR factorization. The old path re-exports the
  same items for one release.
