# Changelog

Public API and user-visible behaviour changes. Release notes for published
versions live in [docs/releases](docs/releases/0.5.0.md).

## Unreleased

### Added

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

### Deprecated

- `rsdeseq2::math::qr`, renamed to `rsdeseq2::math::rank`: the module computes
  matrix ranks and contains no QR factorization. The old path re-exports the
  same items for one release.
