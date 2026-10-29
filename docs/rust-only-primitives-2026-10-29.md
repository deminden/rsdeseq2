# First standalone Rust primitive checkpoint

This is part of the [Rust-only wrapper migration](rust-only-wrapper-migration.md),
not completion of the whole migration. Filtering/BH and normalization primitives
now share standalone Rust calculations. Wald/LRT tails, contrast calculations,
formula-offset composition and optional shrinkage still need migration.

R no longer calculates filtering quantiles, lowess, residual moments, rejection
counts, threshold selection or BH. `evaluate_independent_filtering` serves both
the standalone result-table API and the native R bridge. There are no R numerical
callbacks or frontend-specific arithmetic modes. BH uses one ratio-then-multiply
order throughout the library, matching native R's adjustment order.

The bridge also now passes supplied geometric means directly to the existing
standalone estimator. It previously exponentiated natural-scale inputs, producing
incorrect factors on non-proportional rows. A reproduced example changed from
approximately `0.7663, 1.0217, 1.2772` to the independently checked
`0.5646, 1.8821, 0.9410`. Proportional controls had hidden the defect through final
scale stabilization. The private registered primitive now accepts natural-scale
geometric means, matching both public APIs; its direct-call test uses that input
encoding with the same numerical expectations and tolerance. Discarded R log
conversion and an unused R stabilization implementation were removed.

## Evidence

- A separate executable uses the standalone library without linking or launching
  R. `scripts/check_rust_r_identity.R` compares IEEE-754 bit patterns against the
  installed wrapper: **20 cases**, covering nine filtering cases and eleven
  normalization/metadata cases. These include custom/repeated/unsorted theta,
  tails, missingness, zeros, supplied means, zero/infinite/large means, named
  controls, supplied factors, normalization matrices and weights.
- The checker blocks production calls to R's adjustment/quantile/lowess/mean
  functions and rejects R numerical imports in the bridge. The reference R CI
  gate runs the checker against the checked source-package installation.
- Source package: **Status OK, 1,406 assertions, zero failures/warnings/skips**
  under native R 4.6.2/DESeq2 1.52.0. Its SHA-256 is
  `302cc45d942912a16b8fa12b4f17be9da7ec9f5680d5114f07e4f4742087bc23`;
  all 174 packaged runtime/configuration files match the checkpoint sources.
- Rust: **1,030 passed, one opt-in case ignored**. Core/backend Clippy and workflow
  lint pass. Instrumented boundary tests pass with C ASan/UBSan and Rust ASan;
  the R interpreter, installed dependencies and Rust standard library are not
  instrumented, and leak detection is disabled.
- Frozen GTEx replay: **121/121 per-contrast and 9/9 global gates pass**. Every
  diagnostic and result field except adjusted p-values remains byte-identical
  to v0.3.1. The canonical BH order changes 186,114 adjusted values by at most two
  ULPs (`2.22e-16` absolute); adjustment matches native R exactly when input bits
  are preserved. Frozen oracles and the old output-hash manifest are unchanged.
  This is an intentional numerical unification, not an old-hash acceptance claim.

The checker uses hexadecimal outputs because this environment's R decimal
parser can introduce a one-ULP input change. Bit-preserving transport separates
engine equality from text-import rounding.

## Storage and timing

The initial Rust-only filtering draft retained all candidate columns and sorted
per cutoff. That draft was rejected on performance. The retained implementation
sorts once, reuses the order for each population and stores only rejection counts
and the selected adjustment column. All nine earlier filtering payloads remain
bit-identical after this change.

One clean paired 56,200 × 16 workflow with ten individual result queries takes
9.816 → 8.859 s end to end and 1,354.25 → 1,318.66 MiB peak RSS; all saved arrays
and attributes match at tolerance zero between those two Rust-only versions.
This is not a counted native quantification speed round. The archived R-backed
checkpoint took 8.140 s and 1,158.75 MiB in a separate pass; restoring efficient
batch extraction on the standalone engine remains necessary. An exploratory run
overlapping compilation is excluded from timing evidence.

Raw logs, preserved executables, source packages, conformance payloads and
measurements are under `results/codebase_improvements/rust-only-*` locally.
