# Shared standalone inference probabilities

This checkpoint removes production R calculations of Wald hypothesis statistics
and Normal/t tails, primitive LRT chi-square probabilities, and default Cook's
F cutoffs. R validates/translates arguments and formats returned values; both
interfaces call the same standalone Rust implementation. No R probability
callback, imported R math routine or frontend-specific numerical mode is used.
Package versions remain 0.3.1 and the larger migration remains in progress.

## Implementation

`wald_test_from_estimates` evaluates borrowed estimate/SE vectors without
constructing a fitted model. It calls the same Wald scalar implementation as
fitted coefficient and contrast queries. `lrt_pvalues_from_statistics` reuses
the fitted-LRT chi-square evaluation and caches repeated df. The new native
entry points validate types, dimensions, scalars and storage lengths before
borrowing buffers. They retain no pointers. R retains gene-name alignment and
output formatting; computed LRT p-values preserve aligned gene names.

Missing/nonpositive t df now preserves the selected threshold-hypothesis
statistic and leaves its p-value missing. Previously the Rust scalar path
returned the unthresholded statistic in this case. The new tests check every
supported t alternative, invalid estimates and fitted/primitive bit equality.
Primitive residual df requires a fitted design; callers supplying estimates
alone must provide scalar or per-gene df. Nonfinite inference follows the
existing standalone missing-row policy.

The probability code is organized under `src/probability/`:

- Normal tails use the pure-Rust `libm` erfc implementation, with compensated
  reduction by sqrt(2). The previous approximation differed from R by
  `2.720e-11` at z=1 and failed existing wrapper precision checks.
- Student t tails use a positive incomplete-beta series/continued fraction for
  smaller df and a scaled density integral for large df, avoiding gamma-log
  cancellation and nearly-one beta arguments. The methods follow the
  [incomplete-beta identities](https://dlmf.nist.gov/8.17) and
  [Gauss quadrature](https://dlmf.nist.gov/3.5#v). Quadrature tests independently
  verify polynomial moments through degree 31 and analytic tail integrals.
- LRT df=2 uses its exponential closed form. Large-shape chi-square tails use a
  stable gamma density normalization and integrate outward from the requested
  tail. This corrects a `1.297e-10` error at statistic=df=1,000,000.
- Cook's cutoffs use stable beta-logistic density inversion. The common
  inverse-beta output operation order and the extreme-df chi-square limit
  conventions match the [R qf implementation](https://github.com/wch/r-source/blob/trunk/src/nmath/qf.c).
  Those extreme-df conventions are compatibility approximations, applied
  identically through both interfaces. They are covered around their df
  boundary; a generic mathematically exact F quantile would differ from R there.

`libm` 0.2.16 was already present transitively and is now a direct dependency.
No R runtime is linked or launched by the standalone conformance executable.

## Independent precision evidence

`scripts/references/generate_inference_reference.R` produces the new
`crates/rsdeseq2/tests/data/inference_r462.tsv` using pinned R 4.6.2 statistics.
Inputs and expected results are IEEE-754 hex, avoiding decimal-parser rounding.
The generator and R session are committed. Existing references, release
manifests and precision bounds are preserved.

The panel contains **1,169 scalar cases**: 418 Normal probabilities, 572 t
probabilities, 136 chi-square probabilities and 43 Cook's cutoffs. It includes
branch boundaries, df through 1e300, heavy/subnormal tails, and df near the F
limit-policy boundary. Normal checks use eight machine epsilons in relative
error; the other families use 2e-12 relative error. The absolute floor is only
four smallest subnormal ULPs. The existing tighter wrapper assertions remain
unchanged and pass.

The [scalar precision summary](../data/rust_only_inference_2026-10-29_scalar_precision.tsv)
records the observed errors. Maximum Normal relative error is `3.90e-16`.
Maximum absolute t error is `3.84e-15`; maximum relative t error is `1.13e-13`
on a subnormal probability. Maximum chi-square absolute error is `4.40e-14`.
Maximum Cook's relative error is `1.32e-14`; its largest absolute difference,
`8.37e-11`, occurs at a cutoff near 6,366.

The full eleven-contrast GTEx replay passes **121 per-contrast and nine global
precision gates**, retaining minimum global/target improvements of 597×/664×.
These retain the existing raw Wald comparison scope with Cook's cutoff and
independent filtering disabled; enabled policy has separate Rust/R coverage.
All non-probability result columns are unchanged and all eleven diagnostic
exports are byte-identical to the saved `dee5a4e`/`bfd0de8` outputs. The improved
Normal calculation changes 231,409 p-values by at most `4.940e-11` and 51,986
adjusted p-values by at most `8.902e-11`. These are intentional mathematical
corrections; this checkpoint does not claim that entire result exports match
the old approximation byte for byte.

See the [precision gates](../data/rust_only_inference_2026-10-29_gtex_precision.tsv),
[global gates](../data/rust_only_inference_2026-10-29_gtex_global_precision.tsv) and
[output changes](../data/rust_only_inference_2026-10-29_native_output_changes.tsv).

## Validation and artifacts

- Rust: **1,036 passed**, zero failures, one existing optional full-GTEx test
  ignored in the normal suite. The full replay is independently run above.
- Final R source package: **Status OK**, **1,408 assertions**, no failures,
  warnings or skips. All 236 packaged source files match the working sources,
  allowing only generated DESCRIPTION fields/whitespace.
- **46 standalone/wrapped bit-equality cases** pass, with runtime guards against
  R statistical calls and an imported-R-math symbol check.
- Root/backend Clippy and formatting, Rust 1.97.1, strict documentation, vendor
  synchronization and crate package compilation pass.
- C ASan/UBSan and Rust ASan native-boundary tests pass, including malformed
  primitive buffers and valid numerical controls. R, installed dependencies
  and the Rust standard library are not instrumented; leak checking is off.

The [validation record](../data/rust_only_inference_2026-10-29_validation.json)
binds numerical evidence to the saved executables and checked R source archive.
Detailed logs and immutable artifacts are retained under
`results/codebase_improvements/rust-only-inference-*`.

## Remaining migration

Contrast products/SEs, all-zero decisions, offset/formula arithmetic, bounded
result batches, optional shrinkage and posterior calculations still require
Rust-only implementations. This checkpoint adds no speed or memory claim and
counts no optimization round. Performance/RSS validation will cover the combined
raw-contrast/batch migration; earlier timing evidence remains tied to its
original binaries. The complete acceptance requirements remain in the
[Rust-only migration inventory](rust-only-wrapper-migration.md).
