# Shared standalone contrast-zero rules

This checkpoint moves contrast-specific count-zero decisions and result
overrides out of production R. Numeric design scores, selected-sample zero
checks, global-zero exclusion and Wald/LRT overrides now use standalone Rust.
R aligns metadata and passes borrowed primitive arrays into the existing shared
contrast call. Package versions remain 0.3.1, with the latest-stable dependency
refresh from `9f175f3` retained.

## Shared implementation

`ContrastCounts` borrows validated unsigned or double count matrices in either
storage order. Double inputs must represent finite, nonnegative integer u32
counts. Validation and global-zero detection share one pass, and the resulting
gene flags can be reused for multiple requests. Selected counts are tested for
zero directly; there is no count sum or selected-count matrix allocation.

The numeric and factor selection helpers serve both the existing public native
functions and the borrowed primitive interface. Numeric requests with both signs
select samples using the nonzero coefficient pattern and the existing checked
design-score arithmetic. One-sided requests return false without requiring a
model matrix. The native convenience function retains its early exit before
preparing counts, and invalid requests are checked before a count scan. Factor
requests use the established level/alias rules.

`ContrastZeroPolicy` validates flag/column lengths before mutation. A row is
eligible only when the requested samples are zero and the gene is not globally
zero. Wald sets effect/statistic/p-value to 0/0/1 and retains the SE. LRT changes
only its displayed effect; the stored likelihood-ratio inference is retained.
Missing t df alone does not prevent the established Wald zero override. Primitive
inputs with finite supplied estimates on globally zero rows explicitly test the
exclusion rule rather than relying on missing fitted coefficients to conceal it.

The two identical native scaled-sum implementations are consolidated into a
private math helper without changing the operation order. It serves both contrast
projection and numeric zero selection. Matrix layout changes indexing only.

The R bridge applies zero policy to the owned Rust projection before formatting
outputs. It borrows counts/design buffers and marshals UTF-8 labels, with no R
numerical callback. C validates the context, complete array lengths, dimensions,
enum fields and label buffers before entering Rust. R's remaining helpers handle
identifier alignment and canonicalization. No pointer is retained across calls.

## Verification

- **1,043 Rust tests** pass, with one existing opt-in full-GTEx test excluded
  from the ordinary suite. Backend unit test: one passed. Python support tests:
  45 passed.
- The final R source package reports **Status OK**, **1,450 assertions**, no
  failures, warnings or skips. All **249** packaged source files match the working
  sources, allowing only generated DESCRIPTION fields and DCF whitespace.
- **85 standalone/wrapped bit-equality cases** pass against the checked package.
  Twenty-one new cases cover raw/Normal/t inference, numeric/factor/fractional
  contrasts, factor-reference versus one-sided numeric rules without design
  context, cancellation selecting
  no samples, globally zero exclusion and UTF-8 labels. Targeted runtime guards
  reject production R count sums, matrix products and square roots, alongside
  the existing statistical-function/import guards.
- Independent R calculations in tests verify count selection, covariance SEs,
  directional t probabilities and both override policies. Existing numerical
  references/tolerances remain unchanged. New Rust tests cover storage layouts,
  cached flags across requests, missing context, malformed counts and validation
  of every column before mutation.
- The eleven full GTEx contrasts pass **121 per-contrast and nine global gates**,
  retaining minimum target/global error reductions of **664×/597×**. All **22**
  result/diagnostic exports are byte-identical to `9f175f3`'s saved outputs.
- C ASan/UBSan and Rust ASan boundary checks pass, including malformed zero
  contexts and valid numeric/factor controls. The R interpreter, installed R
  dependencies and Rust standard library are not instrumented; leak detection
  is disabled.
- Root/backend formatting, Clippy, Rust 1.97.1, strict Rust documentation, vendor
  synchronization and crate packaging/compilation verification pass.

The [validation record](../data/rust_only_zero_2026-10-30_validation.json) identifies
the immutable executables and checked source archive. See the
[GTEx gates](../data/rust_only_zero_2026-10-30_gtex_precision.tsv),
[global gates](../data/rust_only_zero_2026-10-30_gtex_global_precision.tsv) and
[export identities](../data/rust_only_zero_2026-10-30_native_output_identity.tsv).
Detailed logs and artifacts remain under
`results/codebase_improvements/rust-only-zero-*`.

## Remaining work

The reusable count context is preparation for bounded result batches; the batch
API and dataset-validation reuse are still pending. Cook's/refit presentation
policy, remaining normalization/offset/formula arithmetic and optional posterior
calculations still need the shared Rust path. In particular, the remaining
global-zero calculation in normal shrinkage is outside this contrast checkpoint.

This checkpoint adds no speed or RSS claim and counts no optimization round.
Combined result/batch timing and memory validation remains required. The
[C01–C23 audit](codebase-improvement-completion-audit-2026-10-29.md) and
[Rust-only migration inventory](rust-only-wrapper-migration.md) remain active.
