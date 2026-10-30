# Shared borrowed Cook's and replacement-row policy

This development checkpoint moves fitted-result Cook's masking and replacement
presentation out of production R. It follows `a7b256f`'s standalone contrast
batches, retaining locfit 0.2.0, L-BFGS-B 0.3.0 and the synchronized dependency
locks. Package versions remain 0.3.1.

## Shared implementation

`CooksPolicy` evaluates borrowed row iterators and optional `LowCountHeuristic`
context. Its output contains p-values and outlier flags only. It creates no
placeholder fitted/result table and performs no BH or independent filtering.
The existing public native `apply_cooks_cutoff*` adapters still recompute BH as
advertised. Stored `ResultRequest` skips the mask-stage BH pass, and R fitted
presentation no longer constructs any intermediate adjusted table. Each path
still computes its own final filtering/BH decision. Native legacy row builders
continue to populate an initial padj column before stored-request masking;
removing that unused initial work remains a later consolidation.

The strict `maxCook > cutoff` comparison, missing maxima, finite-only Cook's
argmax, established `total_cmp` ordering and last-tie selection are retained.
An eligible outlier is spared exactly when at least three raw counts exceed the
count at that sample. The scan can stop after the third larger count. Its result
does not depend on the p-value, alpha or previous requests.

`CountView` borrows native u32 or foreign double counts in either matrix order.
Double storage is validated completely as finite, nonnegative integer u32
values. The contrast count-zero view now uses that same storage/validation;
validation and global-zero detection still share a single pass. No count or
Cook's matrix is copied or transposed for masking.

`ReplacementRows` validates borrowed boolean/integer row flags and aligned
columns. Rows made all-zero by replacement receive effect/SE/stat/p = 0/0/0/1
through the same helper used by native replacement-result merging. Successful
refits and new-zero rows suppress the original maximum Cook's value to zero for
presentation, without changing stored diagnostics. All column lengths are
checked before a standalone column update. Mask evaluation validates before
native table mutation, so a later malformed row cannot leave an updated prefix.

Two-group eligibility now uses an O(samples) time/O(1) scratch partition check
against the second column of a two-column design. It replaces R's quadratic
equality outer products. R only supplies factor codes, declared level counts and
aligned design data. Missing factor observations cannot match a complete finite
design. Label spelling, level order and design-column scaling do not change a
matching partition.

A subsequent source audit of installed DESeq2 1.52.0 `results()` identifies two
remaining parity edges. DESeq2 uses `which.max()`'s first-tie choice and requires
a formula with exactly one variable that is a two-level factor. This checkpoint
preserves the earlier Rust/wrapper last-tie and partition-based behavior. The
new tests establish that legacy rule's migration, not exact DESeq2 equivalence
for tied maxima or formula-origin ambiguity. Independent DESeq2 fixtures and
shared strict origin context are still required before closing that part of C19.
The audit is retained locally as `rust-only-cooks-deseq2-tie-audit.log`.

## Corrected native eligibility

The previous native factor-metadata helper checked only whether observations
contained the numerator and denominator labels. It could enable the special
low-count exception for an additive design with extra coefficients, or for a
design whose partition did not match those labels. The R wrapper required a
matching two-column design. Native helper callers now pass the actual design
and share that stricter eligibility check. An explicit regression has counts
1/10/11/12 and a high Cook's value at the first sample: the matching two-column
design spares the row; adding a third design column must mask it.

Stored `ResultRequest` now infers eligibility from retained formula factors and
a matching design when no override was supplied. Its existing
`low_count_heuristic(bool)` method remains an explicit override. Missing formula
context does not infer eligibility; an eligible request requires retained Cook's
distances. This fixes a default-policy gap with the wrapper rather than changing
fitting arithmetic, probability functions or precision settings.

## R boundary and ownership

The fitted-result wrapper selects cached raw p-values, checks legacy schema
requirements, marshals metadata and attaches Rust output columns. It no longer
assigns statistical zero values, clears maxima, computes factor/design equality
matrices or obtains an intermediate adjusted table just to discard its padj.
The public primitive `applyCooksCutoffRust()` retains its own advertised adjusted
output, using the same mask kernel and one Rust BH call.

C validates complete column/context dimensions, indicator buffers and factor
code descriptors before Rust reads them. A short synchronous Rust preflight
validates flags and reports whether replacement rows require writable column
storage. C then shares untouched effect/SE/stat columns and copies only p-values
when masking; new-zero overrides receive distinct writable copies of all four
columns. The Rust call never borrows shared inputs mutably. Disabled filtering
without new-zero rows shares all four columns. Neither Rust call retains a
pointer or invokes R; allocation/error handling stays outside Rust execution.

Disabled presentation preserves untouched missing-value payloads. Active masking
keeps the previous quiet-NaN encoding for missing p-values. Wrapper probability
and finite/NA Cook's validation is retained, including validation after new-zero
overrides and maximum suppression on refitted rows. The lower-level native
heuristic continues its established finite-only scan for nonfinite Cook's
entries. The two APIs' input-validation contracts remain explicit.

## Verification

- Rust: **1,049 tests passed**, one existing opt-in full-GTEx case excluded from
  the ordinary suite. Backend: one passed. Python support tools: 45 passed.
- The final source package reports **Status OK**, **1,510 assertions**, no
  failures, warnings or skips. All **262** packaged files match working sources,
  allowing only generated DESCRIPTION fields and DCF whitespace.
- **101 standalone/wrapped bit-equality cases** pass against the final checked
  package. Fourteen focused new cases cover mask-only and fitted presentation,
  threshold boundaries, both layouts, refits/new-zero rows, declared levels,
  releveling, mismatched/additive designs, missing labels and disabled masking.
  The previous 87 inference/filter/contrast/batch cases remain in place.
- Independent R reference tests retain the explicit equality-outer rule, finite
  maxima, last ties and three-larger-count comparison. A runtime guard rejects
  production `outer()` calls; a mocked intermediate-BH function must not be
  called. Fixtures vary factor eligibility separately from cutoff modes instead
  of repeating every irrelevant cross-product. Existing DESeq2 numerical tests
  and frozen references are unchanged.
- New Rust tests exercise both storage orders, unsigned/double counts, exact
  cutoffs, tie behavior, missing/infinite entries on the native primitive path,
  masked versus unmasked invalid probabilities, bad counts, malformed flags,
  length validation before mutation, and automatic/explicit stored-request
  context. Original public Cook's and replacement coverage remains.
- Isolated native-boundary tests reject malformed result columns, indicator
  lengths/values, design/count/Cook's dimensions, factor codes and missing
  required buffers. Valid controls prove that both writable-copy and shared-input
  paths preserve their inputs, with masking enabled/disabled and new-zero rows
  present/absent. C ASan/UBSan and Rust ASan pass. The R interpreter, installed R
  dependencies and Rust standard library are not instrumented; leaks are not
  checked.
- All eleven full GTEx contrasts pass **121 per-contrast and nine global gates**,
  with minimum target/global error reductions **664×/597×**. All **22** native
  result/diagnostic exports are byte-identical to `a7b256f`. Those contrasts run
  with Cook's filtering disabled, so the dedicated eligibility/reference tests
  above establish the changed rule's correctness.
- Root/backend formatting, Clippy, Rust 1.97.1, strict Rust documentation,
  read-only vendor checks and crate packaging/compilation verification pass.
  Packaged crate sources, tests and examples also match the working copies.

## One-pass time and memory

All comparisons below use immutable candidate artifacts and the immediately
preceding `a7b256f` checkpoint. Fitting and extraction never run concurrently
within a profiler; compiler, source-check and instrumentation jobs finish before
measurements. Direct identity checks finished 7.62 seconds before the first
measured native process. Neither timing pass was repeated.

| Native tool / checkpoint | Sum over 11 GTEx contrasts | Maximum process RSS |
| --- | ---: | ---: |
| Previous native `a7b256f` | 35.719 s | 302.031 MiB |
| Candidate native | 36.461 s | 302.047 MiB |
| DESeq2, native R 4.6.2 | 189.877 s | 1,877.062 MiB |

Native elapsed increases **2.08%** and summed throughput decreases **2.04%**;
peak native RSS changes by 0.016 MiB. The observed DESeq2/candidate time ratio is
**5.21×**. The first `condition_blocked` pair is **3.838 → 4.630 s** and accounts
for the aggregate increase; all eleven pairs remain in the result. This is an
observed cost with no established cause, not a claim of unchanged native speed.
No 2% optimization round is counted. Follow-up and final combined performance
validation remain required.

All **18** paired R input/fit/output RDS payloads are **byte-identical** and pass
zero-tolerance object comparisons. Their summed process time is
**78.664 → 78.809 s (+0.19%)**; the preceding 17-case control set totals
**73.506 → 73.942 s (+0.59%)**. This mixed workload total is not a general speed
claim. Every process remains below **2,048 MiB**. Maximum R process RSS falls
from **1,137.266 to 1,122.297 MiB** (14.969 MiB).

| R workload, 56,200 genes × 16 samples | Fit baseline / candidate | Ten-result stage baseline / candidate | Stage time ratio | Process time ratio | Peak RSS baseline / candidate |
| --- | ---: | ---: | ---: | ---: | ---: |
| Two-group batch, low-count heuristic active | 1.102 / 1.104 s | 1.772 / 1.509 s | 1.17× | 1.06× | 981.719 / 954.719 MiB |
| Additive-design batch | 2.081 / 2.051 s | 1.669 / 1.637 s | 1.02× | 1.01× | 1,037.219 / 1,035.375 MiB |

Ratios are baseline/candidate. The two-group ten-result stage improves
**14.84%**, its complete process improves **5.64%**, and peak RSS falls
**27 MiB (2.75%)**. The stage includes result-name/request construction, dataset
validation, inference, masking, filtering and R table conversion; it does not
isolate the mask kernel alone. Complete process times additionally include
package startup, fitting and serialization. The earlier batch transport is
available in both packages here; these ratios measure the subsequent policy
change rather than comparing a batch with ten separate calls.

Controls showing costs are kept: the additive ten-individual-call process is
**7.601 → 7.684 s (+1.09%)**, while its ten-result stage improves from
**3.126 → 3.057 s**. The 1,200-gene LRT process is
**2.038 → 2.209 s (+8.36%)**, including fit/result stages of
0.112/0.039 → 0.144/0.042 s. The ashr process is
**2.162 → 2.482 s (+14.79%)**, with its shrink stage 0.130 → 0.134 s. The
largest rlog control is **21.290 → 21.259 s (−0.15%)** in this comparison; its
previous checkpoint's measured cost remains in that separate report. No rlog or
ashr arithmetic was modified, and these single observations do not establish
an explanation for the timing differences. Main README figures remain unchanged.

See the [native summary](data/rust_only_cooks_2026-10-30_native_summary.tsv),
[per-process native times](data/rust_only_cooks_2026-10-30_native_timings.tsv),
[paired R workflows](data/rust_only_cooks_2026-10-30_r_workflows.tsv),
[separate R stages](data/rust_only_cooks_2026-10-30_r_stages.tsv) and
[validation/artifact identities](data/rust_only_cooks_2026-10-30_validation.json).
[Per-contrast gates](data/rust_only_cooks_2026-10-30_gtex_precision.tsv),
[global gates](data/rust_only_cooks_2026-10-30_gtex_global_precision.tsv) and
[unchanged exports](data/rust_only_cooks_2026-10-30_native_output_identity.tsv)
remain independent of the speed measurements.

## Reproduction and remaining scope

Build `cargo build --locked --release --example evaluate_r_interface` and run
`scripts/check_rust_r_identity.R BINARY R_PACKAGE_LIBRARY OUTPUT_DIR` using
`deseq2-r462-native`, with the restored reference/test library in `R_LIBS_USER`.
Use `scripts/check_native_boundary_sanitizers.py` for the recorded instrumentation
scope. Detailed local artifacts use the `rust-only-cooks-*` prefix under
`results/codebase_improvements`.

The paired R profiler retains all previous workflows and adds a full-size
56,200 × 16 two-group ten-result case to exercise borrowed count/Cook's matrices.
It records fit and extraction stages separately and compares complete input,
fitted-state and output payloads. Native timing remains one accepted paired pass
across all eleven full contrasts. No timing run is repeated to replace a result.

The [full migration](rust-only-wrapper-migration.md) remains in progress. Normal
shrinkage's global-zero check, remaining offset/normalization/formula arithmetic,
optional ashr/apeglm and posterior calculations still need the standalone path.
C18's enabled-but-no-refit allocation remains; no native ownership handle is
introduced. The final combined performance audit must also retain previous
observed timing costs rather than declaring every workflow faster.
