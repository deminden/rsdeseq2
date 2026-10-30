# Conditional replacement output and borrowed Cook's factors

This checkpoint addresses C18's enabled-but-no-refit R allocation, together with
two temporary normalization-factor matrices in native Cook's preparation. Its
baseline is `c96e225`. Package versions remain 0.3.1. The wider
[Rust-only migration](rust-only-wrapper-migration.md) remains in progress.

## Output ownership and computation

The R fit adapter previously allocated a full double replacement-count matrix
whenever replacement was enabled. Rust then copied either replaced counts or
original counts into it. The public R adapter retained that matrix only when
`nRefit > 0`, discarding it otherwise.

The fit bridge now allocates the R matrix only after native fitting reports a
retained replacement output. On success, Rust publishes either NULL or a moved,
owned `CountMatrix`; it retains no pointer into R storage. C queries its shape,
allocates and copies the R matrix, then consumes and clears the owner. There is
no persistent R external pointer and no fitting rerun. The public fit still
retains replacement counts under the same `nRefit > 0` condition.

That count includes replacement rows that become all zero, even when no GLM
needs refitting. The new tests explicitly cover that case for Wald and LRT,
alongside the existing mixed refit/new-zero fixture. Original counts, replacement
counts, flags, assays, raw statistics and result metadata preserve their meaning.
The private `.Call` fit output now uses NULL for an unused replacement matrix;
the public R fit already used NULL in that situation.

`fit_abi.h` records the opaque owner, fit output slot, shape/copy and consuming
drop contracts. Before native fitting, C allocates and protects an unwind
continuation while no owner exists. After Rust returns and releases all borrowed
R buffers, copying runs under `R_UnwindProtect`. Cleanup consumes the owner on
normal return or an R error. R allocation and longjmp never occur across a live
Rust fitting/copy frame. Preallocating the continuation also avoids a possible
allocation failure before the unwind context is established; see R's
[cleanup API](https://cran.r-project.org/doc/manuals/r-release/R-exts.html#Condition-handling-and-cleanup-code).

Native Cook's candidate-count preparation now borrows a supplied factor row or
the sample size-factor slice instead of cloning/expanding a full factor matrix.
Recomputing normalized replacement counts likewise borrows supplied factors or
uses the existing direct size-factor normalization operation. Multiplications,
divisions, trimming, validation, replacement eligibility and refit decisions are
unchanged. The two removed native temporaries occur in separate phases; their
sizes must not be added when claiming a peak-memory reduction.

At 56,200 genes × 16 samples, each double matrix has **6.860 MiB** of element
storage. The unused R buffer and two native factor-matrix allocations are
avoided. These are allocation-volume facts, not a claim that process RSS falls
by their sum.

## Lifecycle and numerical checks

The full Rust suite passes **1,065 tests**, plus two backend tests. The owner
test checks dimensions, column-major copy values, malformed copy recovery and
clearing/idempotent drop. One existing opt-in root test remains ignored.

The checked R source package reports **Status OK**, **1,988 assertions**, and
no failures, warnings or skips. The existing DESeq2 replacement parity test is
table-driven over mixed refits and all-zero-only updates for both Wald and LRT;
its original precision tolerances are retained. Boundary child controls cover
disabled, enabled/no-refit, ineligible, mixed and all-zero-only outputs.

A disposable sanitizer build alone enables four injected errors around the
actual fitting/output path: before R allocation, after allocation, native copy
shape rejection and after a successful copy. Every case observes one returned
owner consumed, zero live transport owners, and an exact subsequent recovery fit.
The ordinary checked package contains no injection hook or counters. C ASan/
UBSan and Rust ASan instrumentation is verified. R, installed dependencies and
the Rust standard library are not instrumented; leak detection is disabled.
These controls certify the stated owner lifecycle, not whole-process leak
freedom or an actual process-wide out-of-memory event.

All **189** existing standalone/wrapped identity cases pass. All **121**
per-contrast and **nine** global GTEx precision gates pass, preserving the
664×/597× minimum reductions, and all **22** native result/diagnostic exports
remain byte-identical. Root/backend Clippy, Rust 1.97.1 checks, strict Rust docs,
vendor consistency, crate packaging and all 46 distinct Python support tests
pass. The initially skipped R-enabled Python check is run separately.

The checked R archive matches **290** current source members, allowing only
generated DESCRIPTION fields/DCF whitespace; all **453** packaged crate source,
test and example members match. R allocation profiling is unavailable in the
selected reference build, so no complete R heap trace is claimed. Source-level
conditional allocation, native owner controls and process RSS have separate
evidence scopes.

## First attempt retained

The first candidate deferred the R output but left the native factor temporaries
in place. Its single paired twenty-workflow run preserved every complete R
payload, but the large no-refit case rose **946.344 → 950.547 MiB** (+4.203 MiB)
and **6.482 → 6.794 s**. C18's memory gate was left open. Panel time changed
95.197 → 92.703 s and maximum RSS 1127.516 → 1127.375 MiB; those aggregate values
do not establish success for the target workflow.

The [first-attempt record](data/replacement_buffer_2026-10-30_attempt1_validation.json)
retains that outcome. The unchanged candidate was not remeasured. The second
candidate adds the independent native temporary-factor removal and the
all-zero-only tests before its own single paired measurement.

## Combined candidate measurement

The combined candidate receives one paired pass over all eleven full GTEx
contrasts and one paired pass over all twenty R workflows, each with one native/
BLAS thread. Builds, conformance checks and unrelated fits do not overlap the
measurements. Timed native result hashes and separate untimed diagnostics match
their reviewed exports; all forty-four frozen input/reference files are intact.

| Measurement | Before | Combined candidate | Before / candidate |
| --- | ---: | ---: | ---: |
| Native total elapsed, 11 contrasts | 36.131 s | 36.238 s | 0.9971× |
| Native maximum peak RSS | 301.875 MiB | 302.031 MiB | 0.9995× |
| R total elapsed, 20 workflows | 91.322 s | 91.818 s | 0.9946× |
| R maximum peak RSS | 1127.219 MiB | 1127.172 MiB | 1.0000× |
| Large no-refit peak RSS | 945.391 MiB | 944.375 MiB | 1.0011× |
| Large no-refit elapsed | 6.467 s | 6.489 s | 0.9966× |

The target no-refit peak falls **1.016 MiB**, with elapsed **+0.35%**. C18's
local allocation/lifecycle/RSS gate passes for the combined change. The observed
RSS reduction is smaller than the removed allocation volume; no complete heap
trace or causal attribution to an individual removal is claimed.

Native elapsed rises **0.30%** and R total **0.54%**. Native throughput falls
0.29%, so **no 2% optimization round is counted**. Original DESeq2 takes
192.649 s, giving a **5.316×** DESeq2/combined-native process-time ratio. These
are single-pass observations, not a statistical-equivalence claim; every case
is retained. Each attempt has its own paired baseline, so candidate-only values
from different passes are not treated as an isolated effect estimate.

All **twenty complete R workflow payloads are byte-identical**, with default
strict comparison and no metadata exception. Every case remains below 2 GiB.
The largest time increases are ashr (+0.559 s) and apeglm (+0.431 s); normal
shrinkage instead falls 0.473 s. The large two-group result case rises 4.48 MiB
and 0.149 s. Those costs remain recorded. Panel maximum RSS values can belong
to different cases, and earlier measured memory increases remain open for the
final combined audit.

See the [native timings](data/replacement_buffer_2026-10-30_native_timings.tsv),
[native summary](data/replacement_buffer_2026-10-30_native_summary.tsv),
[R workflows](data/replacement_buffer_2026-10-30_r_workflows.tsv),
[R stages](data/replacement_buffer_2026-10-30_r_stages.tsv),
[complete payload identities](data/replacement_buffer_2026-10-30_r_payload_identity.tsv),
[GTEx gates](data/replacement_buffer_2026-10-30_gtex_precision.tsv),
[native export identities](data/replacement_buffer_2026-10-30_native_output_identity.tsv),
[first-attempt R workflows](data/replacement_buffer_2026-10-30_attempt1_r_workflows.tsv)
and [validation record](data/replacement_buffer_2026-10-30_validation.json).

## Remaining scope

No R numerical callback is added, but production R expression evaluation,
frozen-local VST and optional posterior/s-value calculations still need their
shared standalone paths. Initial redundant native adjustment, final C21
integration and C23's combined performance/memory/source reconciliation remain
open. Earlier measured memory costs remain in that final audit. Main README
claims and published release manifests are unchanged.
