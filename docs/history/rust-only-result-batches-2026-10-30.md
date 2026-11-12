# Bounded result batches on the standalone engine

`resultsList()` reuses one validated fitted-dataset view and one primitive
extraction for several result requests. The native single-call and batch paths
now use the same `PreparedContrastView`. No R numerical callback is added.
The package version remains 0.3.1; this is a local development checkpoint after
`b0855d0`, retaining the stable dependency refresh (`9f175f3`).

## API and shared implementation

```r
requests <- list(
  forward = list(contrast = c("condition", "B", "A")),
  reverse = list(contrast = c("condition", "A", "B"), alpha = 0.05),
  unmasked = list(name = tail(rsdeseq2::resultsNames(dds), 1L), cooksCutoff = FALSE)
)
tables <- rsdeseq2::resultsList(dds, requests, batchSize = 32L)
```

Requests accept the existing result selectors and options. Common `...` arguments
are overridden by each request, including explicit NULL. Outer names and order
are preserved. Dataset, retained Rust fit, primitive Wald/LRT and supported list
inputs keep their existing contextual requirements and result classes/metadata.
Empty requests return an empty list after validating the supplied object.
Cache-only legacy requests do not require reconstruction of absent coefficient
arrays.

`PreparedContrastView` borrows beta, optional SE and covariance matrices. Optional
counts are validated and global-zero flags prepared once per native chunk. It
reuses the established projection, Normal/t inference and zero policy, consuming
one owned result at a time. The R interface borrows column-major storage; native
callers may use row-major storage with the same logical operation order.

The shared C request descriptor carries request options, zero metadata and
caller-owned output arrays. C validates all request buffers and allocates R
outputs before synchronous Rust execution. Rust performs no R calls and retains
no pointers. Single `.Call` names and argument counts remain supported through
a one-request descriptor. A native chunk is limited to 64 requests; the public
default is 32. Returned output storage necessarily grows with the total number
of tables; temporary output arrays are bounded by the chosen chunk size.

R shares selector/metadata adapters and fitted-result finishing between single
and bulk calls. Each request keeps its own Cook's masking, alpha, hypothesis,
filter selection and BH adjustment. No decision is cached across different
p-value populations. Selected fitted results avoid the former preliminary BH
and unfiltered-BH passes before their final Cook's/filtering step.

## Correctness and source checks

- Rust: **1,045 passed**, one existing opt-in full-GTEx case excluded from the
  ordinary suite. Backend: one passed. Python support tools: 45 passed.
- The final R source package reports **Status OK**, **1,488 assertions**, no
  failures, warnings or skips. All **255** packaged source files match the
  working sources, allowing only generated DESCRIPTION fields and DCF whitespace.
- **87 direct standalone/wrapped bit-equality cases** pass against that checked
  package. These include all 85 previous cases plus mixed raw/Normal/t batches
  and a 65-request public call crossing native chunk boundaries. Runtime guards
  reject R statistical/count-reduction/matrix-product calculations in the scoped
  primitive result path, and native import checks reject R probability backends.
- Batch-versus-individual R comparisons are `identical()`, including metadata,
  classes, names, order, scalar/per-gene df, request-specific alpha/filtering,
  Cook's cutoffs, subsets, edited-dataset rejection, replacement/new-zero rows,
  Wald/LRT and legacy cached fits. A spy verifies one dataset validation, one
  primitive extraction and two native calls for four selected requests split
  into chunks of two, alongside a cached request.
- Child-process boundary tests reject malformed lists, short buffers, duplicate
  fields, incompatible contexts, invalid options and 65-request native chunks.
  Positive controls exercise all 64 descriptor positions with mixed projection,
  Normal/t and factor/numeric zero contexts. C ASan/UBSan and Rust ASan pass.
  The interpreter, installed R dependencies and Rust standard library are not
  instrumented; leak detection is disabled.
- All eleven full GTEx contrasts pass **121 per-contrast and nine global gates**.
  All **22** result/diagnostic exports are byte-identical to `b0855d0`.
- Root/backend formatting, Clippy, Rust 1.97.1, strict Rust documentation,
  read-only vendor synchronization and crate packaging/compilation verification
  pass. No frozen reference, precision tolerance, solver setting or package
  version changes.

## One-pass native measurements

The accepted paired pass uses all eleven full GTEx contrasts, one worker and
one BLAS thread. Native total process time includes startup, parsing, full Wald
fitting and TSV output; exact result/diagnostic validation runs separately.
Tool order rotates across contrasts. No timing pass was repeated.

| Tool / checkpoint | Summed time, 11 contrasts | Maximum process RSS |
| --- | ---: | ---: |
| Previous native `b0855d0` | 37.104 s | 302.156 MiB |
| Candidate native | 37.238 s | 302.328 MiB |
| Original DESeq2, R 4.6.2 | 200.934 s | 1,876.984 MiB |

Native elapsed time changes **+0.36%**, and summed throughput **−0.36%**.
The DESeq2/candidate time ratio is **5.40×**. Peak native RSS changes by only
0.172 MiB. This adds no 2% quantification optimization round. No fitting
arithmetic or routing changed, and all 22 exports remain identical.
See the [native summary](../data/rust_only_batch_2026-10-30_native_summary.tsv) and
[per-process measurements](../data/rust_only_batch_2026-10-30_native_timings.tsv).

## Paired R workflow time and memory

All **17** paired complete input/fit/output RDS payloads are **byte-identical**;
the profiler's zero-tolerance object comparison also passes. One pass uses
R 4.6.2 and one worker/BLAS thread per process. Returned tables and retained
fitted-state contents are included in these payloads. Every measured process
stays below the existing **2,048 MiB** budget; the overall maximum changes from
**1,129.781 to 1,133.047 MiB** (3.266 MiB). Main README memory claims are unchanged.

The bulk workload requests ten tables from the same fitted dataset. Its baseline
uses ten individual calls; the candidate uses `resultsList()`. Fit time is
measured separately and is essentially the same for these two paired cases.

| Bulk workload | Fit baseline / candidate | Ten-results baseline / candidate | Extraction time ratio | Whole-process time ratio | Peak RSS baseline / candidate |
| --- | ---: | ---: | ---: | ---: | ---: |
| 1,200 genes × 16 samples | 0.108 / 0.108 s | 0.519 / 0.304 s | 1.71× | 1.09× | 868.406 / 869.750 MiB |
| 56,200 genes × 16 samples | 2.065 / 2.071 s | 2.878 / 1.687 s | 1.71× | 1.20× | 1,090.844 / 1,037.391 MiB |

The larger batch reduces extraction elapsed time by **41.38%** and process peak
RSS by **53.453 MiB (4.90%)**. Ratios above are baseline/candidate; whole-process
time includes package startup, fit, extraction and serialization.

Timing controls are retained explicitly. The unchanged ten-individual-call
workload at 56,200 genes takes **2.882 → 3.003 s** in the extraction stage and
**7.241 → 7.379 s (+1.92%)** for the process; the 1,200-gene extraction remains
0.520 s. Existing 14 control workflows sum to **58.583 → 60.097 s (+2.58%)**.
The largest rlog process is **21.346 → 22.683 s (+6.26%)**, including an rlog
stage of **18.842 → 20.097 s**. These are observed costs in this one pass,
not evidence that every existing workflow became faster. No rlog arithmetic was
changed. This timing cost remains a follow-up in the final performance audit;
no run was repeated to replace it. All 17 workloads together total
**74.968 → 75.100 s (+0.18%)**, a workload-mix figure rather than a general speed
claim. Query-stage gains do not count as a native quantification speed round.

The larger enabled-but-no-refit case has zero refits and peaks at
**941.266 → 940.734 MiB**. Both public fitted states retain zero bytes for
`replaceCounts`/`replaceCooks`. C still allocates its temporary enabled-path
replacement output before Rust reports the refit count. Its double-array payload
is **6.860 MiB** for 56,200 × 16, about 0.73% of this process peak; that is a
storage calculation, not a measured saving. No native ownership handle is added
in this checkpoint, and C18's remaining no-refit allocation is not claimed fixed.
The disabled path's earlier omission remains in place.

See the [paired workflow measurements](../data/rust_only_batch_2026-10-30_r_workflows.tsv),
[separate stage times](../data/rust_only_batch_2026-10-30_r_stages.tsv), and
[validation/source identities](../data/rust_only_batch_2026-10-30_validation.json).
The [GTEx precision gates](../data/rust_only_batch_2026-10-30_gtex_precision.tsv),
[global gates](../data/rust_only_batch_2026-10-30_gtex_global_precision.tsv) and
[unchanged native exports](../data/rust_only_batch_2026-10-30_native_output_identity.tsv)
retain their independent reference bounds.

## Reproduction and evidence

Build the standalone driver with
`cargo build --locked --release --example evaluate_r_interface`, then run
`scripts/validation/check_rust_r_identity.R STANDALONE_BINARY R_PACKAGE_LIBRARY OUTPUT_DIR`
using `deseq2-r462-native` and the restored reference/test library in
`R_LIBS_USER`. Source-package checks use `R CMD check --no-manual
--no-build-vignettes`. Instrumentation uses
`scripts/validation/check_native_boundary_sanitizers.py` with that same R runtime.

The paired R workflow profiler includes both ten individual calls and a bulk
request case at 1,200 and 56,200 genes. For a package predating `resultsList`, the
bulk case executes the same requests individually. It records fit/extraction
stages separately, compares complete fitted/input/output payloads at zero
tolerance, and enforces the existing 2 GiB process-RSS budget. An additional
56,200-gene enabled-but-no-refit workflow assesses the remaining replacement
buffer question. The native profiler retains its accepted-output hash gates and
one paired eleven-contrast timing pass.

Detailed local artifacts are under
`results/codebase_improvements/rust-only-batch-*`. The first direct source build
failed while copying changing compiler intermediates before `.Rbuildignore`
exclusion. The final source build uses a clean snapshot excluding generated
intermediates; its 255-file comparison binds it to the working sources. This
build failure occurred before timing and is preserved in `source-v1/build.log`.

## Remaining scope

C20's bounded transport and validation reuse are implemented. C19's full
presentation policy still needs the shared standalone path: Cook's/refit
presentation, remaining offsets/formula arithmetic, normal-shrinkage zero flags
and posterior calculations remain partly in R. Optional ashr/apeglm still require
standalone implementations. These batches do not claim completion of the
[Rust-only migration](rust-only-wrapper-migration.md) or the full
[C01–C23 audit](codebase-improvement-completion-audit-2026-10-29.md).
