# Shared count, shrinkage-row and offset preparation

This checkpoint moves count validation, normal-shrinkage row decisions and
base-factor/offset composition into standalone Rust. The baseline is `e255eea`;
package versions remain 0.3.1. The broader
[Rust-only migration](rust-only-wrapper-migration.md) is still in progress.

## Shared computation

`CountRows` provides a borrowed count view and cached raw-zero flags in either
storage order. Double-count validation and zero detection happen in one physical
scan. `ContrastCounts` delegates its existing cache to this component; numerical
contrast behavior is unchanged. R count validation calls `CountView` directly,
without allocating R floor/comparison matrices. The fitting adapter uses that
same validator before copying counts into owned native storage; this does not
claim zero-copy fitting.

`NormalShrinkageRows` selects eligible input rows and validates their final and
fitted dispersions. Raw-zero counts and explicit fitted exclusions are distinct
inputs: a nonzero row excluded by weights must also stay out of prior estimation.
Unused missing dispersions on excluded rows remain allowed. Selection preserves
input order, including supplied-result subsets/permutations. R marshals the
selected rows and attaches returned values; it no longer computes row sums,
dispersion eligibility or the selection itself.

`normalization_factors_with_offsets` is shared by native formula workflows and
the R bridge. Effective factors are calculated as **base × exp(offset)**, with
gene/sample factors taking precedence over size factors. Offsets are aligned by
the wrapper, each exponential is evaluated once, and inputs/products must remain
finite and positive. SF estimation still uses the existing Rust implementation.
The kernel supports both borrowed layouts and preserves multiplication order.
Its result is an owned native matrix; the bridge copies that result into R's
column-major output rather than choosing an alternate arithmetic backend.

The native normalization path now validates supplied size factors themselves
when gene/sample factors take precedence. It no longer constructs an unused
`counts / baseSF` matrix. This removes an avoidable full matrix and permits tiny
positive base factors whose effective offset factors are representable. Offset
shape validation remains before size-factor estimation.

## Confirmed correctness fixes

- **Weight-excluded normal shrinkage.** On an independently fitted original
  DESeq2 example, a nonzero gene is excluded by its weights and has missing
  dispersion. Both fits mark it excluded. The old wrapper rejected shrinkage
  because it tested only raw row sums; the new path returns NA like DESeq2.
  Coefficient, expanded-contrast, positive-threshold statistic/p-value/adjustment,
  and supplied-result subset/permutation comparisons pass existing 1e-5 bounds.
- **Unnamed-factor offset provenance.** Effective factors are returned with
  canonical gene/sample names. Previously, unnamed base factors caused saved
  effective factors to fail the identity check on refit, applying offsets twice.
  The reproduced maximum factor difference was **0.336**; it is now **zero**,
  with bit-identical factors/results on refit.
- **Anonymous gene inputs.** Dataset fits previously failed while assigning
  covariance dimensions when gene names were absent. Stable interface row
  identifiers are assigned before preparation; the caller's original dataset
  stays unchanged. Independent DESeq2 numerical comparisons and subset queries
  cover named and anonymous inputs.
- **Explicit cleaned factor aliases.** Without a full model frame, an explicit
  `cell.type` alias resolved a contrast but lost the low-count exception for the
  original `cell type` variable. The adapter now canonicalizes only against the
  unique original variable. Full-frame and R bridge names remain authoritative.

These are demonstrated corrections, not new compatibility rules inferred from
an old wrapper. Their failing baseline probes remain in local evidence. Existing
valid reference inputs, numerical tolerances and release manifests are unchanged.

## Independent evidence and safety

`scripts/generate_row_preparation_reference.R` requires R 4.6.2 and DESeq2 1.52.0,
never loads rsdeseq2 and refuses to overwrite existing references. Supply fresh
output paths to reproduce the frozen files in the Rust/R tests. It records twelve
independent R offset cases and six row-selection cases. Coverage includes SF/NF
precedence, zero offsets, tiny rescued base factors, subnormal products,
overflow/underflow, missing/infinite offsets, raw/fitted exclusions, reordered
rows, invalid dispersions and no eligible genes. Floating inputs/outputs use
exact IEEE hex. Original DESeq2 weighted-shrinkage parity is checked separately
on actual fitted datasets rather than treating the arithmetic fixture as a full
fitting oracle.

All **131** standalone/wrapped conformance cases pass, including the eighteen new
preparation cases; production R numerical functions are guarded against execution
in the direct tests. The full Rust suite passes **1,060** tests, followed by the
final targeted preflight regression. Backend, Clippy, formatting, Rust 1.97.1
minimum-version, documentation, vendor integrity, Python's **46** support tests
and crate package checks pass. The R source package reports **Status OK**,
**1,668 assertions**, no failures, warnings or skips. All **272** source members
match working files, apart from generated DESCRIPTION fields/DCF whitespace;
all **440** packaged crate source/test/example files match.

The new synchronous bridge validates matrix/vector shapes and optional logical
exclusions, protects output allocations and retains no foreign pointer after
return. Malformed-call child tests cover the new entry points alongside valid
count, row-selection and factor-product controls. Fresh C
AddressSanitizer/UndefinedBehaviorSanitizer and Rust AddressSanitizer checks pass.
The R interpreter, installed R packages and Rust standard library are not rebuilt
with instrumentation; leak detection is disabled. This is the recorded bridge
scope, not whole-process leak certification.

All **121** GTEx precision gates and **9** global gates pass. Minimum reported
target/global error reductions remain 664×/597×. All **22** result/diagnostic
exports are byte-identical to `e255eea`. These disabled-Cook's fitting comparisons
complement the specific independent preparation/shrinkage references.

## Timing and memory

One paired native pass covers all eleven full GTEx contrasts. One paired R pass
retains the previous eighteen workflows and adds full-size 56,200 × 16 SF-offset
and gene/sample-factor-offset fits. The same deterministic workflow script runs
both libraries. Compilers, source checks and conformance jobs finish before
measurement. A subsequent external conformance check additionally guarded R's
`exp` and `mean`, with a positive control that detects an R exponential call;
all 131 cases still pass against those same artifacts. Default strict
complete-payload comparison applies; there is no
metadata exception for these named benchmark inputs. Corrected formerly failing
inputs have their separate independent checks above.

| Measurement | `e255eea` | Candidate | Before/candidate |
| --- | ---: | ---: | ---: |
| Native total process time, 11 contrasts | 36.00665 s | 35.43496 s | 1.01613× |
| Native maximum peak RSS | 301.93750 MiB | 302.01563 MiB | 0.99974× |
| R total process time, 20 workflows | 89.57854 s | 88.31486 s | 1.01431× |
| R maximum peak RSS | 1117.79688 MiB | 1127.01563 MiB | 0.99182× |
| R total process time, previous 18 controls | 79.36027 s | 78.47536 s | 1.01128× |

Native throughput improves **1.61%**, below the required 2% optimization-round
threshold; **no round is counted**. Native elapsed time falls 1.59%. Original
DESeq2 takes 188.92884 s, a **5.332×** aggregate DESeq2/Rust time ratio, with
1877.18750 MiB maximum peak RSS. The twenty R workflows take 1.41% less total
time; the previous eighteen controls take 1.12% less. These are single paired
observations rather than significance estimates. The pass was not repeated.

Selected process and stage observations expose both benefits and costs:

| R workflow | Process s before → candidate | Peak RSS MiB before → candidate | Target stage s before → candidate |
| --- | ---: | ---: | ---: |
| SF offsets, 56,200 × 16 | 5.35301 → 4.91086 | 902.65625 → 899.25000 | Fit 2.356 → 2.252 |
| Matrix offsets, 56,200 × 16 | 4.86526 → 4.92864 | 934.75000 → 938.28125 | Fit 2.266 → 2.244 |
| Batched results, 56,200 × 16 | 6.48505 → 6.14992 | 1035.06250 → 1125.03125 | Ten results 1.676 → 1.787 |
| Two-group results, 56,200 × 16 | 4.87337 → 4.76469 | 954.89063 → 1005.40625 | Ten results 1.500 → 1.370 |
| Normal shrinkage, 1,200 × 16 | 2.19003 → 2.15955 | 854.85938 → 862.14063 | Shrinkage 0.115 → 0.115 |
| rlog, 2,048 × 64 | 5.30020 → 5.26563 | 685.89063 → 746.92188 | rlog 3.110 → 3.122 |
| rlog, 4,096 × 96 | 21.73363 → 21.40925 | 742.21875 → 742.93750 | rlog 19.099 → 19.004 |

All twenty complete serialized payloads are **byte-identical**, including fitted
state and metadata. Every workflow remains below the registered 2 GiB RSS budget.
Memory is not uniformly improved: maximum R RSS increases **9.21875 MiB**;
the large batch case increases **89.96875 MiB**, the two-group case **50.51563
MiB**, and the 64-sample rlog control **61.03125 MiB**. Their causes have not
been established. These observed costs stay open for the final combined memory
audit, alongside preceding checkpoints' costs. Native maximum RSS rises 0.07813
MiB. No memory-neutrality claim is made.

The [native timings](data/row_preparation_2026-10-30_native_timings.tsv),
[native comparisons](data/row_preparation_2026-10-30_native_summary.tsv),
[R process metrics](data/row_preparation_2026-10-30_r_workflows.tsv),
[R stage metrics](data/row_preparation_2026-10-30_r_stages.tsv),
[complete payload identities](data/row_preparation_2026-10-30_r_payload_identity.tsv)
and [validation record](data/row_preparation_2026-10-30_validation.json) preserve
the observations and artifact identities. Local immutable artifacts use
`row-preparation-*` under `results/codebase_improvements`. Existing checkpoint
measurements/costs remain preserved, and the main README figures stay unchanged.

## Remaining work

Formula/model-matrix and expanded-contrast arithmetic still execute in R and
require standalone paths. Optional ashr/apeglm and posterior/s-value calculations,
unused initial native row adjustment, enabled-but-no-refit replacement allocation
and the final combined performance/source audit remain open. This checkpoint
does not claim the entire R wrapper is already free of statistical calculations.
