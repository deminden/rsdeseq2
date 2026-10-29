# Codebase improvement execution — 2026-10-29

Implementation follows the [detailed plan](codebase-improvement-plan-2026-10-22.md).
The full plan remains in progress. This document records the first correctness
batch; later execution, delivery, memory and consolidation work is still pending.

## Implemented in the first batch

- Shared C validation checks native buffer types, lengths, matrix dimensions,
  optional inputs, scalar/index ranges and allocation products before dispatch.
  The existing logical size-factor method flag remains supported.
- R validates worker/replacement integer upper bounds before integer conversion.
- Native fits retain unmasked Wald/LRT statistics, including refitted rows and
  the special p-value of one for rows made all-zero by replacement.
- Fitted default/name/contrast results apply Cook's policy before their final
  independent filtering/BH step. Cook's disabling and numeric cutoffs are explicit;
  unused result arguments error. Existing selected/default column layouts remain.
- Shared normalization preparation aligns factors before composition, stores base
  factors separately from effective factors, and applies offsets once on refit.
  Transforms use that preparation, including original exposure normalization for
  blind fitting. Ambiguous older offset fits require original/explicit base input.
- Refitting removes stale package-owned replacement assays and coefficient fields.
- Results, result names and all shrinkage methods obtain a validated fitted view.
  Gene subsets/reordering map cached arrays; changes to counts, samples, design,
  offsets, normalization, weights or checked fitted quantities require refitting.
  Result population adjustments are recomputed after selecting genes.

## Validation recorded

- Complete R tests from the built source package: **1,293 passing assertions,
  zero failures, zero warnings, zero skips**.
- New targeted cases cover child-process malformed native calls, Wald/LRT Cook's
  masking/disabling, direct DESeq2 p-value checks, offset refit idempotence,
  shuffled normalization columns, explicit factor composition, subset BH and
  stale dataset rejection.
- R Rust bridge formatting and all-target Clippy with warnings denied pass.
- The initial source-package check reported one duplicated Rd argument warning;
  its tests passed. The duplicate was corrected and package checks were repeated
  without rerunning the unchanged numerical tests: **Status: OK**. All 157
  runtime source files are hash-identical between those checked packages.
- The preserved CLI executable still has its original hash. No release/version
  bump or publication has been performed.
- Local raw evidence and the preserved pre-change executable/manifest are under
  `results/codebase_improvements/`. These generated files remain outside Git.

No new full-GTEx timing/RSS result is claimed. The primary Rust core and CLI
source have not changed in this batch. R validation overhead and the broader
workflow/transform performance panel still need measurement.

## Remaining work

M1 follow-through includes the larger selector/replacement/threshold matrix,
legacy save/load cases, explicit subset-shrinkage reference cases and broader
transform reuse checks. Existing integration coverage is retained.

Profiling exposed and corrected another native-fit reuse case: rlog no longer
requires a DESeq2 dispersion-function attribute when its fitted trend type is
already retained. Design-aware transforms validate native cached state, including
stored dispFit/base summaries, before reuse. The existing large transform fixture
now checks native-fit rlog against the existing mean-trend DESeq2 reference and
rejects edited counts/dispFit without repeating reference setup.

C09–C10 are implemented: default rank validation scales temporary columns while
preserving original fitting units and explicit absolute-tolerance calls. LRT
entry paths check reduced-model subspace membership before dispersion fitting.
New tests cover mixed units, aliases, alternate nested bases and early invalid
model rejection. The existing Rust suite passed 1,016 tests with one pre-existing
ignored test; the three new cases passed separately. Workspace Clippy passed.
The R suite also passed with the synchronized core. Full GTEx replay passed all
121 per-contrast and nine global gates; all 22 exports matched saved v0.3.1
result/diagnostic files byte for byte. This is precision evidence, not a speed run.

C11 is implemented: a private execution context shares a local pool across builder
clones/stages/refits. Public entry routes install that context, and frozen R rlog
calls use it too. Unspecified workers inherit the caller's context. A test
observes the requested worker count inside fitting in a conflicting ambient pool,
checks clone reuse and rejects zero/unsupported counts. Existing one/four-worker
numerical invariance passes. The workspace passed 1,020 tests with the existing
single ignored case; R tests and core/backend Clippy pass.

C12 is implemented: read-only vendor/version checks derive the manifest from
Cargo metadata, and Rust/R CI reject committed drift. Accepted benchmarks require
a native baseline verified against reviewed output hashes and frozen reference
identity; exploration runs remain explicitly unaccepted. Baseline/binary/manifest
mutation and empty finite coverage fail acceptance. The Python suite passed
37 tests, including drift, version mismatch and benchmark failure cases.

C13 is implemented: tag publication is coordinated by one validation workflow;
Rust, pinned R, full GTEx and MSRV checks precede checked artifact staging.
Publication verifies source/artifact identity, and target builds run numerical
and CLI tests before packaging. The minimum Rust 1.97.1 build passes locally.
Workflow lint passes except for the old linter's unrecognized pre-existing
`macos-15-intel` runner label, which is excluded explicitly. Dependency failure
simulation confirms that Rust/R/parity/MSRV failures block publication. Remote
platform jobs have not been executed by this local implementation.

At the user's request, current comparisons moved to `deseq2-r462-native`:
R 4.6.2 with DESeq2 1.52.0, apeglm 1.34.0 and ashr 2.2-63. All 82 package sources
were retrieved and hashed, their installed versions were verified, and rsdeseq2
was rebuilt against that environment. The complete R suite and source-package
check pass there (Status: OK). Historical reference files remain unchanged.
The new source installer also restored an isolated Rcpp/RcppProgress smoke
environment successfully. The Python suite passes 42 tests.

C14 now has representative stage workloads, preserved executable/R package
baselines and initial native-R workflow/RSS profiles. Their raw output is saved
under `results/codebase_improvements/`; the [memory report](memory-improvements-2026-10-29.md)
records the dimensions, budget and measured baseline. Benchmark smoke execution,
Clippy and MSRV checks pass. An independent Pasilla panel and native R Wald/LRT
references were frozen before memory optimization and are reserved for milestone
validation. No speed/memory improvement is claimed yet.

C15–C17 and the disabled-buffer part of C18 are implemented. Consumed fit arrays
move into their destination, temporary gene outputs are bounded and ordered, and
internal mean/rlog consumers discard global covariance after identical checks.
Thirteen paired R workflows retain exactly the same saved values and attributes.
At 4,096 × 96, rlog peak RSS falls 55.31% with similar elapsed time; all measured
processes stay within 2 GiB. The independent Pasilla Wald/LRT panel passes its
unchanged bounds and all eight exports match the baseline exactly. Rust passes
1,022 tests with one existing ignored case. See the memory report for the full
panel and its limitations. Enabled replacement keeps the existing buffer contract;
a complex opaque result handle is deferred unless a larger measured case earns it.

C22 is implemented as a private output-transaction module with an explicit path
inventory. Existing file identity and filesystem name behavior are checked before
fitting; destination-local temporary files are closed and hashed before any final
file is published. A completion record is published last. Separate path renames
remain non-atomic, and readers must verify the record for coherent generations.
Tests cover hard links/symlinks, case aliases, missing destinations, sidecar
failure, unwritten staging and injected publication failure. Existing CLI numeric
comparisons remain exact. Workspace tests, Clippy, MSRV and 43 Python tests pass.

C19–C21 and C23 remain pending: shared Rust statistical requests, bulk R results,
module consolidation and final release documentation/evidence.
