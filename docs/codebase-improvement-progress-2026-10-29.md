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

C12–C13 remain pending: accuracy-enforced benchmark acceptance, deterministic
vendor checks and combined artifact/publication gates.

C14–C23 remain pending: representative stage/RSS profiling, borrowed ownership,
bounded row storage, selective GLM retention, optional replacement buffers,
shared Rust statistical requests, bulk R results, module consolidation, staged
CLI output and final release documentation/evidence.
