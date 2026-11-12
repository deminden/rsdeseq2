# Codebase review and implementation plan — 2026-10-22

Reviewed commit: `8d4b1e43a68d5aeab0e19536404531c3e2a96333`.

The [detailed implementation plan](codebase-improvement-plan-2026-10-22.md)
expands this review into API decisions, commit boundaries, dependencies,
regression cases and precision/speed/memory acceptance gates.

The next work should prioritize result semantics, normalization provenance,
object validity and execution controls. The review found reproducible workflow
defects with much larger effects than the published numerical error bounds.
There is also a substantial memory opportunity in GLM output retention.

This was a deep review of the numerical pipeline, matrix/design validation,
R/C/Rust bridge, result and transform wrappers, CLI writers, benchmark tools
and release workflows. Small R reproductions were run against the previously
validated local package. The installed `DESeq`, `results` and transform bodies
were checked against the checkout and matched. The current vendored core also
matches the primary core. No implementation was changed and no complete test
suite, timing pass or memory benchmark was rerun for this review.

## 1. Confirmed correctness findings

### F1 — P1: Result selection changes Cook's filtering

Location: [R/results.R](../../r-pkg/rsdeseq2/R/results.R), lines 39–88;
[R/results-fits.R](../../r-pkg/rsdeseq2/R/results-fits.R).

The default fitted-object route starts from the native result table, which
already has Cook's masking. Selecting a coefficient or contrast constructs a
primitive fit and recomputes statistics without carrying the same Cook's
policy through that path.

On the existing 40-gene fixture with `counts[1, 8] <- 2000`, four samples per
group and replacement disabled:

| Request | Gene 1 p-value |
| --- | ---: |
| `results(fit, independentFiltering=FALSE)` | Missing (`NaN`) |
| Same coefficient via `name="conditionB"` | `0.001479817` |
| Default request plus `cooksCutoff=FALSE` | Still missing |

The last request is silently ignored because `results()` accepts `...` without
validating or consuming it. An invented `totallyUnknownArgument=TRUE` was also
accepted. This affects ordinary analysis decisions and calls made by shrinkage
wrappers, not just displayed metadata.

Implementation:

1. Validate result options once, before dispatch. Add explicit supported Cook's
   controls and reject every unused argument with its name.
2. Retain unmasked Wald/LRT statistics separately from presentation results.
   In particular, preserve raw LRT p-values so disabling Cook's filtering can
   recover them.
3. Resolve the selection, apply contrast-specific all-zero behavior, apply the
   appropriate Cook's policy and low-count exception, then perform independent
   filtering/BH. Every fitted-object selector must use this sequence.
4. Keep existing public result shapes where practical. Statistical equality of
   equivalent requests must not depend on matching all optional metadata columns.

Acceptance: default, explicit name, numeric and factor contrasts for equivalent
comparisons agree in all six result fields and missingness under the same
options. Cover Wald, LRT, an actual Cook's outlier, replacement, a low-count
exception, thresholds and changed alpha. Compare against pinned DESeq2 for
supported semantics; explicitly reject unsupported requests.

### F2 — P1: Formula offsets have no stable normalization provenance

Locations: [R/DESeq.R](../../r-pkg/rsdeseq2/R/DESeq.R), lines 82–107 and 271–274;
[R/transforms.R](../../r-pkg/rsdeseq2/R/transforms.R), lines 97–129;
[builder/config.rs](../../crates/rsdeseq2/src/core/implementation/builder/config.rs),
`with_formula_offsets`.

Two failures were reproduced:

* Refitting the returned dataset with the same formula and inputs multiplies
  formula offsets into normalization factors again. With exposures from 0.8
  to 1.2, the second fit's factors divided by the first fit's factors equal
  those exposures. Maximum absolute LFC change was **0.3901497**.
* Supplying the same named normalization factors in reversed sample order
  changes effective normalization and results. Multiplication by formula
  offsets occurs before alignment. Maximum effective-factor difference was
  **0.56** and maximum absolute LFC difference was **0.7733569**.

The transform code has the same multiply-again pattern. Source inspection also
shows different no-normalization-factor semantics: the Rust formula builder
multiplies supplied/estimated size factors by offset multipliers; the R wrapper
creates factors from offsets alone, which then preempt size factors.

Implementation:

1. Introduce one normalization preparation function shared by R fitting and
   transformation. Align every labeled input before arithmetic.
2. Represent base size/normalization factors, formula log offsets and effective
   fitting factors separately. Store provenance when effective factors are
   written back to a dataset.
3. Define explicit-argument precedence, reuse of a previously fitted dataset,
   and the behavior when a user changes its factors or formula. Avoid recovering
   base factors by repeated division of rounded effective factors.
4. Resolve the R/Rust size-factor composition discrepancy with an explicit
   documented contract and matching cross-interface examples.

Acceptance: repeated fitting is idempotent; named input permutations do not
alter results; explicit base factors plus offsets equal manually composed
factors; fit-to-transform reuse applies offsets once; changed formulas cannot
reuse stale offset state. Test non-unit factors and offsets together.

### F3 — P1: Dataset subsetting leaves cached fit state stale

Locations: [R/results.R](../../r-pkg/rsdeseq2/R/results.R), lines 19–37;
[R/lfcShrink.R](../../r-pkg/rsdeseq2/R/lfcShrink.R), stored-fit retrieval and alignment;
[R/DESeq.R](../../r-pkg/rsdeseq2/R/DESeq.R), lines 268–270.

`DESeqDataSet` metadata holds the full `rsdeseq2RustFit`. Ordinary dataset
subsetting updates assays and row metadata but leaves that cached list intact.
The results wrapper uses the cached list without validating it against the
current dataset.

Reproduction: `fit[c(3, 1), ]` has two rows, `gene3,gene1`, but `results()` on it
returned **40 rows starting gene1,gene2,gene3**. Same-size permutations and edits
are particularly concerning because dimension-only checks cannot catch them.
Shrinkage checks some dimensions but also combines current counts with stored
designs, dispersions and other arrays; it needs the same validity policy.

Implementation:

1. Centralize extraction of a validated fitted dataset view for results,
   resultsNames and every shrinkage method.
2. Validate gene/sample identities, order, design, normalization provenance and
   relevant fitted inputs. Reject unsupported mutations explicitly.
3. Support gene subsetting/reordering through one index map applied to every
   gene-shaped field, including covariance and replacement metadata. Recompute
   adjustment/filtering over the requested result population as the reference
   does. A minimal first fix may reject changed objects until this is supported.
4. Handle sample edits separately: a new design or count sample set requires
   refitting. Preserve fast access for an unchanged object and measure the cost
   of validity checks.

Acceptance: row subset and permutation either produce correctly aligned results
or a specific error; same-size count/design/sample modifications cannot silently
reuse incompatible state. Check normal, ashr and apeglm consumers as applicable.

### F4 — P1 hardening: The native boundary does not enforce buffer lengths

Locations: [entrypoint.c](../../r-pkg/rsdeseq2/src/entrypoint.c), lines 47–70 and
168–175; [Rust bridge lib.rs](../../r-pkg/rsdeseq2/src/rust/src/lib.rs), lines 52–79.

Several registered `.Call` handlers check the count matrix but pass other R
buffers to Rust without checking their lengths. For example, normalized-count
handling passes `REAL(size_factors)` and Rust constructs a slice of
`n_samples` elements. A short numeric vector can therefore create an invalid
Rust slice. High-level R validation covers many normal calls, but the native
contract is still unchecked. `catch_unwind` cannot protect against undefined
behavior from a bad buffer.

This is a source-confirmed boundary hazard; deliberately unsafe calls were not
executed during this review.

Implementation: add shared C validators for SEXP types, dimensions, `XLENGTH`,
scalar lengths, integer ranges and checked dimension products before allocations
or pointer extraction. Audit every registered entry point, including optional
buffers. Reject R integers that overflow during `as.integer` conversion before
casting to `size_t`. Longer term, use length-bearing buffer descriptors and a
generated/shared ABI declaration instead of repeated long prototypes.

Acceptance: malformed native inputs produce ordinary R errors in isolated child
processes; valid bridge cases preserve exact outputs. Run appropriate memory
instrumentation on the boundary tests where supported.

### F5 — P2: Rank decisions depend on covariate units

Locations: [math/qr.rs](../../crates/rsdeseq2/src/math/qr.rs), lines 3–5 and 52–62;
[design basics](../../crates/rsdeseq2/src/design/implementation/basics.rs),
`validate_full_rank_with_tolerance`.

The default rank check uses an absolute `1e-10` pivot cutoff. For
`cbind(Intercept=1, x=seq(-3.5,3.5)*1e-12)`, R's `qr()` reports rank 2 while
Rust rejects the fit as rank 1 and labels the nonzero covariate a zero column.

Implementation: separate the public explicitly absolute-tolerance primitive
from a scale-aware default design-validity policy. Establish relative rank
decisions against R on exact aliases, near aliases and mixed-scale columns.
Do not silently rescale the fitted design or change coefficient/prior units.

Acceptance: equivalent unit changes do not cause false zero-column rejection;
true singularities remain rejected; threshold cases have documented behavior.
Keep fitting stability diagnostics separate from structural rank validity.

### F6 — P2: LRT validation does not check model nesting

Location: [core helpers](../../crates/rsdeseq2/src/core/implementation/helpers.rs),
`validate_lrt_pipeline_input`, lines 26–80.

Validation checks sample counts, coefficient counts and each model's rank, but
not whether the reduced design lies in the full design's column space.
A full design `[Intercept, group]` and a one-column alternating `[1,2,...]`
reduced design were accepted: combined rank 3 versus full rank 2. The workflow
reported df=1 and 40 finite p-values. Those are not a valid nested-model LRT.
This finding concerns model validity; no claim is made that DESeq2 rejects
every equivalent arbitrary-matrix request.

Implementation: validate reduced-column membership using the same scale-aware
rank policy as F5 before fitting. Accept a different basis for the same nested
subspace. Make the invalid-request behavior explicit in the API.

Acceptance: nested full/reduced pairs and reparameterized reduced bases pass;
non-nested pairs fail before optimization; valid baseline LRT values stay exact.

## 2. Execution, performance and maintainability

### F7 — P2: Requested threads and execution modes do not control execution

Locations: [builder/config.rs](../../crates/rsdeseq2/src/core/implementation/builder/config.rs),
lines 125–134 and getters;
[IRLS](../../crates/rsdeseq2/src/glm/irls.rs), line 253;
[thread invariance test](../../crates/rsdeseq2/tests/thread_invariance.rs), lines 118–138.

`threads()` stores a field; no fitting code consumes it to create/install a
pool. Parallel loops use the ambient Rayon pool. The R `nproc` argument reaches
this unused setter. The existing thread-invariance test constructs its own
pool, so it verifies determinism but cannot verify that the public setting works.
`ExecutionMode` is likewise stored and read by its getter, without distinct
execution behavior.

Implement a reusable execution context. An explicit worker count must apply to
the whole analysis; an unspecified count can inherit the caller's context.
Reuse the context across stages/refits and avoid nested pool creation or global
environment changes. Document modes as aliases until they have defined distinct
semantics, or reject unavailable modes. Add an observable effective worker count.

Acceptance: exercise the public API in a conflicting ambient pool and verify
the requested count is used inside fitting. Preserve exact outputs at 1 and 4
workers. The published benchmark already sets `RAYON_NUM_THREADS` explicitly;
this finding alone does not invalidate its one-thread measurements.

### F8 — P2: Retained GLM outputs create an avoidable quadratic memory cost

Locations: [IRLS](../../crates/rsdeseq2/src/glm/irls.rs), lines 240–283;
[rlog](../../crates/rsdeseq2/src/transform/rlog.rs), lines 194–206;
[fixed LRT](../../crates/rsdeseq2/src/core/implementation/builder/fixed_lrt.rs),
lines 40–44; [R fit bridge](../../r-pkg/rsdeseq2/src/rust/src/lib/fit.rs), line 130.

Every IRLS fit builds per-gene row outputs, flattens them into global arrays,
and clones beta/covariance arrays to inspect fallback routing. rlog uses
`samples+1` coefficients, retains full covariance through this process, and
then drops the fit when returning the high-level transform. Dispersion mean
fits also construct a full fit when their immediate consumer only needs mu.

One covariance array for 56,200 genes and 61 samples is approximately **1.61
GiB** (`56200 * 62^2 * 8` bytes). For 20,000 genes and 200 samples it is **6.02
GiB**. These are calculated buffer sizes, not measured peak RSS. Clones and
other matrices add to the live allocation. LRT attachments clone reduced-model
fields they can own, and the R replacement path clones the original full fit.

Implement in increasing-risk steps:

1. Borrow arrays for fallback inspection; move reduced-fit fields into their
   destination; avoid cloning the original fit when output ownership permits.
2. Use bounded ordered row batches or direct disjoint output slices to avoid
   retaining two complete representations. Keep arithmetic within each gene and
   output ordering unchanged; preserve deterministic error selection too.
3. Introduce internal output requirements for full inference, mean fitting and
   transforms. Retain all covariance/variance checks needed for identical
   fallback decisions, even when the full covariance is not kept globally.
4. Only then consider omitting proven-unneeded calculations. Do not skip
   covariance computation indiscriminately: it currently influences fallback.
5. Avoid allocating/copying a replacement-count R matrix when replacement is
   disabled or absent, after defining a safe bridge output contract.

Acceptance: exact coefficients, means, fallback rows, iteration/convergence
fields and result exports on existing fixtures; new measured memory scaling
for rlog and larger designs; one paired eleven-contrast timing/RSS pass for
quantification changes. No speedup is claimed by this review.

### F9 — P2: R and Rust maintain separate statistical result engines

Locations: [R/results-statistics.R](../../r-pkg/rsdeseq2/R/results-statistics.R),
[R/results-filtering.R](../../r-pkg/rsdeseq2/R/results-filtering.R),
[glm/wald.rs](../../crates/rsdeseq2/src/glm/wald.rs), and
[independent_filtering.rs](../../crates/rsdeseq2/src/independent_filtering.rs).

The R layer independently implements contrast SEs, Wald alternatives, tails,
BH and lowess-based filtering. Its Wald path loops over genes and creates tail
closures per gene. The Rust fit computes filtering and the default R results
path computes it again. This increases maintenance and wrapper-specific cost.
The existing CLI benchmark does not measure `DESeq -> results -> lfcShrink`
inside an R session.

After fixing F1, introduce a common Rust result request/evaluation API and route
R through a batch bridge. Migrate one primitive at a time, retaining reference
oracles outside production. Keep Bioconductor classes/names/metadata in R.
Measure the end-to-end R workflow and per-stage cost before selecting the next
performance change. Preserve R-tail numerical tolerances during migration.

There are currently 209 public builder methods, and the file split mostly uses
`include!` into shared private namespaces. Keep compatibility forwarding methods
while converging execution on validated request types and a small number of
pipeline functions. Introduce real module boundaries incrementally where this
removes duplicate decisions; a wholesale parser or API rewrite is unnecessary.

## 3. Release and validation improvements

### F10 — P2: Release publication is not gated by combined validation

Locations: [R CI](../../.github/workflows/r.yml),
[crate publication](../../.github/workflows/publish-crates.yml),
[binary release](../../.github/workflows/release-binaries.yml),
[vendor synchronization](../../scripts/release/sync_r_vendor.sh).

Observed workflow gaps:

* R CI synchronizes the vendor before testing, without failing on a resulting
  source diff. A forgotten committed synchronization can be hidden by CI.
* The synchronization script duplicates versions and dependency specifications
  in a hard-coded manifest.
* The crate publish job runs Rust checks but does not depend on R validation.
  Binary publication depends on builds, with no test step on each target.
* Normal Rust CI uses stable Ubuntu only. The advertised minimum compiler and
  shipped macOS/Windows targets lack corresponding routine numerical checks.
* R is pinned, but dependency installation is not a pinned oracle environment.

Implement a read-only vendor/manifest check, derive the vendor manifest from
workspace metadata, and run package checks on the exact staged source artifact.
Use a shared validation workflow that all publication jobs depend on. Add an
MSRV build and a compact numerical smoke suite for shipped targets. Keep a
pinned reference environment separate from a current compatibility lane.

Acceptance: deliberate vendor/version drift blocks CI, an R failure blocks
publication, and the tested source/binary hashes identify the published artifacts.
The current vendor is synchronized; this is a future-regression prevention gap.

### F11 — P2: Broaden workflow coverage and make benchmark success enforce accuracy

Locations: [benchmark_gtex_native_wald.py](../../scripts/benchmarks/benchmark_gtex_native_wald.py),
lines 118–129; [contrast generator](../../scripts/references/generate_gtex_mode_contrasts.py);
[GLM benchmark](../../crates/rsdeseq2/benches/glm_gene_bench.rs).

The speed script records numerical differences but only rejects finite/NA
mismatches unless an optional exact-output directory is supplied. A finite,
incorrect result can finish the standalone benchmark successfully. The separate
precision checker is stronger; make that quality gate part of the benchmark's
accepted-run contract.

The eleven contrast studies vary designs within the same Kidney Cortex cohort.
They are useful regression evidence. Add a compact, independently selected
validation panel spanning small residual df, larger designs, heterogeneous
weights, actual Cook's replacement, local trends and LRT. Prefer one held-out
dataset and targeted synthetic boundary cases over more random splits of the
same cohort. Keep optimization development and final hold-out evaluation distinct.

The existing Criterion `glm_placeholder` benchmark measures a four-observation
likelihood call, and the pipeline benchmark measures normalization only. Replace
those labels/workloads with representative IRLS, dispersion, result extraction
and transform benchmarks. Use diagnostic stage timings and allocation evidence
to choose optimization candidates. Keep one paired full pass per candidate, as
previously requested; broader workflow coverage does not require repeated passes.

Add metamorphic checks for the relationships missed here: equivalent selectors,
labeled permutations, fit reuse and supported subsetting. Keep direct DESeq2
oracles because cross-interface equality alone can reproduce a shared bug.

### Smaller follow-ups

* CLI output paths should be checked for collisions, including input/output
  aliases. Writers currently publish sidecars before the main table, and direct
  file writes can leave partial/mixed artifacts on failure. Stage validated
  outputs and publish a completion manifest; do not claim multiple independent
  renames are a transaction.
* Make CLI dependencies optional for library/R-only builds after measuring build
  cost. Keep source/package synchronization deterministic.
* Replace stale “future”, “placeholder” and earlier-stage API descriptions with
  current scope. Update the R package's placeholder maintainer contact before
  broader distribution.
* Add cancellation/progress at safe orchestration boundaries for long R fits.
  R API calls must stay on the R thread; an interrupt must not longjmp through
  live Rust frames. Design that boundary deliberately.
* Finish the small diagnostic/sidecar assertion improvements identified in the
  earlier test review when changing those areas. Broad further test pruning has
  lower value than the workflow contracts above.

## 4. Implementation sequence and completion gates

Each row is a reviewable change or tightly related pair of changes. Reproduction
cases should become focused regression tests alongside the relevant fix.

| Order | Change | Dependency | Completion gate |
| --- | --- | --- | --- |
| 1 | Reject ignored result arguments; unify Cook's/all-zero/filtering semantics and retain raw statistics | None | F1 examples corrected; pinned Wald/LRT result parity and selector equivalence |
| 2 | Align and compose normalization/offsets once; record provenance | None | F2 idempotence, permutation and explicit-composition checks |
| 3 | Validate cached fitted objects and map supported gene subsets | Normalization provenance from 2 | F3 cases aligned or explicitly rejected; shrinkage consumers checked |
| 4 | Validate every C/Rust boundary buffer and integer conversion | None | Isolated malformed calls error cleanly; valid bridge outputs unchanged |
| 5 | Apply requested workers using one execution context; clarify modes | None | Public API actually honors worker count; 1/4-worker bitwise parity |
| 6 | Establish scale-aware model validity and reject non-nested LRTs | Rank policy before nesting guard | R rank oracles, valid nested bases and invalid pairs; valid pipeline outputs preserved |
| 7 | Gate vendor/package/release identity and enforce benchmark precision | None; complete before release | Intentional drift/failure blocks publication or benchmark acceptance |
| 8 | Remove avoidable fit/covariance copies, then bound row-output retention | Existing fallback invariants | Exact stage/route outputs, measured RSS improvement, no material speed regression |
| 9 | Add output requirements for mean fitting and rlog | 8 | Memory growth bounded by required outputs; transform and fallback parity |
| 10 | Migrate duplicate R result computations to a common batch Rust API | 1 and 4 | All supported R semantics preserved; wrapper timing and precision recorded |

Orders 1–7 are the next correctness/reliability milestone. Memory work should
land separately, followed by statistical-result consolidation. Keep low-level
L-BFGS-B and locfit changes out of this first milestone: the observed defects
arise before or after those algorithms, and changing them would add numerical
risk without fixing these failures.

For each change, preserve the existing frozen references and tolerances. For
refactors and optimizations, require exact baseline outputs on the validated
workflows. For correctness fixes, require the new expected behavior on the
demonstrated failing cases and unchanged behavior elsewhere. Do not require a
known wrong result to remain byte-identical.

Run focused tests first, then the affected integration/reference suites. Before
the next release, run the complete Rust/Python/R package checks on the final
staged artifacts and the existing full GTEx precision checks. For performance
claims, use a single paired eleven-contrast pass with timing and RSS, plus the
specific R/transform workload the change targets. Record actual measurements;
buffer-size calculations and suspected hotspots are not measured gains.

## 5. Minimal reproduction setup

Run from the repository root with the reviewed R package installed and the
pinned DESeq2 reference environment available. The review used R 4.6.1 and
DESeq2 1.52.0. The existing deterministic helper supplies the small count matrix.

```r
library(rsdeseq2)
source("r-pkg/rsdeseq2/tests/testthat/helper-fit-data.R")
ct <- .rsdeseq2_test_counts()
ct[1, 8] <- 2000L
cd <- data.frame(
  condition = factor(rep(c("A", "B"), each = 4)),
  row.names = colnames(ct)
)
dds <- DESeq2::DESeqDataSetFromMatrix(ct, cd, ~condition)
fit <- rsdeseq2::DESeq(dds, fitType="mean", minReplicatesForReplace=Inf)

# F1: same selected coefficient, different missingness.
rsdeseq2::results(fit, independentFiltering=FALSE)$pvalue[1]
rsdeseq2::results(fit, name="conditionB", independentFiltering=FALSE)$pvalue[1]
rsdeseq2::results(fit, cooksCutoff=FALSE, independentFiltering=FALSE)$pvalue[1]
rsdeseq2::results(fit, totallyUnknownArgument=TRUE)

# F3: two input genes versus forty output genes.
rownames(fit[c(3, 1), ])
dim(rsdeseq2::results(fit[c(3, 1), ], independentFiltering=FALSE))

# F2: offsets compound across refits and precede sample alignment.
ct <- .rsdeseq2_test_counts()
cd$exposure <- seq(.8, 1.2, length.out=8)
dds <- DESeq2::DESeqDataSetFromMatrix(ct, cd, ~condition+offset(log(exposure)))
f1 <- rsdeseq2::DESeq(dds, fitType="mean", minReplicatesForReplace=Inf)
f2 <- rsdeseq2::DESeq(f1, fitType="mean", minReplicatesForReplace=Inf)
DESeq2::normalizationFactors(f2)[1, ] / DESeq2::normalizationFactors(f1)[1, ]
max(abs(rsdeseq2::results(f1)$log2FoldChange-rsdeseq2::results(f2)$log2FoldChange))
nf <- matrix(rep(seq(.7, 1.4, length.out=8), each=nrow(ct)),
             nrow=nrow(ct), dimnames=dimnames(ct))
a <- rsdeseq2::DESeq(dds, fitType="mean", normalizationFactors=nf,
                    minReplicatesForReplace=Inf)
b <- rsdeseq2::DESeq(dds, fitType="mean", normalizationFactors=nf[, 8:1],
                    minReplicatesForReplace=Inf)
max(abs(rsdeseq2::results(a)$log2FoldChange-rsdeseq2::results(b)$log2FoldChange))

# F5: validation rejects a full-rank, small-scale covariate.
dm <- cbind(Intercept=1, x=seq(-3.5, 3.5)*1e-12)
qr(dm)$rank
try(rsdeseq2::DESeq(ct, design=dm, fitType="mean", minReplicatesForReplace=Inf))

# F6: non-nested designs still produce LRT p-values.
full <- cbind(Intercept=1, conditionB=rep(c(0, 1), each=4))
reduced <- matrix(rep(c(1, 2), 4), ncol=1)
c(full=qr(full)$rank, reduced=qr(reduced)$rank,
  combined=qr(cbind(full, reduced))$rank)
bad_lrt <- rsdeseq2::DESeq(ct, design=full, reduced=reduced, test="LRT",
                          fitType="mean", minReplicatesForReplace=Inf)
bad_lrt$lrtDf
sum(is.finite(rsdeseq2::results(bad_lrt)$pvalue))
```
