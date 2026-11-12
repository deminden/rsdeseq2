# Detailed codebase improvement plan — 2026-10-22

Status: C01–C23 accepted for the documented development scope; see the
[final validation](dev-validation-2026-11-12.md). The checklist below preserves
the original implementation decisions and acceptance criteria. Publication is
separate and remains on hold.
Base commit: `8d4b1e43a68d5aeab0e19536404531c3e2a96333`.

This expands the findings in the [codebase review](codebase-review-2026-10-22.md)
into implementation decisions, commit boundaries, dependencies and acceptance
criteria. Finding identifiers F1–F11 refer to that review. Names for new types,
helpers, files and commands below are proposals, not existing interfaces.

## 1. Direction and priorities

The first milestone should make normal analysis workflows trustworthy: selecting
a result, refitting a dataset, supplying labeled factors, subsetting genes and
passing input through the native bridge. These have demonstrated failures or a
source-confirmed safety defect. Fixing them has higher immediate value than
changing optimizer tolerances or adding more splits of the same GTEx cohort.

The next performance work should reduce live matrix storage, then remove repeated
result computations between R and Rust. There is concrete source evidence for
these opportunities. Their actual speed and RSS gains still need measurement.

Use four milestones:

| Milestone | Main deliverable | Included work | Exit condition |
| --- | --- | --- | --- |
| M1: trustworthy workflows | Consistent results and reusable fitted objects | C01–C08 | F1–F4 corrected; independent numerical checks and mutation cases pass |
| M2: reliable execution and delivery | Valid model requests, effective thread controls, gated artifacts | C09–C13 | F5–F7 and F10 fixed; benchmark acceptance enforces precision |
| M3: measured memory improvement | Less retained/copy storage with the same fitting decisions | C14–C18 | Exact unaffected outputs, measured RSS improvement, acceptable elapsed time |
| M4: simpler result processing | Shared statistical primitives and a smaller internal execution surface | C19–C23 | R parity preserved; wrapper cost measured; documentation and artifacts agree |

M1 can become a focused correctness release after M2's delivery gates are ready.
It should not wait for every memory or architecture change. Release numbers and
publication are separate decisions after implementation and validation.

### Rules that apply to every milestone

1. Preserve frozen reference files and existing numerical tolerances. Add new
   independent examples when the current references do not cover a failure.
2. Refactors and optimizations preserve outputs, missingness and algorithm routes
   on the unaffected baseline. A demonstrated wrong answer gets a specifically
   documented correction rather than becoming a permanent compatibility rule.
3. Keep fitting arithmetic, thresholds, solver settings and dependency updates in
   separate changes. Avoid mixing them with object handling or memory ownership.
4. Evaluate performance candidates with one paired pass across all eleven GTEx
   contrasts. Additional targeted workloads cover R and transforms, which those
   contrasts do not measure. Repeated whole-study timing passes are not required.
5. A speed optimization counts as a 2% gain only when total quantification
   throughput improves by at least 2% against the immediately preceding accepted
   baseline. Compatible partial gains may be combined before evaluating a round.
6. Record proposed optimizations, outcomes and rejected variants in the existing
   untracked `OPTIMISATION_LOG.md`. Preserve its existing history and keep it
   untracked. Put publishable evidence and reproduction instructions in docs.
7. Do not trade precision for speed. A memory reduction may be accepted without
   reaching the speed-round threshold, but its timing cost must be reported.
8. Main README memory comparisons remain against DESeq2. Development-to-development
   memory comparisons belong in the detailed benchmark documentation.

## 2. C01 — Freeze the implementation and acceptance baseline

### Work

- Record the base commit, tracked working-tree diff, Rust/R package versions,
  dependency lock hashes, compiler, CPU, OS, BLAS and worker settings.
- Preserve the current release executable separately before building candidates.
  Record its hash and that of each tested candidate. Never use a mutable
  `target/release` path as the only identity of a historical baseline.
- Reuse the existing evidence under `results/test_consolidation_followup/`,
  the committed test-consolidation reports and the v0.3.1 output manifest. Record
  which evidence is inherited from earlier validation and which is newly run.
- Preserve the eleven generated GTEx contrasts, their R references, the 44 frozen
  input/reference hashes, and the 22 native result/diagnostic exports.
- Define a small manifest of workflow cases, with the defect, expected behavior,
  reference source, affected interfaces and exact assertions for each case.
- Record current R workflows separately: fit, first results, repeated results,
  multiple contrasts, each shrinkage method, VST and rlog. CLI timings alone
  cannot establish their cost.

The most recent recorded full validation has 1,016 Rust tests, 31 Python tests
and 1,252 R assertions, with the existing slow full GTEx test opt-in. These are
historical evidence counts, not targets to preserve or claims of a new run.
Preserve coverage, not a particular assertion count.

### Baseline classes

| Class | Example | Required treatment |
| --- | --- | --- |
| Unaffected numerical behavior | Existing native Wald exports with Cook's and independent filtering disabled | Exact results and diagnostics on the same platform/build settings |
| Confirmed wrong workflow | Named result bypasses Cook's masking | Correct expected behavior, independently checked against the supported reference semantics |
| Invalid request previously accepted | Non-nested LRT | Specific early error; no fabricated statistical output |
| New supported behavior | Gene subset with valid fit mapping | Explicit contract and independent coverage before enabling |
| Metadata or error-message change | Fit schema version, unsupported argument name | Intentional documented difference without relaxing numerical checks |

### Deliverable and gate

A compact baseline/evidence manifest, with commands pointing to immutable inputs
and binaries. The existing comparator must be invoked with frozen native results
until C12 makes that requirement part of the benchmark interface. No reference
regeneration is part of establishing this baseline.

## 3. Data contracts to establish before broad refactoring

These are internal contracts. Add only fields needed by the changes below and
keep current public constructors and return classes usable.

| Proposed contract | Owns or describes | Important invariant |
| --- | --- | --- |
| `NormalizationState` | Base factors, formula log offsets, effective factors, source and identity metadata | Each offset is composed once after alignment |
| Versioned fitted state | Canonical input identity, raw fit statistics and fit options | Display filtering never destroys raw inference state |
| `ValidatedFitView` | A fit plus validated gene mapping and current-object context | Every array refers to the same genes and fitted samples |
| `ResultRequest` | Selector, test options, Cook's policy, alpha and independent-filtering options | One validation and policy sequence for every entry route |
| `ExecutionContext` | Effective worker policy and reusable pool | An explicit worker request controls the analysis |
| `OutputRequirements` | Which GLM outputs must survive a stage | Retention choices do not change fitting or fallback decisions |

Use existing Rust facilities wherever possible. In particular,
`DeseqFit::results_contrast_request`, its Wald/LRT variants, and the existing
contrast resolution, Wald and filtering helpers already implement substantial
parts of result processing. Extend their options and borrowing model; avoid
creating another independent result engine beside them.

Version fitted state when these contracts become persistent. A version alone is
not a validity check: mutable `DESeqDataSet` content must still be checked against
the fit before use.

## 4. C02 — Harden the R/C/Rust boundary first (F4)

Primary files: [entrypoint.c](../../r-pkg/rsdeseq2/src/entrypoint.c) and
[Rust bridge lib.rs](../../r-pkg/rsdeseq2/src/rust/src/lib.rs).

### Implementation

1. Inventory every registered `.Call` entry and every native buffer argument.
   For each, record SEXP type, permitted dimensions, expected element count,
   whether NULL is allowed, scalar range, index convention and output ownership.
2. Implement shared C validators for numeric/integer vectors, matrices, scalar
   integers, optional buffers and checked dimension multiplication. Check
   `XLENGTH`, not just the presence of a dimension attribute.
3. Validate all inputs before extracting pointers or allocating outputs whose
   size depends on those inputs. Include counts, factors, weights, designs,
   dispersions, priors, frozen transform inputs and shrinkage inputs.
4. Check `n_genes * n_samples`, `n_genes * p * p` and byte-size products against
   all relevant R, C and Rust representable lengths. Do not multiply first and
   check an already-overflowed result.
5. Validate scalar upper bounds in R before `as.integer`, then again at the C
   boundary. Handle NA, Inf, fractional values, negative values and unsupported
   enum codes explicitly. Preserve documented sentinels such as replacement Inf.
6. Audit `PROTECT` lifetimes and cleanup on error. Keep R API calls on the R
   thread and keep R error/interrupt longjmps outside live Rust frames.
7. Document the unsafe Rust slice construction contract beside the helper. A
   panic catcher is not a buffer-length check.

Do not make a wholesale ABI redesign a prerequisite. Length-bearing descriptors
and a shared/generated ABI header can follow when C18/C20 add new bridge calls.

### Regression cases and acceptance

Use a table describing each argument family rather than dozens of nearly
identical tests. Cover a short vector, extra elements where disallowed, wrong
type, inconsistent dimensions, empty scalar, NA scalar, overflow-range scalar,
NULL required buffer and malformed optional buffer. Execute potentially unsafe
`.Call` cases in isolated child R processes with a timeout and explicit expected
exit status. The required outcome is an ordinary R error, never a crash.

Valid wrapper/reference cases must remain exact. Run a supported sanitizer or
memory-check configuration for the boundary suite; report unavailable coverage
explicitly rather than treating a normal run as memory instrumentation.

**Commit boundary:** validators, call-site checks and their tests together.
No solver changes. This is independently useful before the rest of M1.

## 5. C03–C04 — Make result selection and filtering consistent (F1)

Primary files: `R/results.R`, `R/results-fits.R`, `R/results-objects.R`,
`R/results-statistics.R`, `R/results-filtering.R`, and the native fit bridge.

### C03: preserve raw statistics and make options explicit

- Add an explicit `cooksCutoff` argument with supported automatic, disabled and
  numeric-cutoff behavior. Validate options before object-class dispatch.
- Reject every unused named or positional argument in `...`, identifying the
  offending arguments. Audit the primitive constructors for the same issue and
  for early-return paths that may discard explicit overrides.
- Retain unmasked LRT statistics/p-values directly from native LRT state.
  Never reconstruct them from a Cook-masked results table or rounded TSV output.
- Retain enough Cook's and replacement context to implement the policy for any
  supported selector. Clearly distinguish original counts/Cook's data from
  replaced/refitted data, including sample eligibility and low-count exceptions.
- Keep raw coefficient, covariance, degrees-of-freedom and all-zero state
  independent of presentation tables. Wald results can use those primitives.
- Add a fitted-state schema version and a compatibility reader with explicit
  checks for which raw fields are actually available.

For a legacy fit without recoverable raw LRT state, an affected request such as
disabling a previous Cook mask must give a clear refit instruction. Supported
legacy requests may continue only when all required state can be verified.
Do not silently manufacture an unmasked p-value or refit with guessed settings.

Primitive fits constructed without counts/Cook context remain usable. Their
default behavior is documented as no Cook filtering; an explicit automatic or
numeric Cook request requiring missing context errors. An explicit disabled
request remains valid. The request parser must distinguish a missing default
argument from an explicit request.

### C04: route all selectors through one policy pipeline

The common sequence is:

1. Validate the request and obtain a validated fit view.
2. Resolve the default coefficient, name, numeric/list contrast or factor contrast.
3. Calculate the requested effect/SE and unmasked test statistics. For LRT, the
   full-versus-reduced statistic remains the same when only the displayed effect
   changes; reject unsupported combinations of LRT and Wald-specific options.
4. Apply the supported contrast-specific all-zero behavior and invalid-row rules.
5. Apply Cook's masking using the same replacement and low-count rules.
6. Apply independent filtering and BH once for the requested alpha/population.
7. Attach result names, descriptions, metadata and the existing return class.

Keep result-shape compatibility where practical: the default extended table and
the selected-result table need not acquire identical optional columns. Compare
the six statistical fields and their missingness where the statistical request
is equivalent.

Factor contrasts have special all-zero behavior. Do not require arbitrary
factor and numeric requests to match when the reference intentionally treats
them differently. Capture those differences in named test cases.

### Focused regression matrix

| Case | What it must prove |
| --- | --- |
| Two-group Wald, genuine Cook's outlier, replacement disabled | Default, name and equivalent contrast apply the same mask |
| Same fit with `cooksCutoff=FALSE` and with a numeric cutoff | A supported option actually changes the intended rows |
| Eligible replacement fixture | Masking uses the correct replacement/refit policy |
| Low-count Cook's exception | Selection routes preserve the exception |
| Factor levels with both compared groups zero | Reference-specific all-zero behavior survives dispatch |
| Nonzero LFC threshold, directional alternatives and t inference | Effects, SEs, p-values and degrees of freedom use the requested policy |
| LRT with different displayed coefficients | LRT statistics remain raw and selector-independent where intended |
| Two alpha values and filtering on/off | BH/filter population and metadata are recomputed correctly |
| Unknown argument and unavailable context | Actionable early error instead of ignored input |

Use a few deliberate fixtures and table-driven requests. Fit once per fixture
where reuse is the behavior being tested. Compare supported cases with pinned
DESeq2, including exact finite/NA masks. Check the `lfcShrink` routes that call
named results. Metadata assertions should verify meaning, not entire incidental
attribute dumps.

**Exit:** the reproduced default/name/Cook-disabled failure is corrected; raw
LRT recovery works; the unchanged CLI baseline retains its outputs.

## 6. C05–C06 — Make normalization composable and refitting idempotent (F2)

Primary files: `R/DESeq.R`, `R/transforms.R`, normalization helpers and Rust
`builder/config.rs::with_formula_offsets`.

### C05: introduce a shared normalization preparation step

Separate three things: base normalization, the formula log-offset vector and
effective normalization passed to fitting. Align all gene/sample labels first.

For gene g and sample s, the proposed contract is:

```text
base[g,s] = supplied base normalizationFactors[g,s], if present
            otherwise supplied/estimated sizeFactors[s]
effective[g,s] = base[g,s] * exp(formula_log_offset[s])
```

Absent formula offsets mean multiplication by one. Base normalization factors
preempt size factors, preserving the documented factor precedence. Supplied
factors retain their values; composition must not quietly re-center or rescale
the already supplied effective factors.

Use this argument/source precedence:

1. An explicit argument overrides the corresponding object slot. Distinguish
   an omitted argument from an explicitly supplied NULL when clearing a slot.
2. If an object has valid normalization provenance and its current factors equal
   the last factors written by rsdeseq2, reuse the retained base inputs.
3. Changed external factors become new base input after validation. Record the
   new source rather than pretending the previous provenance still applies.
4. Without stored base factors, use the object's supplied factors. Estimate size
   factors only if neither effective base normalization route is supplied.

An important exception is an older rsdeseq2-fitted object with formula offsets
and no normalization provenance. Its stored matrix may already contain the
offsets, so it cannot safely be assumed to be a new base matrix. Require explicit
base factors or reconstruction from the original unfitted input when that
ambiguity cannot be resolved. Include this saved-object migration case in tests.

Store base factors, offsets, source labels, sample/gene identity and the last
effective factors/provenance. Reuse an existing representation where possible;
do not add multiple full dense copies merely to save provenance. Never recover
base factors by repeatedly dividing effective factors by old offsets.

Validate finite log offsets and finite positive effective factors, including
overflow/underflow in `exp` or multiplication. Reject missing, duplicate or
mismatched labels before arithmetic. Align size factors, weights and designs
through the same canonical sample order.

Resolve R and Rust behavior against the explicit-composition contract above.
For DESeq2 parity, supply the manually composed normalization matrix to the
reference; do not assume formula-offset support is identical across interfaces.

### C06: apply the contract to transforms and object write-back

- A second identical `DESeq` call on its returned dataset must reuse the base
  factors and apply the offsets once.
- A changed formula or explicit normalization argument builds new effective
  normalization from the appropriate base, then invalidates prior fit state.
- Separate transform design selection from normalization preparation. The
  proposed blind-transform contract is that blindness changes the fitted design
  and dispersion treatment while retaining supplied normalization/exposure.
  Prepare the original formula's normalization before replacing the fitting
  design with `~1`. Fresh and previously fitted equivalent inputs must agree.
- Verify that blind-transform contract against explicit-factor reference calls
  and document any correction to old behavior. Do not let prior fit history
  decide whether an exposure offset is present.
- For LRT formulas, verify compatible full/reduced normalization offsets.
  Explicitly reject unsupported differing-offset models instead of silently
  dropping an offset from one model.
- Write object assays/metadata consistently. Remove stale package-owned
  replacement assays, flags, obsolete coefficient columns and normalization
  fields when the new fit no longer produces them. Preserve unrelated user data;
  track exactly which columns/assays the package owns.

### Acceptance

Test repeated fit, fresh-versus-refitted transform, reversed named columns,
gene-order permutations, non-unit size factors plus offsets, non-unit base
normalization matrices plus offsets, explicitly composed matrices, changed
offsets and changed factor slots. Include both blind and design-aware transforms.

Expected equality is exact effective factors and unchanged numerical outputs
for equivalent inputs where the same arithmetic is used. Independently check
the corrected behavior against explicit-factor DESeq2 fits. The historical
reproductions with LFC changes of 0.3901497 and 0.7733569 must disappear.

**Commit boundaries:** C05 introduces the normalized input contract and fitting
fix; C06 applies it to transforms and fit-owned state cleanup. Each includes
its own regression cases.

## 7. C07–C08 — Validate cached fits and support gene subsets safely (F3)

Primary consumers: `results`, `resultsNamesRust`, all `lfcShrink` methods and
transform paths that reuse fitted dispersion/trend state.

### C07: make stale state impossible to consume silently

Introduce one fitted-view validator. Record canonical gene/sample identities,
design/full and reduced matrices, normalization provenance, weights, counts and
relevant options at fit time. Reuse existing stored inputs instead of blindly
duplicating them. Normalize integer-versus-double representation when checking
numeric equality: equivalent count values should not fail solely on R type.
Compare the dataset's original counts with original fitted input, keeping
replacement counts in separate state. Define manual edits to package-owned
dispersion/coefficient fields or fitted assays as invalidating the corresponding
cached state; do not allow two contradictory versions of a fitted quantity.

The validator must inspect current mutable content. A stored signature compared
only with another stored signature cannot detect user edits. Use bounded-memory
value comparisons against retained canonical input, with signatures for quick
rejection/diagnostics where useful. Never trust matching dimensions alone.

For an unchanged object, avoid rebuilding formulas or copying every large array
unnecessarily. Still check enough current content to detect same-size edits.
Measure this O(genes × samples) validation cost; a batch-results request can
validate once for multiple contrasts.

Initially reject every changed fitted view with a specific reason and refit
instruction. Cover count changes, sample subset/permutation, design/covariate
changes, factor levels/references, normalization edits and weight changes.
Unrelated metadata such as a gene annotation must not invalidate a fit.

Legacy objects may be admitted only when their available inputs prove the
requested view valid. Otherwise require refitting. No silent automatic refit:
that could use different options, time or normalization than the original fit.

### C08: add supported gene subset/reordering

The final M1 contract supports unique gene subsets and permutations while
requiring the fitted sample universe/order and relevant sample inputs to remain
unchanged. Build one current-gene-to-fitted-row mapping and apply it everywhere.

Maintain a field inventory:

| Field kind | Examples | Mapping behavior |
| --- | --- | --- |
| Per-gene vector | Dispersion, base mean, flags, iterations | Index with gene mapping |
| Gene-by-sample matrix | Counts, weights, means, Cook's, replacement data | Map rows; keep validated sample order |
| Gene-by-coefficient matrix | Beta and SE | Map rows; preserve coefficient identity |
| Gene-by-p-squared matrix | Covariance | Map rows; preserve coefficient layout |
| Sample/design data | Factors, design, model-frame levels | Validate against fitted context |
| Fit-global state | Dispersion prior/trend and fit options | Preserve fitted meaning; do not refit on subset |
| Result-population state | BH and independent-filtering decision | Recompute for the current result population as supported by the reference |

Check current subset counts and normalization against the corresponding original
rows. A matching gene name does not authorize reuse after its count values change.
Reject duplicated/unknown gene identities and ambiguous unlabeled subsets.

For shrinkage, distinguish stored likelihood state from a prior estimated during
the current shrinkage call. Follow the reference's subset semantics for the
latter; do not substitute a slice of previously shrunk full-data results when
that would use a different prior. Apply the same mapping to supplied results
tables and reject ambiguous or mismatched rows.

### Acceptance

Use `fit[c(3,1), ]`, a full reversal, an unchanged object, a one-cell count edit,
same-size sample/design edits, and a subset followed by a count edit. Verify
correct row count/order for results and all shrinkage methods. Check p-adjustment
on the subset against a subset reference, not merely against a slice of the
full-data adjusted p-values.

**Exit:** no consumer mixes current rows with cached arrays from another order.
All unsupported mutations error before statistical processing.

## 8. C09–C10 — Separate model validity from numeric units (F5–F6)

### C09: scale-aware design validation

Primary files: `math/qr.rs`, design validation and their callers.

Keep the existing public function with an explicitly absolute tolerance as an
absolute-tolerance primitive. Introduce a separate default design-validity policy
that evaluates relative scale. Use stable column norms and pivoted rank analysis;
avoid overflow in norm calculation and avoid treating a small but nonzero column
as literally zero.

Before choosing the final tolerance, generate a compact rank-decision table
against the pinned R implementation: exact duplicate columns, intercept aliases,
zero columns, independent columns scaled by powers of ten, and controlled
near-aliases. Include the reproduced `1e-12` covariate case. Freeze the intended
decisions and document genuinely ambiguous numerical-boundary cases.

Scaling a temporary validation matrix is permitted. Scaling the actual fit
design, beta prior, coefficient units or optimization bounds is outside this
change. Structural full rank does not guarantee a well-conditioned fit; expose
conditioning diagnostics separately when useful.

**Acceptance:** unit changes no longer cause false zero-column rejection;
singular designs remain rejected; previously valid model fits retain exact
outputs. Public absolute-tolerance callers retain their documented behavior.

### C10: nested-model validation for LRT

Primary file: `core/implementation/helpers.rs` and every public LRT entry path.

After C09, verify that the reduced design's column space lies inside the full
design's column space under the same scale-aware policy. Prefer a reusable
decomposition/projection check; use augmented-rank checks as an independent
small-case assertion. Equal column names or a smaller column count are not proof
of nesting.

Run this check before dispersion estimation or GLM optimization. Accept nested
models expressed in a different basis or with different column scales. Preserve
the existing residual-df and weighted-design rules for valid models; verify that
both designs refer to the same samples and normalization context. Any additional
gene-specific weighted-rank behavior needs its own reference case rather than an
unexamined change to degrees of freedom.

**Acceptance:** intercept-only reduced models, categorical reduced models and
reparameterized nested bases pass. The reproduced non-nested alternating-vector
case fails early. Valid Wald/LRT output and fallback diagnostics remain exact.

## 9. C11 — Make worker settings effective (F7)

Primary files: builder configuration, pipeline entry points, transform/shrinkage
entry points with worker options, and `tests/thread_invariance.rs`.

### Execution contract

- An explicit positive worker count installs a local Rayon pool for the complete
  analysis. It overrides a conflicting ambient pool for that analysis.
- An unspecified count inherits the caller's context. It does not rewrite global
  environment variables or reset the global Rayon pool.
- Reuse one execution context through dispersion, fitting, replacement/refitting
  and transformation stages. Internal helpers receive/borrow that context.
- Avoid creating pools per stage or per gene. If a caller provides an existing
  compatible context, reuse it. Separate concurrent analyses may own separate
  bounded pools.
- Reject invalid/overflow counts before pool creation. Return pool-construction
  errors as normal API errors.

Make the effective pool size observable inside a representative parallel stage,
through a diagnostic hook or test-only observation. A getter returning the
requested number does not prove enforcement. Pool size is a worker limit, not a
promise that every tiny workload runs simultaneously on all workers.

The current execution-mode variants should retain source compatibility. Document
which modes currently share the strict numerical algorithm. Diagnostic mode can
add observability; do not invent an approximate fast path to justify a name.

### Acceptance

Exercise the public builder/R option with requested one/four workers inside a
conflicting ambient pool, including refits and transforms. Verify actual pool
context in the work, and exact numeric outputs, NA masks, routes and convergence
at one/four workers. Retain the existing determinism test as a separate contract.
Check two independent contexts for interference and clean shutdown on error.

Existing benchmark runs explicitly set `RAYON_NUM_THREADS`; record both requested
and effective workers in future evidence. This fix alone is not a measured gain.

## 10. C12–C13 — Make validation and publication enforce their claims (F10–F11)

### C12a: deterministic vendor and version checks

Primary files: `scripts/release/sync_r_vendor.sh`, manifests and R CI.

1. Derive the vendored package manifest from workspace/package metadata, including
   versions, dependency features, edition and minimum compiler. Eliminate the
   second hand-maintained dependency/version list.
2. Add a read-only check mode. Generate the expected vendor tree in a temporary
   directory, compare it with the committed tree, then fail with a useful diff.
   CI must not synchronize the tested checkout before deciding whether it was
   already correct.
3. Validate Rust crate, R DESCRIPTION, R backend, vendor manifest and lockfile
   consistency. Keep intentional historical benchmark versions exempt.
4. Ensure generation is deterministic and does not depend on a developer's
   absolute paths or an unadvertised scripting-runtime version.

Test stale source, stale manifest, version mismatch, unchanged tree and repeated
generation. A deliberate vendor mismatch must fail the validation job.

### C12b: benchmark acceptance must include precision

Primary files: `scripts/benchmarks/benchmark_gtex_native_wald.py`, the existing comparison
helpers, and Python script tests.

Define two explicit modes:

- **Accepted-run mode:** requires frozen input/reference identity and either an
  exact native-output baseline or a named, fixed numerical contract appropriate
  to an intentional correction. Success means every accuracy gate passed.
- **Exploration mode:** writes measurements and differences but marks the run
  unaccepted. Its output cannot be used as release/optimization acceptance.

Reuse existing comparison/validation functions. Check gene and column coverage,
duplicates, labels/order where required, finite/NA masks, forbidden infinities,
per-field error limits and meaningful nonempty coverage. A table of finite wrong
numbers must fail. So must an all-NA table that vacuously avoids finite errors.
Hash binaries before and after runs and verify frozen input/reference hashes.

Keep the existing full GTEx precision checker authoritative; call its reusable
contract rather than maintain a looser second copy. Add regression cases for
finite-but-wrong output, stale inputs, altered binaries and missing diagnostics.

### C13: one validation dependency graph for publication

Primary files: `.github/workflows/rust.yml`, `r.yml`, `parity.yml`,
`publish-crates.yml` and `release-binaries.yml`.

Create reusable validation jobs/workflows and make every publication path depend
on their successful results for the exact release commit and staged artifacts.

```text
source/vendor/version checks
            |
            +-- Rust tests, formatting, Clippy, minimum-compiler build
            +-- Python comparator/generator tests
            +-- exact R source package build and R CMD check
            +-- target build and compact numerical checks
            +-- required reference/parity evidence
                         |
                  artifact manifest gate
                         |
             crate / R / binary publication
```

Use stable Rust for the normal suite and the advertised minimum Rust version
for a separate locked build. Exercise compact numerical cases on supported
Linux, macOS and Windows environments. For a shipped architecture lacking a
native runner, use a runnable environment where available or explicitly record
the execution gap; cross-compilation alone is not a numerical test.

Split the pinned reference environment from a current compatibility environment.
Pin/reference R, DESeq2, apeglm, ashr and relevant numerical dependencies using
reproducible package sources/locks and checksums, with BLAS/thread/session
information. The currently recorded reference set is R 4.6.1, DESeq2 1.52.0,
apeglm 1.34.0 and ashr 2.2.63. Do not update oracles automatically when an upstream
package changes.

Build the R source tarball once, test that tarball, and retain its hash through
publication. Prepare release-specific README substitutions before packaging,
then freeze the staging tree. For Cargo publication, verify the staged package
contents/manifest against the checked package and publish from that frozen tree;
do not imply `cargo publish` simply uploads an arbitrary supplied archive.
Record any unavoidable regeneration and verify artifact identity at the boundary.

All external publication mutations must wait for the combined gate. Tests running
in a different independent workflow are insufficient unless their exact commit
and successful completion are enforced. Keep release credentials out of PR jobs.

**Acceptance:** intentional R failure, stale vendor and version mismatch each
block the publication dependency graph. Artifacts have a manifest containing
source/staging identity, package/binary hashes, environment and validation result.
Check this in dry-run workflow execution; do not publish a test release merely
to test the graph.

## 11. C14–C18 — Reduce memory without changing fitting decisions (F8)

### C14: measure the live allocations and real stage costs

Replace placeholder benchmark workloads with representative IRLS, dispersion,
result extraction and transform cases. Keep small deterministic stage benchmarks
for development and use complete workflows for acceptance.

Add optional stage diagnostics for elapsed time, dimensions, fallback counts,
requested output retention and effective workers. Keep clocks and profiling data
out of numerical result exports. Default execution should incur no meaningful
profiling allocation or synchronization cost.

Measure:

- Eleven existing native Wald contrasts, with their original settings.
- A moderate-coefficient LRT to expose full/reduced fit ownership.
- R `DESeq -> results` and repeated/multiple-contrast results.
- Replacement disabled, enabled without outliers, and enabled with actual refits.
- rlog with a bounded sample/design-size grid, starting small and increasing only
  after the allocation estimate fits the host's memory budget.

For rlog, one current covariance array costs `genes * (samples+1)^2 * 8` bytes:
about 1.61 GiB at 56,200 genes/61 samples, or 6.02 GiB at 20,000 genes/200 samples.
These are analytical buffer sizes, not peak RSS measurements or expected savings.
Record both allocation estimates and observed peak RSS without conflating them.

Fresh profiling chooses the next hotspot. Existing optimizer/locfit and
cross-product optimizations must be read from the optimization log before
proposing similar work again.

### C15: remove ownership copies first

Primary files: `glm/irls.rs`, `glm/fallback.rs`, `builder/fixed_lrt.rs`, and the
R native fit/replacement bridge.

- Borrow beta/covariance data for fallback inspection instead of cloning full
  arrays to construct a temporary inspection object.
- Give the internal routing helper a borrowed view or validated diagonal/state
  view. Retain the public helper as a forwarding compatibility interface.
- Move reduced-model means, hats and other owned outputs into the destination
  when the source fit is consumed.
- Merge replacement fits with selective moves/updates; avoid cloning the complete
  original fit solely to retain data already owned by the final result.

Audit ownership lifetimes before using `mem::take` or moving fields: diagnostics
and original-versus-replacement semantics must still refer to the right data.
Keep this commit free of arithmetic changes.

**Gate:** exact full arrays and diagnostics on focused/reference cases, followed
by a candidate measurement. Accept a clearly measured RSS reduction with no
material timing regression even if it is below the 2% speed-round threshold.

### C16: bound temporary row-output retention

Start with safe, ordered batches of genes. For each batch, compute per-gene
outputs in parallel and write/move them into the final arrays in canonical order.
Release the batch before starting the next. This is simpler to audit than adding
unsafe concurrent output writes.

Choose a batch bound using output bytes per gene and worker count, then cap it.
Preserve each gene's arithmetic order, initialization, convergence criteria and
fallback options. Do not introduce cross-gene reductions or new thresholds.
If several genes fail, define deterministic reporting by gene order rather than
returning whichever worker fails first.

**Gate:** exact fit/result values and fallback/iteration routes, including error
cases and one/four workers. Record temporary-memory reduction and runtime. Keep
the simpler C15 improvement if batching is slower without sufficient RSS benefit.
Compatible positive partial effects can be evaluated together.

### C17: retain only outputs needed by internal stages

Introduce internal requirements, with the existing public full-inference path
remaining the default:

| Consumer | Required retained output | Potentially avoidable global storage |
| --- | --- | --- |
| Full Wald/LRT inference | Coefficients, covariance and all advertised diagnostics | Only redundant copies |
| Dispersion mean fitting | Fitted means and required convergence/fallback state | Whole-gene covariance and unrelated exported summaries |
| High-level rlog | Transform coefficients/means and required transform state | Global covariance used only for validation/routing |
| Internal diagnostics request | Explicit additional diagnostics | Unrequested diagnostics after their last use |

First change retention, not computation. Compute a gene's covariance/variance
checks where they affect fallback, then discard what the stage does not consume.
Keep any values needed by optimizer fallback locally for as long as required.

The target memory structure for transform stages is their required O(genes ×
samples) / O(genes × coefficients) output plus bounded worker/batch scratch.
Avoid a retained O(genes × coefficients squared) covariance buffer in consumers
that never expose it. Temporary per-gene quadratic work may still be necessary.

Only consider skipping calculations in a later, separately measured change after
proving they cannot affect routes, SE validation, diagnostics or error behavior.

**Gate:** rlog and mean-fit parity, exact fallback identities, iterations and
convergence, measured memory scaling across the bounded grid, and unchanged full
inference outputs. No silent loss of fields from a public full-fit return value.

### C18: avoid unnecessary R replacement buffers

With C02's checks in place, first omit the replacement-count allocation/copy when
replacement is disabled. Pass an explicitly valid absent-buffer representation.
Then handle enabled-but-no-refit runs using an output contract that allocates R
replacement storage only after native computation reports actual replacement.

If a native result handle is needed, define ownership, length query, copy-out,
cleanup/finalizer and error behavior explicitly. Ensure an R allocation failure
cannot leak the handle, and do not let R longjmp through Rust. Include the handle
schema in the ABI declarations shared by both languages.

**Gate:** disabled/no-refit/refit cases preserve counts, assays, flags, diagnostics
and results. Child-process error/cleanup tests and RSS measurements demonstrate
the omitted allocation. A more complex handle design must earn its maintenance
cost; the disabled-path improvement can ship independently.

## 12. C19–C21 — Consolidate statistical processing and internal execution (F9)

Start after C04 has established correct result semantics and C02 has hardened the
bridge. Keep this work separate from the initial correctness release.

### C19: extend the existing Rust results APIs

Create a validated `ResultRequest` and borrowed fitted-state view around the
existing `DeseqFit` contrast methods. Carry Wald options, LRT state, Cook's policy,
alpha and filtering options explicitly. Resolve selectors once and reuse
coefficient/covariance storage.

Build one policy implementation shared by public fitted-state requests and the
R bridge. Avoid the current default-Wald-options assumption in methods that must
serve a nondefault threshold/t request. Keep primitive-fit support by validating
which context exists, rather than fabricating missing fields.

Port one operation at a time: contrast covariance/SE, Wald alternatives/tails,
Cook's policy and filtering/BH. Check against both the current supported R
behavior and pinned DESeq2 at ordinary values and tails. Existing Rust special
functions are not automatically numerically interchangeable with R's.

**Gate:** exact equivalence where arithmetic is unchanged; otherwise all existing
per-gene bounds, missingness and worst-case reference error limits remain intact.
No widened tolerance or degraded tail precision to complete the migration. Keep
an operation on its validated implementation until the shared version passes.

### C20: use a bulk R results bridge

Pass validated primitives plus one or several requests in a single call. Validate
the fitted dataset once for a batch. Borrow covariance and fit arrays across
requests rather than copying them for each contrast. Bound result storage or
process request batches when the requested output itself is large.

R continues to own Bioconductor classes, row/column descriptions, aliases and
metadata assembly. Remove duplicate production R arithmetic only after its
replacement passes the full supported matrix. Keep independent reference code
in tests/generators, not a second production decision path.

Do not reuse BH/filter decisions across different p-value vectors or alpha
values merely because the fit is shared. Cache only validated invariants such as
resolved coefficient indices and unchanged base input.

Measure fit time separately from one-result time, ten-contrast extraction and
shrinkage. Report both stage and end-to-end effects. A faster results-only call
does not establish a 2% improvement in total quantification.

### C21: narrow internal module responsibilities

Converge public builder forwarding methods on a small internal path:

```text
request -> aligned/validated input -> normalization -> dispersion
        -> GLM/refit -> requested statistics -> presentation
```

Use real Rust module boundaries for input validation, fitting, result policy and
execution context where they remove shared mutable assumptions from `include!`
namespaces. Use validated internal request types to prevent impossible option
combinations from reaching fitting.

Retain public builder method signatures/aliases as forwarding adapters in this
cycle. Do not combine a broad API removal, formula-parser rewrite or new generic
framework with numerical consolidation. Move one responsibility at a time,
maintaining numerical coverage and compiler-visible access restrictions.

**Exit:** one authoritative result-policy path, one validated execution path per
workflow, fewer repeated decisions, and measured R wrapper overhead. Reduced
line count by itself is not an acceptance criterion.

## 13. C22 — Make CLI output failure behavior explicit

This follows the P1 correctness work and does not block the first correctness
milestone unless an active output-corruption case is discovered.

- Preflight all main/sidecar paths and input paths. Detect path aliases through
  canonical existing parents, symlinks and existing-file identity; account for
  filesystem case behavior. Reject destructive input/output or sidecar collisions.
- Write all outputs into temporary files on the destination filesystem. Close
  them and validate successful completion before publishing any final path.
- Publish a completion manifest last, containing generation identity, required
  file paths and hashes. Consumers that need coherent generations verify it.
- Document that separate path renames are not one atomic transaction. A failure
  during publishing can leave a generation incomplete; the absent/invalid
  completion marker must expose it. Never claim legacy readers become atomic.
- Consider an optional directory-bundle output for future truly atomic directory
  publication, while retaining current output paths for compatibility.

Test unwritable destinations, path collisions, sidecar write failure and failure
during publication. Assert expected errors and absence of a misleading valid
completion marker; retain useful numeric assertions for successful output.

Cancellation/progress requires a separate safe R-thread orchestration design.
It is a worthwhile follow-up, but is not bundled into numerical or FFI memory
work without an explicit ownership/interrupt design.

## 14. Validation design: add relationships the suite currently misses

Preserve the compact GTEx fixture: 512 genes, all 61 samples and its documented
full-cohort normalization/reference construction. Keep the opt-in full 56,200-gene
test and the eleven full contrast replays. Further shrinking the suite is not
the objective of this work.

### Three complementary forms of evidence

1. **Independent oracle:** pinned DESeq2/R results on the supported request and
   exactly the same input population/settings. This detects bugs shared by our
   R and Rust interfaces.
2. **Behavioral relationship:** equivalent selectors, label permutations, repeated
   fitting, valid subset mapping and explicit-normalization equivalence. These
   catch dispatch/state defects that isolated numerical fixtures missed.
3. **Implementation preservation:** exact baseline arrays, exported results and
   route/iteration diagnostics for a refactor or memory optimization. Final
   p-values alone can conceal a changed optimizer path.

Every new test should name the failure it detects. A useful table row changes
one meaningful condition and asserts a numeric effect or meaningful equality.
Avoid new constructor/header-only smoke tests, repeated full fits for selector
variants, and a Cartesian product of unrelated flags.

### Additional validation panel

Use one compact, independently selected biological cohort plus a small synthetic
boundary panel. Preserve source, seed, sample/gene selection and generation
commands. Suitable initial scope is a few thousand genes and tens of samples;
choose final dimensions for the behavior being exercised, not an arbitrary
number of studies.

| Panel component | Distinct risk covered | Reference/assertion |
| --- | --- | --- |
| Held-out biological cohort | Generalization beyond Kidney Cortex | Full supported numerical result comparison |
| Small residual df and t inference | Tail/df and design assumptions | Pinned per-gene reference values and masks |
| More design coefficients, including an interaction | Covariance, contrast mapping and memory | Contrast SEs, statistics and layout |
| Heterogeneous weights, including rank-failing rows | Weighted preprocessing and skipped rows | Base mean, dispersion, flags and results |
| Cook's outlier with sufficient replicates | Actual replacement/refit behavior | Original/replaced counts, masks and final results |
| Low-count Cook's exception | Policy edge behavior | Expected row-level masking decision |
| Local dispersion trend | Nonparametric route | Trend, dispersion and final statistics |
| Valid/invalid LRT and alternative nested basis | Model nesting and displayed-effect separation | Early error or unchanged test statistic |
| rlog/VST with supplied factors and offsets | Transform normalization and retention | Reference transform values and reuse equivalence |
| Near-singular and mixed-unit designs | Rank validation boundaries | Explicit rank decision and errors |

Keep the biological hold-out out of routine optimization tuning. Evaluate it at
a milestone. If a failure is inspected and tuned against, retain it permanently
as a regression case and select another independent hold-out for later numerical
algorithm work. A random split of the same existing cohort adds a design case;
it does not supply an independent biological hold-out.

### When to run which checks during implementation

- **Each correctness commit:** new focused regressions and the affected R/Rust
  reference tests; use isolated native-boundary tests where relevant.
- **Each core ownership/retention commit:** focused complete-array and routing
  comparisons, then the affected integrations. Synchronize/check the vendor when
  core source changes so the R test uses the intended implementation.
- **Each accepted performance candidate or combined candidate:** one paired
  eleven-contrast run plus its targeted R/transform workload. No full timing run
  is needed for a documentation-only or equivalent module-move commit.
- **Milestone candidate:** complete Rust, Python and R package checks; full GTEx
  precision evidence; independent workflow panel; package/vendor identity.
- **Release candidate:** final staged artifacts pass the same gates, including
  target-specific checks and publication dependencies. Reuse evidence only when
  the validated code/artifact identity is unchanged.

Existing commands to build the implementation runbook around include:

```sh
cargo fmt --all -- --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test --locked --workspace
python3 -m unittest discover -s scripts/tests -p 'test_*.py'
R CMD build r-pkg/rsdeseq2
# Check the exact package filename emitted above, not an ambiguous wildcard.
R CMD check --no-manual --no-build-vignettes rsdeseq2_VERSION.tar.gz
cargo test --locked -p rsdeseq2 --test gtex_real_wald_parity full_gtex_wald -- --ignored --nocapture
```

`VERSION` is a placeholder in this proposed runbook. Also format/lint the separate
R Rust backend. Use an isolated absolute R library path and the pinned reference
environment; require suggested numerical dependencies so checks cannot pass by
silently skipping apeglm/ashr coverage. These are future implementation gates,
not commands run while writing this plan.

## 15. Performance, precision and memory acceptance rules

### Paired measurement protocol

1. Freeze baseline and candidate binaries/packages with hashes.
2. Fix inputs, tool/compiler settings, workers and BLAS threads. Record effective
   settings from inside execution after C11.
3. Run baseline/candidate/DESeq2 in alternating or rotating order by contrast,
   using the existing script with `--repeats 1`. Avoid concurrent benchmark jobs.
4. Collect elapsed time and peak RSS from the same accepted run where possible.
   Normalize OS-specific RSS units explicitly.
5. Compare all output/diagnostic files before accepting measurements. A faster
   wrong output is a failed candidate.
6. For R/transform workloads, use separate child processes per defined workflow
   so a previous allocation does not inflate a later workflow's reported peak.
   Record fit setup separately when reporting repeated-results cost.

For total times across the eleven contrasts:

```text
T_old = sum(previous accepted Rust elapsed times)
T_new = sum(candidate Rust elapsed times)
T_R   = sum(DESeq2 elapsed times)

throughput_gain_percent = 100 * (T_old / T_new - 1)
elapsed_time_reduction_percent = 100 * (1 - T_new / T_old)
speed_ratio_vs_DESeq2 = T_R / T_new
```

A 2% throughput gain means `T_old / T_new >= 1.02`; it is not exactly the same
as 2% less elapsed time. Publish ratios of total times, not the arithmetic mean
of per-contrast ratios, and retain per-contrast results to expose regressions.

For memory, retain per-workload peak RSS and ratios against DESeq2 and the
development baseline. Do not sum peak RSS across sequential contrasts and call
the sum a process peak. A transform's reduced covariance allocation is evidence
about one buffer until a complete-process RSS measurement confirms the effect.

### Decision table

| Outcome | Decision |
| --- | --- |
| Any unexpected numerical/missingness/route regression | Stop; localize the first changed stage; do not accept the candidate |
| At least 2% total throughput gain, precision preserved, memory acceptable | Accept a speed round and update the accepted baseline |
| Smaller compatible positive gain | Keep as a documented partial candidate; combine and measure before counting a round |
| Clear RSS reduction, speed approximately unchanged | Accept as a memory improvement; do not count a 2% speed round |
| More than 2% total slowdown or a conspicuous workload-specific regression | Treat as a blocking investigation, not an unexplained acceptable fluctuation |
| Small negative timing change | Report it; inspect affected stages; do not claim exact speed preservation from one measurement |
| Necessary correctness fix with measurable validation overhead | Minimize it, report it explicitly, and judge correctness separately from speed claims |
| Invalid timing run because of an identified external interruption | Discard that run with the reason recorded; repeat the invalid measurement only |

One paired pass is the chosen operational acceptance procedure. It is not proof
of statistical significance for a small timing difference. Do not introduce
routine repeated passes simply to obtain a favorable number.

Memory budgets for new large-transform cases should be recorded after C14's
baseline and before evaluating candidates. Set a host-safe maximum workload/RSS
and a target allocation reduction per consumer. Avoid arbitrary post-hoc limits
chosen after seeing a candidate. Full inference must preserve its public output;
transform savings must not be achieved by dropping required inference output.

## 16. Commit and dependency ledger

Each row is a reviewable unit. A row may be split mechanically if its tests and
implementation remain understandable, but never mix unrelated numerical changes.
Do not commit a deliberately failing suite as an intermediate milestone.

| Commit | Scope | Prerequisite | Main acceptance evidence |
| --- | --- | --- | --- |
| C01 | Baseline/evidence manifest | Current reviewed commit | Immutable binaries, references and defect cases identified |
| C02 | C/native buffer and scalar checks | C01 | Child-process errors; valid outputs exact |
| C03 | Raw result state, schema and explicit options | C02 | Raw LRT recovery; unsupported arguments rejected |
| C04 | Unified result policy | C03 | Selector/Cook's/threshold/filter reference matrix |
| C05 | Aligned normalization provenance and fit composition | C01 | Explicit-composition and repeated-fit equivalence |
| C06 | Transform normalization and stale fit-owned field cleanup | C05 | Blind/nonblind reuse and cleanup cases |
| C07 | Shared fitted-view validity guard | C03, C05–C06 | Same-size edits detected by every consumer |
| C08 | Gene subset/reordering views | C04, C07 | Correct mapped arrays, subset BH and shrinkage semantics |
| C09 | Scale-aware default rank validation | C01 | Pinned rank table; unchanged valid fits |
| C10 | LRT subspace/nesting guard | C09 | Nested-basis acceptance and early invalid-model error |
| C11 | Effective reusable worker context | C02 | Public setting enforced; one/four-worker exactness |
| C12a | Read-only vendor/version gate | C01 | Deliberate drift fails |
| C12b | Accuracy-enforced benchmark modes | C01 | Finite wrong results cannot be accepted |
| C13 | Combined validation/publication graph | C12a–C12b | R failure blocks publishing; artifact identity recorded |
| C14 | Representative profiling and memory baseline | M1–M2 | Measured stage costs and safe memory grid |
| C15 | Borrow/move instead of full-fit clones | C14 | Exact routes/arrays; measured allocation/RSS change |
| C16 | Bounded ordered output batches | C15 | Exact arithmetic and deterministic error reporting |
| C17 | Internal output retention requirements | C16 | Transform memory scaling with unchanged fallback |
| C18 | Conditional R replacement storage | C02, C06, C15 | Correct lifecycle, replacement output and RSS |
| C19 | Shared Rust result request/evaluation | C04, C08 | Existing and independent statistical precision bounds |
| C20 | Bulk R result bridge | C02, C19 | Wrapper parity, validation reuse and measured cost |
| C21 | Internal modules and forwarding adapters | C11, C19–C20 | Public compatibility and invariant-preserving execution |
| C22 | CLI preflight/staged outputs/completion marker | C12b | Failure behavior and valid numeric exports |
| C23 | Final compatibility, performance and release documentation | Completed release scope | Documentation matches tested artifacts and claims |

The numbering gives a practical execution order. Independent correctness and
delivery changes may move earlier without waiting for unrelated refactors.
Dependencies are conceptual; this plan does not require parallel agents.

### Checkpoints and recovery

- After M1, preserve a known-good correctness checkpoint and its artifacts before
  starting ownership changes.
- After each accepted optimization combination, record its commit/binary hash,
  measurements and exact-output evidence as the next baseline.
- If a candidate fails, retain its hypothesis, diff identity and failure evidence
  in the optimization log. Keep useful independent partial changes separately.
  Do not repeatedly rediscover the same rejected variant.
- Revert or revise the smallest faulty change; avoid rolling back unrelated
  correctness fixes. Do not change reference files to make a regression disappear.
- If a saved-fit schema changes, test round-trip save/load and define supported
  legacy behavior. Do not assume a newer object can be loaded correctly by an
  older package; document that boundary.

## 17. C23 — Documentation and completion criteria

Update [R wrapper documentation](../r-wrapper.md), [compatibility](../compatibility.md),
[development](../development.md), [testing/parity](../testing-parity.md),
[reproducibility](../reproducibility.md) and [release procedures](../release.md) alongside
the relevant implementation, then reconcile them at the release checkpoint.

Document supported result options, primitive-fit limitations, normalization
precedence, formula-offset and blind-transform behavior, supported gene subsets,
fit invalidation, worker semantics and legacy fitted-object handling. Remove
stale future/placeholder descriptions only after checking the actual supported
scope. Correct maintainer metadata before broader package distribution.

Main README updates should contain concise current precision, total speed ratios
and memory comparisons against original DESeq2. Put old-release deltas, detailed
optimization experiments and per-stage development memory costs in docs. Preserve
the reproducible contrast generator and document any new panel generator there.

A milestone/release is ready only when:

- Every finding in its declared scope has a focused regression and clear API
  behavior, including invalid requests and legacy-data limitations.
- Unaffected frozen numerical evidence remains intact, and intentional corrections
  have independent reference evidence with unchanged precision standards.
- The R wrapper, core and vendored source/package agree on the tested code.
- Performance claims identify the exact binary, workload, settings and raw run.
- Memory claims are measured at process level and separate from buffer estimates.
- Required Rust/R/Python/platform and package gates pass for the actual artifact.
- Publication waits for that combined validation; no unrelated untracked files
  or development logs are swept into the release.

## 18. Work deliberately deferred

- L-BFGS-B tolerances, convergence logic and locfit internals: revisit only if
  fresh profiling or stage-by-stage mismatch evidence identifies a remaining
  problem there. The reviewed workflow defects have other causes.
- Another large formula parser rewrite or removal of public builder aliases:
  unnecessary for the proposed correctness and retention contracts.
- More broad test deletion: useful numerical coverage now has greater value than
  reducing counts. Consolidate only when touching a demonstrably repetitive case.
- CLI-only dependency feature separation: measure library/R build cost first,
  then make a small packaging change with the new vendor/artifact gates.
- New sparse/delayed-matrix support, approximate modes and interrupt callbacks:
  separate API/numerical/ownership projects, not incidental additions here.

The first concrete implementation batch is C01–C04: establish the evidence,
harden the native boundary, preserve raw result state, and repair result policy.
Then implement C05–C08 for normalization and fitted-object reuse. This creates a
safe foundation for measured memory work while protecting the precision already
achieved in the released numerical core.
