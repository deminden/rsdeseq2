# Final local plan acceptance

The original [C01–C23 plan](codebase-improvement-plan-2026-10-22.md) now has
**22 checkpoints closed locally**. C19 and C21 are accepted on the requirements
below. C23 remains open until current release artifacts and required platform
validation are reconciled. This is a requirement count, not an estimate of
remaining effort or complete Bioconductor compatibility.

The numerical source is `48917a0`. All 635 paths in its recorded runtime-source
manifest still match. Changes made while gathering this report affect profiling
tools, their acceptance tests/CI wiring, documentation and exported measurements.
No statistical source, reference value, precision tolerance, dependency lockfile
or package version changes to obtain acceptance.

The [machine-readable report](data/final_plan_local_acceptance.json) binds the
local artifacts by SHA-256. Raw process logs, before/after payloads and package
artifacts remain under `results/final-combined-48917a0/` and
`results/native-sparse-sqp/`; these generated directories are not committed.

## C19: shared numerical and behavioral contract

| Requirement | Current implementation | Acceptance evidence |
| --- | --- | --- |
| Validated requests with explicit Wald/LRT, Cook's and filtering policy | `results/request.rs` validates scalar options, checks exactly one stored test and rejects Wald overrides on LRT. Missing required design/Cook's context fails. | Rust result-request, primitive and independent inference tests; 189 actual primitive interface comparisons. |
| Borrowed fitted storage and shared contrast inference | `PreparedContrastView` borrows estimates, covariance, design and counts; `core/fitted_results.rs` and the R inference adapter use the same kernels. | Independent coefficient/list/numeric/factor contrast and tail fixtures; current primitive, fresh-fit and formula-fit identity checks. |
| One statistical masking/adjustment policy | Shared `results/cooks_policy.rs`, `replacement_policy.rs` and `independent_filtering.rs` own masking, replacement-zero values, quantiles, lowess and BH. Each request recomputes its own filtering decisions. | Original policy/filter references; unchanged Rust/R bounds and missingness; later replacement-aware requests and alpha/selector variants. |
| Complete fitting and retained state, not only primitives | The R fit adapter constructs the core builder and uses its shared test dispatch and replacement merge. Later standalone requests explicitly retain replacement plans. | 31 matrix-fit requests: 27 successes, four controlled rejections and 113 later queries; all exported floating-point bits, missing values and acceptance decisions agree. |
| Formula preparation and complete formula fits | Core frame/matrix/expression preparation owns supported transforms, factor coding, offsets and polynomial/term ordering. Custom formula functions reject without execution. | 32 independently generated original-R matrix profiles; 33 actual formula requests: 29 successes, four rejections and 96 later queries, without column remapping. |
| Transforms and posterior methods also run without R calculations | Native saved-local VST, restricted dispersion/likelihood expressions, normal/ashr/apeglm families, controls, mode search, random state and posterior summaries use core modules. | The [25-export supported inventory](native-supported-contract.md) links each family to independent accuracy and actual interface checks. All 29 replay sets retain fields, shapes, acceptance, finite/missing patterns, categorical diagnostics and warnings. |
| Correct sparse SQP compatibility rather than a wider tolerance | Shared SQP now follows the original active-set budget, EM floors, operation order and condition-based direction solve. | [Sparse SQP report](native-sparse-sqp-development.md): nine full and nine primitive original profiles pass the existing `2e-6` bound. The unchanged public reordered-subset test passes its original `1e-5` bound. Ten replay sets intentionally change SQP numbers/diagnostics; nineteen are fully exact. |
| No R numerical fallback | R transports arguments/metadata and assembles objects; C marshals buffers; Rust owns statistical calculations. | The current production-R call guard checks 39 files with eight rejection controls. Actual identity tests block R numerical calls. A separate own C/Rust bridge scan finds no R evaluation/math callback entry points; this lexical check alone is not a behavioral proof. |
| Preserve independent precision limits | Existing references and bounds remain unchanged. Interface identity is checked separately from agreement with original packages. | Full Rust and R suites, 121 per-contrast plus nine global frozen GTEx gates, independent held-out Wald/LRT checks and the full 56,200-gene Wald test pass. All 213 core reference/input files and their four README files are unchanged since `48917a0`. |

The three discarded R-backed experimental commits remain outside the active
branch ancestry. Neither this audit nor historical source reports expand the
documented support boundary to arbitrary R functions, full R grammar, every
Bioconductor mutation convention or unimplemented original methods.

## C21: execution, modules and forwarding

| Requirement | Current implementation and proof |
| --- | --- |
| One validated execution path per supported workflow | Public fitting adapters forward into core builder workflows; normalization/dispersion/GLM/refit and result policy share their core implementations. Actual full matrix/formula requests and retained-state queries check adapter decisions as well as returned tables. |
| Real private module boundaries | `core/execution.rs`, `core/fitted_results.rs`, `core/replacement_fit.rs`, shared result-policy modules and private workflow modules delimit worker context, fitting, retained state and result evaluation. Remaining source splits are not treated as proof merely because line counts fall. |
| Retain public builder signatures and aliases | A balanced lexical scan of public inherent `DeseqBuilder` functions, masking comments/strings, compares against validated 0.4.0 source `d64d2f2`: 210 before and after, zero removed/changed/added signatures. Full workspace/MSRV compilation supplies separate type-check evidence. |
| Reuse worker context safely | `ExecutionContext` shares an `Arc<OnceLock<…>>` pool between stage/refit clones, validates requested bounds and avoids nesting entry into the same pool. Changing worker settings resets the context. Existing concurrent-context, error recovery and one/four-worker numerical tests pass. |
| Bounded bulk results, independent decisions | The native adapter accepts at most 64 requests per chunk and borrows fitted arrays. R validates/extracts dataset primitives once per batch. Filtering, masks and BH remain request-specific. Exact batch/single and selector/error tests pass. |
| Measure fit, one result, ten results and shrinkage separately | Twenty complete child-process workflows retain fit/results/transform/shrinkage stage times and process RSS. A separate six-case experiment times public fitted-list result calls and their outermost native adapters without refitting. Details below make the timing boundary explicit. |

### Measured wrapper costs

The [boundary profile](data/final_plan_r_boundary.tsv) reuses fitted state at
1,200 and 56,200 genes, 16 samples and five coefficients. All six profiled outputs
are `identical` to untraced controls. Trace setup/removal is outside the timers;
nested adapter intervals are not double-counted.

| Genes | Queries | Public call | Total (s) | Native adapter interval (s) | Outside interval (s) |
| ---: | ---: | --- | ---: | ---: | ---: |
| 1,200 | 1 | `results` | 0.004 | 0.001 | 0.003 |
| 1,200 | 10 | individual `results` | 0.045 | 0.006 | 0.039 |
| 1,200 | 10 | `resultsList` | 0.022 | 0.001 | 0.021 |
| 56,200 | 1 | `results` | 0.042 | 0.012 | 0.030 |
| 56,200 | 10 | individual `results` | 0.665 | 0.257 | 0.408 |
| 56,200 | 10 | `resultsList` | 0.177 | 0.118 | 0.059 |

Native intervals include the outer R helper's validation/marshalling and Rust
computation, rather than pure Rust execution. Outside time includes other R
validation/presentation and tracing overhead. Timers have millisecond resolution.
These measurements use fitted-list results, not complete dataset validation;
the twenty workflow profiles cover the latter separately. They establish the
cost boundary and batching behavior, not a new total-quantification speed round.

## C23: current local checks and release limitations

| Gate | Current evidence | Disposition |
| --- | --- | --- |
| Full Rust workspace | 1,154 passing tests, one existing opt-in ignored test; current runtime source matches the tested manifest. | Pass locally |
| Formatting/lint, including separate R backend | Accepted sparse-SQP source has passing workspace/backend format and Clippy logs; no Rust source changed afterward. | Pass locally |
| Exact R source-package check | The checked tarball has 430 source files matching the package tree, allowing only standard DESCRIPTION formatting/build metadata. `R CMD check`: Status OK; 2,867 assertions, no failures/warnings/skips. | Pass locally |
| Python tooling tests with actual R oracle gate | 49 tests pass with the installed checked R package and original ashr available. The new comparison gate rejects wrong moments, missingness, metadata, counts, fitted arrays, unrelated columns and wrong workflow selection. | Pass locally |
| Frozen GTEx | 121/121 per-contrast and 9/9 global gates pass; all 22 native result/diagnostic exports retain the accepted hashes. Full 56,200-gene ignored test passes independently. | Pass locally |
| Biological hold-out | Pasilla Wald/LRT original-reference bounds pass; all eight baseline/candidate table/diagnostic/beta/SE exports are exact. | Pass locally |
| Actual interface identity | Primitive comparisons, full fresh fits and formula fits pass with R calculations blocked, in addition to the accepted posterior/transform replay sets. | Pass locally |
| Paired total time and process memory | Eleven contrasts and twenty R workflows, one paired pass, all observations retained; details below. | Recorded locally with background-load limitation |
| MSRV | `cargo +1.97.1 check --workspace --all-targets --locked` passes. | Pass locally |
| Native boundary sanitizers | C ASan/UBSan plus Rust ASan pass, including four replacement and four formula unwind errors. R/dependencies/Rust standard library are uninstrumented; leak detection is disabled. | Pass for recorded scope |
| Dependencies | Official registry lookup finds all 17 direct workspace dependency declarations at latest non-yanked stable on 2026-11-11, including locfit 0.2.0/L-BFGS-B 0.3.0. No dependency change in this acceptance pass. | Recorded lookup; transitive constraints retain their earlier audit |
| Latest-commit artifacts | Verified crate, source archive, checked R tarball, macOS executable, provenance and checksums must bind to the final documentation/tooling commit. | Required after commit |
| Current platform validation | Earlier released-source CI is historical. Current five-target binary gates and independent Linux regeneration of the new formula/SQP references have not run. | Open |

### Total time and memory

The native panel uses the current CLI and a separately frozen **validated
0.4.0 `d64d2f2` candidate**. That baseline bundle was local and unpublished;
it is not labeled the published release. The R workflow baseline is the
pre-SQP `6e85f68` library. Binaries, source packages, installed libraries, tool
order and reference identities are retained in the hashed local provenance.

- Eleven native contrasts: **43.860 s** current, **45.597 s** baseline,
  **245.675 s** original DESeq2. Original/current total-time ratio: **5.60×**.
  Current versus baseline elapsed difference: **−3.810%**.
- Native maximum whole-process RSS: **301.938 MiB** current,
  **305.813 MiB** baseline and **1,772.109 MiB** original DESeq2.
- Twenty R workflows: **90.745 s** current versus **91.190 s** baseline
  (**−0.487%**). Maximum RSS is **1,134.563 versus 1,133.141 MiB**; all
  cases remain below the predeclared 2 GiB per-process budget.
- Nineteen R payloads are byte-identical. Only ashr posterior mean/SD columns
  intentionally differ in the twentieth. Their errors against original ashr
  are `8.88e-16` / `1.61e-15`, within the existing `2e-6` SQP bound; inputs,
  complete fitted state, metadata and all unrelated output columns remain exact.

Other user processes, a virtualization workload and system indexing were active.
The measurement task itself ran panels sequentially with no concurrent fits or
builds. All observations, load snapshots and individual regressions are retained;
none were discarded or retimed to obtain a favorable total. These data do not
establish a causal 3.8% optimization gain or exact speed preservation. No new
optimization round is counted. The main README compares memory with original
DESeq2 only.

See the [current GTEx comparison](gtex-native-wald-benchmark.md),
[native runs](data/final_plan_native_timings.tsv),
[R runs](data/final_plan_r_workflows.tsv),
[R stage costs](data/final_plan_r_stages.tsv), and
[restricted ashr correction report](data/final_plan_ashr_correction.tsv).

## Reproduce the newly measured boundaries

Use native R 4.6.2 with original ashr 2.2-63 and the checked candidate library.
Choose fresh output paths. `R_SCRIPT`, `BASELINE_LIBRARY` and
`CANDIDATE_LIBRARY` below identify locally installed environments, not
machine-specific repository settings.

```sh
python scripts/profile_r_workflows.py --rscript "$R_SCRIPT" \
  --baseline-library "$BASELINE_LIBRARY" \
  --candidate-library "$CANDIDATE_LIBRARY" \
  --validate-ashr-sqp-correction --output-dir "$FRESH_WORKFLOW_DIR"
"$R_SCRIPT" --vanilla scripts/profile_r_boundary.R "$CANDIDATE_LIBRARY" \
  "$FRESH_BOUNDARY_TSV" \
  "$FRESH_WORKFLOW_DIR/batch_results_1200x16_candidate.rds" \
  "$FRESH_WORKFLOW_DIR/batch_results_56200x16_candidate.rds"
```

Strict exact workflow comparison remains the default. The explicit SQP flag
permits only the two named ashr columns after comparison to the independent
original package; it is not a general numerical tolerance for workflow payloads.

All work remains local. Publication, tagging and pushes await a new user
instruction. Current platform success cannot be inferred from a historical
release or a successful local package check.
