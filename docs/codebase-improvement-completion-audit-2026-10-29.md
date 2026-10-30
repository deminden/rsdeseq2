# Completion audit against C01–C23

The original [review](codebase-review-2026-10-22.md) and
[implementation plan](codebase-improvement-plan-2026-10-22.md) remain the scope.
This audit supplements the [batch progress report](codebase-improvement-progress-2026-10-29.md).
It corrects its earlier suggestion that all acceptance items were complete.
The implementation remains in progress. No tolerance, frozen oracle or package
version has changed to close a gate.

The user subsequently required [one standalone Rust numerical engine](rust-only-wrapper-migration.md)
for both interfaces and removal of the R-backed commits. The three experimental
commits were removed from the active history; independent audit fixes were
replayed. Older R statistical calculations are still present and must be
migrated. Earlier sanitizer/workflow artifacts below certify their recorded
historical source, rather than the not-yet-built replacement engine.

New comparisons use `deseq2-r462-native`: R 4.6.2, DESeq2 1.52.0,
apeglm 1.34.0 and ashr 2.2-63. Historical benchmark results retain their original
checkpoint and executable identity. The measured M4 executable is SHA-256
`813b8552d62760150bd23dad8c22fc6f4dbe365e87af1f585c2aa0dd2968e678`;
the validation JSON previously listed the test-profile executable instead.
That metadata error is corrected without changing the measurements.

| Item | Implementation / evidence | Audit disposition |
| --- | --- | --- |
| C01 | Immutable binaries, 44 reference hashes, 22 native export hashes, independent Pasilla panel, workflow profiles | Preserved. Workflow defect/assertion inventory added below. |
| C02 | Shared C validators and child-process malformed calls | Added actual C ASan/UBSan and Rust ASan build/run. Scope and exclusions below. |
| C03 | Raw Wald/LRT state, explicit options, constructor schema | Fixed constructors ignoring overrides/unused arguments and inherited LRT p-values after stat/df overrides. Explicit NULL now takes precedence. |
| C04 | Selector/Cook's/filter policy tests and native reference comparisons | Existing supported selector matrix retained; no reference changes. |
| C05 | Normalization provenance and aligned composition | Fixed list-input explicit NULL restoring old factors/weights. Explicit dataset weights now write back; clearing removes the assay. |
| C06 | Shared transform preparation and stale-field cleanup | Existing fit/transform numerical coverage retained. Full R suite rerun after normalization fixes. |
| C07 | Validated mutable fitted view | Added checks for maxCooks, betaConv, replacement flags/eligibility/assays and ambiguous legacy offset provenance. |
| C08 | Mapped gene subsets, subset BH, three shrinkage reference routes | Added full reversal and edited-subset rejection; all fitted consumers reject edited counts. |
| C09 | Relative rank validation, explicit absolute-tolerance compatibility | Added frozen 14-case R 4.6.2 table, independent generator and Rust table test; numerical-boundary decisions documented below. |
| C10 | Early reduced-model subspace validation | Existing nested alternate-basis and non-nested rejection tests retained. |
| C11 | Reusable worker context | Added simultaneous independent contexts, error recovery/pool release, observed original/refit/transform contexts and exact one/four-worker outputs. |
| C12a | Read-only vendor/version checker | Existing drift/failure tests retained; synchronization required for source changes. |
| C12b | Accepted/exploration modes and exact result hash gate | Added mandatory reviewed diagnostics and separate untimed result/diagnostic verification; missing/changed diagnostics fail. Accepted-path replay passed for condition_blocked with unchanged references, result and diagnostics. This smoke replay makes no new performance claim. |
| C13 | Coordinated Rust/R/GTEx/MSRV publication dependencies | Existing local graph/failure checks retained. Remote platform jobs have not run locally. |
| C14 | Preserved workflow/stage/RSS profiles and independent panel | Existing measured evidence retained. |
| C15 | Ownership moves | Existing exact comparisons and measured evidence retained. |
| C16 | Bounded ordered row batches | Existing exact comparisons and measured evidence retained. |
| C17 | Internal covariance retention requirements | Existing exact comparisons and measured transform memory reduction retained. |
| C18 | Disabled replacement buffer omitted | Enabled-but-no-refit allocation remains. The batch checkpoint measures 56,200 × 16 at 940.734 MiB with no retained replacement arrays; its avoidable temporary R buffer payload is 6.860 MiB. No native handle is added and no omitted-allocation claim is made for that path. |
| C19 | Borrowed standalone Rust inference, filtering/BH, contrast-zero and Cook/refit presentation | Shared Cook's masks and replacement-row values now pass 101 direct bit cases, a 1,510-assertion source check, instrumented boundary tests and unchanged GTEx exports. Native eligibility now requires a matching two-column design; stored requests infer matching formula context. See the [Cook/refit checkpoint](rust-only-cooks-policy-2026-10-30.md). Exact Cook's tie/formula-origin parity, remaining shrinkage zero flags, offsets/formula and posterior calculations remain open. |
| C20 | Standalone-engine bounded result batches | `resultsList` validates a dataset/extracts primitives once, borrows shared matrices and evaluates chunks of at most 64 requests through `PreparedContrastView`. Single calls use the same kernel/adapters. Exact per-request policies, 87 direct bit cases, source package and instrumentation pass; see the [batch checkpoint](rust-only-result-batches-2026-10-30.md). C19's remaining presentation migration still applies. |
| C21 | Private execution and fitted-result modules, public forwarding APIs | Implemented portions retained; authoritative shared result path awaits C19/C20 completion. |
| C22 | Preflight, staging, hash-verified completion record | Existing collision/failure/numerical tests retained. |
| C23 | Checkpoint-bound docs and source artifacts | Historical artifacts preserved. A final tested source package and documentation reconciliation are required after outstanding implementation. |

The subsequent [stable dependency refresh](dependency-refresh-2026-10-29.md)
uses locfit 0.2.0, L-BFGS-B 0.3.0, statrs 0.19.1 and sha2 0.11.1 across both
Rust dependency graphs. All 15 direct dependencies are latest stable as audited;
upstream-constrained transitive versions are listed explicitly. Numerical exports
and all 14 paired R workflow payloads remain byte-identical to `a6e376a`.
Single-pass native/R times increased 1.12%/1.25%, with essentially flat peak RSS;
no optimization round is counted. This does not complete C19/C20 or the remaining
Rust-only migration.

## Workflow defect and assertion inventory

| Defect / contract | Interfaces | Independent reference | Exact assertion / test |
| --- | --- | --- | --- |
| Display Cook's masking must not destroy raw inference | Fitted default/name/contrast and Wald/LRT | DESeq2 results with enabled/disabled/numeric Cook's cutoff | `test-result-policy.R`: identical missingness, raw p-value recovery, per-route hypothesis results |
| Constructors must honor explicit arguments | Primitive Wald/LRT constructors | Supplied beta/SE; R pnorm/pchisq | `test-primitive-routes.R`: changed LFC and statistic; recomputed LRT df p-value; explicit NULL clearing; named/unnamed unused arguments error |
| Offsets compose once after alignment | Dataset/list fit, refit, transform preparation | Explicit normalization-factor matrix and equivalent offset-free design | `test-normalization-state.R`, `test-fit-integration.R`: exact effective factors/results, repeated-fit idempotence and transform comparisons |
| Explicit NULL clears inherited normalization/weights | Dataset/list fit | Fresh fit with absent factors/weights | `test-normalization-state.R`: exact result equality and absent weights assay |
| Cached state must reflect mutable dataset | Results, resultsNames, normal/ashr/apeglm shrinkage, design-aware transforms | Original cached inputs plus direct reference shrinkage | `test-fitted-view.R`, `test-fit-integration.R`: same-size edits reject, reversal/subset values map correctly, subset BH recomputed |
| Replacement metadata and arrays must be coherent | Fitted Wald/LRT queries | Existing actual replacement fixtures | `test-fit-integration.R`: replaceCounts/replaceCooks edits reject; original numerical reference assertions retained |
| Native buffers must be validated before reads/writes | Registered .Call functions | Known malformed shapes/types/ranges and valid matrix arithmetic | `test-native-boundary.R`: child survives errors, valid normalization and Normal/t tails exact, disabled buffer absent |
| Rank validity should ignore units without accepting aliases | Rust design/fit entry points | Frozen R 4.6.2 default and scaled QR ranks | `model_validity.rs`, `data/rank_r462.tsv`: frozen relative decisions, original matrix bits unchanged, explicit absolute tolerance retained |
| Requested workers control complete workflows | Rust builder, refits, VST/rlog | Simultaneous bounded pools and serial outputs | `execution_control_tests`: observed stage contexts, exact one/four outputs, recovery after error and weak pool release |

## Instrumented native-boundary evidence

`scripts/check_native_boundary_sanitizers.py` builds a disposable package with
C AddressSanitizer/UndefinedBehaviorSanitizer and Rust AddressSanitizer, verifies
the Rust build flags and native sanitizer symbols, then runs the malformed-call
child test and its valid numerical controls. The successful local run is
`results/codebase_improvements/native-boundary-sanitizers-complete/`.
Its `sanitizers.json` reports `passed`; `tests.log` contains the child-test result.
The explicit Rust target avoids instrumenting host compiler proc macros. Direct
interpreter entry preserves the macOS sanitizer environment in the test child.

The R interpreter, installed R dependencies and Rust standard library were not
rebuilt with instrumentation. Leak detection is disabled. This proves the stated
instrumented bridge scope, not whole-process leak freedom or all-platform safety.

## Frozen rank decisions

`scripts/generate_rank_reference.R` requires R 4.6.2/DESeq2 1.52.0 and emits
`crates/rsdeseq2/tests/data/rank_r462.tsv`. It records both R's default QR rank
(`tol=1e-7`) and QR after maximum-absolute column scaling at `tol=1e-10`.
The Rust relative policy matches the latter for all 14 frozen cases. No fitting
design, coefficient unit, optimizer bound or tolerance changed in this audit.

The controlled `1e-9` near-intercept and near-covariate cases have full rank under
the existing relative policy, while R's default QR calls them deficient. At
`1e-12` both policies call them deficient. These deliberate boundary differences
are distinct from the independent `1e-12`-scaled covariate case, which is full
rank in both reference policies. Structural rank does not promise an adequately
conditioned or numerically estimable model.

The [standalone batch checkpoint](rust-only-result-batches-2026-10-30.md) also
records one accepted paired eleven-contrast pass: native elapsed +0.36%, maximum
RSS 302.328 MiB, 5.40× DESeq2 process-time ratio and no optimization round. All
17 paired R workflow payloads are byte-identical and below 2 GiB. Ten-result
batch extraction is 1.71× faster; the existing 14 R controls sum to +2.58% time,
mostly the largest rlog (+6.26%). Those observed costs remain in the final
performance follow-up; the measurements were not repeated or relabeled.

The subsequent [Cook/refit checkpoint](rust-only-cooks-policy-2026-10-30.md)
retains all prior records and adds one paired comparison against `a7b256f`.
Native elapsed is +2.08%, with essentially flat RSS; the first blocked-condition
pair accounts for the aggregate increase. No optimization round is counted.
All 18 R workflow payloads are byte-identical and below 2 GiB; the targeted
two-group ten-result stage improves 14.84% and peak RSS falls 27 MiB. The R
process total changes +0.19%, with slower controls reported explicitly. The
observed native timing cost and final combined performance audit remain open.
