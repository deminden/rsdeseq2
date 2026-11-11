# Completion audit against C01–C23

Current local development routes the default [ashr workflow](shrinkage-component-development.md)
and [negative-binomial apeglm workflow](native-apeglm-development.md) through the
standalone Rust engine. Advanced shrinkage families/methods/options, the remaining
authority audit and final combined checks remain open. The historical ledger
below retains its recorded source scope; C19/C21/C23 are not closed by these
default-method checkpoints. Further development/commits stay local until
the user explicitly authorizes GitHub publication.

The later [general-regression checkpoint](native-regression-development.md)
adds general and `nbinomR` apeglm methods plus restricted custom likelihoods,
with native/R-free output identity and independent original-package posterior
checks. Other advanced methods and the final gates remain open.

The [beta-binomial checkpoint](native-beta-binomial-development.md) adds the
three specialized logit methods, corrects shared BFGS stopping policies and
transports native optimizer warnings. Its original-R warning differences are
recorded explicitly; advanced options and combined completion checks remain open.

The later [random-start checkpoint](native-random-start-development.md) adds
portable native RNG state and both MAP-only random-start methods. Direct/wrapped
posterior, state and error-lifecycle evidence is recorded separately; legacy RNG
modes, advanced options and the final combined gates remain unfinished.

The subsequent [posterior analysis checkpoint](native-posterior-analysis-development.md)
adds shared Laplace/grid intervals, scalar parameter uncertainty, covariance-based
posterior contrasts and grid likelihood overrides. Forty actual native/wrapped
batches agree exactly; original-package comparisons and three explicit source
edge corrections are recorded. Advanced ashr options, warning/method authority
and the final combined C19/C21/C23 gates remain open.

The later [ashr optimizer checkpoint](native-ashr-optimizer-development.md)
adds native ordinary/weighted/variational/compiled-style EM, K=1 controls and
unit/fractional variational priors. Ninety-six complete batches and ninety-seven
optimizer problems agree exactly across native/wrapped interfaces with R math
blocked. Independent categorical tests correct compiled second-probe convergence
and stale/unsafe diagnostics. Remaining likelihood/mode/control/output and final
C19/C21/C23 authority/combined gates stay open.

The subsequent [advanced optimizer checkpoint](native-ashr-advanced-optimizer-development.md)
adds higher-order RRE/MPE cycling and genuine native interior-point fitting.
Fifty-four full and primitive requests preserve numerical/iteration/convergence
bounds and exact wrapped/native identity. Weighted/unweighted 90-decimal KKT
certificates establish convex optima independently; the original Mosek-backed
method was not run. Remaining likelihood/mode/control/output and final combined
C19/C21/C23 gates stay open.

The later [output-selection checkpoint](native-ashr-output-development.md) adds
numeric levels and exact named subsets for the wrapper's posterior fields and
fitted prior. Sixty-seven native/wrapped cases agree exactly; six public profiles
preserve existing precision gates and omitted columns/metadata. The remaining
likelihood/mode/control and full C19/C21/C23 gates stay open.

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

## Subsequent release and migration status

The [0.4.0 release](https://github.com/deminden/rsdeseq2/releases/tag/v0.4.0)
is published at `9324f33`. Its tagged Rust, both R lanes, full parity, MSRV,
checked-source and five target binary checks passed. The original run's final
upload failed because checkout removed downloaded files; `bb7331f` fixes the
ordering and adds guarded recovery. The
[successful recovery](https://github.com/deminden/rsdeseq2/actions/runs/41367071439)
publishes the original checked artifacts without moving the tag. All fifteen
assets, binary/source identities and the crates.io checksum were independently
verified. This closes the platform/publication gate for that released source;
it does not complete the unfinished migration described below.

The development [frozen-local VST follow-up](frozen-local-vst-audit-2026-10-31.md)
now imports saved locfit predictions and translates restricted custom
dispersion expressions into shared Rust. Independent references retain the
existing VST bound; direct complete-output checks agree bit for bit with R
calculations blocked. The historical checkpoint table below retains its
original scope. The subsequent [builder adjustment checkpoint](deferred-result-adjustment-2026-11-10.md)
defers discarded initial BH calculations only at private stages that finalize
row policy, preserving adjusted public constructor/query contracts. Full GTEx
and held-out exports and all twenty paired R payloads remain exact. A combined
native pass measures 1.116% higher throughput, below the 2% counting threshold;
R workflow memory stays below 2 GiB. Optional posterior/s-value calculations,
general method/warning semantics and final combined authority/performance
evidence remain required before C19/C21/C23 can close.

The subsequent [posterior foundation](posterior-primitives-2026-11-10.md) supplies
fixed normal-mixture summaries/tails and ordered s-values with 520 independent
observations and eight sign-rate panels. Native bridge/R-free output bits agree
with R math blocked. Public optional shrinkage routes remain unchanged; prior
fitting, alternate families/likelihoods and the full options/missingness/warning
contract are still open. No unchanged workflow is retimed for unused primitives.

The later [mixture prior fitter](mixture-prior-fitting-2026-11-10.md) supplies
standalone EM/SQP weight fitting, including observation weights and a Dirichlet
penalty. Eight original ashr profiles cover 1,920 observations and complete-fit
bridge/R-free output bits agree. The public optional shrinkage methods still
require their grid/component/likelihood/control and missingness/warning adapters;
this foundation does not close C19/C21/C23 or add a performance claim.

The [normal log-tail support](normal-log-tails-2026-11-10.md) adds direct log
probabilities for the remaining posterior components, with 666 original-R
references and bit-identical native/wrapped output. Existing Wald/public
shrinkage paths and parity bounds remain unchanged. Component/method adoption
and the final audit gates are still open.

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
| C18 | Conditional R replacement output and borrowed native Cook's factors | The [conditional-storage checkpoint](conditional-replacement-storage-2026-10-30.md) completes the local gate: no R replacement matrix for disabled/no-refit paths, moved owned output under R unwind cleanup when needed, and no temporary factor-matrix expansion/clone. Wald/LRT all-zero-only updates retain their counts. Four injected real-fit error paths release owners and recover exactly; all twenty paired R payloads remain identical. Large no-refit RSS falls 1.016 MiB in the combined candidate. The initial R-only attempt's +4.203 MiB outcome is retained. |
| C19 | Borrowed standalone Rust inference, filtering/BH, contrast-zero, Cook/refit, count/shrinkage-row preparation, offset products and formula matrix/frame calculations | The [formula-frame checkpoint](rust-only-formula-frames-2026-10-30.md) removes production R expression, frame/offset and character/logical factor calculations. The [native-parser consolidation](native-parser-shared-evaluation-2026-10-30.md) shares numeric transforms, offsets and matrix products and corrects thirteen independent cases. The previous 189 primitive and 188 expression cases remain exact; the earlier matrix port retains all 94 supported ordered-level references. Native factor/class corrections, warning/method parity, initial adjustment, frozen-local VST and optional posterior calculations remain open. |
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

The [2026-10-30 registry recheck](dependency-recheck-2026-10-30.md) confirms all
15 direct dependencies remain latest stable and updates six newly released
WebAssembly transitive packages in the root lock. Native dependency trees and
the backend lock are unchanged. Locked checks, vendor integrity and crate
packaging verification pass; this inactive update adds no timing claim.

## Version 0.4.0 preparation

The [candidate report](release-0.4.0.md) consolidates subsequent checkpoints into
0.4.0, with new independent factor/class references and the final source and
paired measurement evidence. This version change does not close the unfinished
Rust-only wrapper, general method/warning, C21 authority or C23 combined/platform
audit gates. The earlier statement about unchanged package versions describes the
historical audit above.

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

The [formula-origin correction](cooks-origin-parity-2026-10-30.md) subsequently
resolves the shared kernel's first-tie and original-variable classification
issues using twelve independently generated original-DESeq2 cases. The preceding
implementation disagrees on 20 controlled missingness decisions; the corrected
shared paths match reference p-value/adjustment bits. New fitted origins use
schema 2 with explicit legacy compatibility. The source check passes 1,584
assertions and 113 direct interface cases pass. All 121+9 GTEx gates and 22
unchanged exports pass. Older builder contrast-dependent helper selection,
remaining R arithmetic and the final combined performance audit remain open.

The [native builder correction](builder-cooks-policy-2026-10-30.md) then removes
contrast-dependent Cook dispatch across all 21 builder masking sites. Fifty-six
actual fit/query/context comparisons agree with stored requests, with unchanged
fitted state. Shared metadata preflight preserves alias/shape validation, and
mask-stage BH is omitted before final filtering. The full Rust suite passes
1,055 tests; the R source check remains 1,584 assertions with 113 direct interface
cases. The remaining initial adjustment and Rust-only ports are still required.

The [count/normalization preparation checkpoint](rust-only-preparation-2026-10-30.md)
then shares count validation, raw-zero caching, normal-shrinkage row selection
and offset products with standalone Rust. Independent cases correct weighted
exclusions, unnamed-factor double offsets, anonymous gene identities and an
explicit factor alias. The full Rust suite passes 1,060 tests; the R source check
passes 1,668 assertions and 131 direct interface cases agree. All 121+9 GTEx
gates, 22 native exports and twenty complete paired R payloads are preserved.
One native paired pass improves throughput 1.61%, below the 2% round threshold;
R process total falls 1.41%. Maximum R RSS rises 9.22 MiB, with individual costs
of 89.97 MiB for large batches and 61.03 MiB for the 64-sample rlog retained for
the final combined audit. Eliminating an unused native counts/base-SF matrix
does not complete C18's enabled-but-no-refit replacement-buffer work.

The [formula-matrix checkpoint](rust-only-model-matrix-2026-10-30.md) then replaces
production R model-matrix assembly and expanded list-contrast averaging with
standalone Rust. Forty-eight independent matrices, ten mean-contrast cases and
the exhaustive 2–95 ordered-level sweep preserve exact reference bits; 189 direct
interface cases agree. The source package passes 1,933 assertions and the Rust
suite 1,063 tests. All twenty R workflow payloads and 22 native GTEx exports
remain byte-identical. One native timing pass increases total elapsed 1.26%,
with the blocked random third split retained as the dominant outlier; R time
falls 0.67%, and maximum RSS changes +0.03/−0.53 MiB for native/R. No optimization
round is counted. R expression evaluation and frozen-local VST numerical work
are explicitly still required alongside optional posterior routes and the
previously recorded allocation/performance follow-ups.

The [native formula-transform foundation](native-formula-transforms-2026-10-30.md)
adds standalone finite polynomial/recurrence, scaling and integer-coercion
kernels and uses them in the native formula parser. Its 85 independent R cases
preserve exact values/metadata and invalid-output decisions; 75 of these also
exercise the supported native parser grammar. Integer fractions and range,
large shifted covariates, negative divisors, disabled-centering overflow and
signed-zero behavior receive explicit coverage. The full Rust suite passes
1,065 tests, the checked R package remains at 1,933 assertions, 189 existing
interface identity cases pass, and all 130 previous matrix cases and 22 GTEx
exports retain exact outputs. Production R expression evaluation remains open;
this foundation does not complete C19 or change the other audit dispositions.

Its single paired pass measures native elapsed +0.27% and maximum RSS −2.19 MiB;
the twenty-case R panel measures elapsed −1.05% and maximum RSS +0.84 MiB. All
twenty complete R payloads remain byte-identical and below 2 GiB. Individual
costs, including +8.27 MiB for the large two-group results case, remain recorded.
No optimization round is counted, and the earlier memory/performance follow-up
remains open.

The [conditional replacement-storage checkpoint](conditional-replacement-storage-2026-10-30.md)
then closes C18's local allocation/lifecycle/RSS gate. A temporary owned Rust
output is consumed under R unwind cleanup only when replacement counts must be
retained. Two native factor-matrix materializations are also removed without
changing arithmetic. Independent all-zero-only Wald/LRT cases extend the
existing replacement parity table, and four actual-fit injected errors verify
cleanup and recovery in the instrumented build. The checked source passes
1,988 R assertions, root/backend suites pass 1,065/two tests, 189 direct identity
cases pass, and all 121+9 GTEx gates and 22 exact exports remain intact. The
first attempt's target RSS increase is retained; the combined target pair falls
1.016 MiB with +0.35% elapsed. Native/R totals rise 0.30%/0.54%, all twenty R
payloads are exact and below 2 GiB, and no optimization round is counted. C19,
C21 authority and C23's final combined audit remain open.

The [shared formula-frame checkpoint](rust-only-formula-frames-2026-10-30.md)
subsequently removes production R expression evaluation, missing-row selection,
categorical/response preparation and model-offset summation. One shared Rust
engine serves the actual R bridge and standalone primitive interface. All 188
expression cases (129 preserved plus 59 independent additions) agree exactly,
and 189 earlier interface cases remain exact. The checked source passes 2,445 R
assertions; root/backend suites pass 1,069/three tests. Four injected formula
decode errors release owners and recover under the instrumented boundary tests.
All 121+9 GTEx gates, 22 native exports and twenty complete R payloads remain
unchanged. One paired native pass is essentially flat (elapsed −0.17%); R elapsed
falls 2.13%, dominated by one slower baseline no-refit case, with the other
nineteen cases +0.15%. Maximum native/R RSS is 302.359/1,127.781 MiB; individual
costs and all timings remain recorded. No optimization round is counted.
Legacy parser adoption, warning/method parity, VST/posteriors, initial adjustment
and the final authority/combined audit remain open.

The [native-parser consolidation](native-parser-shared-evaluation-2026-10-30.md)
then routes numeric transforms, coercions, offset addition and matrix products
through the shared engine. Fifty-four actual-parser/R-wrapper cases match
original R after explicit column mapping. A previous-source control disagrees
in thirteen cases (44 matrix and 31 offset entries); those differences are
corrected without changing existing references. Treatment indicators avoid
quadratic contrast storage. The checked source passes 2,445 R assertions; the
root/backend suites pass 1,071/three tests, and fresh sanitizer controls pass.
The prior 189+188 direct cases, all 121+9 GTEx gates, 22 exports and twenty
complete paired R payloads remain exact. Native/R elapsed changes −1.34%/−0.38%,
with essentially flat maxima of 302.063/1,127.594 MiB; individual slower/more
memory-intensive cases remain recorded. No optimization round is counted.
Native factor/class semantics, warning/method parity, VST/posteriors, initial
adjustment, authoritative result routing and final combined/platform gates
remain open.
