# One numerical engine for Rust and R

The user's 2026-10-29 instruction supersedes the earlier permission to retain R
numerical leaves: R is an interface, and the standalone Rust engine must return
identical results to the R-wrapped engine. Production R statistical calculations
and R-specific numerical callbacks are not an acceptable end state.

The local commits `807a107`, `f8a4690` and `bef23cc` were removed from the active
branch. They introduced R-backed Wald probability evaluation, frontend-specific
arithmetic and an R matrix-product backend. Independent input validation,
normalization-state fixes, fitted-state checks, worker tests and benchmark
acceptance fixes were retained. Historical measurements identify their original
artifacts and do not certify a replacement implementation.

## Required architecture

- One implementation of each statistical operation in the standalone Rust
  library. The R bridge calls that same implementation with the same options.
- No callback into R for matrix products, distribution functions, BH, quantiles,
  lowess, moments, offset arithmetic or shrinkage. R cannot choose a different
  arithmetic implementation from the standalone library.
- R may validate/translate arguments, align identifiers, marshal storage, create
  Bioconductor objects and attach descriptions/metadata. Its wrapper must not
  recompute statistics, normalization, thresholds or posterior quantities.
- A native numerical dependency may be used by the Rust engine if it operates
  without R and is used identically through both interfaces. Loading R's math
  runtime or executing R expressions does not satisfy this contract.
- Keep existing supported options and public methods. Preserve independent
  DESeq2 precision bounds and missingness coverage; do not loosen a test to
  conceal a migration regression.

## Migration inventory and order

1. **Establish direct cross-interface evidence.** Compare Rust-only and R-wrapped
   results from identical frozen counts, design, normalization and requests.
   Check every numerical value, missingness, fit diagnostics and filtering
   decision. Keep this separate from the DESeq2 accuracy comparison. Add source
   checks that reject a reintroduced R numerical callback or R statistical path.
2. **Result adjustment and filtering.** Route R table adjustment through the
   standalone BH/filter engine. Move default theta generation, quantiles, lowess,
   residual moments, rejection counts and selection into Rust. Share operation
   order across both interfaces. Sort the filter once per evaluation rather than
   once per quantile, retaining the numerical interpolation contract.
3. **Wald/LRT and contrasts.** Borrow coefficient/SE/covariance primitives for a
   bounded batch. Use the same contrast products, hypothesis policy, probability
   functions, chi-square evaluation and Cook's threshold as Rust-only requests.
   Preserve single-coefficient primitive inputs and contextual errors. Restore
   the batch API using the standalone engine, without R probability callbacks.
4. **Normalization and row policy.** Move geometric-mean conversion, supplied
   size-factor stabilization, factor/offset composition, zero-count contrast
   decisions and Cook's heuristics out of R. Reuse the Rust implementations and
   keep provenance/alignment checks in the wrapper. Audit formula transformations
   separately from object translation so numerical work cannot hide in input
   preparation.
5. **Shrinkage and transforms.** Audit all scale conversion, posterior summaries,
   s-values, prior fitting and threshold processing. Normal shrinkage already
   has a Rust implementation. Existing optional ashr/apeglm routes still call R
   packages: they require standalone implementations with independent posterior
   and tail checks, rather than simply moving their R calls into a native shim.
   Keep existing functionality until the replacement passes its gates.
6. **Reconcile performance, memory and artifacts.** Repeat the appropriate
   paired workflow/holdout checks only after substantive candidate changes.
   Keep the 2 GiB representative-workflow budget and the user's one-pass GTEx
   procedure. A results-only improvement does not count as a 2% total
   quantification round. Test and bind final source artifacts to the accepted
   Rust-only implementation before updating release claims.

## Completion gates

Both interfaces must agree directly, including nondefault thresholds/t tests,
Wald/LRT, replacement/refit, normalization offsets, subsets, batches, transforms
and all supported shrinkage methods. R remains an independent DESeq2 reference
only in tests, fixture generation and comparison tools. The production `.R`
files and bridge must contain no alternate statistical engine or R callback.
Existing independent numerical bounds remain mandatory in addition to equality
between the two interfaces.

The migration is **in progress**. Removing the experimental commits does not
mean the older R calculations listed above have already been removed.

After the 0.4.0 release, the frozen-local VST port imports original saved locfit
knots, residual coefficients and the global polynomial without refitting. Grid
construction, interpolation/extrapolation and integration run in Rust. At the
user's request, restricted custom dispersion expressions also translate into
the same standalone engine; their closures and delayed/active bindings are not
executed. See the [source audit follow-up](frozen-local-vst-audit-2026-10-31.md)
and [wrapper scope](r-wrapper.md). This development work leaves the published
0.4.0 source and binaries intact. The subsequent
[builder adjustment checkpoint](deferred-result-adjustment-2026-11-10.md) removes
discarded initial constructor BH work at private stages that apply final row
policy; public adjusted constructor and convenience-query contracts remain
unchanged. Optional posterior/s-value calculations, general warning/method
semantics and the final authority/performance reconciliation remain open.

The [posterior foundation](posterior-primitives-2026-11-10.md) adds fixed
normal-mixture moments, component probabilities, both tails and ordered s-values
to the standalone engine. Original-package references cover 520 observations
and eight s-value panels; actual bridge and R-free outputs match every bit with
R math blocked. This foundation does not switch public optional methods or
close their prior-fitting, alternate-family, option/missingness and warning
contracts. The inventory in that report preserves their broader migration scope.

The [mixture-weight checkpoint](mixture-prior-fitting-2026-11-10.md) adds a native
EM/SQP fitter over supplied component log likelihoods. Eight independent original
ashr profiles cover 1,920 observations, and complete fits agree bit for bit across
the bridge and R-free process with identical binary inputs. The likelihood offset
and fitting procedure correct a rejected optimizer trial without changing bounds.
Full prior/grid/likelihood/control and missingness/warning adapters still precede
production adoption; the public optional methods and completion gates remain open.

The user's later formula decision disallows custom functions for now; it
supersedes the earlier choice to translate their bodies. The
[built-in expression checkpoint](built-in-formula-expressions-2026-10-30.md)
adds an early custom/masked-call guard and independently checked standalone
expression primitives. That checkpoint left the production model-frame/offset
port open; the guard and evaluator alone did not close C19. The later
[shared formula-frame checkpoint](rust-only-formula-frames-2026-10-30.md)
connects production expression evaluation, row policy, categorical preparation
and offset summation to the standalone engine. R retains symbolic metadata only
on those paths. It adds 188 direct expression identity cases and strict live
whole-frame checks. Legacy parser adoption, general warning/method semantics,
VST/posterior ports and the other completion gates below remain open.

The [native-parser consolidation](native-parser-shared-evaluation-2026-10-30.md)
subsequently routes the compatibility parser's numeric transforms, coercions,
offsets and matrix products through the shared evaluator and matrix plan.
Fifty-four actual-parser/R-wrapper cases agree with original R after explicit
column mapping; a previous-source control disagrees in thirteen cases. Remaining
native factor/class semantics at that checkpoint included ordered coding and
rejection of relevel on ordered factors. The [0.4.0 preparation](release-0.4.0.md)
corrects these through the shared factor/matrix engine, along with unused-level
handling, higher-order contrast policy and product order. Its independent panel
contains 51 matrix/frame cases and nine expanded-factor recoding cases. Expanded
factor recoding for normal shrinkage also moves out of R. General warning/method
semantics, frozen local VST reuse, optional posterior/s-value calculations,
the redundant initial BH pass and final authority/performance audit remain open.
The [frozen-local VST audit](frozen-local-vst-audit-2026-10-31.md) identifies
the saved locfit representation and custom-closure cases that the next port
must address.

The [standalone batch checkpoint](rust-only-result-batches-2026-10-30.md) restores
`resultsList()` using `PreparedContrastView` and bounded native chunks. Both
interfaces share the numerical projection/zero/inference kernels, with 87 direct
bit-equality cases, including 65-request batches. Dataset validation and primitive
extraction happen once per public batch call. Remaining offset, formula and posterior arithmetic above is still required;
the subsequent row-policy checkpoint below removes Cook/refit calculations.

The [first standalone primitive checkpoint](rust-only-primitives-2026-10-29.md)
completes shared filtering/BH and fixes supplied geometric-mean transport, with
20 direct bit-equality cases. The remaining stages above are still required.

The [shared inference checkpoint](rust-only-inference-2026-10-29.md) additionally
moves Wald statistics/tails, primitive LRT probabilities and default Cook's
cutoffs into Rust, with 46 direct bit-equality cases and 1,169 independent
scalar references.

The [shared contrast checkpoint](rust-only-contrast-projection-2026-10-29.md)
additionally moves effect/SE projection into borrowed standalone Rust views and
removes redundant Wald evaluation from native LRT display queries. It passes
64 direct bit-equality cases, a checked R source package and unchanged full GTEx
exports relative to the shared-probability checkpoint.

The [shared contrast-zero checkpoint](rust-only-contrast-zero-2026-10-30.md)
additionally moves count selection, global-zero exclusion and Wald/LRT contrast
zero overrides into Rust, with 85 direct bit-equality cases. Cook's/refit policy,
bounded batches, remaining offsets/formula arithmetic and shrinkage were still
pending at that checkpoint. The later batch checkpoint above restores the bridge;
the other remaining migrations stay open.

The [shared Cook/refit checkpoint](rust-only-cooks-policy-2026-10-30.md) moves
masking, replacement-zero presentation, refit suppression and two-group partition
checks into borrowed standalone Rust. It removes placeholder result construction,
matrix copies and intermediate BH passes from the query bridge. It adds matching
design eligibility to native metadata helpers and automatic formula-context
selection to stored requests, with 101 direct bit-equality cases, 1,510 R
assertions and unchanged full GTEx exports. Normal shrinkage's global-zero check,
remaining offsets/formula arithmetic, optional posterior routes and final
performance reconciliation remain open.

The [formula-origin parity correction](cooks-origin-parity-2026-10-30.md) then
fixes first-tied Cook's maxima and scopes the shared policy by original formula
variables/types, using twelve actual DESeq2 fits/results. New fits persist
validated origin metadata; stored Rust requests and the R presentation adapter
use the same classifier, with 113 direct bit-equality cases. The independently
frozen panel eliminates 20 prior missingness disagreements. Some older native
builder routes still select a helper from the result contrast and need to adopt
the contrast-independent policy. The remaining normalization/formula/shrinkage
ports above continue to be required.

The [builder-policy checkpoint](builder-cooks-policy-2026-10-30.md) subsequently
unifies all 21 native builder masking sites with the original-design classifier.
Intercepts, scaled/mixed contrasts, no-intercept designs and derived/unused frame
columns no longer bypass the low-count exception. Fifty-six actual fit/query/
context combinations match stored requests without changing fitted state. Shared
metadata preflight retains validation; internal masking omits an unused BH pass
before final filtering. Initial row adjustments and the remaining Rust-only
normalization/formula/shrinkage ports remain open.

The [shared preparation checkpoint](rust-only-preparation-2026-10-30.md) then
moves count validation, raw-zero caching, normal-shrinkage row selection and
base-factor × exp(offset) composition into Rust. Both interfaces share the
same kernels, with 131 direct bit-equality cases and a checked source package
passing 1,668 assertions. Independent original-DESeq2 checks also correct
weight-excluded shrinkage and anonymous gene inputs; repeated unnamed-factor
offset fits now preserve effective factors exactly. All twenty paired workflow
payloads and 22 native GTEx exports remain byte-identical on existing valid
inputs. Formula/model-matrix, expanded-contrast and optional posterior/s-value
arithmetic still need standalone paths. Initial native adjustment, C18's
remaining replacement allocation and final combined speed/memory reconciliation
are also still open; the measured individual memory increases remain recorded.

The [formula matrix checkpoint](rust-only-model-matrix-2026-10-30.md) subsequently
moves factor coding, interaction products, matrix assembly and expanded-list
mean weights into shared Rust. Formula matrix construction now uses the same
standalone `ModelMatrixPlan`, with 48 independent matrix references, ten contrast
references, 189 direct identity cases and a checked source package passing 1,933
assertions. All 94 supported ordered-level counts match original R output bits
on the recorded reference environment. Existing GTEx exports and twenty complete
R workflow payloads remain byte-identical. Numerical expression evaluation in
`model.frame`, integration with the older standalone formula parser, frozen-local
VST grids/function evaluation and optional posterior/s-value routes remain open.
The native timing cost is retained for the combined audit; this is no 2% round.

The [native formula-transform foundation](native-formula-transforms-2026-10-30.md)
then corrects integer coercion and replaces uncentered polynomial construction
with the shared polynomial QR operation. Standalone polynomial/recurrence and
matrix-scaling APIs pass 85 independently frozen R cases, including exact
metadata and invalid fitting-output decisions. Native formula parsing uses these
kernels, while production R model-frame expression evaluation still needs to
adopt the shared evaluator. The full previous 130-case matrix sweep and all
22 native GTEx exports remain exact. This numerical foundation does not close
the formula migration or the other remaining items above.

The [conditional replacement-storage checkpoint](conditional-replacement-storage-2026-10-30.md)
subsequently completes C18's local gate: unused R replacement arrays are omitted,
needed arrays are copied from short-lived owned Rust output under R unwind
cleanup, and native Cook's factor temporaries are avoided. Four injected real-fit
errors release owners and recover exactly. All-zero-only Wald/LRT updates retain
their replacement counts; the checked source passes 1,988 assertions, 189 direct
interface cases and twenty exact paired R payloads. The combined target no-refit
pair reduces RSS 1.016 MiB. Its earlier unsuccessful memory attempt remains
recorded. The formula/VST/posterior ports, initial redundant native adjustment
and final combined performance/memory/source reconciliation above remain open.
