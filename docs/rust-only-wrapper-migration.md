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
