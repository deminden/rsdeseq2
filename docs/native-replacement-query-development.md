# Later requests after replacement

`DeseqFit::results_for_request_with_replacement` evaluates a validated request
on the fit and plan returned by `CooksReplacementTestOutput::into_merged_fit`.
Pass the original counts alongside that retained plan. Ordinary
`results_for_request` calls use the same evaluator without replacement context;
existing signatures remain available.

The plan supplies post-refit Cook maxima when refitting occurred. Shared
`ReplacementRows` and `CooksPolicy` restore the established new-zero
effect/SE/stat/p values of 0/0/0/1, suppress repeated outlier masking on refitted
rows, and then apply each request's own independent filtering/BH. Stored fit
matrices and Cook diagnostics are not changed by a query. Row indices are borrowed
and validated for range and strict order, avoiding gene-length flag copies.

Sample eligibility matters: a plan can retain candidate outlier rows even when
no sample is eligible for replacement. Such rows must still undergo ordinary
Cook masking. The plan, R export and standalone query share an eligibility
accessor; the new query path does not treat those candidates as actual refits.
The underlying low-level planner's candidate metadata remains available.

## Missing-value transport

Native optional results represent absent values as `None`. The R interface
supplies its NA marker as a transport primitive; active Cook masking and stored
LRT outputs also retain their established NaN representation. The conformance
driver receives the same marker explicitly and compares raw output bits,
including missing values. It does not replace different missing values with a
single comparison token.

Adjusted-value missing formatting is now a shared Rust helper used by the R
adapter and standalone transport. It preserves the incoming missing payload for
a row retained in the filtering population and otherwise uses the caller's
marker. Population selection, p-values and adjustment decisions are already
calculated by the engine. No R numerical calculation or callback is introduced.

## Checks and limits

The extended fresh-fit checker reuses each fitted model for forward, reverse
and fractional contrasts; automatic, disabled and zero Cook cutoffs; different
alpha/filtering settings; and Wald threshold/t requests. All seven query columns
and all independent-filtering numerical metadata are compared. The local run has
**24 fresh-fit requests and 86 later queries**, with 20 successful fits and four
controlled fitting rejections. Every compared floating-point bit, logical value
and fitting acceptance decision agrees across standalone Rust and the actual R
wrapper. Existing CI invokes this extended checker in both R lanes.

R quantile, lowess, adjustment and distribution functions are blocked during
these calls. Three additional runtime negative controls demonstrate that the
quantile/lowess/adjustment firewall rejects execution. The production lexical
guard passes all 39 R files and its eight negative controls.

All 189 existing primitive interface cases, 78 connected Rust tests, root/example
and backend Clippy, formatting and vendor checks pass. The Rust checks include
independent numerical references, invalid replacement indices and preservation
of distinct missing payloads. All 211 previously committed core reference files
are unchanged. This work adds no fitting arithmetic, solver control, dependency
or tolerance change and makes no new speed/RSS claim.

The [local report](data/native-replacement-query-local.json) binds accepted
sources and artifacts. C19/C21/C23 remain open: factor/formula query context,
supported options/warnings/states and final combined acceptance still need their
full audit. An independent original-DESeq2 probe confirms a separate state
defect: after actual replacement, an unrelated weight-rank-failed row's `allZero`
flag resets to raw-count status in original DESeq2, while the current merged fit
retains its initial skip flag. The p-value remains missing in both. This finding
requires its own fitting-state correction and independent reference coverage;
it is not hidden by cross-interface agreement. No publication is authorized.
