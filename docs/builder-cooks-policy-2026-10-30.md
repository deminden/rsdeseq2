# Design-based Cook's policy across native builder routes

This checkpoint extends the shared original-design policy from `428d1ac` across
all older native builder masking routes. Package versions remain 0.3.1. The
[preceding origin correction](cooks-origin-parity-2026-10-30.md) established the
scientific rule using twelve independent original-DESeq2 cases; the broader
[Rust-only migration](rust-only-wrapper-migration.md) remains in progress.

The subsequent [shared preparation checkpoint](rust-only-preparation-2026-10-30.md)
also corrects the explicit no-frame cleaned factor alias and records its
independent regression alongside count, shrinkage-row and offset preparation.

## Confirmed mismatch and correction

On a six-sample `~ condition` model with a two-level factor and Cook's cutoff 0,
requesting the intercept from a fixed-dispersion Wald builder masked a low-count
row. Querying that same fitted state through `ResultRequest` correctly retained
it. The stored raw statistic, maximum Cook's distance and fitted coefficients
were unchanged: the error came from choosing policy through a contrast matcher.
The failing before-change regression is retained locally.

The same contrast-based selection excluded scaled or mixed numeric contrasts,
no-intercept coefficient names, and model frames containing derived or unused
factors. Some numeric LRT paths always used plain masking. Original DESeq2's
eligibility depends on the original formula variable/type, independently of the
displayed coefficient or contrast.

All **21** native builder masking sites now call one private design-policy
adapter. It shares the original model-frame classifier with stored requests:

- A formula's sole original variable must be a factor with two declared levels.
- Derived/unused frame columns do not change that original-variable decision.
- Matrix origins and originally numeric variables remain ineligible.
- Removed terms and offsets still count as original formula variables.
- Explicit factor metadata remains usable when a primitive caller supplies no
  model frame; a supplied frame takes precedence.

This removes 128 lines of coefficient/contrast and reference-alias inference and
ten repeated dispatch branches. Factor-reference and contrast resolution still
run where needed to calculate/report effect sizes; they no longer choose the
Cook's exception.

The old dispatch also performed incidental metadata validation. That validation
is now explicit in shared preflight: names/aliases and sample shape are checked
before normalization or fitting, including when filtering is disabled. The
stored-query classifier validates its supplied frame too. Invalid metadata is
not silently ignored because the requested contrast happens to bypass a helper.

## Adjustment and numerical contracts

The internal builder adapter applies the shared mask without an intermediate
BH adjustment. Every caller immediately runs final independent filtering/BH.
Public `apply_cooks_cutoff*` helpers retain their adjusted-output contract. Initial
row constructors still populate an unused first adjustment; this checkpoint
does **not** claim all redundant native adjustments have been removed.

No likelihood, dispersion, coefficient, covariance, optimizer, threshold,
fitting arithmetic or frozen reference is changed. Intended output corrections
are limited to p-value masking and its downstream adjustment on affected native
query routes. Disabled-Cook's GTEx exports remain exact. Production `.R` files,
C descriptors and unsafe bridge implementation are unchanged.

## Verification scope

`tests/builder_cooks_policy.rs` compares **56** actual fit/query/context
combinations through three compact tests. Coverage includes fixed-dispersion and
native Wald/LRT, coefficient and numeric queries, intercept/scaled/mixed/
no-intercept displays, transformed and unused factors, matrix/numeric/removed
variable negatives, and no-refit replacement paths. Each table matches a stored
request on the same fit. Each fit also matches the filtering-disabled control;
missing optimizer diagnostics are checked by IEEE payload bits before comparing
remaining fitted state, rather than treating NaN inequality as a fit change.

These are route-wiring comparisons, not a new independent R fitting oracle.
The unchanged twelve-case original-DESeq2 row/origin oracle and existing
independent numerical references provide that separate evidence. Existing
actual-refit and all-zero/weight numerical tests remain in the full suite.

The full Rust suite passes **1,055** tests; the backend test, Clippy, formatting,
Rust 1.97.1 minimum-version builds, documentation, vendor integrity and crate
package verification pass. The checked R source package reports **Status OK**,
**1,584 assertions**, no failures, warnings or skips. Its **266** source members
match working files (generated DESCRIPTION fields/DCF whitespace only), and
**435** packaged crate source/test/example files match. The **113** direct
standalone/wrapped interface cases pass with R numerical functions guarded.
Previous C/Rust boundary sanitizer evidence is inherited for the unchanged
boundary; no fresh instrumentation run is claimed here.

All **121** full GTEx precision gates and **9** global gates pass, with minimum
reported target/global error reductions of 664×/597×. All **22** result/diagnostic
exports are byte-identical to `428d1ac`. No tolerance, reference or release output
manifest is widened or regenerated.

## Timing and memory

One paired native pass covers all eleven full GTEx contrasts against immutable
`428d1ac`, and one paired R pass covers all eighteen workflows. Source checks,
compilers and conformance jobs finish before measurement. Both R packages use
schema 2, so the default strict comparison applies without a metadata exception.
Complete raw RDS payloads and process/stage measurements are retained. Existing
checkpoint costs remain in their original reports and measurements are not
repeated to replace an outlier. Main README figures remain unchanged.

Native summed process time is **36.400 → 36.252 s (−0.41%)**. Total
quantification throughput improves **0.41%**, below the 2% round threshold;
**no optimization round is counted**. DESeq2 totals 190.488 s, giving a
**5.25×** DESeq2/Rust time ratio. Maximum native RSS is
**302.125 → 301.969 MiB (−0.156 MiB)**.

The eighteen R workflows total **78.656 → 79.256 s (+0.76%)**. Maximum RSS
is **1,117.344 → 1,117.141 MiB (−0.203 MiB)**, below the 2,048 MiB budget.
All **eighteen complete raw input/fit/output RDS payloads are byte-identical**,
with zero metadata exceptions.

The full-size two-group ten-result stage is **1.538 → 1.509 s**, with complete
process time **4.967 → 4.957 s**. The full-size additive batch process is
**6.044 → 6.043 s**, while its ten-result stage is 1.654 → 1.664 s. Costs are
retained: full-size individual extraction is **7.461 → 7.783 s (+4.32%)**;
2,048 × 32 rlog is **2.698 → 2.860 s (+6.01%)**, with its rlog stage
0.561 → 0.579 s. The largest rlog process is **21.400 → 21.305 s (−0.44%)**,
with stage 18.904 → 18.786 s. These single observations do not establish a
cause, and no unchanged candidate was rerun to replace a result.

See [native summary](data/builder_cooks_2026-10-30_native_summary.tsv),
[native process times](data/builder_cooks_2026-10-30_native_timings.tsv),
[R workflows](data/builder_cooks_2026-10-30_r_workflows.tsv),
[R stages](data/builder_cooks_2026-10-30_r_stages.tsv),
[complete R payload identity](data/builder_cooks_2026-10-30_r_payload_identity.tsv)
and [validation/artifact identities](data/builder_cooks_2026-10-30_validation.json).
The [GTEx gates](data/builder_cooks_2026-10-30_gtex_precision.tsv),
[global gates](data/builder_cooks_2026-10-30_gtex_global_precision.tsv) and
[unchanged exports](data/builder_cooks_2026-10-30_native_output_identity.tsv)
remain separate precision evidence. Local immutable artifacts use
`builder-cooks-*` under `results/codebase_improvements`.

## Remaining scope

The contrast-dependent builder gap recorded at the preceding checkpoint is
resolved for these masking sites. Initial native row adjustment, normal shrinkage
row preparation, remaining offset/normalization/formula arithmetic, optional
ashr/apeglm/posterior calculations and enabled-but-no-refit allocation still need
completion. The final combined performance and source-artifact audit remains
required before claiming the full development plan or Rust-only migration done.
