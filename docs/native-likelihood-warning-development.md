# Likelihood warnings follow selected operations

The native normal-mixture capture importer preserved recycled component values
but discarded the warning for unequal component-list lengths that are not
multiples. Restricted likelihoods also emitted a full-summary warning sequence
regardless of the requested output subset. Both were behavior gaps in C19/C21;
changing numerical tolerances or optimizers would not address them.

## Shared policy

Capture import now retains one warning-provenance flag when normal-mixture
weight and SD-column lists require nonmultiple recycling. Observation selection
copies that flag alongside the shared parameter handle; it does not copy or
modify parameter arrays. Equal-total-length flattening and exact-multiple
recycling remain warning-free. Component values and their operation order stay
unchanged.

The shared diagnostics module now accepts either restricted expressions or
recycled built-in captures. It follows the original ashr operation order:
initial likelihood fitting, optional flash data, log likelihood/log ratio, then
selected posterior summaries. Numeric output levels and exact named subsets
carry the corresponding diagnostic selections. Duplicate names do not execute a
summary twice. The same selection controls both interfaces.

For the covered recycled likelihood, original ashr emits three warnings for
output level zero, 23 for level one, 74 for level two, 22 for level four and nine
for level five. A mean-only request emits seven; SD-only emits twelve. These
counts describe batch operations, not warnings per observation. Native code
retains their message order and multiplicity without repeating statistical work
just to generate warnings. Ordinary warning-free likelihoods bypass this path.
Likelihood-only mode searches retain warnings for every original trial and the
repeated selected optimum.

`AdaptiveOutputSelection` exposes the additional original diagnostic selectors;
its numerical return fields remain the established posterior mean/SD, s-values
and fitted prior. No public fitting signature, optimizer control, likelihood
kernel, dependency, version or numerical tolerance changes.

## Independent checks

`generate_likelihood_warning_reference.R` runs original ashr 2.2-63 under R 4.6.2.
It records 34 selected-output profiles and two original mode-search traces.
Profiles cover both recycling directions, exact-multiple/equal-total controls,
fixed/fitted priors, missing rows, numeric/named/duplicate output selections and
restricted Laplace expressions with original domain warnings. The restricted
profiles reuse only primitive inputs and expression syntax from an earlier
independent fixture; expected values and warnings are evaluated afresh.

The first before-correction regression requests output level zero and requires
three warnings; the old engine returns none. All 34 selected-output cases now
retain the existing `2e-6` numerical gates and exact warning sequences. The two
mode cases retain the existing mode/posterior/likelihood bounds, exact search
bounds, original evaluation counts and complete warning order. An intermediate
mode-fixture descriptor inherited a null-biased prior label while its original
call used a uniform prior. The descriptor was corrected; all original numerical
outputs, traces and warning arrays stayed unchanged.

Actual native and wrapped requests agree in every exported bit and warning for
all 36 batches with R statistical/search/random routines blocked. Both R CI
lanes run these two conformance checks without requiring cached count fits.
Separate public checks reuse count fits: six ordinary selected-output profiles,
three recycled-likelihood profiles, and six mode profiles retain existing public
precision bounds. Recycled public means/SDs differ from original ashr by at most
`1.12e-16` / `3.89e-16`; output levels one/two/four retain exactly 23/74/22 warnings.
Targeted R package tests also compare warning sequences and selected numbers
against freshly called original ashr.

All 27 accepted replay sets remain unchanged against immutable `a83221d`,
including complete posterior, optimizer, capture, random-state, fresh-matrix and
fresh-formula outputs. All 214 previously committed core reference files remain
unchanged; the warning fixture is additional. Eighteen focused Rust tests, seven
mode/output checks after the trace additions, targeted R tests, both Clippy
scopes, formatting, vendor integrity and the production R numerical-boundary
guard pass. The [local report](data/native-likelihood-warning-local.json) binds
the 628 runtime paths, reference inputs and accepted artifacts.

```sh
Rscript --vanilla scripts/generate_likelihood_warning_reference.R "$FRESH_REFERENCE_JSON"
cargo test --locked -p rsdeseq2 --test ashr_output --test ashr_mode --test location_shrinkage
Rscript --vanilla scripts/check_native_ashr_output.R \
  "$DRIVER" "$R_LIBRARY" "$FRESH_OUTPUT" - crates/rsdeseq2/tests/data/likelihood_warnings_r462.json
Rscript --vanilla scripts/check_native_ashr_mode.R \
  "$DRIVER" "$R_LIBRARY" "$FRESH_MODE_OUTPUT" - crates/rsdeseq2/tests/data/likelihood_warnings_r462.json
```

Pass a cached paired-fit RDS instead of `-` to include the public workflow
comparisons. Count fitting is not repeated for these checks. This correction
closes the demonstrated recycling/output-warning findings, not the full
supported method/option/warning/state inventory or C23 combined precision,
speed/RSS, package/sanitizer/platform and source-bundle requirements. No new
quantification performance claim or optimization round is recorded. Work and
commits stay local; publication remains unauthorized.
