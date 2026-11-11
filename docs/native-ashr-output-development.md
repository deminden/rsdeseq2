# Native adaptive-shrinkage output selection

This unreleased checkpoint implements the original ashr `outputlevel` selection
for fields consumed by `lfcShrink(type="ashr")` in the shared Rust engine.
Numeric levels, including fractional and infinite values, and exact named
subsets select posterior means, posterior SDs, ordinary s-values and the fitted
prior. Level four selects means/SDs/prior; level five selects none of these.
Unknown names and duplicate names have the original inert behavior. Missing
numeric levels and nonscalar numeric levels reject explicitly. NULL, empty
text and logical selections contain none of the consumed fields.

The standalone API is `AdaptiveOutputSelection` with
`adaptive_shrinkage_with_output`. Absent vectors remain absent rather than being
filled or calculated in R. The existing `adaptive_shrinkage` API retains its
complete-output behavior. The real R adapter only presents selected fields and
optional mixture metadata. Result-table construction tolerates omitted columns.
The ashr `data`, `fit_details`, `flash_data`, `loglik`, `logLR` and `post_sampler`
objects are not exposed by this wrapper; naming them does not select posterior
fields. This change does not add a public standalone `ash` object API or an R
sampling callback.

Positive-threshold s-values still use the original estimates and standard errors
with the unpruned fitted prior and the existing normal/alpha-zero query policy.
They require selected `PosteriorMean` and `fitted_g`; ordinary `svalue` need not
be selected. Missing dependencies now give an explicit error rather than passing
NULL objects into subsequent original wrapper calculations. Without any selected
summary, no posterior is assigned or recycled into the query population.
Cumulative sign/null error-rate summaries and their recycling warnings follow
their requested output selection. Other warning/integration-failure authority
remains part of the unfinished full migration audit.

## Evidence

The independent generator runs only original ashr 2.2-63 under R 4.6.2. Its 67
cases cover numeric boundaries, levels four/five, exact and unknown names,
duplicates, empty/missing names, NULL/logical values, invalid numeric shapes,
missing/NaN/zero/infinite-SE observations and positive-threshold queries.
Original posterior values use the existing `2e-6` bound and preserve missing
value distinctions. Threshold dependency rejection is identified as an explicit
wrapper contract, rather than claimed as an original ashr error message.

Twelve focused connected Rust tests pass; the extended output fixture passes
both output tests. The real R bridge and an R-free executable agree on all
selected output bits and acceptance decisions, with R calculations blocked.
Six public 128-gene profiles cover levels one/four, mean-only, SD-only,
posterior-only without a prior, and thresholded level four. Maximum mean, SD and
s-value differences are `6.67e-16`, `4.42e-15` and `2.23e-16`, respectively.
Existing public `1e-3`/`0.02` precision gates remain unchanged.
The R posterior test file passes 61 assertions, including a single-normal
conjugate check for selected means, SDs and ordinary/thresholded s-values.

Eleven complete accepted payload sets replay exactly against the immutable
`dc5f186` executable, including the previous cycling and interior-point requests.
Previously committed numerical oracles remain unchanged. Core/backend Clippy
and vendor integrity pass. [Local evidence](data/native-ashr-output-local.json)
binds the tested sources, fixtures and artifacts.

```sh
Rscript scripts/generate_ashr_output_reference.R "$FRESH_REFERENCE_JSON"
cargo test -p rsdeseq2 --test ashr_output --test adaptive_result \
  --test location_shrinkage --test ashr_optimizers --test advanced_ashr_optimizers
cargo build --release --example evaluate_r_interface
Rscript scripts/check_native_ashr_output.R "$DRIVER" "$R_LIBRARY" \
  "$FRESH_OUTPUT_DIR" "$PAIRED_FITS_RDS"
```

The optional count-fit snapshot used by the checker is the same 1,200-gene paired
workflow snapshot produced by the native random-start validation; only its first
128 genes are used here. It avoids repeating unchanged count fitting.

Estimated modes, remaining likelihood families and SQP controls, the complete
warning/error/authority audit, and final combined precision/speed/memory/
sanitizer/platform/source-package checks remain open. No new runtime dependency,
quantification timing/RSS claim, counted optimization round or publication is
part of this checkpoint.
