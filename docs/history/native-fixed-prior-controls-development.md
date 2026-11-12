# Shared optimizer controls and fixed-prior dispatch

Fixed-prior ashr requests do not fit mixture weights. Original ashr therefore
never inspects their optimizer controls, although it still validates the selected
optimizer name. The native adapters previously parsed controls eagerly. Valid
fixed-prior calls could fail because of an unused control name, container or
value. The standalone conformance adapter also had a separate parser that
accepted some unknown controls during active fits.

## One control policy

`ashr_mixture_fit_options()` now translates primitive controls in the standalone
Rust library. Both complete location/shrinkage routes, the R argument adapter and
the independent primitive conformance driver call that implementation. SQP,
ordinary/weighted/variational/compiled-style EM and native interior-point control
names, numeric bounds and flag rules are maintained in one module. The former
backend SQP parser and duplicate standalone parser are removed.

The factory accepts the portable expression engine's primitive named-list data;
it neither loads R nor calls R functions. Fixed-prior calls validate the
optimizer identifier and skip unused control values. Active calls retain control
validation. The fitting stage still requires a supplied prior when `fixg=TRUE`.
Low-level typed fitting APIs retain their validation and numerical algorithms.

No likelihood, objective, optimizer arithmetic, stopping setting, dependency,
threshold or numerical tolerance changes. Existing direct conformance error text
is retained when adding the core error context would otherwise change its legacy
transport. Actual public methods and result fields remain available.

## Independent evidence

`generate_fixed_prior_control_reference.R` calls original ashr 2.2-63 under
R 4.6.2. Its 26 requests include three unused-control profiles for each of six
optimizer identifiers, six active-fit rejection controls, an unknown optimizer
and a fixed-prior call without a supplied prior. All 18 fixed-prior profiles
succeed; the eight required invalid requests reject. The original environment
reports an informational EM fallback for `mixIP` without REBayes. The native
interior-point implementation remains distinct; no optimizer runs in these
valid fixed-prior profiles.

The first before-correction test fails on `mixSQP` with an unknown but unused
control. An immutable pre-change standalone driver disagrees with the original
acceptance decisions in twelve cases: seven required successes reject and five
active invalid controls are accepted. The shared implementation corrects all
26 decisions. Existing selected posterior and missing-value gates remain `2e-6`;
warning sequences retain their original expectations.

Actual wrapped and R-free requests agree in every exported numerical bit and
acceptance decision. Both R CI lanes now run this 26-case check with R math
blocked. Six public shrinkage profiles reuse cached count fits and compare to
original ashr, with maximum mean/SD/s-value differences of `4.45e-16`, `2.81e-15`
and `2.78e-17` under the unchanged public bounds. Targeted R package tests also
verify successful fixed-prior calls and retained active-control rejection.

All 29 accepted replay sets remain completely unchanged against immutable
`ffd5000`, including the latest matrix/formula fits, later queries and warning
profiles. A preliminary replay detected only an added error-context prefix in
six primitive rejection messages; the conformance adapter retains their prior
text, and the final full replay is exact. Thirty-five focused Rust tests,
targeted R tests, root/backend Clippy, both format checks, vendor integrity and
the 39-file production R numerical guard pass. All 215 previously committed
core reference files remain unchanged. The [local report](../data/native-fixed-prior-controls-local.json)
binds the runtime sources and validation artifacts.

```sh
Rscript --vanilla scripts/references/generate_fixed_prior_control_reference.R "$FRESH_REFERENCE_JSON"
cargo test --locked -p rsdeseq2 --test ashr_output --test ashr_optimizers --test sqp_controls
Rscript --vanilla scripts/validation/check_native_ashr_output.R \
  "$DRIVER" "$R_LIBRARY" "$FRESH_OUTPUT" - crates/rsdeseq2/tests/data/fixed_prior_controls_r462.json
```

Pass a cached paired-fit RDS instead of `-` to include the six public workflow
comparisons. No unchanged count workflow is retimed or refitted for this check.
The broader supported argument/method/warning/state inventory and C23's combined
current-source precision, speed/RSS, package/sanitizer/platform and source-bundle
checks remain open. This checkpoint makes no performance, release or optimization
round claim. Development and commits remain local.
