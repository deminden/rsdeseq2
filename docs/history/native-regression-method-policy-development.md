# Regression method argument policy

This local follow-up extends the [compiled-NB control correction](native-capture-policy-development.md)
across general and beta-binomial posterior methods. It closes another concrete
dispatch gap in the [authority audit](native-numerical-authority-audit.md), while
C19/C21/C23 and final combined validation remain open.

## Shared policy

The standalone Rust `RegressionOptimizer` now retains an unsupported optimizer
until fitting actually needs it. Successful beta-binomial prefits ignore unused
optimizer/bounds controls; genuine CR refinement still rejects an unsupported
optimizer. General and beta-binomial source dispatch also ignores `log.link`.
The R adapter transports these options without evaluating statistics.

One core helper now normalizes scalar bounds for both interfaces and compiled-NB
refinement. Empty vectors and missing/nonfinite endpoints are unbounded; only
the first two entries matter. BFGS ignores bounds. Original apeglm applies bounds
only for the exact `L-BFGS-B` string: its abbreviated `L-B` name still selects
that optimizer inside R, but the surrounding dispatcher discards the bounds.
The new independent profiles preserve that distinction.

No BFGS/L-BFGS-B objective, tolerance, scaling, line search or dependency changes.
Other optimizer families remain unsupported when an actual optimization needs
them. Solver-dependent warning differences on older saturated/random panels
remain documented rather than synthesized to match the original implementation.

## Evidence

The original-only generator uses R 4.6.2 and apeglm 1.34.0. It reuses independent
frozen inputs and defines a balanced interior-count row to isolate valid
beta-binomial prefits. That row was selected solely using original-package
diagnostics. Separate saturated rows exercise actual original refinement.
Tracing original `optimBetabin` records zero calls for valid prefits and positive
calls for refinement cases. No candidate output enters reference generation.

| Check | Result |
| --- | --- |
| Independent method profiles | 21: 17 valid, four controlled rejections |
| Previously rejected valid standalone requests | 15 corrected against immutable `6ca3f27` driver |
| Native/R-wrapped outputs | All 21 batches exact, including random state and warnings; R math blocked |
| Original numerical/missingness/warning-event gates | All profiles pass unchanged bounds; all 17 valid warning-message observations also match |
| Accepted replay | All 24 earlier payload sets unchanged |
| Existing oracle files | All 209 unchanged |
| Focused checks | 13 connected Rust integration tests; 76 R assertions; both Clippy scopes and vendor/format checks pass |

The R tests also pass malformed values for controls that source dispatch ignores:
unused text bounds, an unknown optimizer string and a multi-element `log.link`.
They require finite modes and exact complete-output equality with default calls,
guarding frontend rejection without repeating the independent numerical table.

```sh
Rscript scripts/references/generate_regression_method_policy_reference.R "$FRESH_REFERENCE_JSON"
cargo test -p rsdeseq2 --test regression_method_policy \
  --test compiled_method_policy --test regression_posterior \
  --test heavy_tail --test posterior_analysis
Rscript scripts/validation/check_native_regression_workflow.R \
  target/release/examples/evaluate_r_interface "$R_LIBRARY" "$FRESH_IDENTITY_OUT" \
  "$CACHED_PAIRED_FITS" crates/rsdeseq2/tests/data/regression_method_policy_r462.json identity-only
```

The [local report](../data/native-regression-method-policy-local.json) binds sources,
reference generation, immutable binaries and validation artifacts. Runtime
sources still cover 625 paths. No count fitting, unchanged quantification
retiming/RSS run, counted optimization round or publication occurs here. Remaining
method/state/warning contracts and final shared-route/combined checks stay open.
