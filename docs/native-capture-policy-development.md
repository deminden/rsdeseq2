# Canonical captures and compiled refinement controls

This local checkpoint closes two concrete findings from the
[numerical authority audit](native-numerical-authority-audit.md). The full
C19/C21/C23 contracts and final combined validation remain open. No release,
performance improvement or universal original-R bit parity is claimed.

## Behavior

Canonical `ashr::lik_normalmix` imports now accept stored matrix weights. The
shared Rust engine follows the original constructor's column-major iteration:
when weights and SD captures have equal total lengths, both flatten into one
global mixture, even when their dimensions differ. Otherwise SDs must be a
matrix; component lists recycle in source order. Equal components are not
merged, preserving floating-point addition order. Per-observation SDs retain
their original alignment after missing-row exclusions. These are source
capture rules, rather than per-row mixture-weight semantics.

The compiled negative-binomial methods now interpret controls by method.
`nbinomC` and `nbinomC*` ignore `optim.method`, `bounds` and `log.link`, as the
original MAP-only routines do. `nbinomCR` defers optimizer and bounds use until
prefits disagree, fail to converge, or produce nonpositive posterior variance.
The refinement uses the existing shared Rust BFGS/L-BFGS-B implementation with
the original method's scaling and stopping policy. Its default is BFGS;
previously this path hardwired different L-BFGS-B controls.

For bounded refinement, missing/nonfinite endpoints mean absent bounds; an
empty numeric vector is unbounded, and entries after the first two are ignored.
An unknown optimizer string is harmless if no refinement runs, but rejects
explicitly when required. Other original optimizer families remain unsupported
for actual refinement. An all-zero row with `nbinomCR, log.link=FALSE` rejects,
matching the original failure category with a clearer native error. MAP-only
methods retain their missing-row policy. R imports primitive values and syntax;
all numerical work remains in standalone Rust.

## Validation

The independent generators use R 4.6.2, ashr 2.2-63 and apeglm 1.34.0 in
`deseq2-r462-native`. Candidate outputs do not enter reference generation.
Hard refinement inputs were selected using only original-package calls, and
the fixture records actual original refinement invocation counts. They exercise
both optimizers, tight/empty bounds, deferred unknown controls and rejection.

| Check | Result |
| --- | --- |
| Canonical capture kernels | 8 profiles, including one controlled rejection |
| Complete capture fits | 18 profiles, fixed/fitted priors and missing rows |
| Compiled-method policies | 12 profiles, including two controlled rejections and exact random state |
| Actual native/wrapped identity | All 38 batches agree exactly, with R numerical calls blocked |
| Public workflows on reused count fits | Four matrix-capture and four compiled-method profiles pass unchanged precision/missingness gates |
| Previous accepted outputs | All 21 payload sets unchanged against immutable `f026107` executable |
| Existing frozen oracle files | All 207 files unchanged |
| Focused Rust / R checks | 15 connected integration tests; 115 R assertions; both Clippy scopes and vendor checks pass |
| Production R numerical-call guard | All 39 files pass, with seven rejection controls |

The public matrix profiles differ from original posterior means by at most
`6.67e-16`, SDs by `2.36e-15` and s-values by `1.57e-17`. Compiled public profiles
differ by at most `2.14e-14` in log2 modes, `2.03e-13` in SDs and `2.91e-15` in
s-values. These observed differences do not replace the existing acceptance
bounds. The zero-row rejection and genuine hard refinement cases are covered by
the independent 12-case fixture and actual bridge checks, not by the reused
public count dataset, which contains no all-zero rows.

Reproduction, after building the standalone driver and installing the local R
package into `$R_LIBRARY`:

```sh
Rscript scripts/generate_normalmix_capture_reference.R "$FRESH_CAPTURE_JSON"
Rscript scripts/generate_compiled_method_policy_reference.R "$FRESH_METHOD_JSON"
cargo test -p rsdeseq2 --test normalmix_captures --test compiled_method_policy \
  --test heavy_tail --test regression_posterior --test posterior_analysis
Rscript scripts/check_builtin_kernel_identity.R \
  crates/rsdeseq2/tests/data/normalmix_capture_r462.json \
  target/release/examples/evaluate_r_interface "$R_LIBRARY" "$FRESH_KERNEL_OUT"
Rscript scripts/check_location_shrinkage_identity.R \
  target/release/examples/evaluate_r_interface "$R_LIBRARY" "$FRESH_LOCATION_OUT" \
  crates/rsdeseq2/tests/data/normalmix_capture_r462.json
Rscript scripts/check_native_regression_workflow.R \
  target/release/examples/evaluate_r_interface "$R_LIBRARY" "$FRESH_METHOD_OUT" \
  "$CACHED_PAIRED_FITS" crates/rsdeseq2/tests/data/compiled_method_policy_r462.json identity-only
Rscript scripts/check_native_capture_policy.R \
  "$R_LIBRARY" "$FRESH_PUBLIC_OUT" "$CACHED_PAIRED_FITS"
```

The [local report](data/native-capture-policy-local.json) binds runtime sources,
references, generators, executables and artifacts. The regression reference
assertions are shared without changing their numerical gates. Count fits are
reused; unchanged quantification paths are not retimed at this step. No new
dependency or counted optimization round is introduced. Full warning
multiplicity, remaining option/supported-state contracts, cross-route authority
and the final combined speed/RSS/source-package checks still need completion.

The [likelihood warning follow-up](native-likelihood-warning-development.md) now
preserves the nonmultiple-recycling warning sequences omitted at this checkpoint,
including selected-output and mode-search multiplicity. The original capture
values and all previously accepted payloads stay unchanged. Other supported-state
and final combined gates remain open.
