# Large-shape numerical corrections and wrapper guard

This unreleased checkpoint corrects independent large-count/log-F probability
and density failures in the shared Rust engine. The affected kernels do not run
an optimizer. No `rcompat-lbfgsb`, optimizer control, dependency lock or existing
precision bound changes here. The full development plan remains active.

## Observed defect and correction

At Poisson/binomial counts of `1e8`–`1e10`, direct subtraction of large log-gamma
terms lost density precision, and fixed series/fraction iteration limits could
return NaN probabilities and moments. The accepted `32a765c` executable fails
299 value checks in the new expanded panel. This includes finite original-R
probabilities/moments becoming NaN, not merely last-bit disagreement.

Large-shape gamma and beta kernels now integrate the exact distributions in
scaled coordinates. Stable `log1pmx` and existing Stirling corrections normalize
the density without subtracting large log gammas. Beta integration uses a
centered logit coordinate, reflecting unequal shapes to avoid cancellation;
ordinary/log probabilities retain their existing caller contracts. Small/unequal
shape normalization uses a stable gamma-ratio expression. Large log-F densities
share the same normalized equations. These are changes of integration coordinates,
not normal approximations. The equations follow
[DLMF incomplete gamma definitions](https://dlmf.nist.gov/8.2),
[incomplete beta definitions](https://dlmf.nist.gov/8.17) and
[Stirling normalization](https://dlmf.nist.gov/5.11).

The implementation reuses the existing native Gauss–Legendre tail integrator;
its visibility is expanded without changing its arithmetic. The new routine is
isolated in `shrinkage/large_shape.rs`. Existing ordinary-shape paths retain their
arithmetic. The independent symmetric-beta half-mass check extends to shape
`1e14`. No R routine or numerical runtime is called by the standalone executable.

## Evidence and reproduction

The original-only R generator supplies 41 kernel profiles and 287 rows: Poisson
identity/log, binomial identity/logit, log-F and its infinite-df limits, the new
branch boundary, up to twelve standardized tail units, and unequal beta shapes.
An independent Python generator evaluates density equations with mpmath at 100
decimal digits. It records 90-digit values and the input reference hash; candidate
outputs do not enter either generator.

| Check | Result |
| --- | --- |
| Original kernel log CDF/density bound | Unchanged `2e-8 + 2e-12 * abs(reference)` |
| Original kernel moment bound | Unchanged `2e-6 + 2e-12 * abs(reference)` |
| Accepted-baseline failing value checks | 299 |
| Corrected failing value checks | 0 |
| Maximum original log-CDF difference | `1.034364e-8` |
| Maximum original log-density difference | `2.395666e-10` |
| Maximum independent 100-digit density difference | `1.064394e-10` |
| Maximum original first/second-moment difference | `2.219336e-13` / `1.110224e-13` |
| Actual R-wrapped/standalone kernels | All 41 batches agree bit for bit with R math blocked |
| Earlier accepted replay | Twenty complete payload sets unchanged against `32a765c` |

Focused checks include 32 Rust integration tests, native quadrature and bounded
optimizer controls/gradient/hard GTEx replays, both Clippy scopes, vendor integrity
and 111 R assertions. The six frozen hard GTEx optimizer cases still satisfy
existing objective (`5e-6`) and coefficient (`2e-3`) comparison bounds. Those
checks establish bounded agreement, not universal bitwise optimizer parity.

```sh
Rscript scripts/references/generate_large_shape_likelihood_reference.R "$FRESH_REFERENCE_JSON"
python scripts/references/generate_large_shape_density_truth.py "$FRESH_REFERENCE_JSON" "$FRESH_TRUTH_JSON"
cargo test -p rsdeseq2 --test large_shape_likelihood
Rscript scripts/validation/check_builtin_kernel_identity.R \
  crates/rsdeseq2/tests/data/large_shape_likelihood_r462.json \
  target/release/examples/evaluate_r_interface "$R_LIBRARY" "$FRESH_IDENTITY_OUT"
Rscript scripts/validation/audit_r_numerical_boundary.R r-pkg/rsdeseq2/R "$FRESH_AUDIT_JSON"
```

The [local report](../data/native-large-shape-local.json) binds sources, generators,
references, executables and comparison artifacts. Quantification fitting is not
rerun for unused optional kernels. No new timing/RSS claim, optimization round,
release or publication is made.

## Wrapper guard and remaining authority work

The unused `.rsdeseq2_low_count_cooks_spares_row` R statistical helper is removed;
active Cook's decisions already use the shared native policy. The lexical guard
parses all 39 production R files and rejects direct statistical/dynamic calls,
including numerical functions passed as common apply callbacks. Seven negative
controls verify rejection. Quoted formula/likelihood templates are inspected as
syntax data, not executed. The R workflow now schedules this guard; its exact command passes locally.
Remote CI has not been run.
Its scope is explicitly lexical: runtime firewalls and full interface tests
remain required, and reviewed metadata getters remain permitted.

The [authority audit](native-numerical-authority-audit.md) lists confirmed pending
method/capture policies. C19/C21/C23 remain open until those contracts and the
combined precision/speed/memory/sanitizer/platform/source-bundle gates are closed.
