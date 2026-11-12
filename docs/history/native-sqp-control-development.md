# Native SQP numerical controls

This unreleased checkpoint completes numerical translation of the original
mixsqp control names in the local ashr adapter. `MixtureFitOptions::sqp` exposes
the added native policies through `MixtureSqpOptions`. The shared engine handles
active-set convergence, small search directions, sufficient decrease, step-size
reduction/minimum, Hessian regularization growth, automatic active-set budgets
and explicitly requested truncated SVD. The default native path and accepted
full-matrix arithmetic remain unchanged.

Backtracking is bounded to 10,000 probes and stops if multiplication cannot
reduce the step. Near-unit reduction factors therefore report nonconvergence
instead of spending an effectively unbounded time in one SQP iteration.
Augmented matrix dimensions are checked before reserving spectral workspace.

| R control | Shared numerical policy |
| --- | --- |
| `convtol.activeset` | Bound-multiplier convergence tolerance |
| `zero.threshold.searchdir` | Stop the quadratic subproblem at small directions |
| `suffdecr.linesearch` | Armijo sufficient decrease |
| `stepsizereduce`, `minstepsize` | Backtracking factor and minimum step |
| `identity.contrib.increase` | Growth of the Hessian diagonal correction |
| `maxiter.activeset=NULL` | Resolve the budget to `min(20, components+1)` |
| `tol.svd` | Optional spectral truncation of the augmented likelihood |
| `normalize.rows` | Validated; likelihoods are already row-normalized upstream |

The existing `eps`, `numiter.em`, `maxiter.sqp`, `maxiter.activeset`,
`convtol.sqp` and `zero.threshold.solution` controls remain supported. The R
bridge and R-free conformance transport share the same core policies. The bridge
control translator is separated from argument/family handling. Informational
`verbose` input remains accepted; original console-progress and complete
warning/error formatting are still part of the authority audit.

## Spectral fitting and source policy

Absent `tol.svd`, or a zero value, retains the accepted full-matrix route. With a
positive explicit threshold, native nalgebra SVD decomposes the row-normalized
likelihood augmented by nonzero Dirichlet pseudo-observations. At least two
singular directions are retained. Rank discovery preserves original doubling
and the fallback used when the partial solver requests more than half the
smaller matrix dimension. Four or fewer components use the full matrix.
Reconstructed negative entries adjust the stabilizing offset, as in original
mixsqp. Final likelihood and Dirichlet-penalty diagnostics use the original,
unapproximated observations.

The native decomposition is deterministic. An independent LAPACK SVD reference
is recorded alongside, rather than replacing, unmodified original IRLBA results.
Native spectral fitting does not consume R random state; the existing posterior
random-start state protocol remains separate and unchanged.
Both references pass the same `2e-6` numerical bounds. A rejected native trial
incorrectly added the stabilizing offset to low-rank EM updates. Independent
original one-step probes show that EM omits this offset even though the SQP
objective includes it; the implementation now preserves that distinction.

The numerical implementation is independent Rust and nalgebra code, without
copying/linking the R/C++ optimizer. Control meanings follow the
[authors' reference](https://stephenslab.github.io/mixsqp/reference/mixsqp.html).
SVD reconstruction currently feeds the dense native objective/Hessian; this
feature does not claim original mixsqp's compact-factor performance benefits.
No runtime dependency is added.

## Evidence and reproduction

Original references use R 4.6.2, ashr 2.2-63, mixsqp 0.3-54 and IRLBA 2.4.1.
Thirty-three full profiles and thirty-three primitive problems cover every
numerical control, combinations, finite iteration budgets, automatic limits,
asymmetric grids, Dirichlet priors and six invalid controls. Twenty-seven valid
profiles cover 1,296 observations. Genuine rank-seven/rank-six truncation cases
exercise likelihood-only and augmented-prior decomposition; small matrices also
exercise the original fallback. Original weights, posterior means/SDs/s-values,
likelihoods and penalized likelihoods retain `2e-6` gates. Convergence decisions
match. An analytical first-EM state proves the direction cutoff takes effect;
an independently constructed rank-two matrix verifies spectral reconstruction.

Twenty connected Rust tests, the spectral unit test and 70 R posterior
assertions pass. Both complete batch results and raw optimizer diagnostics agree
bit for bit across the real R bridge and separate R-free executable for all 33
requests, including rejected inputs. R statistical/SVD functions are blocked
with a positive control. Six public 256-gene workflows retain the existing
`1e-3`/`0.02` gates; maximum mean/SD/s-value errors are `8.89e-16`, `7.25e-15`
and `1.23e-17`.

Thirteen complete accepted payload sets replay exactly against immutable
`5865f7e`, including estimated modes, selected outputs and earlier posterior and
regression cases. Existing frozen oracles and tolerances are untouched.
Core/backend Clippy and vendor integrity pass. [Local evidence](../data/native-sqp-controls-local.json)
binds source, references and artifacts.

```sh
Rscript scripts/references/generate_sqp_control_reference.R "$FRESH_REFERENCE_JSON"
cargo test -p rsdeseq2 --test sqp_controls --test mixture_fit --test ashr_mode \
  --test ashr_output --test ashr_optimizers --test advanced_ashr_optimizers --test location_shrinkage
cargo test -p rsdeseq2 --lib rank_two_likelihood
cargo build --release --example evaluate_r_interface
Rscript scripts/validation/check_location_shrinkage_identity.R "$DRIVER" "$R_LIBRARY" \
  "$FRESH_BATCH_OUT" "$FRESH_REFERENCE_JSON"
Rscript scripts/validation/check_native_ashr_optimizers.R "$DRIVER" "$R_LIBRARY" \
  "$FRESH_OPTIMIZER_OUT" "$PAIRED_FITS_RDS" "$FRESH_REFERENCE_JSON"
```

The paired-fit snapshot is the existing native random-start workflow output;
count fitting is not repeated. Remaining built-in likelihood families, complete
method/warning/error/authority and final combined precision/speed/memory/
sanitizer/platform/source-bundle gates remain open. No quantification timing/RSS
claim or counted optimization round is added. All work remains local.
