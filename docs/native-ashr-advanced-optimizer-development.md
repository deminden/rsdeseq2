# Local cycling and interior-point mixture fitting

This unreleased checkpoint extends the [native EM implementation](native-ashr-optimizer-development.md)
with cycling reduced-rank extrapolation (RRE), minimal-polynomial extrapolation
(MPE), squared/unsquared polynomials, and a separate native interior-point solver.
It does not complete the remaining likelihood/mode/control/output migration or
the final C19/C21/C23 gates. All development and artifacts remain local.

## Shared numerical engine

`SquaremOptions::polynomial_order` and `polynomial_method` select higher-order
cycling. The default order remains one and preserves the accepted scalar path.
R imports `control=list(K=2,method="mpe",square=TRUE)`. RRE solves the
normalized Gram residual problem; MPE solves a least-squares annihilation problem.
Squared extrapolation convolves coefficients over the stored iterate sequence.
Stabilization, objective safeguards, fallback probes, initialization and budget
ordering follow the original method contract. Orders exceeding the retained
component count fail explicitly. Workspace is bounded by the requested order
and component count; dimension arithmetic is checked before allocation.

The implementation uses native nalgebra QR and the polynomial equations in
[Roland, Varadhan and Frangakis (2007)](https://doi.org/10.1007/s11075-007-9094-2).
No GPL SQUAREM/ashr implementation is copied or linked. An affine contraction
with a known fixed point independently checks both polynomial methods and both
squaring policies. The compiled-style method retains its original scalar policy
when K is supplied. Its numeric method control is validated; a character method
fails just as in original ashr, even when the compiled calculation ignores the
numeric method selection.

`MixtureOptimizer::InteriorPoint(InteriorPointOptions)` minimizes the actual
unregularized convex mixture objective using a positive log-barrier path and
scaled Newton equations. Positive-boundary steps, a decreasing-barrier schedule
and Armijo search preserve feasibility. Scaling by current weights avoids
inverse-square weight conditioning. The default SQP optimizer is not called by
this route. Background on barrier/Newton methods is available in
[Boyd and Vandenberghe's Convex Optimization](https://web.stanford.edu/~boyd/cvxbook/).

Native `optmethod="mixIP"` runs this solver directly, independently of optional
R solver packages. `control=list(rtol=...)` specifies an upper bound on the
normalized convex dual gap. A mass-scaled accuracy target can tighten it to
protect likelihood precision. The native API also bounds Newton/centering
iterations. Exhaustion remains visible as nonconvergence; objective evaluations
and iterations are reported separately. Returned-point residuals are refreshed
on exhaustion. Zero observation weights with positive prior mass remain valid;
Dirichlet parameters below one require variational fitting.

## Independent comparisons and their scope

The original-package environment is R 4.6.2, ashr 2.2-63 and SQUAREM 2026.1 in
`deseq2-r462-native`. Fifty-four full requests and fifty-four captured primitive
problems cover orders two/three, both polynomials, squared/unsquared variants,
ordinary/weighted/variational/compiled-style methods, low budgets, tight and
monotone controls, and weighted/unweighted interior fitting. Six original
compiled-method type errors remain explicit rejections. Valid full requests
cover 2,304 observation rows.

Existing component/weight/concentration/summary/likelihood bounds remain `2e-6`.
Original cycling convergence decisions and iteration counts match. Nonconvergence
warnings match; the previously documented original compiled diagnostic bounds
warning is omitted by safe native code. Fourteen focused Rust integration tests
and the independent affine test pass, together with core/backend Clippy and
vendor integrity.

**The original Mosek-backed mixIP method was not executed.** REBayes is absent
from the reference environment. Interior-point comparison therefore uses original
unregularized mixSQP fits and a separate 90-decimal mathematical oracle, with
that provenance explicit in both fixtures and reports. The mathematical generator
reads only frozen original likelihood/input bits, enumerates active faces, solves
the likelihood score equations, and verifies all active and inactive KKT
conditions. Concavity certifies the global maximum independently of the native
solver and original SQP. Weighted/unweighted results satisfy `2e-7` weight and
likelihood bounds against those certificates. Categorical/Dirichlet solutions,
prior-only fitting, invalid priors/tolerances and exhausted budgets add independent
coverage. This evidence does not claim direct parity with Mosek's implementation
or its stopping diagnostics.

Actual R/native and separate R-free processes agree bit for bit for every full
request and primitive optimizer result/error, with R statistical calculations
blocked. Five public 256-gene workflows cover RRE, MPE, weighted/variational
cycling and interior fitting. Maximum mean/SD/s-value errors are `8.81e-8`,
`2.65e-8`, `7.12e-10`; the interior profile's mean error is `6.16e-11`. Public
`1e-3`/`0.02` gates are unchanged. Direct original-Mosek parity remains unmeasured.

Immutable `d033995` replay preserves nine accepted payload sets exactly,
including all 97 previous EM/diagnostic problems and 96 previous complete EM
batches, default location/query results, general/NB/beta-binomial/random-start
and posterior-analysis requests. All previously committed numerical oracles
remain untouched. [Local evidence](data/native-ashr-advanced-optimizer-local.json)
binds source, fixtures and artifacts.

## Reproduction and remaining scope

```sh
Rscript scripts/generate_ashr_advanced_optimizer_reference.R "$FRESH_REFERENCE_JSON"
python scripts/generate_interior_mixture_truth.py "$FRESH_REFERENCE_JSON" "$FRESH_TRUTH_JSON"
cargo test -p rsdeseq2 --test advanced_ashr_optimizers --test ashr_optimizers --test location_shrinkage --test mixture_fit --test source_convergence --test mixture_grid
cargo test -p rsdeseq2 --lib polynomial_extrapolation_solves
cargo build --release --example evaluate_r_interface
Rscript scripts/check_location_shrinkage_identity.R "$DRIVER" "$R_LIBRARY" "$FRESH_BATCH_OUT" "$FRESH_REFERENCE_JSON"
Rscript scripts/check_native_ashr_optimizers.R "$DRIVER" "$R_LIBRARY" "$FRESH_OPTIMIZER_OUT" "$PAIRED_FITS_RDS" "$FRESH_REFERENCE_JSON"
```

The mathematical generator requires mpmath (recorded version 1.4.1), a development
oracle dependency only. No runtime dependency is added. Fresh paths preserve
older evidence, and the optional count-fit snapshot avoids repeating count fits.

Remaining work includes estimated modes, other likelihood families, remaining
SQP controls/output subsets, the complete warning/error/authority contract and
final combined precision/speed/memory/sanitizer/platform/source-bundle checks.
No unchanged quantification timing or RSS is retested here, no new performance
claim or counted optimization round is made, and no publication is authorized.
