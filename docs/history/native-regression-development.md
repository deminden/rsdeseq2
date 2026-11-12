# Local general-regression shrinkage development

This unreleased checkpoint extends the shared Rust posterior engine with apeglm
`general` and analytic `nbinomR` methods. It follows the accepted
[default negative-binomial checkpoint](native-apeglm-development.md).
BFGS and bounded L-BFGS-B, numerical derivatives, Laplace uncertainty,
empirical prior adaptation, scale conversion and s-values run in Rust. The R
wrapper imports inputs and constructs result objects. No R likelihood function,
optimizer or statistical routine is invoked by candidate calculations.

## Imported likelihoods

General likelihood functions have five arguments in the order
`(y, x, beta, param, offset)` and return one log-density per observation.
They may contain one pure expression using arithmetic, elementary transforms,
comparisons, `ifelse`, `pmin`, `pmax`, `%*%`, `as.numeric`, `as.double`, `lgamma`,
`dnorm`, `dnbinom`, `dbinom` and `dt`. Density arguments support positional,
exact and partial-name binding; the density log flag is scalar. Supported names may be explicitly qualified by
their original `base` or `stats` namespace. Matrix products and vector recycling
are evaluated by the shared Rust engine.

Captured values must be stored primitive numeric/logical scalars. The importer
rejects arbitrary calls, assignments, reductions, captured arrays, masked
functions, active bindings and unevaluated promises. It validates even unused
branches without calling them; supported `ifelse` branches are evaluated lazily.
The canonical original `apeglm::logLikNB` is recognized by its namespace and
syntax, then represented as a native negative-binomial likelihood. Formula
custom functions remain disallowed.

General fitting supports central/noncentral Student-t priors, Gaussian nuisance
priors and the original empty/duplicate nuisance-index behavior. `nbinomR`
retains its analytic gradient, zero-centered Cauchy prior, initialization and
log-link-dependent zero-row policy. Missing `param` is rejected, including for
custom expressions that do not use it, matching original apeglm's preprocessing.
A numeric parameter containing missing values may be passed to an expression
that does not evaluate it.

## Evidence and precision

The original reference uses `deseq2-r462-native`: R 4.6.2, DESeq2 1.52.0 and
apeglm 1.34.0. Twenty-one posterior profiles cover Gaussian, algebraic Gaussian,
Student-t, Poisson and custom/built-in negative-binomial likelihoods, weights,
offsets, bounds, alternate priors, intercept adaptation, nuisance indices,
numeric matrix conversion, unused missing parameters and invalid inputs.
Complete coefficient modes/SDs and sign/threshold summaries retain the existing
public bounds: `0.01` log2 units for mode/SD, `0.02` for sign quantities and
`0.05` for threshold quantities. Prior scale/variance checks remain tight.

Twenty density-expression panels cover large sizes/counts, binomial and Student-t
distributions, missing values, infinities and lazy branches. Nine observations
exceed a new leaf gate of `2e-10 + 2e-12 * abs(reference)` when compared with
original R. Independent 90-decimal log-gamma PMFs and integration of the defining
Normal/chi-square mixture show that these are reference cancellation/tail errors.
The same tight leaf gate passes against that independent truth. Original R
observations remain stored separately, and no existing posterior bound or frozen
reference was relaxed. The largest discrepancies are about `1.96e-9` log units
for negative-binomial size `1e8`, and `0.01167` log units at `x=-8`, `df=100`,
`ncp=0.7`. R emitted two tail-precision warnings during reference generation.
R documents cancellation and accuracy limitations in its
[noncentral Student-t implementation](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/TDist.html).

The actual R/native bridge and separate R-free executable agree bit for bit on
selected posterior columns, s-values, adapted priors, nuisance indices and
convergence status for all valid profiles; both reject the two invalid profiles.
The candidate runs with R statistical routines blocked, including matrix products,
densities, gamma functions, optimization and original apeglm routines. Thirty-five
import-policy assertions pass. Public general, `nbinomR` and translated custom
negative-binomial workflows pass against original DESeq2 on 1,200 genes. The
[machine-readable report](../data/native-regression-local.json) records errors and
source/artifact identities.

The default apeglm regression check still passes the thirteen frozen profiles,
a normalization-factor profile and the public default/sign/threshold workflow.
No quantification timing or memory claim is added for unchanged fitting paths.

The subsequent [beta-binomial checkpoint](native-beta-binomial-development.md)
adds the specialized logit family and corrects the common BFGS stopping policy.
All original general-regression references remain frozen and pass; the updated
1,200-gene general/`nbinomR` errors are recorded in that checkpoint. The original
report above retains its own source scope.

## Reproduction and remaining scope

```sh
Rscript scripts/references/generate_regression_posterior_reference.R "$FRESH_R_REFERENCE"
python3 scripts/references/generate_regression_density_truth.py "$FRESH_R_REFERENCE" "$FRESH_TRUTH"
cargo test -p rsdeseq2 --test regression_posterior --test heavy_tail
cargo build --release --example evaluate_r_interface
Rscript scripts/validation/check_native_regression_workflow.R "$DRIVER" "$R_LIBRARY" "$FRESH_OUTPUT" "$PAIRED_FITS_RDS"
```

The optional truth generator requires mpmath; it is a development dependency,
not a runtime dependency. It imports no candidate implementation. Fresh output
paths preserve earlier evidence. The checker reuses the accepted Rust count fit
and fits only the newly changed posterior methods against independent references.

The later [random-start checkpoint](native-random-start-development.md) adds the
two MAP-only random-start methods using portable native RNG state. Grid intervals, parameter uncertainty,
contrasts, other optimizer policies, and the complete warning/control surface
remain unfinished. Advanced ashr work, C21's authority audit and C23's final
combined precision/speed/memory/sanitizer/platform/source-bundle checks also
remain open. This checkpoint is local, does not change package versions and
must not be treated as release completion or publication authorization.
