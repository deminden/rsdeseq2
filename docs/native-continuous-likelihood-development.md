# Native continuous observation likelihoods

This unreleased development adds canonical `ashr::lik_normalmix` and
`ashr::lik_logF` observation likelihoods to the standalone Rust location engine.
The actual R wrapper calls that same implementation. R supplies stored primitive
parameters and object metadata; it does not evaluate likelihood functions,
probabilities, integrals or posterior summaries. This checkpoint does not close
C19/C21/C23 or authorize publication.

## Inputs and numerical policy

Normal-mixture noise accepts shared nonnegative weights with shared noise SDs or
one SD row per complete observation. Log-F accepts positive scalar or per-row
degrees of freedom, including infinite values. Parameter storage is shared;
selecting a posterior row clones an `Arc`, rather than the complete parameter
matrix. Missing-row exclusion during prior fitting retains original row indices.
Matrix storage order is explicit in the standalone API and bridge.

These likelihoods use uniform prior components, including point masses. Original
ashr rejects normal prior components with these likelihoods, and so does the
native engine. The original normal-mixture helper's unnormalized noise-weight
moment arithmetic is preserved; its moments do not reweight noise components
by their conditional interval masses. The original log-F helper supplies no
second moment, so its posterior SD remains missing rather than being invented.
Original log-F missing-observation rejection is retained.

Rust evaluates normal-mixture probabilities/moments and log-F beta/gamma
probabilities, log densities and truncated first moments. Positive gamma tails
use series/continued fractions; infinite integration bounds use native
Gauss–Kronrod quadrature. These are standalone mathematical implementations,
with no R runtime linkage or copied R statistical routine. Mathematical context:
[DLMF incomplete-gamma expansions](https://dlmf.nist.gov/8.9) and the
[Gauss–Kronrod quadrature reference](https://netlib.org/quadpack/dqk15i.f).

Ordinary CDF arithmetic matters here: the original moment helper subtracts
near-one CDFs, amplifying their last-bit differences. Native ordinary beta
complements and small-argument gamma complements retain that operation order
instead of reconstructing every ordinary CDF from a log CDF. A public infinite-df
case initially differed by about 2.46 in posterior mean; its corrected error is
below `5e-10`. Independent hand-selected near-one intervals and full batches
now guard that cancellation. Source integration failures, density overflow and
midpoint fallbacks remain observable rather than being silently replaced with
a different conditional-mean policy. Both-infinite-df probability endpoints
also retain the original half-mass value at zero.

The wrapper recognizes canonical constructor syntax and reads already-stored
primitive captures. It never calls the constructors or their numerical closures.
Modified closures, active bindings and delayed captures do not execute. Other
restricted custom likelihood expressions keep their existing translation path.
See the [original ashr interface](https://stephens999.r-universe.dev/ashr/doc/manual.html)
for its likelihood/prior conventions.

## Independent evidence

The generator `scripts/generate_continuous_likelihood_reference.R` runs only
original ashr 2.2-63 under R 4.6.2. Its committed fixture contains eleven kernel
profiles (99 observation rows) and 21 complete batches: 18 valid batches with
108 rows and three controlled rejections. Coverage includes row parameters,
unnormalized and zero-SD mixtures, finite/infinite bounds, small/large/infinite
degrees of freedom, missingness and cancellation. Existing reference files are
not regenerated or changed.

| Check | Bound or outcome |
| --- | --- |
| Original kernel log CDF/density | `2e-8 + 2e-12 * abs(reference)` |
| Original kernel first/second moments | `2e-6 + 2e-12 * abs(reference)`; original NA/NaN/infinities retained |
| Original complete mean/SD/s-value/CDF/log likelihood | absolute `2e-6` |
| Actual wrapped versus standalone requests | Exact output bits for eleven kernel and 21 complete batches |
| Eight cached public 128-gene profiles | Seven valid; one matching rejection; posterior bound `1e-3`, s-value bound `0.02` |
| Maximum valid public mean difference | `7.714315e-7` |
| Previously accepted replay | Fifteen complete payload sets unchanged against `3b85219` |
| Focused checks | 24 Rust integration tests, two numerical unit tests and 83 R assertions |

Public coverage includes matrix-valued noise SDs, row-specific log-F degrees of
freedom, infinite degrees, positive thresholds and estimated modes. The original
second-infinite-df public profile constructs an invalid likelihood matrix; both
engines reject it, with their separately recorded error messages. Positive-threshold
s-values preserve DESeq2's normal-likelihood query policy even when the fitting
likelihood is overridden. Reference s-value missingness uses the original helper.
All candidate R computations run with statistical entry points blocked.

The [local report](data/native-continuous-likelihood-local.json) binds source,
executables, fixture and validation artifacts. Development is on the local
branch only. No new runtime dependency, quantification speed/RSS claim or counted
optimization round is introduced. Cached count fits are reused; unchanged count
quantification is not repeated for this optional-likelihood change.

## Reproduction and remaining work

```sh
Rscript scripts/generate_continuous_likelihood_reference.R "$FRESH_REFERENCE_JSON"
cargo test -p rsdeseq2 --test continuous_likelihood
cargo build --release -p rsdeseq2 --example evaluate_r_interface
Rscript scripts/check_location_shrinkage_identity.R \
  target/release/examples/evaluate_r_interface "$R_LIBRARY" "$FRESH_IDENTITY_OUT" \
  crates/rsdeseq2/tests/data/continuous_likelihood_r462.json
Rscript scripts/check_native_continuous_likelihoods.R \
  target/release/examples/evaluate_r_interface "$R_LIBRARY" "$FRESH_WORKFLOW_OUT" \
  "$CACHED_PAIRED_FITS_RDS"
```

Poisson/binomial likelihoods, the complete method/warning/error/authority audit
and final combined precision, speed, memory, sanitizer, platform and source-bundle
gates remain open. The published release and historical measurements retain
their original source identities.
