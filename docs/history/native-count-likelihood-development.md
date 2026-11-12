# Native Poisson and binomial observation likelihoods

This unreleased checkpoint adds Poisson identity/log and binomial identity/logit
likelihoods to the shared Rust engine. Canonical R likelihood objects supply
stored primitive counts, scales, trials and link names. The wrapper evaluates
neither their closures nor their statistical helpers. It does not complete the
full C19/C21/C23 authority and release audit, and publication remains unauthorized.

## Shared numerical and input policy

`BuiltinLocationLikelihood::poisson` and `::binomial` accept nonnegative finite
observations and scalar or complete-row parameter vectors. Poisson scales must
be finite and positive; binomial observations cannot exceed aligned trial counts.
Explicit `PoissonLink` and `BinomialLink` enums select the four supported links.
Parameter vectors remain in shared storage; row selection does not clone them.

Rust evaluates the gamma/beta log probabilities, densities and truncated first
and second moments. The original zero-mass fallback selects the denser endpoint,
choosing the right endpoint on a tie. Tiny adjusted shapes for log/logit links
are computed and checked in Rust. Log/logit fits retain original moments on the
intensity/probability scale, rather than returning moments of log intensity or
log odds. The [original ashr manual](https://stephens999.r-universe.dev/ashr/doc/manual.html)
describes the likelihood conventions; numerical comparisons here use installed
ashr 2.2-63, not the newer online manual's package version.

Uniform-prior convolution likelihoods select their CDF subtraction direction
across the complete component/observation matrix. An observation-wise absolute
difference would change valid/undefined rows. This global decision is shared by
prior fitting, raw/pruned summaries, mode profiling, CDF queries and the direct
fixed-mixture posterior API. Existing continuous paths retain their arithmetic
and do not allocate an additional row-index vector for this feature.

The native adapter also implements count-specific defaults: omitted modes are
estimated, omitted null weights become one, identity-link support constrains
new priors, and automatic grids use the source count-derived observations/SEs.
Poisson mode intervals use scaled counts (or their log-adjusted values);
binomial identity uses its constrained range, and logit uses adjusted log odds.
The first entry of a numeric mode pair expands Poisson identity bounds; other
count links retain their original range policy. Explicit native search intervals
remain explicit. Supplied fixed priors retain their original endpoints.

Original normal-prior rejection and count missing-input rejection are retained.
The installed original argument check accepts mixed dummy observations when at
least one estimate is zero and at least one SE is one; native requests preserve
that behavior. It is not silently tightened to an all-zero/all-one check. Two
independent mixed-observation profiles guard the resulting global-orientation
and NaN policies. Complete warning/error ordering and general authority remain
part of the final audit.

The R importer verifies canonical function syntax against the installed
constructor's two branches and reads only forced primitive captures. It checks
stored metadata alignment, and Rust verifies the adjusted log/logit captures.
Modified closures, active bindings, delayed captures, inconsistent metadata and
tampered adjusted counts do not execute as R calculations. No R statistical
routine was copied into the engine, and no R numerical runtime is linked by the
standalone executable. There are no added runtime dependencies.

## Independent validation

`generate_count_likelihood_reference.R` creates twelve kernel profiles (108 rows)
and 24 complete profiles: twenty valid profiles (120 rows) and four rejections.
They cover scalar/row parameters, zero and large counts, both links, finite and
infinite bounds, fixed/fitted/FDR priors, missingness, and batch orientation.
`generate_count_mode_reference.R` independently records twelve original mode
searches, including explicit/automatic grids and expanded seeds. Their bounds
match bit for bit and evaluation counts match the original search trace.

| Evidence | Acceptance |
| --- | --- |
| Original kernel log CDF/density | `2e-8 + 2e-12 * abs(reference)` |
| Original kernel moments | `2e-6 + 2e-12 * abs(reference)` |
| Complete posterior mean/SD/s-value/CDF/log likelihood | absolute `2e-6`, original missingness retained |
| Original mode and mode-fit numerical fields | absolute `2e-6`; exact bounds and evaluation count |
| R-wrapped versus standalone requests | Exact bits for 12 kernel, 24 complete and 12 mode batches |
| Sixteen public cached 128-gene profiles | Fixed/default modes, automatic grids and default FDR across all four links; unchanged `1e-3` posterior/`0.02` s-value bounds |
| Independent truths | Exponential gamma and uniform beta kernels/moments; direct posterior batch orientation |
| Earlier accepted outputs | Seventeen replay sets unchanged against `543a57e` |

The shared test checker now serves continuous and count fixtures without
repeating the same numerical assertions. Existing committed oracles are
unchanged. Actual candidate R calls run with statistical functions and scalar
optimizers blocked. The [local report](../data/native-count-likelihood-local.json)
binds the source, reference files, executables and validation artifacts.

```sh
Rscript scripts/references/generate_count_likelihood_reference.R "$FRESH_REFERENCE_JSON"
Rscript scripts/references/generate_count_mode_reference.R "$FRESH_MODE_REFERENCE_JSON"
cargo test -p rsdeseq2 --test count_likelihood --test continuous_likelihood
cargo build --release -p rsdeseq2 --example evaluate_r_interface
Rscript scripts/validation/check_location_shrinkage_identity.R \
  target/release/examples/evaluate_r_interface "$R_LIBRARY" "$FRESH_IDENTITY_OUT" \
  crates/rsdeseq2/tests/data/count_likelihood_r462.json
Rscript scripts/validation/check_native_count_likelihoods.R \
  target/release/examples/evaluate_r_interface "$R_LIBRARY" "$FRESH_WORKFLOW_OUT" \
  "$CACHED_PAIRED_FITS_RDS"
```

No unchanged quantification timing/RSS run or optimization round is claimed.
Cached count fits are reused. Remaining work is the complete option, parameter,
method/warning/error and authoritative-routing audit, followed by combined
precision, speed, memory, sanitizer, platform and source-bundle checks.
