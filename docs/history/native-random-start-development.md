# Local native random-start shrinkage development

This unreleased checkpoint adds MAP-only `nbinomC*` and `betabinC*` using shared
Rust random-number generation. It follows the
[beta-binomial checkpoint](native-beta-binomial-development.md). R imports stored
state and presents results; it never generates a start vector or calls an R
probability, random-generation or fitting routine.

## State and fitting contract

`PosteriorRng` supports Mersenne-Twister and L'Ecuyer-CMRG with Inversion normals.
`from_seed` initializes the R default generator; `from_seed_with_kind` initializes
either portable family entirely in Rust. `from_r_seed` imports a serialized R
state and `r_seed` exports it. Both sampling-kind headers are preserved, including
the older rounding header, although these posterior methods do not sample indices.
Imported parallel-stream states use the same generator equations.

The implementations use the published MT19937 and combined-recursive equations,
plus documented R initialization/state conventions. Normal inversion starts with
the existing statrs quantile and refines it against the shared compensated normal
tail. No GPL generator or qnorm implementation is copied or linked, and no new
runtime dependency is added. See the original
[Mersenne Twister paper](https://doi.org/10.1145/272991.272995),
[L'Ecuyer paper](https://doi.org/10.1287/opre.47.1.159) and
[R state documentation](https://stat.ethz.ch/R-manual/R-devel/library/base/html/Random.html).

Each fit draws one `N(0, 0.5^2)` coordinate per design coefficient, then reuses
that single vector for every gene. An entirely zero negative-binomial batch
still consumes its batch vector. Random starts remain MAP-only and use the
existing family-specific optimization policy. Prior adaptation, offsets,
normalization and result scaling remain in the common engine.

The standalone APIs `fit_negative_binomial_posterior_with_rng` and
`fit_regression_posterior_with_rng` borrow mutable native state. The latter uses
`RegressionMethod::BetaBinomialCStar`. After a fitting error following a draw,
the caller still owns the advanced state. Early validation failures consume no
draws. The R wrapper preserves this lifecycle by writing back returned state
before reporting a fitting error. An unseeded request is initialized in Rust.
Use an explicit seed or imported state for reproducibility.

Every signed 32-bit state-word pattern is preserved, including R's integer-NA
representation of an unsigned word. Such words are data, not missing observations.
The native bridge restores their bits on input and returns the same bits on
output. Active/delayed state bindings, nonprimitive storage, unwritable bindings,
malformed states and unsupported RNG modes fail explicitly before generation.
Other uniform/normal modes, including user callbacks and Box-Muller, are not
implemented. R's serialized uniform state does not contain Box-Muller's cached
normal value, so importing it alone would not establish equivalent state.

## Independent evidence

Original references use R 4.6.2 and apeglm 1.34.0 in `deseq2-r462-native`.
Thirty generator profiles cover initial seeding, both uniform families, warmed
states crossing MT boundaries, legacy sampling headers, imported parallel
streams, zero-word endpoint handling, integer-NA words, and CMRG seeding-range
rejection boundaries. All 1,920 uniforms and serialized transitions match R
exactly. Normal draws pass a strict 16-ULP component gate; the observed maximum
is recorded in the machine-readable report. This comparison to original R is
separate from the exact native/wrapped identity contract.

The posterior panel has 85 negative-binomial and 79 beta-binomial profiles:
1,270 fitted/input rows across three seeds and both portable families, plus two
explicit integer-NA state cases. It reuses frozen original-reference inputs,
preserves the existing posterior bounds and retains the original numerical
values. Repeated state vectors are pooled in the fixture; expansion changes no
input bits or draw ordering. The first 162 profiles remain unchanged when the
two bridge edge cases are appended.

Actual R/native and separate R-free processes compare selected mode/SD columns,
s-values, adapted priors, nuisance indices, convergence, warnings and final RNG
state bit for bit. Candidate calls run with R random generation, seeding,
probability, matrix and optimizer routines blocked. Public checks use 256 genes
and both methods/families, compare directly with original apeglm and verify
`.Random.seed` exactly. Additional checks cover post-draw errors, early rejection
and native initialization of an unseeded session. Import tests verify that state
inspection does not execute active or delayed bindings.

The [local report](../data/native-random-start-local.json) binds source/artifact
identities, exact state evidence and public precision errors. Original optimizer
warning counts can differ near floating-point plateaus; actual native warnings
remain exposed and identical across the native/wrapped interfaces. These
observations remain part of the final authority audit.

## Reproduction and remaining scope

```sh
Rscript scripts/references/generate_posterior_rng_reference.R "$FRESH_RNG_REFERENCE"
Rscript scripts/references/generate_random_start_posterior_reference.R "$FRESH_POSTERIOR_REFERENCE"
cargo test -p rsdeseq2 --test posterior_rng --test regression_posterior --test heavy_tail
cargo build --release --example evaluate_r_interface
Rscript scripts/validation/check_native_regression_workflow.R "$DRIVER" "$R_LIBRARY" "$FRESH_OUTPUT" "$PAIRED_FITS_RDS" "$FRESH_POSTERIOR_REFERENCE"
```

The generators load only original R/apeglm. Fresh paths protect earlier evidence.
The checker reuses accepted count fits and saved original posterior results when
provided. Deterministic methods remain covered by their unchanged references.

The later [posterior analysis checkpoint](native-posterior-analysis-development.md)
adds grid intervals, scalar parameter uncertainty and posterior contrasts, while
preserving every accepted random draw/state and posterior output. Other optimizer/
RNG modes, advanced ashr methods/options, the final warning/authority audit and combined
precision/speed/memory/sanitizer/platform/source-bundle checks remain unfinished.
No unchanged quantification timing or RSS is remeasured here and no optimization
round is claimed. All work stays local; this is not a release-completion claim.
