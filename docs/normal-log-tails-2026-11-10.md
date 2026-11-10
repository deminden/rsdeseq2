# Standalone normal log tails

The posterior migration now has standalone normal log-CDF and log-upper-tail
primitives in `rsdeseq2::math`. They support likelihood calculations for the
remaining uniform/truncated-normal components when ordinary probabilities
underflow. The existing Wald probability path and public shrinkage routes are
unchanged.

For ordinary positive arguments, the new function takes the logarithm of the
existing standalone normal tail. Negative arguments use `log1p` to retain small
log probabilities near one. Beyond eight standard deviations, it evaluates the
normal Mills ratio with a 64-level Laplace continued fraction, then combines it
with the normal log density. This is a direct Rust implementation of the
[continued-fraction identity](https://dlmf.nist.gov/7.9.E1). Infinite arguments
have their limiting values, NaN propagates, and no R runtime or callback is used.

The [independent generator](../scripts/generate_normal_log_reference.R) calls
only original R 4.6.2 `stats::pnorm(log.p=TRUE)`. Its
[frozen panel](../crates/rsdeseq2/tests/data/normal_log_r462.json) has **666
arguments** covering ordinary values, both sides of the eight-sigma transition,
probability-underflow boundaries, extreme finite log tails through roughly
`1.4e154`, infinities and NaN. Local regeneration is byte-identical. The
[Rust test](../crates/rsdeseq2/tests/normal_log.rs) requires each finite result to
be within 16 ULPs of the independent reference, including subnormal values.
The maximum observed difference is **3 ULPs** for each function in this panel;
this is a measured panel result, not a claim about every binary64 argument.

The [direct interface check](../scripts/check_normal_log_identity.R) transports
actual input bits and requires every output bit to agree between the R-free
executable and the native R bridge. It blocks R `log`, `log1p`, `exp`, `expm1`
and `pnorm` during bridge calls and verifies its guard with a positive control.
All 666 arguments pass for both log tails. Both R CI lanes run this check.
The conformance driver shares numeric/bit input parsing across the new tail,
prior-fitting and posterior panels; their earlier complete-fit checks still pass.

The [validation record](data/normal_log_validation_2026-11-10.json) binds the
source, generator/fixture and checked artifacts. Local validation passes 1,084
Rust tests, three backend tests, 2,667 R assertions with source check status OK,
47 Python tests with one platform skip, both Clippy/MSRV 1.97.1 graphs and vendor
integrity. All 329 archived R/C/Rust files match the checked source. No existing
parity bound or oracle changed. No unchanged public workflow is retimed and no
new speed, memory or optimization-round claim is made.

Uniform/truncated-normal components, other likelihoods, the complete ashr
prior/grid/options/missingness/warning adapter and apeglm still precede public
adoption. This checkpoint does not close the
[Rust-only migration](rust-only-wrapper-migration.md) or its final combined
authority, performance and artifact gates.

```sh
Rscript scripts/generate_normal_log_reference.R /tmp/normal-log.json
cmp /tmp/normal-log.json crates/rsdeseq2/tests/data/normal_log_r462.json
cargo test --locked -p rsdeseq2 --test normal_log
cargo build --release --locked -p rsdeseq2 --example evaluate_r_interface
Rscript scripts/check_normal_log_identity.R \
  target/release/examples/evaluate_r_interface R_CHECK_LIBRARY /tmp/normal-log-identity
```

Use the pinned native R environment for regeneration and a checked candidate
library for the identity comparison. Keep the output directory fresh.
