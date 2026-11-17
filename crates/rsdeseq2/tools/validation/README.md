# Internal numerical validation tools

These Rust programs compare fixtures or implement the standalone side of the R
interface identity checks. They are not user-facing analysis examples.

Cargo target names are explicit in the crate manifest and remain stable:

```sh
cargo build --release -p rsdeseq2 --example evaluate_r_interface
```

The existing Cargo example target type keeps these programs out of the installed
CLI set. Its executable remains at `target/release/examples/evaluate_r_interface`.
Shared protocol modules are used by integration tests through explicit paths.

## Golden output harness

`golden` hashes every output field of a fixed workflow panel so structural
refactoring can prove it changed no bits. It runs the 4 × 4 toy reference,
`gtex_auto_dispersion`, `pasilla_holdout` and the compact GTEx panel through
ratio/poscounts size factors; Wald and LRT fits with parametric, local and mean
trends; observation weights; normalization factors; Cook's replacement; VST and
rlog; and normal, ashr and apeglm shrinkage, each at 1 and 4 worker threads.

```sh
cargo run --release --locked -p rsdeseq2 --example golden -- --check
cargo run --release --locked -p rsdeseq2 --example golden -- --check --case pasilla_holdout
GOLDEN_VERBOSE=1 cargo run --release --locked -p rsdeseq2 --example golden -- --check
```

Each manifest line is `input, workflow, threads, field, sha256`. Values are
hashed from raw `f64` bit patterns, so signed zeros and NaN payloads (R's
`NA_real_`) count. Fields hashed from `Debug` text carry a `~debug` suffix; this
is exact for all non-NaN values. A workflow that fails records its error text.

Manifests live in `golden/<profile>.<arch>-<os>.tsv`. The numeric profile is
`r-linux-x86_64` on Linux x86_64 and `portable` elsewhere, matching the current
compile-time platform arithmetic. Only an explicitly NUMERIC change may rewrite a
manifest (`--write`), and its pull request must show the per-field differences.
A missing manifest fails `--check`; generate it on that target first.
