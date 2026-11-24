# Internal numerical validation tools

These Rust programs compare fixtures or implement the standalone side of the R
interface identity checks. They are not user-facing analysis examples.

Cargo target names are explicit in the crate manifest and remain stable:

```sh
cargo build --release -p rsdeseq2 --example evaluate_r_interface
```

Run the standalone probe at `target/release/examples/evaluate_r_interface`.
It selects the same host R arithmetic as the wrapper by default. Set
`RSDESEQ2_NUMERIC_PROFILE` for both processes to compare a different profile.
Integration tests also use the shared protocol modules.

## Golden output harness

Run this harness from the repository checkout after downloading the Pasilla
validation fixture.

`golden` hashes every output field of a fixed workflow panel so structural
refactoring can prove it changed no bits. It runs the 4 × 4 toy reference,
`gtex_auto_dispersion`, `pasilla_holdout` and the compact GTEx panel through
ratio/poscounts size factors; Wald and LRT fits with parametric, local and mean
trends; observation weights; normalization factors; Cook's replacement; VST and
rlog; and normal, ashr and apeglm shrinkage, each at 1 and 4 worker threads.

```sh
python3 scripts/references/fetch_pasilla_validation.py
cargo run --release --locked -p rsdeseq2 --example golden -- --check
cargo run --release --locked -p rsdeseq2 --example golden -- --check --case pasilla_holdout
GOLDEN_VERBOSE=1 cargo run --release --locked -p rsdeseq2 --example golden -- --check
```

Each manifest line is `input, workflow, threads, field, sha256`. Values are
hashed from raw `f64` bit patterns, so signed zeros and NaN payloads (R's
`NA_real_`) count. Fields hashed from `Debug` text carry a `~debug` suffix; this
is exact for all non-NaN values. A workflow that fails records its error text.

The portable profile uses `golden/portable.tsv` on every target. Use
`RSDESEQ2_NUMERIC_PROFILE=portable` to select it explicitly. Host R profiles
use `golden/<profile>.<arch>-<os>.tsv`. A missing manifest fails `--check`.
Only an explicitly numerical change can rewrite an existing manifest
(`--write`). Record the per-field differences and their causes.
