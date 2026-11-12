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
