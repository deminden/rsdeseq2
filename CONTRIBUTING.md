# Contributing

Statistical computation lives in the Rust crate `crates/rsdeseq2`; the R
package and CLI translate arguments and assemble objects. See
[development](docs/development.md) for repository layout, validation commands
and performance rules, and [testing and parity](docs/testing-parity.md) for
the R and GTEx checks.

## Two kinds of change

Every pull request is labelled exactly one of:

- **IDENTICAL**: numerical output must not change by a single bit. Prove it
  with the golden harness:
  `cargo run --release --locked -p rsdeseq2 --example golden -- --check`.
- **NUMERIC**: output may move within stated bounds. Report per-field maximum
  absolute and log10 differences against the previous commit and the R
  references, and significance flips at padj 0.05 and 0.1. Never relax a
  tolerance or rewrite an existing R reference to make a test pass.

Do not mix the two in one pull request, and keep file moves and logic changes
in separate commits. Aim for under about 600 changed lines unless the change
is mechanical.

## Gates

```sh
cargo fmt --all -- --check
cargo clippy --locked --workspace --all-targets -- -D warnings -A missing_docs
cargo test --locked --workspace
cargo run --release --locked -p rsdeseq2 --example golden -- --check
python3 scripts/validation/check_docs.py
python3 -m unittest discover -s scripts/tests -p 'test_*.py'
python3 scripts/validation/check_module_cycles.py
scripts/release/sync_r_vendor.sh --check
```

After changing the core, run `scripts/release/sync_r_vendor.sh` and commit the
regenerated vendor tree separately.

## Code rules

- No `unwrap`/`expect` in library code (denied outside tests); return
  `DeseqError`, using `errors::invariant` for states earlier validation rules
  out.
- No `unsafe` outside the x87 kernels and the R bridge; every unsafe block has
  a `SAFETY:` comment.
- New dependencies need a stated reason in the pull request.
- Record public API changes in [CHANGELOG.md](CHANGELOG.md); removals need a
  deprecation release first.
- Performance claims need a paired eleven-contrast run and peak RSS
  (`scripts/benchmarks/paired_gtex_rust_timing.py`).
