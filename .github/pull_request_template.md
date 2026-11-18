## Task

Hardening task number and label: **IDENTICAL** / **NUMERIC** (choose one; never both).

## What changed

<!-- Five lines or fewer. -->

## Proof

- IDENTICAL: `cargo run --release --locked -p rsdeseq2 --example golden -- --check`
  reports every field identical (1 and 4 threads).
- NUMERIC: per-field maximum absolute and log10 differences against the
  previous commit and against the R references, and significance flips at
  padj 0.05 and 0.1 (any flip needs maintainer sign-off).

## Performance (Phase 5 tasks)

Paired eleven-contrast timing (`scripts/benchmarks/paired_gtex_rust_timing.py`),
throughput ratio against the previous commit, and peak RSS.

## Public API

Changes and the `CHANGELOG.md` entry, or "none".

## Review findings that did not hold

<!-- Optional. -->

## Gates

- [ ] `cargo fmt --all -- --check`
- [ ] `cargo clippy --locked --workspace --all-targets -- -D warnings -A missing_docs`
- [ ] `cargo test --locked --workspace`
- [ ] golden harness `--check`
- [ ] `python3 scripts/validation/check_docs.py`
- [ ] `python3 -m unittest discover -s scripts/tests -p 'test_*.py'`
- [ ] `python3 scripts/validation/check_module_cycles.py`
- [ ] `scripts/release/sync_r_vendor.sh --check`
