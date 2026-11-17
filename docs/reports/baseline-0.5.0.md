# Hardening baseline: 0.5.0

Measured on 2026-11-17 from source commit `ede9157` (the 0.5.0 tree) before any
hardening change. Later phases report against these numbers in the
[hardening progress report](hardening-progress.md).

## Environment

| Item | Value |
| --- | --- |
| Machine | Cloud VM, 2 vCPU Intel Xeon @ 2.10 GHz, 7 GiB RAM |
| OS | Ubuntu 24.04.5, glibc 2.39, Linux x86_64 |
| Toolchain | rustc/cargo 1.97.0 (see note below) |
| Numeric profile | `r-linux-x86_64` (compile-time Linux x86_64 arithmetic) |

The workspace declares `rust-version = "1.97.1"`. That toolchain could not be
downloaded in the measuring environment, so all commands used 1.97.0 with
`--ignore-rust-version`. Rust does not change floating-point semantics between
patch releases, and the golden manifest was produced with the same toolchain
that checks it here; CI regenerates nothing and checks with its own stable
toolchain.

## Tests and build

| Measure | Value |
| --- | --- |
| Tests passed / failed / ignored | 1,159 / 0 / 1 |
| Test executables | 78 (76 integration binaries, library, CLI) plus doctests |
| `cargo build --release --locked --timings` (clean) | 90.2 s wall, 155.6 s user, peak RSS 0.89 GiB |
| — of which the `rsdeseq2` library unit | 44.5 s |
| `cargo test --workspace --no-run` (clean debug profile) | 117.3 s wall |
| `cargo test --workspace` (built) | 35.5 s wall |
| Golden harness (`--check`, release) | 7,750 fields, about 11 s |

Commands were run in this order on an otherwise idle machine:

```sh
cargo build --release --locked --timings
cargo test --locked --workspace --no-run
cargo test --locked --workspace
```

## Eleven-contrast GTEx timing

Rust-only result path, one Rayon worker, three repetitions, with
`scripts/benchmarks/paired_gtex_rust_timing.py` (the
`benchmark_gtex_native_wald.py` command without R). Input hashes matched
`docs/data/gtex_dispersion_reference_sha256.tsv`. Binary SHA-256
`cbe90d2ac4d364d03a2a8f6634bc25b4c52538c98d6714c243100e79499ce8b1`.

| Contrast | Median s | Repetitions s | Peak RSS MiB |
| --- | ---: | --- | ---: |
| condition_blocked | 10.62 | 10.46, 10.80, 10.62 | 253.6 |
| sex_female_vs_male | 7.54 | 7.45, 7.54, 7.54 | 246.6 |
| male_age4_vs_age3 | 6.33 | 6.60, 6.29, 6.33 | 213.5 |
| blocked_random_rep01 | 11.08 | 11.15, 11.08, 10.90 | 253.5 |
| blocked_random_rep02 | 10.80 | 10.73, 10.80, 10.95 | 253.5 |
| blocked_random_rep03 | 12.05 | 11.98, 12.05, 12.18 | 253.2 |
| blocked_random_rep04 | 12.13 | 11.78, 12.30, 12.13 | 253.5 |
| blocked_random_rep05 | 12.03 | 12.03, 12.14, 11.82 | 253.3 |
| unblocked_random_rep01 | 7.68 | 7.68, 7.55, 7.78 | 246.5 |
| unblocked_random_rep02 | 7.78 | 7.69, 7.78, 7.99 | 246.6 |
| unblocked_random_rep03 | 7.67 | 7.67, 7.65, 7.67 | 246.4 |
| **Sum of medians / max** | **105.71** | | **253.6** |

Repetition spread reaches 4% on single contrasts, so only paired passes with
the 2% summed-throughput rule from [development](../development.md#numerical-and-performance-changes)
can establish a speed-up. These absolute times are specific to this VM.

## Precision gates on this platform

`scripts/validation/check_gtex_dispersion_improvement.py` with the 0.5.0
binary passes all 9 global-maximum gates and 113 of 121 per-contrast gates on
Linux x86_64. The eight failures are non-target "no worse than 1.05 × before"
gates whose committed `before` values were not recorded on this platform:

| Contrast | Field | Committed before | 0.5.0 here |
| --- | --- | ---: | ---: |
| condition_blocked | result dispersion | 2.24e-05 | 4.17e-05 |
| condition_blocked | result lfcSE | 1.38e-06 | 1.05e-05 |
| condition_blocked | result stat | 1.19e-05 | 7.91e-05 |
| condition_blocked | result pvalue | 4.35e-07 | 1.63e-05 |
| condition_blocked | stage dispMAP | 1.38e-05 | 4.17e-05 |
| condition_blocked | stage dispersion | 2.24e-05 | 4.17e-05 |
| blocked_random_rep02 | result stat | 1.56e-05 | 2.32e-05 |
| blocked_random_rep02 | result pvalue | 6.82e-06 | 1.37e-05 |

NUMERIC tasks must therefore compare against the previous commit measured on
the same platform as well as against the R references, and must not treat
these eight pre-existing gaps as regressions they caused. The checker also needs
the historical `rust_fit_diagnostics.tsv` files for its stage-route summary; the
comparison tables above are written before that step.

## Structural metrics

From `python3 scripts/validation/hardening_metrics.py` (static counts; see the
script for exact rules):

| Metric | 0.5.0 |
| --- | ---: |
| Production `unwrap` / `expect` / `unreachable!` | 79 / 18 / 22 (119) |
| `include!` sites | 80 |
| `pub fn` | 882 |
| `pub fn fit*` builder methods | 142 |
| Mutual top-level module pairs | 11 |
| Integration-test binaries | 76 |
| Platform `cfg` sites (`linux` + `x86_64`) | 31 |
| R entry points / not routed through `run_backend` | 34 / 7 |
| `unsafe` tokens in the R bridge | 287 |
| Repeated-block share of production lines (10-line windows) | 10.8% |

The repeated-block share uses a different detector from the review's 5.8%
exact-clone figure; compare it only with itself over time.
