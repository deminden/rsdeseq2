# Hardening progress

One row per completed phase, measured as in the [0.5.0 baseline](baseline-0.5.0.md).
Times come from the same 2-vCPU Linux x86_64 environment unless a row says
otherwise; structural counts come from
`python3 scripts/validation/hardening_metrics.py`.

| Phase | Commit | Production unwrap + expect + unreachable! | Repeated-block share | `pub fn` | `include!` | Test binaries | Release build s | Test run s | Eleven-contrast s | Peak RSS MiB |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Baseline 0.5.0 | `ede9157` | 119 | 10.8% | 882 | 80 | 76 | 90.2 | 35.5 | 105.7 | 253.6 |
| Phase 1 | `hardening/1.6-features` | 119 | 10.3% | 881 | 80 | 76 | 115.3 | 33.5 | 103.0 (paired: 0.5.0 109.2) | 253.6 |

## Targets

| Metric | 0.5.0 | Target |
| --- | ---: | --- |
| Production unwrap + expect + unreachable! | 119 | 0 outside documented invariants |
| `include!` sites | 80 | 0 |
| `pub fn fit*` on the builder | 142 | 1 `run` plus deprecated wrappers, then 1 |
| Mutual top-level module pairs | 11 | 0 |
| Integration-test binaries | 76 | 2 |
| Platform `cfg` sites that change numerics | 31 | 0 (runtime profile) |
| R entry points not routed through `run_backend` | 7 of 34 | 0 |
| Significance flips against DESeq2, gated contrasts | not yet measured (Phase 6.4) | 0 |

## Phase 0 notes

- The golden manifest `crates/rsdeseq2/tools/validation/golden/r-linux-x86_64.x86_64-linux.tsv`
  holds 7,750 field hashes; every field is identical at 1 and 4 threads.
- Workflows that legitimately fail are pinned by their error text: parametric
  trends on the 9-gene and 4-gene inputs, and size-factor estimation on the
  9-gene input (which supplies full-cohort factors for every other workflow).
- Manifests for the `portable` profile (macOS, Windows, Linux aarch64) still
  need to be generated on those hosts with `--write`; until then `--check`
  fails there by design.
- Hardening lint warnings at baseline (workspace `[lints]`, warn only): 461 in
  the library and CLI, 997 across all targets. See the
  [lint baseline](baseline-0.5.0.md#lint-baseline).

## Phase 1 notes

- Every task was IDENTICAL; the golden manifest is unchanged since Phase 0 and
  passes at 1 and 4 threads after each task.
- Eleven-contrast outputs are byte-identical to 0.5.0. One paired pass (2
  repetitions, 1 thread) gave a 1.061 throughput ratio over 0.5.0. Unpaired
  sums drift by about 6% between sessions on this VM, so only the paired ratio
  is meaningful.
- The release build is slower (115 s against 90 s) because of thin LTO with one
  codegen unit (task 1.5); the binary is 16% smaller.
- Panic-site counts are unchanged by design: Phase 1 is deduplication. Phase 2
  removes them.
- The R bridge now builds without the CLI dependencies (task 1.6).

