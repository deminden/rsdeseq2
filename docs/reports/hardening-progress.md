# Hardening progress

One row per completed phase, measured as in the [0.5.0 baseline](baseline-0.5.0.md).
Times come from the same 2-vCPU Linux x86_64 environment unless a row says
otherwise; structural counts come from
`python3 scripts/validation/hardening_metrics.py`.

| Phase | Commit | Production unwrap + expect + unreachable! | Repeated-block share | `pub fn` | `include!` | Test binaries | Release build s | Test run s | Eleven-contrast s | Peak RSS MiB |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Baseline 0.5.0 | `ede9157` | 119 | 10.8% | 882 | 80 | 76 | 90.2 | 35.5 | 105.7 | 253.6 |
| Phase 1 | `hardening/1.6-features` | 119 | 10.3% | 881 | 80 | 76 | 115.3 | 33.5 | 103.0 (paired: 0.5.0 109.2) | 253.6 |
| Phases 2–4, 7, 8.1–8.2 (in progress) | `b60bd26` | 12 | 9.8% | 889 | 0 | 1 | not re-measured | 26 (one binary; 37 before 7.1) | see 5.1 (paired 1.119× vs Phase 1 at 2 threads) | 253.7 |

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
| Significance flips against DESeq2, gated contrasts | 0 on the GTEx gate (measured from Phase 6.4) | 0 |

## Phase 0 notes

- The golden manifest `crates/rsdeseq2/tools/validation/golden/r-linux-x86_64.x86_64-linux.tsv`
  holds 7,750 field hashes; every field is identical at 1 and 4 threads.
- Workflows that legitimately fail are pinned by their error text: parametric
  trends on the 9-gene and 4-gene inputs, and size-factor estimation on the
  9-gene input (which supplies full-cohort factors for every other workflow).
- Portable now uses a shared manifest. The golden harness checks 8,032
  fields against it on macOS and Linux x86_64. The host R profiles have
  separate manifests.
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


## Phases 2–8 notes (to `b60bd26`)

- Production unwrap/expect/unreachable!: 119 → 12. The remaining 12 are
  documented invariants with function-level allows (2.3). The library denies
  unwrap, expect and undocumented unsafe (2.4); the R bridge now denies
  undocumented unsafe blocks too, and all 208 carry SAFETY comments (8.2).
- `include!` sites: 80 → 0 in the core crate (3.1, one parent per PR; rustdoc
  item list unchanged for each). Mutual module pairs: 11 → 2, both documented
  and checked in CI (3.2).
- Integration-test binaries: 76 → 1 (7.1). Cold test build 66 s → 38 s;
  test executables 7.1 GB → 0.42 GB.
- R entry points without a panic guard: 7 of 34 → 0 (8.1).
- `builder_fit_methods` is still 142: the `run` consolidation (4.2) waits for
  the design note.
- 6.10 (NUMERIC) changed only the toy 4×4 parametric golden cases, from an
  error to the local-trend fallback DESeq2 performs. Golden manifest now has
  8,032 fields.
- Found by new tests: unknown formula calls such as `~ q(age)` recursed until
  the stack overflowed (fuzzing, fixed in 7.3a); DESeq2's parametric-to-local
  fallback and local-trend row rule were missing (edge-case probing, 6.10).

## Mutation testing sample (7.4)

`cargo mutants` on `numeric/validate.rs` and `multiple_testing.rs` (lib and
integration tests), 48 mutants: 45 caught, 3 missed. All 3 missed mutants
were in the BH helpers used by independent filtering:
`rejections_for_filter` returning a constant, `<` → `<=` at the alpha
boundary, and `||` → `&&` in `bh_adjust_f64` validation. Unit tests now kill
all three. The weekly `coverage and mutation` workflow covers `numeric/`,
`multiple_testing`, `normalization` and `independent_filtering`; line
coverage needs llvm-tools, which this environment cannot install, so it is
measured in CI only.

## Development after the bundled checkpoint

The 17 November review lists the improvement targets. The implementation
records through `db04ff0` cover the `run` migration, runtime numeric
profiles, rank and transform parity, fuzzing, reference breadth, and R
vendor generation. Compatibility wrappers remain through the planned
deprecation release. The scaled rank reference generator was restored from
its original commit and reproduces the checked-in fixture byte for byte.

| Area | Remaining work |
| --- | --- |
| Independent references | Run the informational R 4.5 and Bioconductor 3.22 lane in hosted CI |
| Architecture | Review the [candidate 1.0 public surface](../design/run-entry-point.md#candidate-public-surface-for-10) before narrowing exports |
| Validation | Run the five hosted platform panels, including pinned full-table restoration from a shallow checkout |
| Release | Align package versions and verify the release artifacts when preparing 0.6.0 |

IRLS now writes directly into preallocated output matrices. Its covariance
diagnostic excludes skipped genes. The two remaining circular dependencies
were removed while keeping existing `glm::nb` and `independent_filtering`
imports available. The transaction fuzzer covers staging, paths that alias,
incomplete outputs and failure at each rename. The posterior fit separates
offset validation, optimization and Laplace covariance construction. The
in-place QR trial changed 2,236 macOS golden fields. A scratch-buffer QR trial
matched both golden profiles and full GTEx output but achieved only 0.914× the
existing throughput. The branch retains the existing QR solver.

The Parity workflow restores the two large GTEx tables from a pinned Git
commit and verifies their SHA-256 hashes against the compact fixture manifest.
Routine tests use the tracked compact panel.

### Local verification and benchmark

- Rust workspace tests: 1,220 passed, including the explicitly run
  56,200-gene GTEx replay that the default suite skips. Both macOS and portable
  golden checks matched 8,032 fields. The installed R package passed 2,874
  expectations without skips on R 4.6.2 and DESeq2 1.52.0.
- On macOS arm64, one paired four-worker run over eleven full GTEx
  contrasts took 19.42 s versus 20.91 s for the baseline from `09fe127`.
  The separate executables used the portable profile and produced
  byte-identical result TSVs. Throughput increased 7.65% and peak RSS was
  258.23 versus 261.67 MiB. This is a single pass.
- The configured macOS arm64 release panel passed 115 tests. The installed
  R package check above ran on the same numerical implementation before
  the posterior preparation/Laplace extraction.
- The synthetic Cook's replacement/rlog prior probe ran at 20,000 × 100 and
  60,000 × 800 with 1, 2, 4 and 8 workers. At 60,000 × 800, the eight-worker
  run took 0.142 s for replacement and 1.416 s for rlog prior estimation.
  Peak process RSS was 3,218 MiB. All runs returned the same candidate count
  and rlog prior bits for each matrix. These are single measurements.
