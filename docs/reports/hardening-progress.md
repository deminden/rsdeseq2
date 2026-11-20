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

| Area | Status |
| --- | --- |
| Architecture | The [public API review](../design/run-entry-point.md#reviewed-public-surface-for-10) identifies the workflow, fitted-result and standalone-stage entry points. The 0.7 migration handles alias removal. |
| Release | Versions and lockfiles agree at 0.6.0. All four Rust archives built from their packaged contents and regenerated byte for byte. The R source archive installed in both reference environments. |
| Local validation | macOS arm64 release panel, both golden profiles, full 56,200-gene GTEx replay and shallow-clone fixture restoration passed. |
| Earlier Bioconductor | R 4.5.3 with DESeq2 1.50.2 passed all 2,874 expectations with zero skips after the cached binding fix in `583c950`. R 4.6.2 with DESeq2 1.52.0 passed the same checks. |
| Other targets | See the [local platform results](#local-platform-validation). The hosted 0.6.0 matrix is pending. |

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


## Local platform validation

The dispersion refactor at `3626502` borrows the selected counts, offsets and
weights during mean fitting. Weighted dispersion updates keep the R
compact-position weight rule. The optimized CLI and its `1522aac` baseline
produced byte-identical result and fit-diagnostic TSVs for eleven full GTEx
contrasts under the portable profile.

The R wrapper selects the system math library on supported hosts. Linux arm64
uses binary128 scalar reductions and the OpenBLAS ARMV8 dot-product order.
Each release target runs the host R numerical panel and the portable golden
manifest. Linux arm64 also runs its exact formula panel and its host R
manifest.

| Target | Local results |
| --- | --- |
| macOS arm64 | 114 release-panel tests passed. The workspace suite, both 8,032-field golden checks and the separate 56,200-gene GTEx replay passed. |
| Linux arm64 | 1,221 workspace tests, 129 release-panel tests, both 8,032-field golden checks and both full GTEx replay tests passed in Ubuntu 24.04. |
| Linux x86_64 | 114 release-panel tests, both golden manifests and both full GTEx replay tests passed in Ubuntu 24.04 under QEMU 10.0.13. |
| macOS x86_64 | The locked workspace and all targets cross-compiled on macOS arm64. The runtime panel is pending. |
| Windows x86_64 | The locked workspace libraries and binaries cross-checked on macOS arm64. The test build requires a Windows SDK and the runtime panel is pending. |

Original Linux arm64 R 4.6.2 with ICU 74.2 and OpenBLAS 0.3.26 independently
regenerated all seven formula panels exactly. The installed R package and
standalone Rust agreed with those references on 191 expression cases, 54
native formulas, and 60 standard and expanded factor cases. The independent
uniform posterior checks passed their original bounds. The R formula-related
testthat subset ran 849 expectations with zero failures. Its ordered-factor
fitting test requires DESeq2, which was absent from that Linux arm64 R library.

The macOS host R profile passed all 121 per-contrast and nine global GTEx
precision gates with the frozen references. A separate portable run passed
107 of 121 per-contrast gates and all nine global gates. The portable
significance check reported zero flips at adjusted p-value thresholds 0.05
and 0.1, and all top-100, top-500 and top-2,000 sets matched. The eight genes
at the failed maxima had differences in the fitted dispersion trend. Three
also had different MAP iteration counts. The beta iteration counts matched
for all eight.

The initial paired trials had variable wall times under concurrent load.
`41d7509` restored contiguous validation and likelihood evaluation for
full IRLS fits while retaining borrowed rows for dispersion refits. In one
paired four-worker pass against `1522aac`, the final build took 31.26 s
versus 33.00 s across eleven contrasts. That is 1.056× throughput, with
byte-identical result TSVs. Maximum peak RSS was 264.0 MiB versus 263.5
MiB. A separate `time -l` check on the blocked condition contrast measured
264.34 MiB for the candidate and 262.12 MiB for the baseline. Both
8,032-field golden checks and the full 56,200-gene replay passed after this
change. The final host R results and diagnostics matched the previous build
byte for byte across all eleven contrasts.

The macOS R 4.5.3 and 4.6.2 packages each passed 2,874 testthat expectations
with zero failures and zero skips, and the SIGINT check confirmed that the
interrupted R session could fit again. The interface probe matched 33
complete formula fits and their later queries after selecting the same
host profile as the R wrapper. Linux arm64 passed the C ASan/UBSan and Rust
ASan check, including all four replacement and four formula unwind cases.
Locked minimum-Rust tests and the dependency audit passed locally.

Four Rust archives built from packaged contents and regenerated byte for byte
at `41d7509`. Its R source archive passed all 2,874 expectations on each of
R 4.5.3 and R 4.6.2 with zero skips. Both package checks recorded two
vendoring notes. R 4.6.2 also recorded the compiled-symbol warning.
The Python validation suite passed all 71 tests with the installed R package
and full GTEx tables available.


### Linux arm64 profile verification

The `dc2e0f9` source passed 1,221 workspace tests on macOS arm64 as well as
Linux arm64. Linux x86_64 passed its 114 release-panel tests under QEMU
10.0.13. All 8,032 portable hashes matched on those three targets, and each
host profile matched all 8,032 fields in its own manifest.

On macOS arm64, the rebuilt binary passed all 121 per-contrast and nine
global GTEx precision gates under `r-host`. The rebuilt R source archive
installed on R 4.5.3 with DESeq2 1.50.2 and on R 4.6.2 with DESeq2 1.52.0.
Each installation passed 2,874 expectations with zero skips. All four Rust
archives built offline from their packaged contents. The Python suite
passed 71 tests with the rebuilt R package and full GTEx tables available.
