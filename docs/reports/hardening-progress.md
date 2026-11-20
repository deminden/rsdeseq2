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
`multiple_testing`, `normalization` and `independent_filtering`. At that
checkpoint, line coverage was measured in CI.

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
| macOS x86_64 | The locked workspace cross-compiled, and the unit and integration test executables linked on macOS arm64. The runtime panel is pending. |
| Windows x86_64 | The locked workspace libraries and binaries cross-checked on macOS arm64. The test build requires a Windows SDK and the runtime panel is pending. |

Original Linux arm64 R 4.6.2 with ICU 74.2 and OpenBLAS 0.3.26 independently
regenerated all seven formula panels exactly. The installed R package and
standalone Rust agreed with those references on 191 expression cases, 54
native formulas, and 60 standard and expanded factor cases. The independent
uniform posterior checks passed their original bounds. After restoring all
82 pinned package versions, the Linux arm64 R package passed the full suite
of 2,874 expectations with zero failures and zero skips.

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


The local instrumented suite measured 86.66% line coverage for the rsdeseq2
crate with nightly Rust and cargo-llvm-cov 0.9.1. All nine matrix and checked
arithmetic tests selected by CI also passed under Miri with the portable
profile.

The integer `sum()` follow-up at `2ef9ff1` passed 1,222 workspace tests on
macOS arm64 and Linux arm64. The formula panel also passed on emulated
Linux x86_64. Six integer-promotion and missing-value cases matched the
original R result types and bits through both arm64 R wrappers. Linux
arm64 passed 130 release-panel tests. macOS arm64 and Linux arm64 each
passed the two 56,200-gene GTEx replay tests and matched their host and
portable golden hashes.

The archives built from `fa2b7b1` include that follow-up. All four Rust
archives built offline and reproduced byte for byte. The R archive passed
2,874 expectations with zero failures and zero skips on macOS R 4.5.3,
macOS R 4.6.2 and Linux arm64 R 4.6.2. The R 4.6.2 checks recorded the
compiled-symbol warning and two vendoring notes. The R 4.5.3 check
recorded both vendoring notes.

The rebuilt macOS CLI passed all 130 GTEx precision gates under `r-host`.
The portable run passed 107 per-contrast and all nine global gates.


### Portable dispersion precision

The portable likelihood now evaluates a centered log-gamma difference.
The likelihood and its derivatives use `log1p` for the mean-dispersion
term. An 80-digit decimal replay covered 200 individual likelihood values. The
median error decreased from 12.5 to 1 binary64 ulp. For count 1, mean
0.01 and log dispersion −25, the result moved from
`-0.009992687141410671` to the rounded reference
`-0.010000000000138185`.

On the eleven frozen GTEx contrasts, the portable profile passed 116 of
121 per-contrast gates and all nine global gates. The five remaining
comparisons cover the blocked condition p-value and fitted dispersion,
the first blocked random split's fitted dispersion, and the third and
fourth blocked random splits' log2 fold changes. The significance check
found zero flips and full top-set overlap. The host R profile passed all
130 gates, with the result and diagnostic TSVs matching the previous
archived run byte for byte. The new portable manifest records 2,094
changed field hashes. All 8,032 fields matched across macOS arm64, Linux
arm64 and emulated Linux x86_64.

For the blocked condition's remaining p-value maximum
(`ENSG00000168056.15`), the 80-digit posterior replay found a positive
MAP objective gain of about `1.24e-10` at iteration eight. Both binary64
paths measured a decrease of about `2.33e-10` at that step. The underlying
objective magnitudes were about `6.45e5`. Original DESeq2 refits using
the saved portable gene-wise estimates reproduced both remaining Rust
trend outputs within `7e-12`. In the fourth blocked split, a fixed-input
refit of `ENSG00000217228.1` still differed by `2.11e-6` between Rust
profiles with the same offsets and dispersion, after both fits reached
the L-BFGS-B fallback.

One paired four-worker run of the eleven portable contrasts took 34.84 s
with the new likelihood and 35.64 s with the prior build from `c6caa1b`.
Peak RSS was 264.27 MiB and 263.06 MiB, respectively.

The `adc4027` archives built offline and all four Cargo archives
reproduced byte for byte. The R archive passed all 2,874 expectations
with zero skips on Linux arm64 R 4.6.2 and on macOS R 4.6.2 and R 4.5.3.
Both R 4.6.2 checks recorded the existing compiled-symbol warning and
vendoring notes. The R 4.5.3 check recorded the vendoring notes.

### Portable GLM fallback finite differences

The portable fallback evaluates the per-sample negative-binomial loss
change across the bounded L-BFGS-B stencil. The loss change cancels
terms that depend only on the count and dispersion. The optimizer keeps
the R finite-difference step and uses the direct objective difference
for wide stencils or extreme means. An 80-digit decimal calculation
checks the one-sided stencil at a coefficient bound and the two-sided
stencil in a weighted design with a four-billion-count row.

The `96ae0c7` portable run passed 117 of 121 per-contrast gates and all
nine global gates with the original frozen references and tolerances.
The maximum fold-change error decreased from `1.50e-6` to `1.13e-6`
in blocked random split 03, clearing the gate. In blocked random split
04 it decreased from `5.63e-6` to `2.60e-6`. The remaining failures are
that fold-change gate, the blocked condition's p-value and fitted
dispersion gates, and blocked split 01's fitted dispersion gate. There
were zero significance flips and complete top-set overlap across the
eleven contrasts. The host profile passed all 130 gates. All 22 host
result and diagnostic TSVs matched `adc4027` byte for byte.

The change updates 620 of 8,032 portable field hashes. The Cargo
archives built offline and reproduced byte for byte. The R source
archive passed all 2,874 expectations with zero skips under macOS
R 4.5.3 and 4.6.2 and Linux arm64 R 4.6.2. The checks recorded the
same compiled-symbol warning and vendoring notes as `adc4027`.

In `e8b8969`, compensated summation of the signed loss changes brought
the portable release build to 118 of 121 per-contrast gates and all
nine global gates. Blocked split 04 cleared its fold-change gate at
`1.01e-6`, and the split 03 maximum decreased further to `1.07e-6`.
The remaining comparisons are the blocked condition's p-value and
fitted dispersion and blocked split 01's fitted dispersion. The
significance panel found zero flips and complete top-set overlap. The
host release build passed all 130 gates and matched the 22 archived
host TSVs byte for byte. All 8,032 portable field hashes matched on
macOS arm64, Linux arm64 and QEMU Linux x86_64. The workspace tests and
both full 56,200-gene replay tests passed on macOS and Linux arm64.
All four `e8b8969` Cargo archives built offline and reproduced byte for
byte. Each of the three local R checks again passed 2,874 expectations
with zero skips and recorded the same package-check notes and warning
profile as `96ae0c7`.

One paired four-worker run over the eleven contrasts took 17.12 s with
`e8b8969` and 16.98 s with `adc4027`. Peak RSS was 258.06 MiB and 259.47
MiB, respectively. The paired result files matched the precision runs.

Compensating the dispersion likelihood reduction introduced six new
failures across a separate four-contrast trial. A centered derivative
reduced the median error from about 418 million binary64 ulps to one
ulp on 200 cases against an 80-digit digamma recurrence. At count
and mean 1 and log dispersion −25, it recovered the negative derivative
`-6.943971932417719e-12` from the existing positive result `2.7800684085535026e-5`.
Applying that derivative at shape greater than 10 increased the
blocked-condition p-value error. Cutoffs at 1,024 and 65,536 left the
three gaps in the two blocked contrasts. The dispersion trials were
kept in local experiment worktrees.

An 80-digit replay of the three gene-wise fits with the largest
influence on the remaining trend differences found positive gains at
iteration 12 in both Rust profiles. They ranged from `1.50e-10` to
`5.75e-10`. The host profile rounded the two blocked split 01 gains
negative. The portable profile rounded the blocked condition gain
negative. In each case the rejected step changed the final gene-wise
estimate supplied to the trend fit.

### Portable MAP posterior gain recovery

When the portable MAP line search rejects a short step close to the
rounded Armijo bound, it uses two-point Gauss–Legendre quadrature of
the posterior score to estimate the gain. If that gain is close to the
resolution of the posterior and the score changes sign, it checks
whether a shorter step also meets the Armijo bound. A high-precision
test checks a case in which the gain of either step is smaller than
one binary64 ulp of the posterior.

The `39af7c0` portable release build passed 119 of 121 per-contrast
gates and all nine global gates using the frozen references. The
blocked-condition maximum p-value error fell from `7.35e-7` to
`4.35e-7`, clearing its gate. The two fitted-dispersion comparisons in
the blocked condition and blocked split 01 remain open. The
significance check found zero flips and full top-set overlap. All 130
host gates passed, and the 22 result and diagnostic TSVs matched
`e8b8969` byte for byte.

The change updates 506 of the 8,032 portable field hashes. All portable
hashes agree on macOS arm64, Linux arm64 and QEMU Linux x86_64. The host
manifests also passed on those three environments. The Rust workspace
and both full GTEx replay tests passed on macOS and Linux arm64. The
macOS arm64 release panel passed 114 tests. Rust 1.97.1 and Clippy
passed. The Intel macOS tests cross-linked, and the Windows MSVC
libraries and binaries passed type checking. Their native runtime
checks remain pending.

All four Cargo archives built offline and reproduced byte for byte.
The R archive passed all 2,874 expectations with zero skips on macOS
R 4.5.3 and 4.6.2 and Linux arm64 R 4.6.2. Both R 4.6.2 checks recorded
the existing compiled-symbol warning and vendoring notes. The R 4.5.3
check recorded the vendoring notes.

In one paired four-worker run, `39af7c0` took 18.07 s and `e8b8969`
took 17.14 s. Peak RSS was 259.55 MiB and 259.22 MiB, respectively.

Applying a stronger gain check in the gene-wise stage cleared the
blocked-condition fitted-dispersion gate on trial. It introduced two
other fitted-dispersion failures across the eleven-contrast panel.
Each new failure followed the changed estimate of one gene. In both
profiles, 80-digit replays found positive gains of about `1.6e-10`
at the rejected full steps. Both profiles had previously accepted a
shorter step for those rows. Centered scores and midpoint trials still
left regressions on the frozen panel, so the MAP-specific change is
the retained implementation.

### Reusable IRLS QR workspace

The default IRLS QR solver reuses its weighted design, response,
reflectors, basis and coefficient buffers across iterations. It keeps
the previous reduction order. The triangular solve preserves the LU
signed-zero arithmetic and uses the original LU path for nonfinite
entries and subnormal pivots. A separate allocation probe measured
zero allocation or reallocation calls across 64 warmed 61 × 5 solves.
The previous solver made 576 calls for the same inputs.

The `3ef6e9c` release binary matched `39af7c0` byte for byte across all
44 portable and host result and diagnostic TSVs. The portable run
passed 119 of 121 per-contrast and all nine global gates. The host run
passed all 130 gates. All 8,032 fields in each portable and host golden
panel matched on macOS arm64, Linux arm64 and emulated Linux x86_64.
The workspace tests and both full GTEx replay tests passed on macOS and
Linux arm64. The Intel macOS tests cross-linked, and Windows MSVC
libraries and binaries type-checked.

The two-repeat, four-worker benchmark over the eleven GTEx contrasts
took 21.55 s with reusable QR and 22.24 s with `39af7c0`. Throughput
increased 3.2%. Peak RSS was 259.81 MiB and 260.83 MiB, respectively.
The result files matched byte for byte. The earlier run overlapped
the release build and background indexing.

All four Cargo archives built offline and reproduced byte for byte.
The R source archive passed 2,874 expectations with zero skips in each
of the three local R environments and recorded the same R-version
package-check warning and notes as `39af7c0`.

At this checkpoint the opt-in normal-equations solver used an owned
LU factorization. A detached trial reused its workspace and matched the
existing coefficients across 64 solves
with zero warmed allocation calls. The existing solver made 320 calls
in that probe.

### Reusable IRLS normal-equations workspace

In `42075c0`, the opt-in solver reuses its crossproduct matrix, response
and active-column buffer across iterations. Row pivots permute the
response in place. A regression test matched the previous pivoted LU
bits across a successful solve, a singular solve and a successful
retry with the same buffers. The integrated allocation probe measured
zero calls across 64 warmed solves. The previous solver made 320 calls
on the same inputs.

The public normal-equations fit matched the independent DESeq2
high-dispersion covariance fixture in both profiles on macOS arm64,
Linux arm64 and emulated Linux x86_64. The two-profile GTEx validation
matched all 44 archived result and diagnostic TSVs from `3ef6e9c` byte
for byte. All four Cargo archives built offline and reproduced byte
for byte. The R archive passed 2,874 expectations with zero skips in
each local R environment, with the existing package-check notes and
R 4.6.2 compiled-symbol warning.

Four formulations of the small-shape dispersion likelihood grouped
the Stirling log terms before cancellation. The most accurate trial
reduced the median error from four binary64 ulps to one and the maximum
from 1,692 to 67 across 200 cases checked at 80-digit precision. Every
formulation introduced new failures when replayed through the two
blocked fits. These formulations stay in the experiment worktree.

### IRLS output retention

In `1d54903`, inference keeps the fitted means, leverage and covariance
outputs. Dispersion mean fitting returns its fitted means and rlog
transforms return their fitted coefficients. These paths omit the
leverage matrix and share reusable buffers for the covariance checks
that determine optimizer fallback. Extreme inputs still run the same
element-wise overflow checks. The public rlog fit functions keep their
full diagnostics.

All 44 GTEx result and diagnostic TSVs match `42075c0` byte for byte
across both profiles. The host profile passed all 130 gates. Portable
passed 119 of 121 per-contrast gates and all nine global gates. The
workspace and both golden manifests passed on macOS and Linux arm64.
The full macOS GTEx replay and all 47 portable transform tests passed.
A host-reference formula test reproduced the same polynomial discrepancy
when forced to run in the portable profile on the earlier executable.
The targeted IRLS tests and portable golden panel passed on emulated
Linux x86_64. Intel macOS tests linked and Windows MSVC libraries and
binaries type-checked.

Over two paired four-worker runs of all eleven contrasts, summed median
time was 15.80 s versus 15.93 s. Median RSS was 256.61 versus 258.21 MiB
and peak RSS was 259.97 versus 259.92 MiB. The 0.8% throughput increase
is below the 2% threshold for a GTEx optimization round. Seven paired
rlog runs per fixture produced identical output. Median RSS fell from
12.33 to 11.83 MiB for 4,096 genes and 16 samples, and from 12.19 to
11.45 MiB for 2,048 genes and 32 samples. Median throughput increased
5.4% on the first fixture and decreased 0.1% on the second. The first
fixture had one slower startup that reached 13.38 MiB, above the
baseline's 12.61 MiB peak.

All four Cargo archives built offline and reproduced byte for byte.
All three local R package checks passed 2,874 expectations with zero
skips. R 4.5.3 reported the existing two package notes. R 4.6.2 also
reported the existing compiled-symbol warning.

A 60-digit Newton step for the two remaining portable dispersion fits,
with the final Gamma regression rows fixed, increased the maximum
deviations from the independent DESeq2 references.

### Shared offsets and narrow QR projection

In `b8c5811`, weighted and expanded beta-prior fits pass shared size
factors through to IRLS. The intercept shortcut normalizes one row at
a time. Frozen rlog builds its gene-specific factor matrix directly
from the input offsets and intercepts. Focused integration checks
passed 90 host-profile and 34 portable-profile tests.

The detached one-column QR projection matched both 8,032-field golden
manifests and all eleven GTEx result TSVs. It improved throughput by
7.55% and 4.08% in two paired runs, with peak RSS changes of −8.11 and
+3.69 MiB. A 32-sample rlog fixture was 4.6% slower with the one-column
projection. In `b9b550e`, designs with up to four coefficients use the
one-column projection and wider designs use the batched projection.
The integrated version passed 31 IRLS unit tests, 72 host-profile
integration tests and 10 portable-profile integration tests. The two
commits together matched all 44 GTEx result and diagnostic files from
`1d54903` and both macOS golden manifests.

### Preliminary dispersion mean replay

An independent R 4.6.1 and DESeq2 1.52.0 replay reproduced the frozen
gene-wise estimates for the rows with the largest trend influence.
For `ENSG00000268892.1`, Rust also reproduced the dispersion within
`1.1e-14` when it received DESeq2's preliminary fitted means. The
largest difference in those means was `4.21e-4` in the blocked
condition and `5.20e-3` in blocked split 01. The R-compatible optimizer
reduced both differences below `7.2e-8` in the fixed-input replay.

A trial used the compatible optimizer for preliminary fits and checked
short MAP steps against the integrated score. The trial cleared both
open fitted-dispersion comparisons. Eight other portable comparisons
and seven host comparisons failed the full panel. All nine global gates
passed in each profile. Replaying two of the borderline convergence
decisions in DESeq2 with the trial's own inputs also changed the stopping
steps from the frozen output. The trial and replay inputs were archived
locally, and the optimizer policy was restored.

### Direct Householder projection for wider designs

In `6ab4479`, IRLS applies the stored reflectors directly to the weighted
response when the design has more than four coefficients. The triangular
solve uses the same reusable scratch. With four workers and seven paired
repetitions per fixture, rlog throughput rose 21.7% for 4,096 genes × 16
samples and 32.5% for 2,048 genes × 32 samples. Peak RSS changed from
11.80 to 12.25 MiB and from 11.56 to 11.45 MiB, respectively. The largest
change between builds was `2.05e-14`. On both fixtures, maximum absolute
error against DESeq2 1.52.0 remained below `6.3e-11`.

Exactly 16 rlog fields changed in each 8,032-field manifest. The full
portable manifest matched on macOS and Linux arm64. The toy and
automatic-dispersion portable cases also matched on emulated Linux
x86_64, where the full host manifest was regenerated. Targeted IRLS,
covariance and rlog tests passed. All 44 GTEx result and diagnostic
files again matched byte for byte across both profiles. A detached trial
that used the direct projection for narrow designs introduced six
failures across the six portable blocked contrasts, so those designs
retain their earlier projection order.

The combined macOS R 4.6.2 package check passed 2,874 expectations with
zero skips. It recorded the existing compiled-symbol warning and two
vendoring notes.

### Mutation workflow verification

The scheduled job verifies a nonempty mutation selection and propagates
baseline and infrastructure failures. Local discovery found mutants in
all six matrix entries. The real MAP sample caught a mutation of the
upper dispersion bound and identified an equivalent mutation that
discarded an unused grid-start override. The override was removed, and
both pre-change golden manifests matched after that cleanup.
