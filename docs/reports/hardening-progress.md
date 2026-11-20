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

### CLI arguments and matrix storage

In `1f1e817`, Wald and LRT use typed argument groups and share contrast
resolution, path staging and sidecar writing. The supplied beta-prior
routes also consume the grouped normalization and reporting inputs.
All 52 focused CLI tests passed. Both subcommands retained their flag
order and conflict diagnostics. `28fe321` preserves the command-specific
descriptions pinned by the help snapshots. A portable-profile replay
matched the result table and five sidecars byte for byte for each
subcommand.

`34b7908` adds checked reservations to generated and expanded matrices.
The shared allocator reports bytes requested when a reservation fails.
All four public rlog design and prior constructors return a typed
dimension error for an overflowing sample count. The affected transform,
Cook's distance, matrix and GLM tests passed.

In `9b83380`, fitted models share immutable design matrices and coefficient
names when cloned. Weighted summaries and expanded coefficients reuse
their per-row buffers. The focused design, weighted summary and GLM
tests passed. Both 8,032-field macOS golden manifests matched across
the three commits.

The workspace follows the later crate-boundary decision in `e7ad144`.
The three reusable `rcompat` crates are separate. The CLI remains in the
`rsdeseq2` package so its installed binary and Rust entry point ship
together.

### Portable QR and MAP convergence

In `b921f0a`, the portable profile uses direct Householder projection for
all design widths. For short MAP proposals, an integrated score checks
either rounded Armijo decision when the gain is close to the bound.
The portable optimizer also continues when the rounded positive gain
falls within the roundoff band of the existing convergence tolerance.
An R 4.6.1 replay confirmed that the reference continued past iteration
eight for `ENSG00000175899.14`, with a gain of `1.0063e-6`.

The release build passed all 121 per-contrast and nine global gates
in both profiles. Portable fitted-dispersion maxima fell from
`4.98e-7` to `3.02e-7` for the blocked condition and from `7.83e-7` to
`6.76e-7` for blocked split 01. The significance check found zero flips
at adjusted p-value cutoffs of 0.05 and 0.1. Every top-set check matched.
The host profile retains its narrow-design projection order after that
profile's direct-projection trial failed nine comparisons. All 22 host
result and diagnostic files matched the previous milestone byte for byte.

The portable manifest updates 2,228 of 8,032 fields. It matches across
macOS arm64, Linux arm64 and emulated Linux x86_64. The macOS and Linux
arm64 host manifests also matched their existing hashes.

Two paired four-worker runs of the eleven GTEx contrasts measured
16.856 s versus 17.265 s in summed median time. Throughput rose 2.43%
and peak RSS fell from 259.44 to 259.19 MiB. The timed outputs matched
the candidate precision run byte for byte.

The batched macOS checks passed 1,231 workspace tests, both separate
56,200-gene GTEx tests and 80 portable-profile integration tests. Clippy
and Rust 1.97.1 checks passed. Intel macOS test executables linked and
Windows MSVC libraries and binaries type-checked. The macOS R 4.6.2
package passed 2,874 expectations with zero skips and recorded the
existing compiled-symbol warning and two vendoring notes.

### Reduced-output covariance and coefficient buffers

In `42a0b8e`, dispersion mean fitting and transform-only rlog omit
the second covariance matrix product when the diagonal is provably
positive and all entries can be scaled safely. Uncertain cases use
the full calculation and keep the existing overflow and optimizer
fallback behavior. Both fits also omit the standard-error and
optimizer-start matrices. Public inference fits keep their diagnostics.

Each numerical profile passed 59 focused IRLS, covariance and rlog
tests, matched all 8,032 golden fields, and passed all 121 per-contrast
and nine global GTEx precision gates. All 44 GTEx result and diagnostic
TSVs matched the preceding checkpoint byte for byte. Clippy passed for
the library.

Two paired four-worker runs across eleven GTEx contrasts measured
15.452 s versus 16.274 s in summed median time. Throughput increased
5.32%, and peak RSS was 259.31 versus 257.97 MiB. Seven paired rlog
runs per fixture produced byte-identical outputs. Median throughput
was 0.6% lower for 4,096 genes and 16 samples and 3.2% higher for
2,048 genes and 32 samples. Median RSS changed by −0.03 MiB and
+0.08 MiB, respectively.

In `14b2700`, IRLS keeps the transposed design and its factorization
for the fit, and reuses the response and initial coefficients across
rows. The 32 focused IRLS tests and both rank-deficient tests passed,
as did Clippy and both 8,032-field golden checks.

### Batched local validation and reference tooling

The `915ae25` release build passed both sets of 130 GTEx precision
gates. All 44 result and diagnostic TSVs matched the preceding precision
checkpoint byte for byte. Two paired runs against `42a0b8e` measured
15.724 versus 15.687 s in summed median time and 259.52 versus 259.55
MiB in peak RSS. The initial workspace reuse produced no separate
throughput gain on GTEx. The seven paired rlog runs per fixture matched
byte for byte. The 16-sample case was 3.7% slower and the 32-sample
case varied by 0.1%.

In `ae83ba9`, the integration suite uses one support module for
reference readers, conversions and shared design fixtures. The macOS
workspace batch passed 1,232 unit and integration tests and one doctest.
Both 56,200-gene replay tests passed separately. Clippy and the locked
Rust 1.97.1 check passed. Python tooling tests passed, including the
separately enabled full GTEx fixture check.

The R source archive from `915ae25` passed 2,874 expectations with zero
skips on R 4.5.3 and 4.6.2. The checks retained the two vendoring notes.
R 4.6.2 also retained the compiled-symbol warning.

In `1de7cd5`, the standalone shrinkage crate carries its own four
unit-test fixtures. Two keep their original bytes and two reproduce
only the independent analytic and grid leaves used by the unit tests.
All four unpacked Rust archives built offline and passed 217 library
tests with the workspace's host R profile. An exploratory portable
run differed from the R student-component normalization reference
at the existing 16-ulp bounds. The archive check selects the same
host R profile as ordinary fixture tests.

The `just generators` entry point lists reference generators, output
paths and replay checks. The inventory test checks that each top-level
reference generator is registered.

### Portable Student-t uniform posteriors

In `8133dcf`, portable Student-t calculations use `libm` Gamma values
and the unit-shape beta identity. Interval probabilities use the outer
tail, and narrow intervals use local integration on the original
scale. Same-tail moments subtract small tail integrals to avoid losing
the component variance after rounding. The host R profiles retain the
original calculation and its reference.

The portable test checks 630 observations across 30 prior and
degrees-of-freedom combinations. Its independent reference evaluates
likelihoods, weights, conditional moments and CDFs with 100 working
decimal digits at the existing test tolerances. For the symmetric
10-degree case at estimate −15 and standard error 0.1, the posterior
SD is 1.56674834635824 against 1.56674834635824 from the independent
reference. The original ashr fixture records 1.59286532675430.

The focused Student-t integration test passed on macOS arm64, Linux
arm64 and emulated Linux x86_64. The macOS and Linux arm64 host R
bounded-posterior panels passed. Both macOS golden manifests matched
all 8,032 fields. Release targets now run the portable Student-t replay
alongside the existing portable golden and host R checks.

At `8133dcf`, the portable Gaussian panel reported the same 19 mismatches
as `a85ca59`, all in one ultra-narrow uniform fixture. For its estimate
−40 and standard error 8, the 90-digit midpoint log density is
−15.49838007488450867. The averaging error at the fixture width is
below `1e-24`. The portable replay and the original ashr fixture both
differ from that value by more than `1e-3`.

### Portable Gaussian uniform posteriors and release checks

In `96d2361`, Gaussian likelihoods and CDFs use differences in the outer
tail and share narrow-interval integration with Student-t. The portable
profile computes far-tail variances around the nearer endpoint and
combines mixture variances around component means. For the uniform
interval [−1, 1] at estimate 40 and standard error 0.01, the portable
posterior mean is 0.999997435897773 and its SD is 2.56410205836264e-6
within 1e-12 of the independent reference. The host R profile retains
its original ashr result.

The mathematical Gaussian reference checks 650 observations across ten
priors with unchanged tolerances, including the ultra-narrow case and
four far-tail standard deviations exposed by the independent comparison.
All bounded-posterior tests passed in both profiles on macOS arm64,
Linux arm64 and emulated Linux x86_64. Each release target now runs the
full bounded-posterior panel in both profiles.

The macOS checkpoint passed 1,233 Rust unit and integration tests and one
doctest, plus both separate 56,200-gene replays. The native release-profile
check passed 115 host R tests and all three portable bounded-posterior
tests. Both golden manifests matched their 8,032 fields. Clippy and the
Rust 1.97.1 checks passed. Intel macOS test executables linked, and Windows
MSVC libraries and binaries type-checked.

The release executable passed 130 GTEx precision gates in each numeric
profile. All 44 result and diagnostic TSVs matched the preceding checked
release. The four staged Rust archives reproduced byte for byte and
passed 217 library tests offline. The R source archive passed 2,874
expectations with zero skips on both R 4.5.3 and R 4.6.2. Both checks
retained the two vendoring notes, and R 4.6.2 also retained the
compiled-symbol warning.

### Structured design and GLM diagnostics

Rank validation reports the computed rank, column count and zero-column
labels. LRT nesting failures expose the full and combined ranks. Ambiguous
coefficient resolution retains the requested names and matching columns in
model-matrix order. All-zero GLM errors retain both fit and source rows,
including selected fits. Intercept arithmetic errors retain the rejected
value and the gene and sample when available. These cases retain their
existing messages and coarse error kinds.

Two targeted integration runs passed 150 and 100 checks respectively.
Workspace Clippy passed for all targets. The selected-row unit test passed
and verifies that fit and source indices can differ. The Rust assertions
in the affected validation paths now inspect structured fields. The
remaining text-only diagnostics and their Rust assertions are still under
migration.

Negative-binomial likelihoods now report the stage, rejected value and
available row and sample coordinates. Weighted term and partial-sum
overflows identify the particular sample. Cox–Reid column-sum overflows
report the actual sum and its gene, sample and design column. The messages
and coarse error kinds are preserved. Ordered-factor contrast tests now
check the structured prerequisite and canonical factor name through both
Wald and LRT entry points.

The targeted batches passed 88 integration tests and one unit test.

Formula and model-frame lookups now report unresolved names and
collisions by scope. The error includes the matching labels for a
collision and the factor-level role or observed sample index when
available. Model-frame and builder character contrasts share level
resolution. The tests verify exact-name priority, alias collision order
and retained diagnostics across Wald and LRT routes. Existing error
messages and coarse kinds are preserved.

The focused design, contrast and LRT run passed 210 tests. The final
assertion updates passed a 16-test Wald and formula run. Workspace Clippy
passed for all targets after the implementation changes.

Dispersion trends and variance priors now retain rejected arithmetic
values and the associated gene or coefficient when known. Failed
parametric predictions inside the fit loop map filtered rows back to
their original input indices. The no-eligible-genes error records
`min_disp` and the number of rows considered. A single-gene beta prior
now rejects variance overflow for a retained coefficient and continues
to use the configured wide variance for an intercept.

The targeted dispersion reference, trend, VST and beta prior batches
passed 137 integration tests and 12 unit tests.

Posterior work errors now distinguish cancellation, the requested
evaluation count exceeding its budget, and arithmetic overflow in grid
sizing. Invalid controls retain their values. Summary shape errors name
the array and axis, and coefficient mismatches identify both coefficients.
Regression fixture checks match variants for invalid likelihood shapes,
missing parameters, unsupported optimizers and missing compiled prefits.
The messages and coarse error kinds are preserved.

All 15 targeted posterior and regression tests passed, including the
exact grid budget boundary and cancellation during active work. Workspace
Clippy passed for all targets.

Optimizer validation now exposes the rejected control or array input
through structured variants. Numeric errors identify the observation and
component for likelihood cells, and the relevant index for weights and
Dirichlet parameters. Compiled EM parser errors reach Rust callers as
typed errors and retain the original diagnostic in the fixture transport.
The 19 targeted optimizer and mixture tests and workspace Clippy passed.

Formula validation now identifies the resource that exceeded its bound,
the observed count or lower bound, and the configured limit. Unsupported
functions retain their names and the required namespace when known.
Invalid wire lengths expose the declared count and remaining bytes, and
model-matrix response errors distinguish an invalid index from a factor.
All 150 targeted formula, design and model-matrix tests and workspace
Clippy passed.

Analysis preflight now identifies incompatible settings by variant.
Contradictory contrast metadata and linear-mu replacement are rejected
before formula evaluation. A beta-prior coefficient outside the standard
design now returns its bounds error. Stored LRT, replacement, sample-level
and beta-binomial validation expose the rejected settings, row or value.
The integration assertions on these paths now match variants and fields.

The 226-test workflow batch passed. The final preflight update passed
the six run and beta-prior tests, and workspace Clippy passed. The
run-equivalence generator reproduced the test source byte for byte
after formatting.

Checked dispersion operations now expose the failed value and the
available gene, coefficient or reduction index. The residual diagnostic
recovers the original gene index after selection. Invalid p-values,
unwritten output files and mismatched LRT offsets retain their relevant
indices, paths and values. Native error assertions now match structured
fields, and the CLI stderr tests continue to check the user-facing text.
Six targeted numeric tests, 25 application unit tests, 91 integration
tests and workspace Clippy passed.

### Size-factor estimation

The median pass filters nonfinite log geometric means once and reserves one
reusable sample buffer. This preserves control-gene ordering and multiplicity.
On macOS arm64, five alternating runs per executable compared the change
with `37f8c8f`. Size-factor output hashes matched in both numeric profiles
for all synthetic cases at 20,000 × 6, 20,000 × 100 and 60,000 × 800.

In the sparse ratio trials, 80% of the genes had a zero in one sample.
Median runtime at 60,000 × 800 fell from 142.7 to 98.0 ms in the portable
profile and from 128.3 to 86.5 ms in the host R profile. Dense ratio and
poscounts timings were within 4% at the two larger sizes. Row-wise
experiments with blocks of eight and 32 samples took 6–24% longer on
the largest matrix and required more ratio scratch space. The median
pass keeps one sample buffer and uses linear-time selection.

The batched macOS verification after the error migration and normalization
change passed 1,238 Rust tests and one doctest, plus both separate 56,200-gene
GTEx replays. Workspace Clippy passed. The optimized golden harness matched
8,032 fields in each of the host R and portable profiles at one and four
workers.

### Linear dispersion projection and IRLS setup

Linear dispersion fits reuse the normalized projection for rough dispersion
estimates and move that buffer into the count-scale means. The projection
borrows the design transpose and reserves its matrix products through the
shared allocation helper. Initial IRLS design and crossproduct buffers use
the same helper. Matrix products keep their previous grouping and column-major
storage.

The 166 focused tests and workspace Clippy passed. A seven-coefficient
regression checks every fitted mean bit against the original multiplication
on nalgebra's wide matrix path. Both the macOS host R and portable golden
manifests matched all 8,032 fields at one and four workers.

### Solver workspace allocation

IRLS reserves its per-gene status buffer and worker scratch through the
shared allocator. QR, normal equations and covariance calculations propagate
allocation errors separately from singular or nonfinite calculations.
Worker scratch is initialized at the first fitted row and reused.
Cox-Reid uses checked crossproduct and derivative buffers. If its optional
design cache cannot reserve memory, it uses the existing streaming path.
The IRLS optimizer penalty reuses the gradient and Hessian buffers.

The 164 focused Rust tests and workspace Clippy passed. The oversized
workspace regression covers QR, normal equations and covariance and checks
that each workspace can complete a valid solve after the error. Both
numeric profiles matched their 8,032-field golden manifest.

The covariance calculation now uses its existing product matrix for pivoted
LU on wider designs. The four-coefficient case uses nalgebra's original
closed form with a stack copy. The rare QR path for irregular pivots also
solves in its existing triangular buffer. Its lower triangle is cleared
before reuse. Tests compared covariance bits with nalgebra at four, five
and nine coefficients, including pivot swaps and reuse after a singular fit.
The QR comparison includes a subnormal pivot and a successful retry. Each
profile passed 85 focused tests, including the independent covariance
reference, and matched all 8,032 golden fields. Workspace Clippy passed.

### Mixture optimizer allocation

The mixture fitter normalizes each likelihood row in its final buffer.
SQP allocates its first Hessian when needed and reuses it across iterations.
Both SQP policies check the allocation of their active Hessians and
factorization inputs. The interior-point path scales the Hessian in place.
The optional low-rank approximation builds the augmented likelihood
directly in column-major storage and reserves its reconstructed output
after rank selection.

Each numeric profile passed three mixture unit tests and 15 targeted
integration tests, including the SQP and advanced optimizer references.
Workspace Clippy passed and both golden manifests matched all 8,032 fields.

### Posterior and polynomial workspaces

The negative-binomial, beta-binomial and general posterior fits share the
finite-difference Hessian calculation and reuse the perturbed coefficient
buffer. Laplace inversion reserves its wide output through the shared
allocator. Compiled beta-binomial fits reuse the first inverse when the
refinement check accepts it. BFGS updates the inverse metric in place and
reuses the metric and vectors across iterations. Polynomial extrapolation
returns allocation errors through the fitter and keeps its existing recovery
from a singular coefficient solve.

Each profile passed 22 focused tests and workspace Clippy passed. The BFGS
unit test compares the in-place rank-two update against the previous matrix
expression bit for bit.

### Batched allocation checkpoint

The macOS arm64 workspace at `6206eee` passed 1,242 unit and integration
tests and one doctest. Both separate 56,200-gene GTEx replays passed.
Both profiles matched the 8,032-field golden manifests at one and four
workers. The locked Rust 1.97.1 all-targets check passed. Intel macOS test
executables cross-linked, and Windows MSVC libraries and binaries
type-checked. Those platforms still need native runtime verification.

The optimized CLI passed all 130 precision gates in each numeric profile
over eleven GTEx contrasts. All 44 result and diagnostic TSVs matched
the outputs from `96d2361` byte for byte.

The self-contained R source archive passed 2,874 expectations with zero
failures and zero skips on both R 4.5.3 with DESeq2 1.50.2 and R 4.6.2 with
DESeq2 1.52.0. Both checks retained the hidden-file and portable-filename
notes from vendoring. R 4.6.2 also retained the compiled-symbol warning.

All four staged Cargo archives reproduced byte for byte. They built
offline from the packaged sources and passed 221 library tests in total.

### Typed posterior fit stages

Posterior batches select one row likelihood and pass its MAP result into the
Laplace stage when uncertainty is requested. A separate summary step writes
false-sign and threshold rates and applies grid results when available.
Primary and opposite-start warnings retain their output ordering, and
random starts retain one shared draw per batch.

The 23 focused posterior tests passed in each profile, all 11 shrinkage
unit tests passed under the host R profile, and workspace Clippy passed.
Both golden manifests matched all 8,032 fields.

### Labelled TSV parsing

Count, design, normalization-factor and observation-weight readers share a
typed matrix parser. Matrix and two-column readers reuse one CSV record per
input. The reader policies preserve the existing label checks and error
locations. A regression checks the line of a ragged count record and the
flattened cell index of a nonfinite design value after a quoted label.

The 55 I/O unit tests, 47 TSV and CLI integration tests, and workspace
Clippy passed.

### Minimum supported Rust version

The workspace and R source package now require Rust 1.96.0. This is the first
release with the stable `core::range::Range` type used by the public matrix
index accessors. The workspace, R backend, CI jobs and installation docs
declare the same version.

On macOS arm64, Rust 1.96.0 passed the locked offline all-targets check and
221 workspace library tests. Both numeric profiles matched their 8,032-field
golden manifests. The vendored source archive installed with Rust 1.96.0 on
R 4.5.3 and passed 189 native bridge and fit-object expectations with zero
failures or skips.

### Platform runtime checkpoint

The release panel now includes the posterior RNG, regression and TSV tests.
It runs the host R golden check on every target with a committed host
manifest, as well as the portable golden check on every release target.

The `6eeb479` source passed the release-profile checks below.

| Target | Host R tests | Portable bounded-posterior tests | Golden fields in each profile |
| --- | ---: | ---: | ---: |
| macOS arm64 | 126 | 3 | 8,032 |
| Ubuntu 24.04 arm64 | 142 | 3 | 8,032 |
| Ubuntu 24.04 x86_64 under QEMU 10.0.13 | 115 in the initial panel, 11 in the added group | 3 | 8,032 |

The added group also passed under the portable profile on Linux x86_64.
Both separate 56,200-gene replays passed under each numeric profile on
both Linux targets. The fixture hashes were checked before running those
replays. Intel macOS test executables linked, and Windows MSVC libraries
and binaries type-checked. Native runtime results for those two targets
and the hosted release run remain open.

### macOS R shared-library linking

The macOS package exposes the R initialization symbol and strips unused code.
The `.Call` registration table retains all 30 routines with the same arities.
A paired arm64 build with R 4.6.2 and Rust 1.96.0 reduced the shared library
from 9,151,320 to 5,934,904 bytes. R's compiled-code checker stopped reporting
the two `exit` imports and still reported the Rust runtime `abort` import.

The stripped build passed all 2,874 R expectations with zero failures
or skips. A SIGINT interrupted a posterior grid, and the same R process
then completed a second native fit.

### R backend release profile

The R backend builds in a separate Cargo workspace. Its release profile now
uses thin LTO, one codegen unit and panic unwinding. In a paired macOS arm64
build with R 4.6.2 and Rust 1.96.0, the stripped library fell from 5,934,904
to 5,091,368 bytes. The registration export and the Rust runtime `abort`
diagnostic were present in both builds.

The new build passed 2,874 expectations with zero failures or skips, and
a SIGINT cancelled a posterior calculation before the same R process
completed a second fit. A single paired workflow pass matched all four
payloads exactly. For the synthetic 56,200-gene × 16-sample batch-results
case, elapsed time fell from 11.00 to 7.88 s and peak RSS rose from 975 to
1,010 MiB. The ashr, apeglm and rlog cases differed by less than 2% in
elapsed time and less than 1 MiB in peak RSS.

### R bridge buffer validation

Shared Rust slice helpers now check alignment and the maximum byte span
before borrowing nonempty C buffers. Count conversion, dense matrix
conversion, expanded contrast coefficients and apeglm contrast matrices
use the same fallible allocator as the core.

The locked offline R build used Rust 1.96.0. Both focused Rust boundary
tests and all 361 selected expectations on R 4.6.2 passed without
failures or skips. All 40 posterior analysis batches matched standalone
Rust output exactly. A SIGINT cancelled the posterior calculation, and
the same R process completed a second fit.

### R bridge output initialization

The Rust bridge borrows fresh R output storage as slices of
`MaybeUninit<T>` and writes values directly into those slots. This
avoids a separate initialization pass over large matrices. The formula
owner slot and the duplicated columns used in result presentation use
the checked helper for already initialized storage.

On macOS arm64, all five Rust bridge tests passed with Rust 1.96.0, and
the replacement-count copy into uninitialized storage passed under Miri.
R 4.6.2 passed all 2,874 package expectations without failures or skips.
The 40 posterior batches still matched the standalone Rust output
exactly, and a second native fit succeeded after SIGINT cancellation.

### Sample offsets in R contrast shrinkage

Expanded-model normal shrinkage passes sample size factors directly to
the GLM. This removes the genes × samples copy and lets IRLS compute
the log offsets once per sample. Per-cell normalization factors retain
precedence when supplied.

All 72 selected R 4.6.2 expectations passed. A 128-gene test matched
all baseline output bytes with shared, weighted, and varying per-cell
offsets. In one paired synthetic run with 20,000 genes and 100 samples,
both builds produced byte-identical output in 0.21 s, and peak RSS
fell from 194.5 to 178.9 MiB.

### Direct normalized-count output in R

The R wrapper borrows the original column-major counts and per-cell
factors. The core streams checked divisions into the R result matrix.
Invalid factors and division overflow still report the first failing
cell in logical row order.

All 24 selected Rust tests and 88 R 4.6.2 expectations passed. Paired
synthetic runs with 60,000 genes and 128 samples produced byte-identical
results. Peak RSS fell from 269.8 to 182.0 MiB with sample factors and
from 387.0 to 240.6 MiB with per-cell factors. The smaller probe also
matched all four error messages exactly.

### Borrowed base metadata in R

The R base mean and variance path borrows counts, factors and weights.
It reuses one sample-sized row, computes the mean once, and writes the
mean, variance and raw all-zero flag directly into the result buffers.
Normalization and weight errors retain their precedence over failures
in the row statistics.

All 25 selected Rust tests and 264 R 4.6.2 expectations passed. Two
paired 60,000-gene × 128-sample cases produced byte-identical output.
Peak RSS fell from 212.2 to 182.0 MiB with sample factors and from
388.1 to 299.3 MiB with per-cell factors and weights. The small-case
comparison matched single-sample variance, raw-zero flags and all
four error messages exactly.

### Borrowed size-factor estimation in R

The estimator reads validated counts directly from R's column-major
matrix. The default call skips building a control-index array in C.
Explicit controls retain their order and multiplicity. The standalone
Rust path retains its direct row-slice access.

All 14 selected Rust tests and 260 R 4.6.2 expectations passed. Dense
ratio, sparse ratio and poscounts runs at 60,000 genes by 128 samples
produced byte-identical R output, with about 30 MiB less peak RSS per
run. Standalone Rust checksums matched in both numeric profiles. The
60,000-gene × 800-sample timings differed by at most 3.4%. An initial
5.3% difference for the medium portable ratio case was 0.1% over 12
additional alternating runs.

### Direct formula-offset output in R

Formula offset composition shares factor validation and offset
multipliers between the owned Rust matrix and the direct R output.
Sample-specific products are computed once per sample. The R path
writes the effective factors in column order and preserves the first
logical row error when multiple products fail.

The focused Rust tests and 277 R 4.6.2 expectations passed. Paired
60,000-gene × 128-sample runs matched byte for byte with sample factors
and per-cell factors, and both used about 58 MiB less peak RSS.
The small-case comparison matched all five diagnostics exactly.

### Batched R bridge and normalization verification

On macOS arm64, source commit `b6da6e8` passed 1,247 Rust unit and
integration tests and one doctest. Both separately selected GTEx tests
passed for all 56,200 genes in each numeric profile, using the
SHA-256-verified full fixtures. Both golden manifests matched all
8,032 fields at one and four workers. Workspace formatting, Clippy and
the locked all-targets check with Rust 1.96.0 passed.

One source archive built from the committed files installed offline
with Rust 1.96.0 in R 4.5.3 with DESeq2 1.50.2 and R 4.6.2 with
DESeq2 1.52.0. Its SHA-256 was
`e8ce1d3d23a719b3120ce3378311c87209d9ef12ef448eeb3c653e568873b435`.
Each installation passed all 2,874 expectations without skips or
failures and completed a fresh native fit after SIGINT cancellation.
Both shared libraries exported only the R initialization symbol.

The archive build still warns about long filenames in vendored
dependencies. The hosted release matrix and native Intel macOS and
Windows runtime results remain open.

### Portable R source archive

The source package stores registry dependencies and their Cargo checksums
in an XZ archive, as recommended by the [CRAN Rust packaging guidance](https://cran.r-project.org/web/packages/using_rust.html).
Configure extracts the archive and builds with the bundled offline
configuration. The Windows configure entry point generates the offline
Rtools makefile through the same extraction step.

The outer archive fell from 15,897,994 to 8,616,840 bytes. None of its
540 paths exceeded 100 bytes or had hidden components. All 4,308
previously shipped registry files matched byte for byte, and three
registry metadata files survived packaging that the old build had
omitted. Linux arm64 unpacked the archive and resolved the locked
dependencies offline. The Windows configure script was exercised
under a macOS shell. Native Windows compilation remains open.

On macOS arm64, R 4.6.2 passed the hidden-file and portable-filename
checks and retained one compiled-code warning for the Rust runtime's
`abort` symbol. Examples and tests were disabled for that package
check. The same archive installed offline on R 4.5.3 and passed all
80 selected native-boundary and row-preparation expectations.

### Windows x64 R library cross-link

The Windows package build selects Rust's UCRT target and uses LLVM-MinGW
for C compilation and linkage. Rtools provides the build utilities. The
linker resolves the Windows imports required by the Rust static library
and links the unwind library statically. The DLL exports
`R_init_rsdeseq2` for routine registration.

The generated Windows makefile built the offline Rust archive with
Rust 1.96.0 and linked all 14 C translation units against CRAN's R 4.5.3
headers and DLL on Linux arm64. The check used R's own Windows
`Makeconf` and `winshlib.mk`. The x86-64 DLL imported R and UCRT and
had no `libunwind.dll` or `msvcrt.dll` import. Native Windows
installation and runtime checks remain open.

### Native Windows R workflow

The R workflow installs the Linux reference lane's checked source archive
on Windows x64 with Cargo offline. The smoke script checks DLL
registration, normalization, recovery after a native error and agreement
between one- and two-worker fits. The Windows job is required by the
reusable R workflow used in release validation.

The smoke script passed on the macOS arm64 R 4.5.3 and 4.6.2
installations. Workflow lint passed and the x64 LLVM-MinGW driver path
was verified in the pinned ZIP. The hosted Windows run remains open.

### Mutation workflow verification

The scheduled job verifies a nonempty mutation selection and propagates
baseline and infrastructure failures. Local discovery found mutants in
all six matrix entries. The real MAP sample caught a mutation of the
upper dispersion bound and identified an equivalent mutation that
discarded an unused grid-start override. The override was removed, and
both pre-change golden manifests matched after that cleanup.

### Expanded beta-prior formula offsets

The standalone expanded beta-prior result helpers compute effective sample
factors once and pass them as shared GLM offsets. This reduces offset
storage from genes × samples to samples for these formula fits. Supplied
cell factors still use a matrix, with checked allocation for the composed
factors. Coefficient and contrast results share the fit and replacement
implementation. The Cook's step reuses the normalized counts for both
distance calculation and replacement planning.

All 53 focused results tests passed under both macOS numeric profiles.
The formula cases compare coefficient and contrast fits against explicit
cell factors, including nonuniform observation weights and an actual
Cook's replacement refit.

### R installation toolchain and build jobs

The R installer checks Cargo and rustc against the backend manifest and
reports both versions before compilation. It can select a toolchain in the
personal Cargo directory when the tools on PATH are missing or too old.
Vendor generation checks the Rust minimum in the workspace, the backend
manifest and the R package metadata. Archived installs use two Cargo jobs
and begin compiling C after the Rust build finishes. The archived Cargo
recipe clears the external jobserver settings so Cargo honors its job limit.
[CRAN guidance for Rust packages](https://cran.r-project.org/web/packages/using_rust.html)

On macOS arm64, the source archive installed with an empty Cargo home,
R 4.6.2, Rust 1.96.0 and parallel Make. The install log showed the requested
job limit and build order, and the native install smoke test passed. Local
shell checks exercised both generated recipes, the personal toolchain
fallback, missing and outdated tools, and metadata version mismatches.

### Borrowed rows in fixed-dispersion inference

The Wald, LRT and expanded beta-prior workflows borrow counts, offsets,
dispersion values and weights from the original dataset. Mixed all-zero
or weight-failed inputs pass the same selected row indices to the fit and
likelihood calculation. Fits with every row active return the original
output buffers directly. The intercept shortcut reuses its working-weight
buffer across genes.

The 49 GLM unit tests and 268 pipeline and results integration tests passed
under each macOS numeric profile. Both optimized golden checks matched all
8,032 fields at one and four workers. The host R profile also passed the
three standalone beta-prior workflow cases and the compact fixed-dispersion
GTEx reference; the compact GTEx reference passed under the portable profile.

The additional `many_samples` beta-prior comparison exceeded the
`1e-9` threshold under the portable profile. The maximum relative gap
was `3.007726454745714e-7` in the p-value for `gene0140`, with zero
missingness differences or significance flips. The same command on
`dfe9e1c` reproduced the exact gap.

### Portable MAP stopping on the beta-prior reference

For `gene0140`, the portable MAP line search accepted an estimated gain of
`1.93e-11` after its rounded objective had declined. One objective ulp
was `2.91e-11`. The recovered step kept the same score sign and shifted
the fitted dispersion to `0.16218280510633146`. In R 4.6.2 with DESeq2
1.52.0, the trace stayed at `0.16218277997255204` from iterations 9 to 34.
Iteration 35 moved it to `0.16218277997255501`.

After an accepted step, the portable profile retains that point when a
later recovered gain falls below both the convergence tolerance and one
objective ulp and the original proposal keeps the score sign. Proposals
that cross zero can still use the existing shorter-step recovery. The
500-sample beta-prior panel's maximum relative error fell from
`3.008e-7` to `1.484e-12`. Release validation now runs the frozen
beta-prior cases under the portable profile.

The portable CLI passed all 121 per-contrast and nine global precision
gates on SHA-verified GTEx references. The eleven contrasts had zero
significance flips at adjusted p-value cutoffs 0.05 and 0.1 and complete
top-set overlap. The portable golden manifest updates 86 fields from
the compact GTEx Wald and LRT workflows with normalization factors.
All 8,032 host R fields matched the existing macOS manifest.
On macOS, 77 focused tests passed in each profile, along with Clippy
and the locked Rust 1.96.0 library check.
The revised portable manifest and all three beta-prior tests also passed
on Linux arm64 and Linux x86_64 under QEMU 10.0.13.

### Sample-only offsets in formula workflows

Formula runs with sample-specific offsets compose the effective factors
once and fit with a shared sample vector. The returned fit and any Cook's
refit materialize the normalization-factor matrix for public metadata.
Dispersion starts retain the column-mean summation order of the equivalent
dense normalization matrix. Cloned builders share supplied cell factors
and observation weights.

A paired macOS arm64 fixed-dispersion Wald probe used 20,000 genes, 100
samples, and two workers. Across three runs per version, median peak RSS
fell from 114.39 MiB at `09b8923` to 83.86 MiB, a 26.7% reduction. The
normalization-factor and coefficient checksums matched on all three pairs.

The 174 focused tests passed under both numeric profiles. The formula
comparison checks early size factors, linear gene-wise dispersions, and
expanded beta-prior results against explicit cell factors. Both 8,032-field
golden checks matched. Clippy, format, and the module-cycle check passed.

### Batched local 0.6.0 release validation

Source `02f2841` passed 1,248 workspace unit and integration tests and
one doctest on macOS arm64. Both explicitly selected 56,200-gene GTEx
tests passed under each numeric profile after checking the fixture
hashes. Both 8,032-field golden manifests matched. All-targets Clippy,
the locked Rust 1.96.0 all-targets check, and formatting passed.

The four staged Cargo archives regenerated byte for byte. All four
built from the packaged sources offline and passed 221 library tests.
The R source archive (SHA-256
`9f7cf1ad4ff597f4e01990f96c3a628a05117b0c19185cf5e3142bd424d05880`)
built offline with Rust 1.96.0 in empty Cargo homes. R 4.5.3 with
DESeq2 1.50.2 and R 4.6.2 with DESeq2 1.52.0 each passed 2,874
expectations in 150 tests with zero skips and failures. The native
smoke checks and real SIGINT recovery passed in both installations.
Each shared library exported only `R_init_rsdeseq2`. R 4.6.2 passed
package examples and tests. The Rust runtime's `abort` symbol produced
the single compiled-code warning.

The compact end-to-end GTEx gate now detects a positive probability
rounding to zero and rejects non-finite values. Release validation
runs the gate under the portable profile. Both focused tests passed
under both macOS profiles and under the portable profile on Linux
arm64 and Linux x86_64 under QEMU, with zero significance flips at
both cutoffs in each pipeline run. Native Intel macOS and Windows
runtime validation remains open.

### IRLS optimizer fallback scratch

The bounded fallback reuses its analytic gradient and Hessian across
genes. Within one gene, the objective and gradient callbacks reuse
an evaluation at identical coefficient bits. The portable finite
difference stencil, inference means and working weights reuse storage
across rows. Unfloored means go directly into the output row.

One paired macOS arm64 allocation probe forced fallback on all 384 rows
of a 61-sample, four-coefficient design. The portable profile made
71,028 allocation or reallocation calls versus 84,097 at `8617b33`.
The host R profile made 73,362 versus 84,751. Numerical output and
optimizer-diagnostic checksums matched for both pairs.

The 104 focused tests passed under each profile, including the
GTEx optimizer reference and a scratch reuse check across genes.
Both 8,032-field golden manifests matched. Clippy and formatting passed.

### Borrowed rows and frozen offsets in rlog

Estimated and frozen rlog pass active gene indices into IRLS and borrow
the original counts and factors. Prior estimation reuses its normalized
counts buffer for log-fold changes. Frozen rlog validates the scaled
factors in row order and composes them in a reusable sample-sized worker
buffer during fitting.

A paired portable-profile probe on macOS arm64 used 8,192 genes, 16
samples, per-cell factors, two workers, and an all-zero row every five
genes. Across three runs, median peak RSS fell from 13.98 to 11.09 MiB
for estimated rlog and from 10.64 to 8.03 MiB for frozen rlog against
`4c4aa82`. Transformed values, intercepts and prior checksums matched on
all pairs.

Each profile passed 140 focused GLM, normalization and transform tests.
Both 8,032-field golden manifests matched. Clippy, formatting and the
module dependency check passed. The separate exact-R formula-transform
panel passed in the host R profile. Its degree-10 polynomial bit
comparison also fails under the portable profile at the baseline commit.

### Weighted rlog prior storage

Prior estimation stores one weight per gene. The shared weighted quantile
merges tied values and accumulates normalized weights in the sorted pair
buffer. Stable sorting preserves tie order, and the cumulative search
retains the existing cutoff and interpolation. The regression covers
tied absolute values and a fractional position.

A 20,000-gene, 100-sample probe on macOS arm64 used four workers and
three paired runs per profile. Median peak RSS fell from 126.89 to
81.34 MiB for the portable profile and from 126.94 to 81.34 MiB for
the host R profile against `7ee8515`. Prior variance bits matched
on every pair. Each profile passed 78 focused tests and its 8,032-field
golden check. A single workspace run for the combined rlog changes
passed 1,251 unit and integration tests plus one doctest. The default
workspace run skipped the two full GTEx replays previously run
separately.

### Native rlog smoke and full GTEx probability replay

The installed-package smoke test now calls estimated and frozen rlog
with per-cell factors and an interior all-zero gene. It compares the
active results against a separate compact fit, using two workers for
the masked fit and one for the compact fit. It also checks the zero
row and the returned intercepts and prior.

The R source at `c293eec` produced archive SHA-256
`d8e2fbacc40feb9cce46a9f24b48077600c2dfa344d0dc9498e45ca3138ebdbe`.
R 4.6.2 on macOS arm64 installed it with Rust 1.96.0 and an empty
offline Cargo home. The expanded smoke test passed, as did 113
expectations across eight transform tests against DESeq2 1.52.0 with
zero skips. The same archive cross-linked against Windows R 4.5.3
with Rust 1.96.0. The x64 UCRT DLL exported only `R_init_rsdeseq2`.
The Intel macOS test executables cross-linked. Native Windows and
Intel macOS runtime results remain open.

The fixed-dispersion compact and full GTEx panels now also require a
maximum absolute log10 error of `1e-5` for p-values and adjusted
p-values. The release workflow replays the full fixture under both
numerical profiles with one release test build. The frozen reference
hashes and compact selection were verified before replay. The full
normalization and Wald tests passed under both profiles on macOS arm64,
Linux arm64 and Linux x86_64 under QEMU. On all three targets, maximum
log10 error was `4.142e-7` under the host R profile and `7.167e-8`
under the portable profile for both probabilities. No significance
decisions flipped at either cutoff. The compact checks passed under
both profiles on macOS. Integration Clippy, the Rust 1.96.0 all-targets
check, workflow lint, formatting and documentation checks passed.

### Observation-weight rank workspace and diagnostics

Observation-weight preprocessing reuses one design matrix and two
retained-index buffers across genes. Both the weighted design and
each thresholded Cox–Reid design use the same elimination routine
and tolerance as before. The normalized weight matrix and scratch
buffers reserve capacity with checked allocation.

Weight input and weighted base-metadata multiplication errors
expose the zero-based gene and sample and the rejected value through
`CalculationLocation::GeneSample`. Their displayed messages retain
the established indices. The rank regression alternates restricted
and full thresholded designs to exercise buffer reuse.

A paired macOS arm64 probe preprocessed 20,000 genes and 100 samples,
with three design columns and a mix of usable, all-zero and rank-deficient
weight rows. Allocation and reallocation calls fell from 198,385 at
`a8dc71b` to five. Median time over three runs fell from 64.2 to
52.7 ms. The normalized-weight and rank-flag checksum matched on all
three pairs.

Each profile passed 132 focused integration tests and the 8,032-field
golden check. The final rank and weight regression run passed nine
tests in each profile. The numeric rank unit check, Clippy, formatting
and the module dependency check also passed.

### IRLS likelihood scratch and diagnostic retention

IRLS now reuses the negative-binomial gamma-term buffer across rows
in the same Rayon task. Means-only and transform fits skip the
inference diagnostic vectors and evaluate the final likelihoods in
fit order without collecting them. Full inference stores replacement
likelihoods only for rows the optimizer refits and writes those
values back to their original fit positions.

A paired macOS arm64 probe used 8,192 genes, eight samples, two workers
and the portable profile. Across three runs, median allocation and
reallocation calls fell from 16,582 to 8,393 for the rlog transform
and from 24,789 to 16,605 for the rlog inference fit against `b624785`.
Transform, likelihood and standard-error checksums matched on all
pairs. Median runtime was within approximately 5% on both workloads.

Each numeric profile passed 180 focused IRLS, likelihood, dispersion
and transform tests and matched all 8,032 golden fields. Integration
Clippy, formatting, documentation and module dependency checks passed.

### Batched observation-weight and IRLS release checks

The `8a0aeac` source passed 1,251 Rust unit and integration tests plus
one doctest. Both compact end-to-end GTEx tests and both full
56,200-gene fixture tests passed under each numeric profile, with
zero significance flips at either cutoff. Workspace Clippy and the
locked Rust 1.96.0 all-targets check passed.

The exact R source archive (SHA-256
`1dfa86df2dd1b56c85bf506d3e726c2deb149225892eddae438a9a3eb8e593e5`)
installed offline with Rust 1.96.0 and fresh Cargo homes. R 4.5.3
with DESeq2 1.50.2 and R 4.6.2 with DESeq2 1.52.0 each passed 2,874
expectations across 150 tests with zero skips and failures. The
expanded native smoke and actual SIGINT recovery passed in both R
installations. Their only exported routine was `R_init_rsdeseq2`.
Both binaries still import the Rust runtime's `abort` symbol.

The four staged Cargo archives regenerated byte for byte. The verifier
now extracts regular files and directories beneath the expected crate
name, using the same path checks on Python versions without the
`tarfile` extraction filter. Tests cover executable files, traversal,
foreign roots and links. The final extraction reproduced all 818 files
from the four archives byte for byte. All four archives built offline
and passed 222 library tests under the system Python 3.9 verifier.

The remaining wide-covariance allocations in the previous rlog probe
came from matrix multiplication packing. A detached macOS arm64
trial avoided packing for six to twelve coefficients and matched the
fit checksums. The first version cut the eight-sample transform from
8,393 to 203 median allocation calls, with runtime about 3% higher.
With fixed-size columns, the median runtime improved on the
five- and six-sample probes and varied by less than 3% on the eight-
and eleven-sample probes. The accepted release checks above cover
the nalgebra implementation.

### Cook's robust dispersion and post-refit masking

Cook's distance groups samples once from borrowed design rows. The
trimmed variance pass processes independent genes in parallel with
one sort buffer per Rayon task. It traverses cells and samples in the
original order when checking for arithmetic errors. Replacement uses
the same scratch allocation pattern for trimmed means. Post-refit
`maxCooks` substitutes zero during the scan and avoids a full copy of
the Cook's distance matrix.

A paired macOS arm64 probe used 20,000 genes, 100 samples, four
replicated design cells and the portable profile. With four workers,
the median over three runs for robust dispersion fell from 44.6 to
7.8 ms against `67370a4`. Cook's distance, replacement and masked
maxima together fell from 55.9 to 17.1 ms. Their allocation and
reallocation calls fell from 440,453 to 192, and median peak RSS
fell from 107.14 to 91.58 MiB. Checksums matched under one and four
workers for every pair.

Each numeric profile passed 107 focused Cook's, replacement,
thread-invariance and end-to-end GTEx checks, and both 8,032-field
golden manifests matched. The regression covers the first failing
gene, masked and unmasked NaNs, all-sample masking and signed zero
in model-matrix cells. Integration Clippy, formatting, documentation
and the module dependency check passed.

### Shared counts in Cook's replacement

Count matrices share immutable count buffers and axis labels when
cloned. Cook's replacement builds the requested candidate counts and
outlier flags, then copies original counts only if an eligible
replacement changes a value. When a copy is needed, workers apply the
changed rows in parallel. Saturated designs and replacement decisions
that preserve all counts reuse the input buffer.

The portable profile's 60,000-gene, 800-sample no-change replacement
probe used four workers. Across five pairs against `091f786`, median
peak RSS fell from 1,337.02 to 1,147.05 MiB and median time fell from
152.2 to 140.8 ms. A separate 20,000-gene, 100-sample stress probe
marked every cell as an outlier. Median time there was 6.41 ms versus
6.08 ms and peak RSS was 57.78 versus 60.16 MiB. Replacement counts,
candidate counts and refit-flag checksums matched on all pairs.

The 60,000 × 800 robust-dispersion probe used four replicated design
cells. Median times with 1, 2, 4 and 8 workers were 0.968, 0.467, 0.383
and 0.198 s. All worker counts and paired runs matched the original
`67370a4` dispersions. At eight workers, that source took 1.289 s.
Current median peak RSS at eight workers was 552.62 MiB.

Both profiles passed 183 focused tests and matched all 8,032 golden
fields. The regression checks shared backing storage when flagged
outliers preserve the counts and a separate count buffer when an
eligible replacement changes a value. Integration Clippy, formatting,
the module dependency check and documentation checks passed.

### Reused normalized counts in Cook's refit plans

The builder and expanded beta-prior workflows pass their owned
normalized counts into refit planning. The planner renormalizes only
cells whose raw counts changed. The public planner continues to
compute replacement normalization from the counts and offsets, so
caller-supplied normalized counts still act only as the seed for
replacement predictions.

A paired 20,000-gene, 100-sample fixed-dispersion Wald run used four
workers and the portable profile. Across nine pairs against
`f34d0ab`, median peak RSS fell from 111.69 to 96.36 MiB. Median
runtime was 180.0 versus 177.8 ms. Normalized counts, weighted
metadata in the targeted regression, and final probabilities
matched their reference. The nine timed pairs matched on all
normalized-count, base-metadata and result checksums.

Each profile passed 211 focused checks and matched all 8,032 golden
fields. The regression exercises both offset layouts with weighted
base metadata and checks that the owned matrix retains its buffer.
It also verifies that the public planner accepts a seed matrix
whose values differ from the normalized counts it returns.
All-targets Clippy, formatting and the module dependency check passed.

### Batched Cook's release verification

The release CLI at `c9aa13d` includes all three Cook's changes.
The paired baseline is `67370a4`, the commit before those changes.
Both were built with Rust 1.99.0 on macOS arm64. The eleven-contrast
gate used the frozen reference hashes and four Rayon threads per fit.

| Arithmetic profile | Dispersion gates | Minimum target reduction | Minimum global reduction | Result and diagnostic files matching baseline |
| --- | ---: | ---: | ---: | ---: |
| Portable | 130/130 | 236× | 414× | 22/22 |
| R host | 130/130 | 664× | 597× | 22/22 |

Both profiles had zero adjusted p-value decision flips at 0.05 and
0.1 against DESeq2. The top 100, 500 and 2,000 genes ranked by
p-value matched in every contrast.

One alternating-order pair per contrast with four threads and the
portable profile measured 14.657 s for the baseline and 13.681 s for
the candidate. The 1.0714× throughput exceeded the 1.02× gate, and
all result-only TSV hashes matched. Maximum peak RSS was 258.91 MiB
for the baseline and 259.11 MiB for the candidate.

The 3,745-gene compact end-to-end gate passed in both profiles through
`DeseqBuilder::run`. It now checks sample and gene alignment, all 61
estimated size factors, and all gene-level replacement flags against
DESeq2, including the 50 genes flagged for replacement. The existing
default-filtering, missingness, probability and decision bounds passed.

The generated R source archive at `c9aa13d` has SHA-256
`5a1f8180e9037877fbd2e02c883199da70b47e7b22abb635af6154b6f1aa9c6f`.
It installed offline into fresh R libraries using empty Cargo homes
and Rust 1.96.0. The R 4.5.3 and R 4.6.2 lanes each passed 150 tests
and 2,874 expectations with zero skips or failures. Both passed the
installed-package native smoke test.

## Review plan verification for 0.6.0

| Priority | Implementation and evidence |
| --- | --- |
| 1. Golden harness | The [golden harness](../development.md) covers 8,032 fields. Portable hashes have replayed on macOS arm64, Linux arm64 and emulated Linux x86_64. |
| 2. Per-gene parallelism | MAP, Cook's distance, robust dispersion and replacement, and the rlog prior run in parallel with ordered errors and results. See the [thread-scaling probe](#local-verification-and-benchmark). |
| 3. Borrowed offsets | Sample factors stay in a sample vector through fitting. Selected dispersion, Wald, LRT, beta-prior and rlog fits borrow source rows. Where the API returns gene-by-sample factors, the matrix is created for that output. |
| 4. Numerical profiles | The crate and CLI default to portable arithmetic. The R bridge selects its host profile. The [platform workflow](../../.github/workflows/platform-checks.yml) checks the release panel and portable golden hashes on five targets. Native Intel macOS and Windows x64 runtime results are still open. |
| 5. One workflow entry point | The builder, CLI and R bridge use `DeseqBuilder::run(&Analysis)`. The [API migration](../design/run-entry-point.md) keeps deprecated wrappers in 0.6 and schedules removal in 0.7. |
| 6. Real-data probability gates | The [compact GTEx end-to-end gate](../reproducibility.md#compact-and-full-gene-replay) estimates normalization and dispersions with default filtering. It checks Cook's replacement, log-scale probability, nonfinite values, underflow and significance flips. The [full replay](#native-rlog-smoke-and-full-gtex-probability-replay) checks all 56,200 genes under both profiles, including log-scale error for p-values and adjusted p-values. |
| 7. Shared checked arithmetic | `rcompat-numeric` owns the checked arithmetic, matrix allocation and validation primitives. The [summation rules](../summation-rules.md) specify the reduction used for each R primitive. |
| 8. Structured errors | The non-exhaustive error enum carries domains, dimensions, indices and rejected values. The [diagnostics work](#structured-design-and-glm-diagnostics) covers the early input stages, GLM, contrasts and probability calculations. [Weight errors](#observation-weight-rank-workspace-and-diagnostics) expose both gene and sample. |
| 9. Rank and initial coefficients | `aaea914` implemented R's rank rule and `1977237` aligned rank-deficient starts. The independently generated rank reference includes scaled near-collinear designs. |
| 10. Modules and dependency cycles | Core and R bridge use modules. The module checker reports zero mutual top-level dependency pairs. It runs in Rust CI. |
| 11. Cargo features and workspace | The four Cargo packages separate numerical, formula, shrinkage and DESeq workflow code. Formula, shrinkage and CLI features are checked independently. The CLI ships from the workflow crate and the R native backend has its own package manifest. |
| 12. Release build and IRLS buffers | Release uses thin LTO, one codegen unit and unwinding. IRLS reuses worker buffers, writes outputs directly and skips unrequested covariance, hat and [diagnostic retention](#irls-likelihood-scratch-and-diagnostic-retention). The [QR work](#direct-householder-projection-for-wider-designs) records the measured range used by the direct solver. |
| 13. Test layout, fuzzing and fixtures | One integration binary shares a support module. The suite has property tests and seven fuzz targets, including parsing and transactional publication. The large GTEx tables are restored by hash. Reference breadth covers small, sparse and 500-plus-sample designs, and a previous Bioconductor release. |
| 14. CI and toolchain | Pinned workflows cover locked builds, caching, Rust 1.96.0, dependency audit, semver, coverage, mutation tests, Miri and R sanitizers. The [batched release checks](#batched-observation-weight-and-irls-release-checks) record offline builds with the minimum toolchain. |
| 15. R bridge | Native entries are guarded and CI denies undocumented Rust unsafe blocks. The offline package passed both R reference environments and SIGINT recovery. R 4.6.2 reported one compiled-code warning for the Rust runtime's `abort` import. |
| 16. CLI and R vendor | The CLI shares typed argument groups and [labelled TSV readers](#labelled-tsv-parsing). The R release builds from an offline vendor archive generated from the locked dependencies. [Contributor documentation](../../CONTRIBUTING.md) links to the checks and generator inventory. |

The [quirk register](../deseq2-quirks.md) links intentional compatibility
behavior to its pinning tests. Formula evaluation, ashr and apeglm
retain their experimental scope. The [release procedure](../release.md)
tracks the remaining hosted checks and publishing separately from local
verification.
