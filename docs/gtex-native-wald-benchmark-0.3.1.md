# GTEx native Wald parity, speed and memory: 0.3.1

The patch release preserves the validated numerical behavior from `9acb44b`
and integrates twenty accepted optimisation rounds. Every contrast includes
all 56,200 GTEx v8 Kidney Cortex genes. The native pipeline estimates size
factors, gene-wise dispersions, the parametric trend, prior/MAP dispersions,
and Wald results.

## Contrasts and settings

The full sample set contains 61 samples: 52 male and nine female. The male-age
comparison uses 52 samples, with 30 in age bin 4 and 22 in age bin 3.

| Contrast | Design | Samples | Runs |
| --- | --- | ---: | ---: |
| Original blocked condition split | `~ perm_block + condition` | 61 | 1 |
| Male vs female | `~ sex` (`sexmale` tested) | 61 | 1 |
| Male age bin 4 vs 3 | `~ age_bin` | 52 | 1 |
| Stratified blocked random splits | `~ perm_block + group` | 61 | 5 |
| Balanced unblocked random splits | `~ group` | 61 | 3 |

The generator preserves the original condition assignments for the original
contrast. Blocked random splits preserve the original A/B totals within each
sex/age permutation block. Unblocked random splits assign 30 B and 31 A samples.
All random assignments use seed `20261014`. The directory
`sex_female_vs_male` retains its original name; its tested `sexmale`
coefficient is male relative to female.

References use R 4.6.1 and DESeq2 1.52.0, `fitType="parametric"`,
`betaPrior=FALSE`, and `minReplicatesForReplace=Inf`. Cook's result filtering
and independent filtering are disabled in both tools. No final dispersions
are supplied to Rust. The public R wrapper's replacement feature has separate integration coverage;
it is disabled for this matched numerical comparison.

## Precision

Errors are absolute differences from the frozen DESeq2 reference, over matching
finite genes. All contrasts retain complete gene coverage and the reference's
finite/NA pattern. Base means are checked at floating-point precision.

| Contrast | Finite LFC rows | Max LFC error | Max SE error | Max p-value error | Max adjusted p-value error |
| --- | ---: | ---: | ---: | ---: | ---: |
| Original blocked | 48,683 | `7.84e-07` | `1.38e-06` | `4.35e-07` | `1.14e-06` |
| Male vs female | 48,683 | `1.09e-09` | `1.35e-05` | `2.98e-05` | `7.31e-11` |
| Male age 4 vs 3 | 48,088 | `1.47e-07` | `2.34e-05` | `2.16e-05` | `8.25e-07` |
| Blocked random 01 | 48,683 | `4.23e-06` | `2.75e-05` | `3.74e-05` | `4.09e-06` |
| Blocked random 02 | 48,683 | `8.87e-08` | `4.67e-06` | `6.82e-06` | `6.94e-09` |
| Blocked random 03 | 48,683 | `1.31e-06` | `4.77e-06` | `9.78e-06` | `3.48e-09` |
| Blocked random 04 | 48,683 | `1.97e-06` | `7.23e-06` | `9.63e-06` | `1.2e-09` |
| Blocked random 05 | 48,683 | `7.15e-06` | `8.59e-06` | `1.25e-05` | `1.26e-05` |
| Unblocked random 01 | 48,683 | `5.28e-09` | `1.2e-05` | `1.72e-05` | `4.44e-09` |
| Unblocked random 02 | 48,683 | `1.51e-09` | `1.37e-06` | `1.66e-06` | `2.11e-06` |
| Unblocked random 03 | 48,683 | `1.25e-09` | `3.1e-06` | `2.53e-06` | `2.87e-05` |

[Machine-readable precision results](data/gtex_native_wald_parity_0.3.1.tsv)
include LFC median, p99 and maximum errors, plus maximum SE, statistic,
p-value, adjusted p-value and final-dispersion errors. The
[dispersion correction report](gtex-dispersion-improvement.md) contains the
stage comparisons, remaining dispersion tails, and all 121 per-contrast and
nine global maximum-error checks. The 0.3.1 release's 22 result/diagnostic files
match the validated 0.3.0 / `9acb44b` outputs exactly; their
[SHA-256 comparison](data/release_0.3.1_gtex_output_sha256.tsv) records that check.

## Speed

Measured on 2026-10-14 on Apple M5/macOS, with one fitting/BLAS thread and one
process run per tool and contrast. The complete suite took **41.824 s**
with Rust 0.3.1, **165.443 s** with the published Rust 0.3.0 binary, and
**227.954 s** with DESeq2 1.52.0. Ratios of these summed times are
**3.96×** for 0.3.0 / 0.3.1 and **5.45×** for DESeq2 / Rust.
The direct release comparison corresponds to **74.7% less time**.
Individual DESeq2 / Rust ratios range from **3.56× to 7.12×**.

| Contrast | Rust 0.3.1 (s) | DESeq2 (s) | DESeq2 / Rust | Rust 0.3.0 (s) | 0.3.0 / 0.3.1 |
| --- | ---: | ---: | ---: | ---: | ---: |
| Original blocked | 4.01 | 21.62 | 5.40× | 17.64 | 4.40× |
| Male vs female | 2.98 | 16.42 | 5.51× | 11.40 | 3.82× |
| Male age 4 vs 3 | 2.30 | 15.31 | 6.66× | 10.32 | 4.49× |
| Blocked random 01 | 3.93 | 21.87 | 5.56× | 17.45 | 4.43× |
| Blocked random 02 | 6.12 | 21.83 | 3.56× | 17.59 | 2.87× |
| Blocked random 03 | 3.98 | 21.91 | 5.51× | 17.81 | 4.48× |
| Blocked random 04 | 3.96 | 21.81 | 5.51× | 17.06 | 4.31× |
| Blocked random 05 | 4.75 | 21.91 | 4.61× | 17.39 | 3.66× |
| Unblocked random 01 | 3.89 | 24.65 | 6.33× | 15.57 | 4.00× |
| Unblocked random 02 | 3.38 | 24.06 | 7.12× | 11.67 | 3.45× |
| Unblocked random 03 | 2.54 | 16.58 | 6.54× | 11.56 | 4.56× |

All processes run sequentially. Tool order rotates across contrasts to balance
which tool runs first. Timings include startup, count/design parsing, full
fitting and result-table serialization. Reference comparison and separate fit-
diagnostic validation are outside these timed processes. Inputs, result flags,
thread counts and frozen references match across tools. R startup and package
loading are included. These single-run measurements describe this machine and
workflow; they do not estimate repeated-run uncertainty.

## Memory

Peak whole-process resident memory is recorded from the same timed processes.
Rust 0.3.1 uses **258.22–348.86 MiB**.
DESeq2 / Rust peak-RSS ratios are **5.20–6.46×**, equivalent to
**80.8–84.5% lower RSS for Rust**.
Relative to 0.3.0, per-contrast Rust RSS changes range from
**-0.7% to +11.2%**. The speed improvements therefore
include a modest memory increase on some contrasts.

| Contrast | Rust 0.3.1 peak RSS (MiB) | DESeq2 peak RSS (MiB) | DESeq2 / Rust | Rust 0.3.0 peak RSS (MiB) |
| --- | ---: | ---: | ---: | ---: |
| Original blocked | 338.75 | 1824.30 | 5.39× | 338.56 |
| Male vs female | 294.00 | 1796.38 | 6.11× | 270.47 |
| Male age 4 vs 3 | 258.22 | 1668.02 | 6.46× | 258.97 |
| Blocked random 01 | 335.69 | 1821.34 | 5.43× | 313.77 |
| Blocked random 02 | 342.55 | 1871.52 | 5.46× | 338.59 |
| Blocked random 03 | 335.62 | 1785.88 | 5.32× | 314.59 |
| Blocked random 04 | 335.62 | 1808.97 | 5.39× | 337.72 |
| Blocked random 05 | 348.86 | 1812.73 | 5.20× | 313.77 |
| Unblocked random 01 | 295.36 | 1701.00 | 5.76× | 275.41 |
| Unblocked random 02 | 296.08 | 1758.36 | 5.94× | 270.41 |
| Unblocked random 03 | 292.72 | 1757.48 | 6.00× | 294.81 |

RSS includes the interpreter/libraries, input matrices and native allocations
throughout the entire process lifetime. MiB = 2^20 bytes. It is the native
macOS maximum resident set size, not an R object-size or incremental allocation
estimate. Both tools use one fitting/BLAS thread.

## Evidence and optimisation history

Committed evidence:

- [Per-contrast time, memory and direct release ratios](data/gtex_native_wald_speed_memory_0.3.1.tsv)
- [All 33 raw process measurements and result hashes](data/gtex_native_wald_runs_0.3.1.tsv)
- [Native resource reports](data/gtex_native_wald_resources_0.3.1.txt)
- [Hardware, settings and binary hashes](data/gtex_native_wald_provenance_0.3.1.json)
- [Published baseline archive and binary provenance](data/gtex_native_wald_baseline_provenance_0.3.1.json)
- [R session](data/gtex_native_wald_R_session_0.3.1.txt)
- [Input/reference hashes](data/gtex_native_wald_input_sha256_0.3.1.tsv)
- [Timed parity against DESeq2](data/gtex_native_wald_timed_parity_0.3.1.tsv)
- [All 121 precision gates](data/gtex_release_0.3.1_comparison.tsv) and [nine global gates](data/gtex_release_0.3.1_global_maxima.tsv)
- [Result/diagnostic byte equality](data/release_0.3.1_gtex_output_sha256.tsv)
- [Twenty accepted optimisation rounds](data/quant_speed_rounds_0.3.1.tsv)

Each optimisation round passed a 2% total-throughput threshold in one paired
pass across eleven contrasts, with byte-identical results and diagnostics.
The product of those per-round ratios is 4.3407×. Host timing changed during
the optimisation work; that cumulative estimate has a different provenance
from the direct release ratio above. Optimisation timings also include fit-
diagnostic output. [The 0.3.0 report](gtex-native-wald-benchmark-0.3.0.md)
retains its original measurements.

## Reproduce

The [contrast generator](../scripts/generate_gtex_mode_contrasts.py) is committed.
It generates the original, biological, blocked random and unblocked random
contrasts from the documented local GTEx source with seed `20261014`.
For fixture generation and source-data paths, see the
[dispersion correction report](gtex-dispersion-improvement.md#reproduction-and-regression-coverage).
Keep frozen references unchanged for release replay:

```sh
cargo build --release -p rsdeseq2 --locked
python3 scripts/check_gtex_dispersion_improvement.py \
  --workers 1 --output-dir results/benchmarks/release_0.3.1/gtex
python3 scripts/benchmark_gtex_native_wald.py \
  --binary target/release/rsdeseq2 \
  --baseline-binary /path/to/rsdeseq2-0.3.0 --baseline-label 0.3.0 \
  --validated-results-root results/benchmarks/release_0.3.1/gtex \
  --rscript /path/to/deseq2-r461/bin/Rscript \
  --repeats 1 --rust-threads 1 --blas-threads 1 \
  --output-dir results/benchmarks/release_0.3.1/speed_memory
```

Use the platform CLI archive from the published v0.3.0 release for the baseline.
Omit `--baseline-binary` for a two-tool Rust/DESeq2 comparison. The committed
[Python benchmark](../scripts/benchmark_gtex_native_wald.py) and
[R helper](../scripts/benchmark_gtex_native_wald.R) save all raw timings, peak
RSS, result hashes, R/hardware provenance and frozen-input checksums.
