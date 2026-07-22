# Benchmarks

The benchmark suite compares supported CLI commands with the corresponding
DESeq2 operations. It reports CLI-stage timings, memory, and
numeric parity for each checked stage instead of treating a full `DESeq()` run
as one opaque benchmark.

## What Is Measured

The speed/RSS benchmark runner measures:

- `rsdeseq2 size-factors` versus `DESeq2::estimateSizeFactorsForMatrix()`,
- `rsdeseq2 base-mean` versus DESeq2 size factors plus normalized row means,
- elapsed wall-clock time from `/usr/bin/time -v`,
- maximum resident set size from `/usr/bin/time -v`,
- max absolute output difference between Rust and DESeq2 for each run,
- per-group medians, min/max ranges, and median absolute deviation for elapsed
  time and peak RSS.

These are process-level measurements. They include command startup, TSV
parsing, output writing, R startup, and DESeq2 package loading. Use
`cargo bench -p rsdeseq2` when you need Rust inner-loop microbenchmarks.

## Run

Synthetic moderate run:

```bash
python3 scripts/benchmark_speed_memory.py \
  --genes 10000,50000 \
  --samples 16 \
  --repeats 2 \
  --rscript /path/to/Rscript-from-an-R-4.6-DESeq2-1.52-environment \
  --output results/benchmarks/speed_memory_r461_compensated_2026-07-22.tsv
```

Real count matrix run:

```bash
scripts/benchmark_rsdeseq2.sh \
  --rscript /path/to/Rscript \
  --counts-file /path/to/raw_counts.tsv.gz \
  --repeats 3 \
  --output results/benchmarks/real_speed_memory.tsv
```

The count table must be tab-delimited with gene IDs in the first column. Plain
TSV and gzip-compressed TSV inputs are both supported.

## Moderate Process Benchmark

On 2026-07-22, a two-repeat process-level run used R 4.6.1 and DESeq2 1.52.0.
The dimensions were chosen to be large enough to include real parsing and
matrix work, but small enough to rerun during ordinary development.

Median process-level results:

| operation | genes | samples | `rsdeseq2` elapsed | DESeq2 elapsed | speedup | `rsdeseq2` RSS | DESeq2 RSS | RSS ratio | max abs diff |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `size-factors` | 10,000 | 16 | 0.010 s | 3.865 s | 386.5x | 6.01 MiB | 601.64 MiB | 100.05x | `4.885e-15` |
| `base-mean` | 10,000 | 16 | 0.030 s | 3.890 s | 129.67x | 7.06 MiB | 602.05 MiB | 85.29x | `5.684e-13` |
| `size-factors` | 50,000 | 16 | 0.080 s | 4.380 s | 54.75x | 11.63 MiB | 633.15 MiB | 54.43x | `4.663e-15` |
| `base-mean` | 50,000 | 16 | 0.140 s | 4.625 s | 33.04x | 16.54 MiB | 638.58 MiB | 38.62x | `4.547e-12` |

Elapsed values come from `/usr/bin/time`, whose output is quantized to 0.01 s,
and each median contains two runs. The speed and RSS ratios therefore describe
this process-level measurement rather than a precise inner-loop estimate.

The raw and summary files for this run are:

- `results/benchmarks/speed_memory_r461_compensated_2026-07-22.tsv`
- `results/benchmarks/speed_memory_r461_compensated_2026-07-22_summary.tsv`

## Historical Real Matrix Run

On 2026-05-27, a three-repeat run used a real GTEx muscle count matrix with
73,321 genes and 818 samples. The runner materializes compressed inputs before
running either CLI, so both tools read the same uncompressed count table.

Observed medians, with min-max elapsed ranges in parentheses:

| operation | tool | elapsed | max RSS | max absolute difference |
| --- | --- | ---: | ---: | ---: |
| `size-factors` | `rsdeseq2` | 3.48 s (3.46-3.49) | 237 MiB | `3.15e-14` |
| `size-factors` | DESeq2/R | 24.67 s (24.59-24.94) | 2.03 GiB | |
| `base-mean` | `rsdeseq2` | 4.07 s (3.95-4.25) | 695 MiB | `5.47e-09` |
| `base-mean` | DESeq2/R | 25.88 s (25.76-26.56) | 2.47 GiB | |

The primitive CLI speedups were 7.1x for size factors and 6.4x for base means.
Peak RSS was about 8.7x lower for size factors and 3.6x lower for base means.

## Publication Normalization Parity

The publication-style normalization sweep compares Rust primitive outputs with
offline DESeq2 outputs generated ahead of time from the same count matrices.
It does not call R to produce expected values during comparison.

| output | coverage | median elapsed | max RSS | max abs diff | max rel diff | mismatches |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `size-factors` | 17 tissues, 8,731 samples | 3.401 s | 237.9 MiB | `2.132e-14` | `6.501e-15` | 0 |
| `normalized-counts` | 17 tissues, 612,699,575 cells | 25.231 s | 693.9 MiB | `1.937e-7` | `9.887e-15` | 0 |
| `base-mean` | 17 tissues, 1,185,270 genes | 3.571 s | 694.7 MiB | `6.519e-9` | `8.214e-15` | 0 |

The sweep launched six tissue tasks concurrently. Its per-task elapsed values
include shared CPU and filesystem contention and are included as absolute run
measurements, not as release-to-release speed evidence. The unrounded sweep
summary is stored in
[`docs/data/normalization_r461_summary.tsv`](data/normalization_r461_summary.tsv).

To isolate the cost of compensated accumulation, two alternating sequential
heart-tissue replays used identical R 4.6.1 references and inputs:

| operation | v0.2.4 elapsed median (range) | v0.2.5 elapsed median (range) | v0.2.4 RSS median | v0.2.5 RSS median |
| --- | ---: | ---: | ---: | ---: |
| `size-factors` | 2.090 s (2.073–2.107) | 2.108 s (2.101–2.115) | 138,908 KiB | 138,928 KiB |
| `normalized-counts` | 14.477 s (14.339–14.616) | 14.317 s (14.306–14.328) | 396,044 KiB | 395,912 KiB |
| `base-mean` | 2.099 s (2.080–2.118) | 2.107 s (2.100–2.114) | 396,822 KiB | 396,814 KiB |

The absolute medians differ by +0.87%, -1.11%, and +0.39%, while peak-RSS
medians differ by +20 KiB, -132 KiB, and -8 KiB. Two runs do not establish a
material runtime or memory change. Precision did improve: v0.2.5 heart
size-factor maximum absolute and relative errors are `1.199e-14` and
`4.839e-15`; the v0.2.4 values were `2.265e-14` and `1.305e-14`.
The unrounded comparison is stored in
[`docs/data/normalization_heart_release_comparison_r461.tsv`](data/normalization_heart_release_comparison_r461.tsv).

## Wald Parity Checks

The Wald parity check compares `rsdeseq2 wald` result tables with offline
DESeq2 result tables. To reproduce a check, prepare:

- the raw count table used by DESeq2,
- the design/sample metadata used by the contrast,
- the DESeq2 size-factor mode or normalized-count source used by the contrast,
- the DESeq2 result table with `gene`, `baseMean`, `log2FoldChange`, `lfcSE`,
  `stat`, `pvalue`, and `padj` columns.

Command shape:

```bash
python3 scripts/real_data_parity.py \
  --study-root /path/to/study-inputs-and-reference-outputs \
  --binary target/release/rsdeseq2 \
  --contrast-size-factors estimate \
  --contrast tissue:null_type:replicate_id \
  --output results/benchmarks/wald_parity.tsv \
  --diagnostics-output results/benchmarks/wald_parity_diagnostics.tsv \
  --diagnostics-limit 20
```

The validation set contains publication-style full-blocked contrasts covering:

- a high-error adjusted-p-value tail,
- propagation through Cook's replacement and refitting,
- a p-value and missingness edge case.

| contrast class | rows | elapsed | peak RSS | max p-value abs diff | max adjusted-p abs diff | mismatch counters |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| heart full-blocked permutation, rep01 | 69,045 | 128.625 s | 2,878.9 MiB | `5.700e-4` | `8.748e-4` | 0 |
| heart full-blocked permutation, rep02 | 69,045 | 141.314 s | 2,905.9 MiB | `1.647e-4` | `1.494e-10` | 0 |
| testis full-blocked permutation, rep01 | 71,625 | 123.999 s | 2,822.3 MiB | `1.829e-5` | `1.146e-10` | 0 |
| pancreas full-blocked permutation, rep01 | 68,542 | 110.582 s | 2,649.7 MiB | `2.198e-4` | `2.958e-6` | 0 |

Together these checks cover 278,257 result rows with zero missing-row or
finite/NA-pattern mismatches across result columns. The remaining differences
are finite numeric drift in the fitted model path; the mismatch counter does
not impose a tolerance on finite differences.
The unrounded measurements are stored in
[`docs/data/wald_r461_contrast_summary.tsv`](data/wald_r461_contrast_summary.tsv).

DESeq2 preserves pre-optimizer IRLS hat diagonals after L-BFGS-B fallback, and
those diagonals feed Cook's replacement decisions. rsdeseq2 preserves the same
fallback hats while taking the optimizer's beta, standard error, fitted mean,
log-likelihood, and convergence diagnostics.

### L-BFGS-B precision

The `rcompat-lbfgsb` 0.2.1 dependency was compared with 0.1.6 by replaying the
same 512 bounded negative-binomial objective cases against saved R 4.6.1
`optim(..., method="L-BFGS-B")` outputs. The reference was generated on
x86_64 Linux with OpenBLAS 0.3.32, one BLAS thread, and exact 17-digit input
serialization.

| dependency | endpoint + objective within scan tolerance | objective within scan tolerance | exact endpoint + objective | exact function/gradient counts |
| --- | ---: | ---: | ---: | ---: |
| `rcompat-lbfgsb` 0.1.6 | 493/512 | 507/512 | 0/512 | 311/512 |
| `rcompat-lbfgsb` 0.2.1 | 512/512 | 512/512 | 512/512 | 512/512 |

The practical scan thresholds are maximum parameter error `5e-3`, absolute
objective error `1e-5`, or relative objective error `1e-8`. The exact run sets
all three endpoint/objective tolerances to zero. The fixture generator is
`scripts/generate_lbfgsb_synthetic_stress_fixtures.R`; it refuses R versions
other than 4.6.1. These figures characterize the recorded fixture environment and do
not imply that floating-point-sensitive callbacks are bit-exact across other
platforms or math-library configurations. The aggregate counts are stored in
[`docs/data/lbfgsb_synthetic_stress_summary.tsv`](data/lbfgsb_synthetic_stress_summary.tsv).

#### Effect on end-to-end real-data precision

A controlled replay held the rsdeseq2 source, release profile, input, design,
and saved DESeq2 reference constant and changed only
`rcompat-lbfgsb` 0.1.6 versus 0.2.1. The target was the 65,580-gene kidney
full-blocked permutation contrast with 78 samples and estimated size factors.

| result column | 0.1.6 finite difference | 0.2.1 finite difference | 0.2.1 analytic gradient |
| --- | ---: | ---: | ---: |
| `log2FoldChange` | `3.11e-14` / `3.79e-12` / `7.70e-04` | `3.15e-14` / `3.79e-12` / `7.70e-04` | `2.46e-14` / `3.03e-12` / `7.70e-04` |
| `lfcSE` | `1.74e-12` / `1.88e-10` / `1.95e-06` | `1.77e-12` / `1.88e-10` / `5.10e-06` | `1.38e-12` / `1.50e-10` / `5.10e-06` |
| `stat` | `5.27e-12` / `3.70e-11` / `1.25e-03` | `5.33e-12` / `3.70e-11` / `1.25e-03` | `4.18e-12` / `2.96e-11` / `1.25e-03` |
| `pvalue` | `2.76e-12` / `4.30e-11` / `6.50e-05` | `2.78e-12` / `4.30e-11` / `6.50e-05` | `2.37e-12` / `4.38e-11` / `6.50e-05` |
| `padj` | `0` / `8.08e-11` / `4.50e-05` | `0` / `8.06e-11` / `7.29e-05` | `0` / `8.19e-11` / `7.29e-05` |

The dependency-only change was not a material end-to-end improvement. In the
separate analytic-gradient comparison, median/p99 errors changed from
`3.15e-14` / `3.79e-12` to `2.46e-14` / `3.03e-12` for LFC, from
`1.77e-12` / `1.88e-10` to `1.38e-12` / `1.50e-10` for SE, and from
`5.33e-12` / `3.70e-11` to `4.18e-12` / `2.96e-11` for the Wald statistic
(about 22%/20% lower). Maximum errors were unchanged.
Most genes do not use the fallback, and full
Wald differences also include dispersion estimation, IRLS, covariance, Cook's
handling, and BH propagation. The measurements are stored in
[`docs/data/lbfgsb_real_data_precision.tsv`](data/lbfgsb_real_data_precision.tsv).

Only a small fraction of rows used L-BFGS-B:

| diagnostic | result |
| --- | ---: |
| Kidney rows using L-BFGS-B | 26/65,580 (`0.040%`) |
| L-BFGS-B rows with Cook's replacement/refit optimization | 9 |
| LFCs changed by 0.2.1 among those rows | 18: 10 closer to DESeq2, 8 farther |
| Rows using L-BFGS-B across eight real-data contrasts | 305/535,178 (`0.057%`) |

For eight optimizer rows in the saved fixture that did not undergo replacement, the
median/max Rust-versus-DESeq2 dispersion differences were
`4.77e-08` / `1.42e-06`, while the corresponding optimizer beta-target
differences were `9.82e-05` / `7.71e-04`. Fixed-input optimizer parity therefore
does not eliminate end-to-end drift: flat or weakly identified objectives can
amplify small upstream dispersion differences, and callback arithmetic or
backup initialization can alter the optimization trajectory.

Three whole-workflow timing runs had overlapping ranges: 0.1.6 median
`26.79 s` (range `25.87–28.09 s`) and 0.2.1 median `25.94 s` (range
`25.76–28.51 s`). The observed 3.2% median reduction is within the run-to-run
variation and does not demonstrate a runtime gain. The optimizer-use counts,
input differences, and timings are stored in
[`docs/data/lbfgsb_real_data_route_summary.tsv`](data/lbfgsb_real_data_route_summary.tsv).

A second controlled three-run comparison isolated the beta-gradient change.
Finite differences had median `28.58 s` (range `28.46–28.66 s`); the analytic
gradient had median `28.36 s` (range `28.18–28.77 s`). The 0.8% median change
and overlapping ranges do not demonstrate a material whole-workflow speedup.

#### R 4.6.1 callback arithmetic and stability check

Three trajectory-sensitive details determine compatibility on high-error rows:
R's matrix product uses fused multiply-add in the recorded reference environment,
DESeq2 passes natural-log QR backup coefficients directly to the log2
optimizer, and `optim()` obtains finite-difference gradients from R's
negative-binomial callback values. The final/refit fallback independently
evaluates the same probability identities. Gene-wise dispersion uses the
reduced analytic objective and gradient; the callback-compatible route is
limited to final/refit fitting.

Final/refit fallback runs both endpoints and applies a stability check that
does not use reference data. It substitutes the analytic endpoint only when
the analytic solver succeeds, the compatible objective exceeds the analytic
objective by more than
`10 * factr * f64::EPSILON * max(abs(compatible_objective), abs(analytic_objective), 1)`,
and at least one coefficient differs by more than ten finite-difference steps.
This preserves compatible callback trajectories while rejecting materially
divergent fits in flat nuisance directions.

#### Compensated normalization and high-error benchmark

The versioned high-error set contains the 100 highest-error v0.2.4 rows from a
69,045-gene validation contrast, with each row tied to its highest-error result
column. The v0.2.4 binary uses the v0.2.4 Rust sources, and the saved reference
was produced by R 4.6.1 with DESeq2 1.52.0. The scorer requires every frozen
gene to be present, so the diagnostics file must cover the complete contrast;
it never estimates a missing row from a truncated-table cutoff.

The input fixture is stored in
[`docs/data/wald_frozen_worst100_r461.tsv`](data/wald_frozen_worst100_r461.tsv).
Generate diagnostics with `--diagnostics-limit 69045`, then run the scorer
below. The measured v0.2.5 median absolute error is `6.064e-10`, the mean is
`1.261e-4`, and the maximum is `1.526e-3`.
The unrounded aggregate measurements are stored in
[`docs/data/wald_frozen_worst100_r461_summary.tsv`](data/wald_frozen_worst100_r461_summary.tsv).

```bash
python3 scripts/score_frozen_worst_genes.py \
  --fixture docs/data/wald_frozen_worst100_r461.tsv \
  --diagnostics results/benchmarks/frozen_worst100_diagnostics.tsv \
  --report-only
```

| 100-row measure | v0.2.4 | v0.2.5 |
| --- | ---: | ---: |
| median absolute error | `1.464e-4` | `6.064e-10` |
| mean absolute error | `3.794e-4` | `1.261e-4` |
| maximum absolute error | `3.094e-3` | `1.526e-3` |
| rows improved | — | 89/100 |
| rows improved at least 10x | — | 78/100 |

The absolute v0.2.4 measurements were `1.464e-4` median, `3.794e-4` mean, and
`3.094e-3` maximum. Using unrounded values, v0.2.5 is 241,402x lower at the
median, 3.009x lower at the mean, and 50.68% lower at the maximum. Compensated
per-gene log accumulation closely reproduces R's extended-precision row means
and prevents a boundary dispersion estimate from contaminating the fitted
trend. The maximum describes the complete set: 89/100 rows improved and 11/100
did not.

The v0.2.5 validation run measured 128.625 s and 2,947,976 KiB peak RSS. The
single v0.2.4 run measured 120.719 s and 2,945,168 KiB. The absolute differences
were +7.906 s and +2,808 KiB (+6.55% and +0.095%). These one-run values do not
establish a release-to-release runtime or memory change; in particular, they
are not evidence of a whole-workflow speed gain.

Three additional full-blocked contrasts completed in 141.314 seconds (heart),
123.999 seconds (testis), and 110.582 seconds (pancreas). Their measured maximum
Wald-statistic errors were `1.768e-3`, `5.615e-5`, and `2.803e-4`, respectively.

## Interpret Carefully

Read the speedups above as primitive CLI speedups, not full-workflow DESeq2
speedups. The DESeq2 reference is run as a separate R process because that is
the reproducible command-line comparison users can repeat. For end-to-end Wald
runtime, compare contrast-level runs with the same counts, design, size-factor
mode, replacement settings, and independent-filtering settings.

If DESeq2 is not installed in the selected R environment, DESeq2 benchmark rows
fail clearly in the raw output rather than being substituted by any fallback.
