# Benchmarks

The benchmark suite measures validated command surfaces with apples-to-apples
reference comparisons against DESeq2. It reports CLI-stage timings, memory, and
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
  --output results/benchmarks/speed_memory_r460_2026-07-20.tsv
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

## Current Moderate Run

On 2026-07-20, a two-repeat process-level run used R 4.6.0 and DESeq2 1.52.0.
The dimensions were chosen to be large enough to include real parsing and
matrix work, but small enough to rerun during ordinary development.

Median process-level results:

| operation | genes | samples | `rsdeseq2` elapsed | DESeq2 elapsed | speedup | `rsdeseq2` RSS | DESeq2 RSS | RSS ratio | max abs diff |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `size-factors` | 10,000 | 16 | 0.010 s | 4.06 s | 406x | 6.0 MiB | 601.8 MiB | 100x | `5.33e-15` |
| `base-mean` | 10,000 | 16 | 0.020 s | 3.77 s | 188x | 6.9 MiB | 602.0 MiB | 87x | `5.68e-13` |
| `size-factors` | 50,000 | 16 | 0.080 s | 4.27 s | 53x | 12.0 MiB | 633.0 MiB | 53x | `5.33e-15` |
| `base-mean` | 50,000 | 16 | 0.135 s | 4.45 s | 33x | 16.6 MiB | 638.6 MiB | 38x | `4.55e-12` |

The raw and summary files for this run are:

- `results/benchmarks/speed_memory_r460_2026-07-20.tsv`
- `results/benchmarks/speed_memory_r460_2026-07-20_summary.tsv`

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
| `size-factors` | 17 tissues, 8,731 samples | 0.872 s | 238.6 MiB | `4.53e-14` | `2.64e-14` | 0 |
| `normalized-counts` | 17 tissues, 612,699,575 cells | 3.772 s | 694.4 MiB | `1.94e-07` | `9.89e-15` | 0 |
| `base-mean` | 17 tissues, 1,185,270 genes | 0.764 s | 694.9 MiB | `6.52e-09` | `8.21e-15` | 0 |

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

The post-fix check deliberately reran a small but informative set of
publication-style full-blocked contrasts instead of a fresh full sweep:

- the previous worst aggregate adjusted-p tail,
- the diagnosed Cook's replacement cascade,
- the earlier one-row p-value/NA anomaly.

| contrast class | rows | elapsed | peak RSS | max p-value abs diff | max adjusted-p abs diff | mismatch counters |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| blood full-blocked permutation | 69,817 | 1,875.14 s | 6,006.6 MiB | `4.13e-03` | `3.21e-03` | 0 |
| muscle full-blocked permutation | 70,671 | 872.53 s | 5,982.8 MiB | `2.18e-03` | `2.43e-03` | 0 |
| pancreas full-blocked permutation | 68,542 | 105.89 s | 2,677.6 MiB | `1.37e-03` | `8.32e-04` | 0 |

Together these checks cover 209,030 result rows with zero missing-row or
finite/NA-pattern mismatches across result columns. The remaining differences
are finite numeric drift in the fitted model path; the mismatch counter does
not impose a tolerance on finite differences.

The diagnosed issue was not an `rcompat-lbfgsb` optimizer bug. DESeq2 preserves
the pre-optimizer IRLS hat diagonals after L-BFGS-B fallback, and those hat
diagonals feed Cook's replacement decisions. Rust now preserves the same
pre-optimizer hats for fallback rows while still taking the optimizer's beta,
standard error, fitted mean, log-likelihood, and convergence diagnostics.

### L-BFGS-B precision refresh

The `rcompat-lbfgsb` 0.2.0 dependency was compared with 0.1.6 by replaying the
same 512 bounded negative-binomial objective cases against saved R 4.6.0
`optim(..., method="L-BFGS-B")` outputs. The oracle was generated on
x86_64 Linux with OpenBLAS 0.3.32, one BLAS thread, and exact 17-digit input
serialization.

| dependency | endpoint + objective within scan tolerance | exact endpoint + objective | exact function/gradient counts |
| --- | ---: | ---: | ---: |
| `rcompat-lbfgsb` 0.1.6 | 495/512 | 0/512 | 317/512 |
| `rcompat-lbfgsb` 0.2.0 | 512/512 | 512/512 | 512/512 |

The practical scan thresholds are maximum parameter error `5e-3`, absolute
objective error `1e-5`, or relative objective error `1e-8`. The exact run sets
all three endpoint/objective tolerances to zero. The fixture generator is
`scripts/generate_lbfgsb_synthetic_stress_fixtures.R`; it refuses R versions
other than 4.6.0. These figures characterize the pinned fixture stack and do
not imply that floating-point-sensitive callbacks are bit-exact across other
platforms or math-library configurations.

#### Effect on end-to-end real-data precision

A controlled replay held the current rsdeseq2 source, release profile, input,
design, and saved DESeq2 reference constant and changed only
`rcompat-lbfgsb` 0.1.6 versus 0.2.0. The target was the 65,580-gene kidney
full-blocked permutation contrast with 78 samples and estimated size factors.

| result column | 0.1.6 finite difference | 0.2.0 finite difference | 0.2.0 analytic gradient |
| --- | ---: | ---: | ---: |
| `log2FoldChange` | `3.11e-14` / `3.79e-12` / `7.70e-04` | `3.15e-14` / `3.79e-12` / `7.70e-04` | `2.46e-14` / `3.03e-12` / `7.70e-04` |
| `lfcSE` | `1.74e-12` / `1.88e-10` / `1.95e-06` | `1.77e-12` / `1.88e-10` / `5.10e-06` | `1.38e-12` / `1.50e-10` / `5.10e-06` |
| `stat` | `5.27e-12` / `3.70e-11` / `1.25e-03` | `5.33e-12` / `3.70e-11` / `1.25e-03` | `4.18e-12` / `2.96e-11` / `1.25e-03` |
| `pvalue` | `2.76e-12` / `4.30e-11` / `6.50e-05` | `2.78e-12` / `4.30e-11` / `6.50e-05` | `2.37e-12` / `4.38e-11` / `6.50e-05` |
| `padj` | `0` / `8.08e-11` / `4.50e-05` | `0` / `8.06e-11` / `7.29e-05` | `0` / `8.19e-11` / `7.29e-05` |

The dependency-only change is not a material end-to-end improvement. The
production analytic-gradient path is: relative to 0.2.0 finite differences,
median/p99 errors fell by 22%/20% for LFC, 22%/20% for SE, and 22%/20% for the
Wald statistic. Dominant maxima remain controlled by other pipeline stages.
Most genes do not use the fallback, and full
Wald differences also include dispersion estimation, IRLS, covariance, Cook's
handling, and BH propagation. The tracked comparison is
[`docs/data/lbfgsb_real_data_precision.tsv`](data/lbfgsb_real_data_precision.tsv).

The full per-gene diagnostics make the dilution explicit:

| diagnostic | result |
| --- | ---: |
| Kidney rows routed through L-BFGS-B | 26/65,580 (`0.040%`) |
| Routed rows with Cook's replacement/refit optimization | 9 |
| Routed LFCs changed by 0.2.0 | 18: 10 closer to DESeq2, 8 farther |
| Eight-contrast hard-real bundle routing | 305/535,178 (`0.057%`) |

For the eight fixture-covered kidney optimizer rows that were not subsequently
replaced, the median/max Rust-versus-DESeq2 dispersion differences were only
`4.77e-08` / `1.42e-06`, but the corresponding optimizer beta-target
differences were `9.82e-05` / `7.71e-04`. Thus the remaining optimizer-labeled
tail is primarily an input-parity problem: 0.2.0 reproduces R exactly when the
objective inputs are frozen, but a flat or weakly identified real objective can
amplify a tiny upstream dispersion difference into a visibly different beta.

Three whole-workflow timing runs had overlapping ranges: 0.1.6 median
`26.79 s` (range `25.87–28.09 s`) and 0.2.0 median `25.94 s` (range
`25.76–28.51 s`). The observed 3.2% median reduction is descriptive noise-scale
evidence, not a demonstrated runtime gain. The tracked route, input-drift, and
timing summary is
[`docs/data/lbfgsb_real_data_route_summary.tsv`](data/lbfgsb_real_data_route_summary.tsv).

A second controlled three-run comparison isolated the beta-gradient change.
Finite differences had median `28.58 s` (range `28.46–28.66 s`); the analytic
gradient had median `28.36 s` (range `28.18–28.77 s`). The 0.8% median change
and overlapping ranges do not demonstrate a material whole-workflow speedup.

## Interpret Carefully

Read the speedups above as primitive CLI speedups, not full-workflow DESeq2
speedups. The DESeq2 reference is run as a separate R process because that is
the reproducible command-line comparison users can repeat. For end-to-end Wald
runtime, compare contrast-level runs with the same counts, design, size-factor
mode, replacement settings, and independent-filtering settings.

If DESeq2 is not installed in the selected R environment, DESeq2 benchmark rows
fail clearly in the raw output rather than being substituted by any fallback.
