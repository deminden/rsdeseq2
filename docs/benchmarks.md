# Benchmarks

The current benchmark suite measures workflow surfaces with apples-to-apples
reference comparisons against original DESeq2. It reports CLI-stage timings
rather than one monolithic `DESeq()` call so each validated stage has its own
speed, memory, and numeric-parity record.

## What Is Measured

The speed/RAM benchmark runner measures:

- `rsdeseq2 size-factors` versus `DESeq2::estimateSizeFactorsForMatrix()`,
- `rsdeseq2 base-mean` versus DESeq2 size factors plus normalized row means,
- elapsed wall-clock time from `/usr/bin/time -v`,
- maximum resident set size from `/usr/bin/time -v`,
- max absolute output difference between Rust and DESeq2 for each run,
- per-group medians, min/max ranges, and median absolute deviation for elapsed
  time and peak RSS.

This is a process-level benchmark. It includes CLI/R startup, TSV parsing, and
output writing. That makes it useful for end-user command behavior, but it is
not a pure inner-loop microbenchmark. Use `cargo bench -p rsdeseq2` for Rust
microbenchmarks.

## Run

```bash
scripts/benchmark_rsdeseq2.sh \
  --rscript /home/den/miniforge3/envs/rnaseq451/bin/Rscript \
  --genes 1000,10000 \
  --samples 8,16 \
  --repeats 3
```

Outputs:

- `results/benchmarks/speed_memory.tsv`
- `results/benchmarks/speed_memory_summary.tsv`

## Latest Local Run

On 2026-05-24, a three-repeat run against DESeq2 1.46.0 in the local
`rnaseq451` R environment measured the current primitive CLI paths on synthetic
matrices with 1,000 or 10,000 genes and 8 or 16 samples.

Observed medians:

- `rsdeseq2`: 0.0021-0.0091 seconds, 3.25-5.5 MiB maximum RSS.
- DESeq2/R reference process: 3.41-3.66 seconds by median elapsed time,
  661-673 MiB maximum RSS.
- Max absolute output difference versus DESeq2: at most `3.41e-12`.

The full table is in `results/benchmarks/speed_memory_current_summary.tsv`.

## Real Dataset Run

On 2026-05-27, a three-repeat run used a real publication-data
`muscle_raw_counts.tsv.gz` matrix with 73,321 genes and 818 samples. The
benchmark runner materializes compressed `.tsv.gz` count inputs into its
temporary directory before running either CLI, so both tools read the same
uncompressed count table.

Observed medians, with min-max elapsed ranges in parentheses:

| operation | tool | elapsed | max RSS | max absolute difference |
| --- | --- | ---: | ---: | ---: |
| `size-factors` | `rsdeseq2` | 3.48 s (3.46-3.49) | 237 MiB | `3.15e-14` |
| `size-factors` | DESeq2/R | 24.67 s (24.59-24.94) | 2.03 GiB | |
| `base-mean` | `rsdeseq2` | 4.07 s (3.95-4.25) | 695 MiB | `5.47e-09` |
| `base-mean` | DESeq2/R | 25.88 s (25.76-26.56) | 2.47 GiB | |

The README reports this real-data table because both operations have matching
reference outputs. The resulting primitive CLI speedups were 7.1x for size
factors and 6.4x for base means. Peak RSS was about 8.7x lower for size
factors and 3.6x lower for base means. The full table is in
`results/benchmarks/real_muscle_speed_memory_2026-05-27_summary.tsv`.

For a quick smoke run:

```bash
scripts/benchmark_rsdeseq2.sh \
  --rscript /home/den/miniforge3/envs/rnaseq451/bin/Rscript \
  --genes 1000 \
  --samples 8 \
  --repeats 1
```

To benchmark an existing real count matrix instead of synthetic counts, pass a
plain or gzip-compressed tab-delimited count table with gene IDs in the first
column:

```bash
scripts/benchmark_rsdeseq2.sh \
  --rscript /home/den/miniforge3/envs/rnaseq451/bin/Rscript \
  --counts-file /path/to/raw_counts.tsv.gz \
  --repeats 3 \
  --output results/benchmarks/real_speed_memory.tsv
```

## Interpret Carefully

Rust speedups in this benchmark should be read as primitive CLI speedups, not
as full-workflow DESeq2 speedups. DESeq2 package loading and R process startup
are included, because users pay that cost when running the reference command as
a process. The summary file reports medians, min/max ranges, and median
absolute deviations so repeated runs are less sensitive to one noisy sample
while still showing spread.

If `DESeq2` is not installed in the selected R environment, DESeq2 rows fail
clearly in the raw output rather than being substituted by any fallback.

## Real-Data Parity Sweep

On 2026-05-27, `scripts/real_data_parity.py` compared the Rust CLI with
offline DESeq2 outputs from a fresh publication-data study. The script does not
call R; it treats saved DESeq2 outputs as fixtures, derives full-tissue size
factors from the reference normalized-count matrices, and compares only outputs
where the current Rust CLI has matched inputs. The latest run used the current
release binary and completed with zero swaps for every command. The full fresh
driver run took 319.78 seconds wall time, with 710,656 KiB peak RSS and zero
swaps.

Command shape:

```bash
python3 scripts/real_data_parity.py \
  --study-root /path/to/decor_method_study \
  --binary target/release/rsdeseq2 \
  --tissue kidney \
  --tissue liver \
  --tissue pancreas \
  --tissue heart \
  --tissue muscle \
  --contrast kidney:full_blocked_permutation:1 \
  --contrast-size-factors estimate \
  --output results/benchmarks/real_data_parity_2026-05-27_fresh.tsv
```

Add `--diagnostics-output <path>` to the real-data parity script when working
on numerical parity. For contrast runs this also asks the CLI for Cook's
replacement row metadata plus original-fit and replacement-refit diagnostic
sidecars, then joins the relevant `dispGeneEst`, `dispFit`, `dispMAP`,
`dispersion`, dispersion/beta iteration and convergence fields, deviance,
Cook's summaries, Rust fallback-optimizer iteration, start/final objective,
and projected-gradient fields, beta, and beta standard-error values onto the
largest result-table differences.

Primitive parity results at floating-point parity:

| output | coverage | median elapsed | max RSS | max abs diff | max rel diff | mismatches |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `size-factors` | 5 tissues, 1,998 samples | 1.55 s | 237 MiB | `2.62e-14` | `1.99e-14` | 0 |
| `normalized-counts` | 5 tissues, 138,321,118 cells | 7.03 s | 693 MiB | `1.19e-07` | `9.74e-15` | 0 |
| `base-mean` | 5 tissues, 341,286 genes | 1.64 s | 694 MiB | `4.66e-09` | `6.73e-15` | 0 |

The same run also reconstructed a kidney Wald result for the
`condition_B_vs_A` contrast with design `~ perm_block + condition`, using
split-estimated size factors. The CLI Wald path applies Cook's
outlier replacement/refit before final Cook's masking and independent filtering,
matching the saved reference result shape closely:

| metric | mean abs diff | median abs diff | p99 abs diff | max abs diff | mismatches |
| --- | ---: | ---: | ---: | ---: | ---: |
| `baseMean` | `1.13e-12` | `8.88e-16` | `6.82e-12` | `6.52e-09` | 0 |
| `log2FoldChange` | `2.17e-08` | `3.77e-14` | `3.33e-12` | `7.70e-04` | 0 |
| `lfcSE` | `1.57e-10` | `2.33e-12` | `1.66e-10` | `8.26e-07` | 0 |
| `stat` | `3.19e-08` | `6.07e-12` | `3.44e-11` | `1.25e-03` | 0 |
| `pvalue` | `3.64e-09` | `3.03e-12` | `4.20e-11` | `6.50e-05` | 0 |
| `padj` | `2.12e-08` | `0` | `7.87e-11` | `4.50e-05` | 0 |

The current focused Wald run used
`results/benchmarks/real_data_parity_non_lbfgsb_start_probe.tsv`, covered
65,580 genes and 78 retained samples, took 151.0 s, reached 610 MiB peak RSS,
and reported zero swaps. The largest previous non-optimizer standard-error
tail was caused by pre-clamping the MAP dispersion starting value to
`maxDisp`; preserving starts above `maxDisp` and clamping only final stored
dispersions reduced `lfcSE_max_abs` from `3.27e-04` to `8.26e-07` and
`lfcSE_mean_abs` from `3.06e-08` to `1.57e-10` on this contrast. Remaining
maximum differences are now dominated by other hard rows, including beta
fallback/statistic tails.
