# Benchmarks

The current benchmark suite measures only primitives that have an
apples-to-apples comparison with original DESeq2. It does not benchmark full
`DESeq()` because full end-to-end parity is not implemented yet.

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

Primitive parity results:

| output | coverage | median elapsed | max RSS | harshest max abs diff | harshest max rel diff | mismatches |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `size-factors` | 5 tissues, 1,998 samples | 1.55 s | 237 MiB | `2.62e-14` | `1.99e-14` | 0 |
| `normalized-counts` | 5 tissues, 138,321,118 cells | 7.03 s | 693 MiB | `1.19e-07` | `9.74e-15` | 0 |
| `base-mean` | 5 tissues, 341,286 genes | 1.64 s | 694 MiB | `4.66e-09` | `6.73e-15` | 0 |

The same run also reconstructed one full-blocked real contrast with
split-estimated size factors and a numeric `perm_block + condition` design.
The CLI Wald path now applies the implemented Cook's outlier replacement/refit
stage before final Cook's masking and independent filtering, matching the saved
reference result shape much more closely:

| output | contrast coverage | status |
| --- | ---: | --- |
| `wald_results` | 65,580 genes, 78 retained samples | Missingness matches the saved reference for baseMean, log2 fold change, lfcSE, Wald statistic, p-value, and adjusted p-value when size factors are estimated on the retained split. Median abs diffs are `1.04e-13` for log2 fold change, `5.70e-12` for lfcSE, `1.75e-11` for Wald statistic, `6.49e-12` for p-value, and `0` for adjusted p-value. P99 abs diffs are `1.31e-11`, `6.51e-10`, `1.29e-10`, `5.31e-11`, and `9.60e-11`, respectively. The harshest max abs diffs are `5.67e-04` for log2 fold change, `3.27e-04` for lfcSE, `9.21e-04` for Wald statistic, `4.79e-05` for p-value, and `4.50e-05` for adjusted p-value. Runtime was 128.0 s with 610 MiB peak RSS and zero swaps in the latest focused rerun. |

That contrast is now useful as a hard regression target for the next numerical
work: after aligning replacement/refit dispersion-function reuse with DESeq2,
the largest remaining real-contrast differences are split between two tails.
The largest log2-fold-change and Wald-statistic rows are non-replaced
dispersion-outlier rows that route through the pure-Rust L-BFGS-B beta
optimizer (`betaIter=100`). In the latest focused run the hardest row used 79
fallback/polish iterations and reduced its selected-coefficient objective from
`690.9302` to `481.4264`, leaving an analytic projected-gradient norm of
`1.29e-06`. The L-BFGS-B path is followed by a bounded projected-gradient polish
for rough exits; without that polish, a small number of real fallback rows
stopped with much wider log2-fold-change tails. The largest standard-error rows
converge quickly but have very large MAP dispersions. The benchmark harness
uses the split-level size-factor path that matches the saved contrast.
