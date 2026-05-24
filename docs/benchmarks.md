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

On 2026-05-24, a five-repeat run used a real `decor_method_study`
`muscle_raw_counts.tsv` matrix with 56,937 genes and 881 samples. The source
file stores some integer counts in scientific notation, so the CLI count reader
accepts integer-valued numeric fields such as `1e+05`.

Observed medians, with min-max elapsed ranges in parentheses:

| operation | tool | elapsed | max RSS | max absolute difference |
| --- | --- | ---: | ---: | ---: |
| `size-factors` | `rsdeseq2` | 1.15 s (1.14-1.28) | 199 MiB | `3.86e-14` |
| `size-factors` | DESeq2/R | 26.71 s (24.87-27.32) | 1.90 GiB | |
| `base-mean` | `rsdeseq2` | 1.38 s (1.33-1.44) | 581 MiB | `4.47e-07` |
| `base-mean` | DESeq2/R | 27.55 s (25.58-28.59) | 2.28 GiB | |

The README reports this real-data table because both operations have matching
reference outputs. The resulting primitive CLI speedups were 23.2x for size
factors and 20.0x for base means. Peak RSS was about 9.8x lower for size
factors and 4.0x lower for base means. The full table is in
`results/benchmarks/real_muscle_speed_memory_current_summary.tsv`.

For a quick smoke run:

```bash
scripts/benchmark_rsdeseq2.sh \
  --rscript /home/den/miniforge3/envs/rnaseq451/bin/Rscript \
  --genes 1000 \
  --samples 8 \
  --repeats 1
```

To benchmark an existing real count matrix instead of synthetic counts, pass a
tab-delimited count table with gene IDs in the first column:

```bash
scripts/benchmark_rsdeseq2.sh \
  --rscript /home/den/miniforge3/envs/rnaseq451/bin/Rscript \
  --counts-file /path/to/raw_counts.tsv \
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
