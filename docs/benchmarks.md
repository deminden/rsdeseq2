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
- max absolute output difference between Rust and DESeq2 for each run.

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

On 2026-05-23, a three-repeat run against DESeq2 1.46.0 in the local
`rnaseq451` R environment measured the current primitive CLI paths on synthetic
matrices with 1,000 or 10,000 genes and 8 or 16 samples.

Observed medians:

- `rsdeseq2`: 0.0075-0.02 seconds, 3.3-5.6 MiB maximum RSS.
- DESeq2/R reference process: 3.54-3.79 seconds, 676-689 MiB maximum RSS.
- Max absolute output difference versus DESeq2: at most `3.5e-12`.

The full table is in `results/benchmarks/speed_memory_summary.tsv`.

## Real Dataset Run

On 2026-05-23, a three-repeat run used a real `decor_method_study`
`muscle_raw_counts.tsv` matrix with 56,937 genes and 881 samples. The source
file stores some integer counts in scientific notation, so the CLI count reader
accepts integer-valued numeric fields such as `1e+05`.

Observed medians:

| operation | tool | elapsed | max RSS | max absolute difference |
| --- | --- | ---: | ---: | ---: |
| `size-factors` | `rsdeseq2` | 1.26 s | 199 MiB | `3.86e-14` |
| `size-factors` | DESeq2/R | 26.93 s | 1.90 GiB | |
| `base-mean` | `rsdeseq2` | 1.63 s | 581 MiB | `4.47e-07` |
| `base-mean` | DESeq2/R | 27.21 s | 2.28 GiB | |

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
a process. The summary file reports medians so repeated runs are less sensitive
to one noisy sample.

If `DESeq2` is not installed in the selected R environment, DESeq2 rows fail
clearly in the raw output rather than being substituted by any fallback.
