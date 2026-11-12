# Benchmark procedure

The [GTEx comparison](reports/gtex-native-wald-benchmark.md) contains the measured
precision, process-time ratios and peak memory against original DESeq2 for the
eleven full Wald contrasts. The main README reports the same comparison.

## Full quantification

Use the same counts, design, normalization, test, trend and filtering/replacement
settings for both tools. Restore the pinned R environment and generate the panel
as described in [reproducibility](reproducibility.md). It includes biological
contrasts and fixed-seed blocked/unblocked random splits, retaining all genes.

The matched runner is `scripts/benchmarks/benchmark_gtex_native_wald.py`.
Accepted mode requires reviewed native result and diagnostic exports, their SHA
manifest, and unchanged input/reference hashes. It verifies exports with separate
untimed native fits before recording result-only timing. Exploration mode reports
measurements without accepting them as optimisation/release evidence.

Time and memory cover the entire process: startup, libraries, input parsing,
fitting and output serialization. Peak RSS is per process, not a sum of peaks.
R's interpreter and loaded libraries remain included. Use one fitting/BLAS worker
for the documented comparison. Record source, executable, package, compiler,
platform, settings and input/reference hashes.

A time ratio is **DESeq2 seconds / Rust seconds**. The aggregate ratio divides
the sums over the same contrasts; it does not average per-contrast ratios.
Group ratios divide group medians. RSS ratios similarly divide DESeq2 peak RSS
by Rust peak RSS. Keep units and scopes explicit.

For optimisation acceptance, use one paired pass across the full panel, alternate
execution order and keep all observations. Count a round only at 2% higher summed
throughput with unchanged numerical/output gates and acceptable memory. Combine
partial gains and record attempts in the ignored local optimisation journal.
Small changes observed under uncontrolled external load are not established
performance gains. Report diagnostic-export timing separately from result-only timing.

## Primitive and inner-loop measurements

The process runner also compares size-factor and base-mean CLI operations with
original DESeq2 operations. These primitive ratios do not describe whole-workflow
quantification speed.

```sh
python3 scripts/benchmarks/benchmark_speed_memory.py \
  --genes 10000,50000 --samples 16 --repeats 2 \
  --rscript "$CONDA_PREFIX/bin/Rscript" \
  --output results/benchmarks/primitive-comparison.tsv
scripts/benchmarks/benchmark_rsdeseq2.sh \
  --rscript "$CONDA_PREFIX/bin/Rscript" \
  --counts-file /path/to/raw_counts.tsv.gz --repeats 3 \
  --output results/benchmarks/real-primitives.tsv
cargo bench --locked -p rsdeseq2
```

Count tables have gene IDs in the first column; the real-input runner accepts
plain or gzip-compressed TSV. The runners select GNU `/usr/bin/time -v` or native macOS `time -l`
for process resource measurements.

## R workflows

`scripts/benchmarks/profile_r_workflows.py` compares complete fitted/shrinkage/
transformation payloads, process time, peak RSS and stage costs for two installed
R packages. Strict payload equality is the default. Explicit numerical migration
modes are scoped to reviewed corrections; they are not a general tolerance escape.
Use its `--help` for library and output arguments. Wrapper boundary costs and
Rust-only stages must not be compared with mismatched process scopes.

Keep raw case measurements and provenance with the current report. Do not
present a saved binary's timings as a new measurement of another source version.
