# GTEx native Wald parity and speed: 0.3.0

This comparison uses all 56,200 GTEx v8 Kidney Cortex genes in each of eleven
contrasts. The native workflow estimates size factors, gene-wise dispersions,
the parametric trend, prior/MAP dispersions, and Wald results. The numerical
correction is the validated implementation from `9acb44b`.

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
are supplied to Rust. The public R wrapper's replacement feature is tested
separately; it is disabled for this matched numerical comparison.

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
| UnBlocked random 01 | 48,683 | `5.28e-09` | `1.2e-05` | `1.72e-05` | `4.44e-09` |
| UnBlocked random 02 | 48,683 | `1.51e-09` | `1.37e-06` | `1.66e-06` | `2.11e-06` |
| UnBlocked random 03 | 48,683 | `1.25e-09` | `3.1e-06` | `2.53e-06` | `2.87e-05` |

[Machine-readable precision results](data/gtex_native_wald_parity_0.3.0.tsv)
include LFC median, p99 and maximum errors, plus maximum SE, statistic,
p-value, adjusted p-value and final-dispersion errors. The
[dispersion correction report](gtex-dispersion-improvement.md) contains the
stage comparisons, remaining dispersion tails, and all 121 per-contrast and
nine global maximum-error checks. The release's 22 result/diagnostic files
match the saved `9acb44b` outputs exactly; their
[SHA-256 comparison](data/release_0.3.0_gtex_output_sha256.tsv) records that check.

## Speed

Measured on 2026-10-14. Rust is **1.24–1.52× faster** in these single-process
runs. A ratio above one means less elapsed time for Rust.

| Contrast | Rust (s) | DESeq2 (s) | Time ratio: DESeq2 / Rust |
| --- | ---: | ---: | ---: |
| Original blocked | 16.84 | 21.55 | 1.28× |
| Male vs female | 11.13 | 16.46 | 1.48× |
| Male age 4 vs 3 | 10.09 | 15.36 | 1.52× |
| Blocked random 01 | 17.42 | 21.90 | 1.26× |
| Blocked random 02 | 17.16 | 21.22 | 1.24× |
| Blocked random 03 | 17.05 | 21.51 | 1.26× |
| Blocked random 04 | 16.99 | 21.59 | 1.27× |
| Blocked random 05 | 17.14 | 21.64 | 1.26× |
| UnBlocked random 01 | 11.53 | 16.41 | 1.42× |
| UnBlocked random 02 | 11.68 | 16.33 | 1.40× |
| UnBlocked random 03 | 11.46 | 16.62 | 1.45× |

The recorded run uses an Apple M5 on macOS 27.0.1, with Rust fitting limited to
one Rayon thread and DESeq2 running without BiocParallel, with one BLAS thread.
Fits run sequentially, with no simultaneous benchmark fits. Each table entry
is one process run per tool and contrast; these are single-run timings, not
replicate medians or confidence intervals. The runner supports additional
repetitions and reports their medians and ranges.

Wall-clock time includes process startup, input parsing, full native Wald
fitting, and result TSV serialization for both tools. R package loading is
included. Compilation, contrast preparation, reference generation and the
post-run parity comparison are excluded. Speedup is DESeq2 seconds divided by
Rust seconds. This is a comparison of these CLI/R-script workflows on this
machine, rather than an inner-loop or universal speed claim.

[Per-contrast timing data](data/gtex_native_wald_speed_0.3.0.tsv),
[raw process runs](data/gtex_native_wald_speed_runs_0.3.0.tsv),
[hardware/thread provenance](data/gtex_native_wald_speed_provenance_0.3.0.json),
[R session](data/gtex_native_wald_speed_R_session_0.3.0.txt), and
[input/reference checksums](data/gtex_native_wald_speed_input_sha256_0.3.0.tsv)
are committed. The [timed-run parity table](data/gtex_native_wald_timed_parity_0.3.0.tsv)
records both tools against the frozen reference. Rust's timed result files
match the validated outputs exactly. R under the explicit one-thread BLAS
profile differs from the original frozen R run by at most `1.72e-9` in LFC
and `9.00e-7` in p-value; the
precision table above consistently uses the frozen oracle.

## Reproduce the contrasts and references

The generator and reference helper are committed in the repository:

* [`generate_gtex_mode_contrasts.py`](../scripts/generate_gtex_mode_contrasts.py)
  creates the biological designs and deterministic random splits from the
  checked-in complete Kidney Cortex count matrix and sample metadata.
* [`generate_gtex_mode_contrast_references.R`](../scripts/generate_gtex_mode_contrast_references.R)
  generates the DESeq2 results, dispersion diagnostics and size factors offline.

From the repository root, with the pinned R/DESeq2 environment available:

```sh
python3 scripts/generate_gtex_mode_contrasts.py \
  --rscript /path/to/Rscript
cargo build --release -p rsdeseq2 --locked
python3 scripts/check_gtex_dispersion_improvement.py \
  --output-dir results/benchmarks/release_0.3.0/gtex
```

Use `--counts-only` to generate designs and counts without calling R. The
generator writes to ignored `results/benchmarks/gtex_mode_contrasts` by
default. To prepare a separate reference bundle, supply `--output-dir PATH`
and use that path as `--fixture-root` for the checker and benchmark. Full
source URLs, matrix checksum and extraction instructions are in the
[checked-in fixture notes](../crates/rsdeseq2/tests/data/gtex_v8_real_wald/README.md).

The checker validates frozen input/reference hashes as well as numerical
errors. Recreating a reference with a different R/BLAS build can change its
last floating-point digits and stopping decisions; inspect those differences
before treating that environment as the recorded oracle.

## Reproduce the speed comparison

```sh
python3 scripts/benchmark_gtex_native_wald.py \
  --rscript /path/to/Rscript \
  --rust-threads 1 --blas-threads 1 --repeats 1 \
  --output-dir results/benchmarks/release_0.3.0/speed
```

[`benchmark_gtex_native_wald.py`](../scripts/benchmark_gtex_native_wald.py) and
its [R helper](../scripts/benchmark_gtex_native_wald.R) are committed. They
write raw timings, per-contrast summaries, numeric comparisons to the frozen
reference, R session information, thread settings, hardware details and input,
reference and binary hashes. The benchmark leaves all frozen references
unchanged. Use `--repeats 3` for repeated measurements; odd/even repetitions
alternate tool order. `--contrasts NAME ...` selects a named subset.
