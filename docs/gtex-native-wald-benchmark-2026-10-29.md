# Historical GTEx native Wald comparison — 2026-10-29

This historical development measurement used checkpoint `807a107`, whose
R-backed changes have since been removed from the active branch. Its executable
SHA-256 is `813b8552d62760150bd23dad8c22fc6f4dbe365e87af1f585c2aa0dd2968e678`.
This is development-checkout evidence; package versions remain 0.3.1 and this work
has not created a new release.

## Input, environment and scope

All eleven contrasts contain 56,200 GTEx v8 Kidney Cortex genes. The original
blocked split, sex comparison, five stratified blocked random splits and three
balanced unblocked random splits use 61 samples. The male-age comparison uses
52 samples (30 in bin 4 and 22 in bin 3). `sex_female_vs_male` is the historical
directory name; the tested `sexmale` effect is male relative to female.

The committed generator uses seed 20261014 and preserves original condition
assignments/within-block totals. Native R 4.6.2 / DESeq2 1.52.0 is used for the
current process comparison, in `deseq2-r462-native`, with OpenBLAS/LAPACK 3.12.0,
clang 23.1.3 and R's `-O3 -mcpu=native` configuration. The Rust executable uses
Cargo's release profile. Apple M5/macOS 27.2, one fitting/BLAS worker, one process
run per tool/contrast, rotating tool order. These are distinct studies, not
repeated estimates with confidence intervals.

Both tools estimate size factors, gene dispersions, parametric trend and prior/MAP
dispersions. Neither gets final dispersions supplied. Cook's masking/replacement
and independent filtering are disabled for the matched native quantification.
The R wrapper's replacement, shrinkage and result metadata are tested separately.

## Precision

Current columns below compare directly with the new native R results. All
finite/NA masks and gene coverage match. Maxima are over each group's contrasts.

| Contrast group | Contrasts | Max LFC error | Max SE error | Max p-value error | Max adjusted p-value error |
| --- | ---: | ---: | ---: | ---: | ---: |
| Original blocked | 1 | `7.84e-07` | `1.38e-06` | `3.73e-07` | `1.14e-06` |
| Male vs female | 1 | `1.41e-09` | `1.35e-05` | `2.98e-05` | `7.31e-11` |
| Male age 4 vs 3 | 1 | `1.47e-07` | `2.34e-05` | `2.16e-05` | `8.25e-07` |
| Blocked random splits | 5 | `7.15e-06` | `2.75e-05` | `3.74e-05` | `1.26e-05` |
| Unblocked random splits | 3 | `5.28e-09` | `1.20e-05` | `1.72e-05` | `2.87e-05` |

Worst errors across all 618,200 gene/contrast rows are baseMean 5.59e-9,
dispersion 1.24e-4, LFC 7.15e-6, SE 2.75e-5, statistic 1.66e-4, p-value 3.74e-5
and adjusted p-value 2.87e-5. These are reference differences, not drift from the
previous native fit: all 22 result/diagnostic exports are byte-identical to the
reviewed v0.3.1 outputs. All 121 per-contrast and nine global dispersion gates
pass. The separate full 56,200-gene Rust oracle also passes. Historical R 4.6.1
references and their 44 input/reference hashes remain unchanged.

## Speed and memory

| Contrast group | Rust (s) | DESeq2 (s) | DESeq2 / Rust |
| --- | ---: | ---: | ---: |
| Original blocked | 3.82 | 20.24 | 5.30× |
| Male vs female | 2.61 | 14.73 | 5.63× |
| Male age 4 vs 3 | 2.25 | 13.77 | 6.11× |
| Blocked random splits | 3.98 | 20.12 | 5.05× |
| Unblocked random splits | 2.54 | 15.13 | 5.95× |

| Contrast group | Rust peak RSS (MiB) | DESeq2 peak RSS (MiB) | DESeq2 / Rust |
| --- | ---: | ---: | ---: |
| Original blocked | 301.95 | 1763.56 | 5.84× |
| Male vs female | 240.95 | 1690.38 | 7.02× |
| Male age 4 vs 3 | 211.42 | 1600.14 | 7.57× |
| Blocked random splits | 300.12 | 1757.20 | 5.85× |
| Unblocked random splits | 240.94 | 1634.42 | 6.78× |

Summed process time is 36.225 seconds Rust and 195.035 seconds DESeq2: **5.38×**.
The preserved pre-memory checkpoint `ff863e2` takes 36.279 seconds in the same
pass, a 0.15% throughput difference. This does not qualify as a 2% speed round;
it supports unchanged native throughput in this single paired pass. The final
scope includes output staging, hashing and completion-record publication.

Native peak RSS is 211.42–301.95 MiB, versus DESeq2 1,600.14–1,763.56 MiB.
Whole-process RSS includes input/interpreter/libraries; it is not just fitted
buffer storage. Groups show medians across distinct contrasts. Internal memory
and earlier R workflow changes are recorded in the [memory report](memory-improvements-2026-10-29.md).
These figures do not validate the pending [Rust-only wrapper migration](rust-only-wrapper-migration.md).

## Committed evidence

- [Per-contrast summary](data/current_gtex_summary_2026-10-29.tsv) and
  [individual timings](data/current_gtex_timings_2026-10-29.tsv).
- [Direct native-R precision](data/current_gtex_direct_parity_2026-10-29.tsv),
  [frozen-reference comparison](data/current_gtex_parity_2026-10-29.tsv),
  [22 export identities](data/current_gtex_export_identity_2026-10-29.tsv).
- [Accepted run](data/current_gtex_acceptance_2026-10-29.json),
  [binary/hardware provenance](data/current_gtex_provenance_2026-10-29.json),
  [R session](data/current_gtex_R-session_2026-10-29.txt) and
  [locked comparison environment](reproducibility.md#current-comparison-environment).

## Reproduction

The [Python contrast generator](../scripts/generate_gtex_mode_contrasts.py),
[offline reference generator](../scripts/generate_gtex_mode_contrast_references.R)
and [process benchmark](../scripts/benchmark_gtex_native_wald.py) are committed.
With Python 3.10 or newer, restore the pinned R environment and generate a separate
exploratory panel; do not overwrite historical frozen references:

```sh
export R_LIBS_USER="$PWD/results/r-reference-library"
python scripts/install_r_reference.py --rscript "$CONDA_PREFIX/bin/Rscript"
python scripts/generate_gtex_mode_contrasts.py \
  --rscript "$CONDA_PREFIX/bin/Rscript" \
  --output-dir results/benchmarks/reproduced_gtex_r462
cargo build --release --locked -p rsdeseq2
python scripts/benchmark_gtex_native_wald.py --exploration --repeats 1 \
  --rscript "$CONDA_PREFIX/bin/Rscript" \
  --fixture-root results/benchmarks/reproduced_gtex_r462 \
  --output-dir results/benchmarks/reproduced_gtex_comparison
```

Use native R 4.6.2; `CONDA_PREFIX` above points to its environment. Newly generated
references may differ in last bits and remain exploration evidence until reviewed.
Accepted optimization evidence instead requires the preserved reviewed native
result and diagnostic directory, its committed SHA manifest and unchanged frozen
references. Accepted runs verify result/diagnostic exports in separate untimed
native fits before recording the result-only timing scope. A missing or changed
diagnostic export fails acceptance; these verification fits are excluded from
the reported fit process time and RSS:

```sh
python scripts/benchmark_gtex_native_wald.py --repeats 1 \
  --rscript "$CONDA_PREFIX/bin/Rscript" \
  --validated-results-root results/benchmarks/release_0.3.1/gtex \
  --output-dir results/benchmarks/current-accepted
```

For prior published measurements, see the [0.3.1 record](gtex-native-wald-benchmark-0.3.1.md)
and [0.3.0 record](gtex-native-wald-benchmark-0.3.0.md).
