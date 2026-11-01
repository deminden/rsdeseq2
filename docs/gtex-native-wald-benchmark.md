# GTEx native Wald comparison — 0.4.0

Final measurements use native R 4.6.2 / DESeq2 1.52.0 in
`deseq2-r462-native`, Apple M5/macOS 27.2, one Rust/BLAS worker and one paired
pass over eleven full contrasts. All 56,200 genes are retained in each contrast.
The tested CLI SHA-256 is
`ea197ba12ba7c1ae2391408f60fd1215fe8d3ab5754db9cd4a18743d9dbbea72`.
The directory named `sex_female_vs_male` tests male versus female.

## Precision

Both tools estimate normalization and gene/trend/MAP dispersions. Cook's
masking/replacement and independent filtering are disabled for the matched
native quantification; the R workflow panel checks those paths separately.
These differences compare the actual timed DESeq2 outputs with the candidate,
including matching finite/NA patterns and gene coverage.

| Contrast group | Contrasts | Max LFC error | Max SE error | Max p-value error | Max adjusted p-value error |
| --- | ---: | ---: | ---: | ---: | ---: |
| Original blocked | 1 | `7.84e-07` | `1.38e-06` | `3.73e-07` | `1.14e-06` |
| Male vs female | 1 | `1.41e-09` | `1.35e-05` | `2.98e-05` | `3.42e-11` |
| Male age 4 vs 3 | 1 | `1.47e-07` | `2.34e-05` | `2.16e-05` | `8.25e-07` |
| Blocked random splits | 5 | `7.15e-06` | `2.75e-05` | `3.74e-05` | `1.26e-05` |
| Unblocked random splits | 3 | `5.28e-09` | `1.2e-05` | `1.72e-05` | `2.87e-05` |

All 121 per-contrast and nine global frozen-reference gates pass; all 22 native
result and diagnostic exports are byte-identical to `fada553`. Earlier frozen
reference files, tolerances and release manifests remain unchanged. The
independent Pasilla Wald/LRT checks also pass with all eight exports exact.

## Observed process time and memory

| Contrast group | Rust (s) | DESeq2 (s) | DESeq2 / Rust |
| --- | ---: | ---: | ---: |
| Original blocked | 3.75 | 19.08 | 5.09× |
| Male vs female | 2.58 | 14.49 | 5.61× |
| Male age 4 vs 3 | 2.23 | 13.46 | 6.03× |
| Blocked random splits | 3.90 | 19.60 | 5.02× |
| Unblocked random splits | 2.54 | 14.71 | 5.78× |

| Contrast group | Rust peak RSS (MiB) | DESeq2 peak RSS (MiB) | DESeq2 / Rust |
| --- | ---: | ---: | ---: |
| Original blocked | 302.06 | 1763.86 | 5.84× |
| Male vs female | 240.94 | 1689.47 | 7.01× |
| Male age 4 vs 3 | 213.02 | 1600.69 | 7.51× |
| Blocked random splits | 301.84 | 1757.19 | 5.82× |
| Unblocked random splits | 240.89 | 1632.77 | 6.78× |

Summed time is **36.27 seconds for Rust versus 189.25 seconds for DESeq2**, an observed **5.22×** ratio.
Groups report medians over distinct contrasts. Peak whole-process memory includes
interpreter/libraries; no outlier was removed and no unchanged case was retimed
within the final paired pass. These are process observations, not estimates of
run-to-run variability.

The final candidate takes 36.273 s versus the previous candidate's 36.105 s
(0.465% difference), with 302.063 versus 302.031 MiB maximum RSS. Twenty paired
R workflows take 88.356 versus 88.397 s, with 1,133.750 versus 1,131.188 MiB
maximum RSS; all 20 complete payloads are byte-identical and below 2 GiB.
The main README compares memory against original DESeq2 only.

## Evidence

- [Final per-contrast summary](data/release_0.4.0_final_native_summary.tsv),
  [individual runs](data/release_0.4.0_final_native_timings.tsv),
  [acceptance](data/release_0.4.0_final_native_acceptance.json), and
  [provenance](data/release_0.4.0_final_native_provenance.json).
- [Final R workflow runs and payload hashes](data/release_0.4.0_final_r_workflows.tsv).
- [Frozen precision gates](data/release_0.4.0_gtex_precision.tsv),
  [global gates](data/release_0.4.0_gtex_global_precision.tsv), and
  [22 export identities](data/release_0.4.0_gtex_output_sha256.tsv).
- [Linux exact formula reference](linux_numerical_parity.md) and
  [release scope](release-0.4.0.md).
- [Earlier candidate observation](gtex-native-wald-benchmark-0.4.0-checkpoint.md)
  retains its original measurements and artifact identities.

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
  --validated-results-root results/benchmarks/release_0.4.0/gtex \
  --validated-output-hashes docs/data/release_0.4.0_gtex_output_sha256.tsv \
  --output-dir results/benchmarks/current-accepted
```

For historical records, see [2026-10-29](gtex-native-wald-benchmark-2026-10-29.md),
[0.3.1](gtex-native-wald-benchmark-0.3.1.md) and
[0.3.0](gtex-native-wald-benchmark-0.3.0.md).
