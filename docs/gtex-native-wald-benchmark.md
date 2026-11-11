# GTEx native Wald comparison — current local development

Measurements at `48917a0` use native R 4.6.2 / DESeq2 1.52.0 in
`deseq2-r462-native`, Apple M5/macOS 27.2, one Rust/BLAS worker and one paired
pass over eleven full contrasts. All 56,200 genes are retained in each contrast.
The tested CLI SHA-256 is
`0387fdccecf88c7a6bbd4375aa95f9ad8e54b877a17bce28c682d7c4cd692f11`.
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
result and diagnostic exports match the unchanged accepted 0.4.0 manifest. Earlier frozen
reference files, tolerances and release manifests remain unchanged. The
independent Pasilla Wald/LRT checks also pass with all eight exports exact.

## Observed process time and memory

| Contrast group | Rust (s) | DESeq2 (s) | DESeq2 / Rust |
| --- | ---: | ---: | ---: |
| Original blocked | 7.80 | 39.32 | 5.04× |
| Male vs female | 3.72 | 29.70 | 7.98× |
| Male age 4 vs 3 | 2.73 | 15.78 | 5.79× |
| Blocked random splits | 4.37 | 22.83 | 5.22× |
| Unblocked random splits | 2.59 | 15.16 | 5.84× |

| Contrast group | Rust peak RSS (MiB) | DESeq2 peak RSS (MiB) | DESeq2 / Rust |
| --- | ---: | ---: | ---: |
| Original blocked | 297.42 | 1764.75 | 5.93× |
| Male vs female | 241.64 | 1691.23 | 7.00× |
| Male age 4 vs 3 | 213.00 | 1579.69 | 7.42× |
| Blocked random splits | 300.47 | 1756.91 | 5.85× |
| Unblocked random splits | 242.73 | 1632.47 | 6.73× |

Summed time is **43.86 seconds for Rust versus 245.68 seconds for DESeq2**, an observed **5.60×** ratio.
Groups report medians over distinct contrasts; group ratios divide those medians.
Peak whole-process memory includes
interpreter/libraries; no outlier was removed and no unchanged case was retimed
within the final paired pass. These are process observations, not estimates of
run-to-run variability.

The current candidate takes 43.860 s versus the validated local `d64d2f2`
candidate's 45.597 s (−3.810% elapsed), with 301.938 versus 305.813 MiB maximum
RSS. This baseline is an unpublished local bundle, not the published release.
Twenty paired R workflows take 90.745 versus 91.190 s, with 1,134.563 versus
1,133.141 MiB maximum RSS, all below the predeclared 2 GiB budget. Nineteen
payloads are byte-identical; the two corrected ashr posterior columns in the
twentieth pass the independent original oracle with all other fields exact.
The R baseline is the pre-SQP `6e85f68` library. The main README compares
memory against original DESeq2 only.

Other user calculations, virtualization and system indexing were active during
the run. Load snapshots are retained with the raw logs. The measurement task
ran its panels sequentially without concurrent fits/builds; no observation was
discarded or retimed. Small development speed differences are not established
by this loaded single pass, and no optimization round is counted.

## Evidence

- [Current per-contrast summary](data/final_plan_native_summary.tsv),
  [individual runs](data/final_plan_native_timings.tsv),
  [group aggregation](data/final_plan_native_groups.tsv), and
  [source/acceptance/provenance hashes](data/final_plan_local_acceptance.json).
- [Current timed original-R precision](data/final_plan_native_precision.tsv),
  [R workflow runs and payload hashes](data/final_plan_r_workflows.tsv),
  [stage times](data/final_plan_r_stages.tsv) and
  [measured wrapper boundary costs](data/final_plan_r_boundary.tsv).
- [Frozen precision gates](data/release_0.4.0_gtex_precision.tsv),
  [global gates](data/release_0.4.0_gtex_global_precision.tsv), and
  [22 export identities](data/release_0.4.0_gtex_output_sha256.tsv).
- [Final local acceptance and limitations](final-plan-local-acceptance.md).
  The new formula/SQP independent Linux references and current platform jobs
  remain open. [Earlier Linux reference evidence](linux_numerical_parity.md)
  and [0.4.0 release scope](release-0.4.0.md) retain their historical scope.
- The preceding release comparison remains in the
  [release report](release-0.4.0.md) and its original
  [summary](data/release_0.4.0_final_native_summary.tsv),
  [runs](data/release_0.4.0_final_native_timings.tsv) and
  [provenance](data/release_0.4.0_final_native_provenance.json).
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
