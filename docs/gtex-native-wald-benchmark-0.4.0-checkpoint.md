# Historical GTEx native Wald comparison — 2026-10-30 candidate

Measured on 2026-10-30 using native R 4.6.2 / DESeq2 1.52.0 in
`deseq2-r462-native`, Apple M5/macOS 27.2, one Rust/BLAS worker and one paired
pass over eleven full contrasts. The tested CLI is SHA-256
`c41f9357e61ef612ec53bc6f45ae930f05e4782afb30d4710e3b2f2294376dfa`. All 56,200 genes are retained:
618,200 gene/contrast rows. Ten designs have 61 samples; the male-age design
has 52. The historical `sex_female_vs_male` directory tests male versus female.

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
| Original blocked | 3.92 | 20.53 | 5.23× |
| Male vs female | 2.63 | 15.01 | 5.71× |
| Male age 4 vs 3 | 2.29 | 14.19 | 6.20× |
| Blocked random splits | 4.26 | 26.32 | 6.17× |
| Unblocked random splits | 2.60 | 16.05 | 6.17× |

| Contrast group | Rust peak RSS (MiB) | DESeq2 peak RSS (MiB) | DESeq2 / Rust |
| --- | ---: | ---: | ---: |
| Original blocked | 301.95 | 1682.03 | 5.57× |
| Male vs female | 242.91 | 1689.48 | 6.96× |
| Male age 4 vs 3 | 211.38 | 1600.11 | 7.57× |
| Blocked random splits | 301.98 | 1735.72 | 5.75× |
| Unblocked random splits | 242.89 | 1640.25 | 6.75× |

Summed time is 40.669818 s for Rust and 242.484798 s for DESeq2 (5.9623×).
The preceding `fada553` checkpoint takes 38.339707 s;
the candidate is 6.078% slower in this pass.
Its maximum RSS is 307.671875 MiB versus the preceding checkpoint's
303.031250 MiB. No optimization round is counted.

**Performance conclusions are provisional.** During blocked random splits,
a virtual-machine process used about 66% CPU and another Rust process about
7%; machine load averages were 8.65/6.93/5.80. This task ran no builds or other
fits alongside either measurement panel. All observations are retained. For
example, blocked random split 3 takes 4.365 s before and 5.280 s after, while
DESeq2 takes 45.731 s; split 4 takes 5.230 s before and 5.851 s after. These
observations do not establish a causal regression, isolated improvement or
unchanged speed. The complete combined performance audit remains open.

Group time/RSS entries are medians across distinct contrasts. Ratios divide
those displayed group medians. Whole-process RSS includes the interpreter,
libraries and input; it does not measure only fitted arrays. Timing covers
process startup, input parsing, fitting and result TSV export. Separate untimed
fits verify diagnostic and result hashes before timing. The acceptance record
certifies frozen inputs and exact outputs; it does not remove the background-load
limitation on performance interpretation.

## R workflows

One paired pass over twenty workflows compares the checked 0.4.0 source package
with `fada553` in separate libraries. Complete raw RDS payloads match byte for
byte in all twenty cases, with no metadata exception. This covers normalization
offsets, repeated/batched queries, Wald/LRT, replacement/no-refit, VST/rlog and
normal/ashr/apeglm shrinkage. R versions and dependency packages are shared.

Total process time is 112.764132 s before and
111.325130 s after (-1.276%).
Maximum RSS is 1140.062500 MiB before and
1137.765625 MiB after; every case remains below 2 GiB.
The largest rlog takes 28.033 s before and 31.461 s after, with 741.375 and
745.938 MiB respectively. All per-case costs and gains remain in the record;
background activity also limits this panel's performance interpretation.

## Evidence

- [Validation/source identities](data/release_0.4.0_validation.json),
  [native acceptance](data/release_0.4.0_native_acceptance.json) and
  [native provenance](data/release_0.4.0_native_provenance.json).
- [Per-contrast summary](data/release_0.4.0_native_summary.tsv),
  [individual runs](data/release_0.4.0_native_timings.tsv) and
  [background-load observation](data/release_0.4.0_background_load.json).
- [Direct R 4.6.2 precision](data/release_0.4.0_direct_r462_precision.tsv),
  [frozen precision gates](data/release_0.4.0_gtex_precision.tsv),
  [global gates](data/release_0.4.0_gtex_global_precision.tsv),
  [22 export identities](data/release_0.4.0_gtex_output_sha256.tsv).
- [R workflows](data/release_0.4.0_r_workflows.tsv),
  [stage timings](data/release_0.4.0_r_stages.tsv),
  [20 exact payloads](data/release_0.4.0_r_payload_identity.tsv) and
  [R session](data/release_0.4.0_R-session.txt).
- [Independent holdout](data/release_0.4.0_pasilla_validation.json) and
  [0.4.0 release scope](release-0.4.0.md).

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
