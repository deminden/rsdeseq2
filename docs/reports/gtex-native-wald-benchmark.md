# GTEx native Wald comparison

The recorded comparison uses native R 4.6.2 / DESeq2 1.52.0 in
`deseq2-r462-native`, Apple M5/macOS 27.2, one Rust/BLAS worker and one paired
pass over eleven full contrasts. Each contrast retains all 56,200 genes.
The measured CLI SHA-256 is
`447ff5b2df1d2f55b94d4940c2c6507cb5063c859e4d6e3c667bae3c8982a780`.
These are measurements of that binary, not a newly timed 0.5.0 build.
The current implementation separately passes [0.5.0 CI](../testing-parity.md#ci-status).
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
result and diagnostic exports match the reviewed acceptance manifest. Frozen
reference files, tolerances and export hashes remain unchanged. The
independent Pasilla Wald/LRT checks also pass with all eight exports exact.

## Observed process time and memory

| Contrast group | Rust (s) | DESeq2 (s) | DESeq2 / Rust |
| --- | ---: | ---: | ---: |
| Original blocked | 5.96 | 31.03 | 5.20× |
| Male vs female | 2.75 | 17.69 | 6.42× |
| Male age 4 vs 3 | 2.38 | 14.54 | 6.10× |
| Blocked random splits | 4.39 | 23.71 | 5.40× |
| Unblocked random splits | 3.43 | 19.63 | 5.72× |

| Contrast group | Rust peak RSS (MiB) | DESeq2 peak RSS (MiB) | DESeq2 / Rust |
| --- | ---: | ---: | ---: |
| Original blocked | 304.62 | 1571.23 | 5.16× |
| Male vs female | 243.05 | 1689.53 | 6.95× |
| Male age 4 vs 3 | 213.61 | 1600.56 | 7.49× |
| Blocked random splits | 304.41 | 1681.58 | 5.52× |
| Unblocked random splits | 244.09 | 1643.09 | 6.73× |

Summed time is **43.33 seconds for Rust versus 242.12 seconds for DESeq2**, an observed **5.59×** ratio.
Groups report medians over distinct contrasts; group ratios divide those medians.
Peak whole-process memory includes
interpreter/libraries; no outlier was removed and no unchanged case was retimed
within the final paired pass. These are process observations, not estimates of
run-to-run variability.

Other user calculations, virtualization and system indexing were active during
the measurement. Fits/builds within the panel ran sequentially, and no run was
discarded or retimed. Small timing changes are not established gains under this
uncontrolled external load.

## Evidence

- [Per-contrast summary](../data/native_sqp_warning_native_summary.tsv),
  [individual runs](../data/native_sqp_warning_native_timings.tsv),
  [group aggregation](../data/native_sqp_warning_native_groups.tsv), and
  [source, acceptance and provenance hashes](../data/gtex_comparison_provenance.json).
- [Timed original-DESeq2 precision](../data/native_sqp_warning_native_precision.tsv).
- [Frozen precision gates](../data/release_0.4.0_gtex_precision.tsv),
  [global gates](../data/release_0.4.0_gtex_global_precision.tsv), and
  [reviewed native export hashes](../data/release_0.4.0_gtex_output_sha256.tsv).

The detailed R workflow procedure also checks complete payloads, elapsed time and
RSS; its [raw observations](../data/review_fixes_2026-11-12_R.tsv) and
[source/acceptance record](../data/r_workflow_provenance.json) retain their
own timing scope. They are not part of the Rust-versus-DESeq2 result-only ratio.

## Reproduction

The contrast generator, independent original-R fitting helper and process
benchmark runner are committed. [Reproducibility](../reproducibility.md) describes
biological/random designs, the fixed seed, fresh-reference exploration and the
frozen acceptance gate. [Benchmarks](../benchmarks.md) defines time/RSS scope,
ratios and required provenance.
