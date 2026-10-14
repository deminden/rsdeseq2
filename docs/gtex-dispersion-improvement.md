# GTEx dispersion parity correction

The target was a 100-fold reduction in maximum absolute errors against the
same 11 GTEx contrasts. All nine relevant global maxima fell by at least
596.7-fold. Within the five affected unblocked contrasts, every targeted
per-contrast maximum fell by at least 664-fold. Six blocked contrasts retained
their previous error limits, and every numeric result field kept the
reference's finite/NA pattern. The complete comparison contains 121 passing
per-contrast checks and nine passing global-maximum checks over 11
56,200-gene contrasts.

## Cause and correction

DESeq2 automatically uses linear fitted means when an unweighted design has
as many distinct sample rows as coefficients. Two-group models satisfy this
condition. The previous CLI always used GLM fitted means for dispersion
estimation. This changed gene-wise dispersion estimates, their parametric
trend, the prior variance, and MAP dispersions.

The automatic estimator now chooses the same branch as DESeq2. It uses the
existing dispersion likelihood, Cox-Reid correction, likelihood-improvement
check, and fallback grid. The explicit GLM estimator remains available for
forced `linearMu=FALSE` calculations. The automatic path is used by the CLI,
top-level Rust workflows, and R backend. Supplied observation weights select
GLM means, including uniform weights.

No L-BFGS-B, locfit, trend, or MAP optimizer changes were needed to achieve the
measured reduction. Correcting the upstream fitted means removes the broad
trend/MAP discrepancies. In all five affected contrasts, both the trend's
eligible gene set and its final robust residual screen match the reference;
the MAP grid-fallback and dispersion-outlier sets also match. The fitted trend
coefficients differ by at most `9.9e-10`, and reconstructed robust prior
variances differ by at most `1.4e-9`.

## Maximum errors

These are the maximum absolute differences over all 11 contrasts, including
the blocked controls. The largest before and after differences can belong to
different genes or contrasts. Several remaining global maxima come from the
blocked controls, whose fitting branch already matched DESeq2.

| Field | Before maximum | After maximum | Global reduction |
| --- | ---: | ---: | ---: |
| Gene-wise dispersion | `5.20e1` | `1.24e-4` | `4.20e5` times |
| Fitted dispersion trend | `2.41` | `7.44e-7` | `3.25e6` times |
| MAP dispersion | `2.85` | `1.06e-4` | `2.69e4` times |
| Final dispersion | `2.85` | `1.24e-4` | `2.30e4` times |
| Log2 fold change | `1.80e-2` | `7.15e-6` | `2.52e3` times |
| LFC standard error | `1.17e-1` | `2.75e-5` | `4.25e3` times |
| Wald statistic | `9.89e-2` | `1.66e-4` | `5.97e2` times |
| P-value | `3.01e-2` | `3.74e-5` | `8.03e2` times |
| Adjusted p-value | `5.91e-2` | `2.87e-5` | `2.06e3` times |

The five affected unblocked models improve more strongly for LFC: their
maximum falls from `1.80e-2` to `1.47e-7`. The smallest per-contrast reduction
across the targeted fields is 664-fold. The global maximum values are saved in
[data/gtex_dispersion_global_maxima.tsv](data/gtex_dispersion_global_maxima.tsv).

Base means already agreed at floating-point precision and were checked for
preservation rather than a 100-fold reduction. The remaining errors are finite
numeric tails; this evidence covers the specified Kidney Cortex data and
models, with Cook's filtering/replacement and independent filtering disabled.

The full per-contrast values are in
[data/gtex_dispersion_improvement.tsv](data/gtex_dispersion_improvement.tsv).
[data/gtex_dispersion_stage_routes.tsv](data/gtex_dispersion_stage_routes.tsv)
records trend coefficients, gene selection, derived prior variances, and
fallback/outlier routes. Input and reference hashes are recorded in
[data/gtex_dispersion_reference_sha256.tsv](data/gtex_dispersion_reference_sha256.tsv).

## Reproduction and regression coverage

The 11 models include the original blocked split, five new stratified random
splits, three balanced unblocked random splits, a male/female comparison, and
age bin 4 versus 3 among males. The seed is `20261014`. Reference outputs use R
4.6.1 and DESeq2 1.52.0 with parametric dispersions, no beta prior, and no
replacement. These commands use the complete checked-in Kidney Cortex count
fixture and generate references offline:

```sh
python3 scripts/generate_gtex_mode_contrasts.py --rscript /path/to/Rscript
cargo build --release -p rsdeseq2
python3 scripts/check_gtex_dispersion_improvement.py
```

The checker uses the saved baseline maxima and frozen reference hashes. It
checks a 100-fold improvement in every affected maximum, controls for blocked
fits and base means, and rejects gene coverage or finite/NA mismatches. The
recorded hashes require the original numeric references; differences in R/BLAS
arithmetic on another environment may produce different reference hashes.

The compact checked-in regression preserves eight formerly high-error genes
and an all-zero row, using full-dataset size factors. It independently checks
DESeq2's automatic, forced GLM, and weighted mean-fitting branches, including
fitted means before dispersion optimization. Top-level and normalization-factor
checks ensure the automatic choice reaches the public APIs. Runtime Rust
calculation does not call DESeq2.

```sh
cargo test -p rsdeseq2 --test auto_dispersion
cargo test --workspace
```
