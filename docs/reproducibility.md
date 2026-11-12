# Reproducibility

## Reference environment

Use `deseq2-r462-native`: R 4.6.2, DESeq2 1.52.0, apeglm 1.34.0 and ashr
2.2-63. The [environment lock](../config/r-reference/r-reference-environment-2026-10-29.json)
records compiler/BLAS settings; the [source lock](../config/r-reference/r-reference-sources-2026-10-29.tsv)
pins all 82 non-base R packages by version, URL and SHA-256.

```sh
conda activate deseq2-r462-native
export R_LIBS_USER="$PWD/results/r-reference-library"
python3 scripts/environment/install_r_reference.py --rscript "$CONDA_PREFIX/bin/Rscript"
python3 scripts/environment/install_r_reference.py --check-installed \
  --library "$R_LIBS_USER" --rscript "$CONDA_PREFIX/bin/Rscript"
```

Use the selected environment's Rscript explicitly. Installation verifies source
hashes and dependency order without changing the locks. Record compiler, BLAS,
platform and thread settings for measurements on another machine. Fixtures with
an R 4.6.1 generation requirement retain that requirement; do not regenerate
those expectations under another R version and call them equivalent.

## Generate biological and random GTEx contrasts

The committed [Python generator](../scripts/references/generate_gtex_mode_contrasts.py)
uses the full GTEx v8 Kidney Cortex fixture: 56,200 genes and 61 samples.
It produces eleven designs:

- The original blocked condition contrast.
- Male versus female and male age-bin 4 versus 3. The directory
  `sex_female_vs_male` encodes male versus female despite its name.
- Five random splits preserving condition-group sizes within the existing blocks.
- Three unblocked random splits preserving the overall condition-group sizes.

The random seed is `20261014`. Design rows and count columns retain the same
sample order. Full-cohort designs reuse all gene rows; the age comparison selects
male samples. The [R helper](../scripts/references/generate_gtex_mode_contrast_references.R)
independently fits original DESeq2; it never imports candidate results.

Generate into a fresh output directory:

```sh
python3 scripts/references/generate_gtex_mode_contrasts.py \
  --rscript "$CONDA_PREFIX/bin/Rscript" \
  --output-dir results/benchmarks/reproduced-gtex
cargo build --release --locked -p rsdeseq2 --bin rsdeseq2
python3 scripts/benchmarks/benchmark_gtex_native_wald.py --exploration --repeats 1 \
  --rscript "$CONDA_PREFIX/bin/Rscript" \
  --fixture-root results/benchmarks/reproduced-gtex \
  --output-dir results/benchmarks/reproduced-comparison
```

`--counts-only` on the Python generator reproduces inputs/designs without fitting
new R references. Fresh references can differ in their final bits across recorded
environments. Exploration mode reports those differences; it does not silently
replace the accepted frozen reference set. See the [measured comparison](reports/gtex-native-wald-benchmark.md).

## Eleven-contrast acceptance

Restore the reviewed full contrast/reference bundle under
`results/benchmarks/gtex_mode_contrasts/`, or pass its location with `--fixture-root`.
The checker requires all eleven designs, result and diagnostic references, and
exact hashes from [the reference manifest](data/gtex_dispersion_reference_sha256.tsv).
A missing or changed reference fails rather than triggering regeneration.

```sh
cargo build --release --locked -p rsdeseq2
python3 scripts/validation/check_gtex_dispersion_improvement.py \
  --workers 1 --output-dir results/checks/gtex
```

This runs 121 per-contrast and nine global precision gates, including finite/NA
patterns and gene coverage. The committed gate baseline is
[data/gtex_dispersion_improvement.tsv](data/gtex_dispersion_improvement.tsv).
Keep its values unchanged during candidate validation.

Accepted timing additionally requires the reviewed native result/diagnostic
exports. Pass their directory explicitly:

```sh
python3 scripts/benchmarks/benchmark_gtex_native_wald.py --repeats 1 \
  --rscript "$CONDA_PREFIX/bin/Rscript" \
  --validated-results-root /path/to/reviewed-native-results \
  --validated-output-hashes docs/data/release_0.4.0_gtex_output_sha256.tsv \
  --output-dir results/benchmarks/accepted-comparison
```

That manifest's filename identifies the reviewed acceptance set; validation
preserves its hashes. The runner verifies native result and diagnostic exports
with untimed fits before recording result-only timing. See [benchmark scope](benchmarks.md).

## Compact and full-gene replay

All-gene counts, sample metadata, fixed dispersions and original DESeq2 outputs
are checked in under
[`tests/data/gtex_v8_real_wald/`](../crates/rsdeseq2/tests/data/gtex_v8_real_wald/README.md).
Routine tests use a 512-gene selection plus full normalization; the opt-in replay
fits all genes. Replaying requires no R installation or data download.

```sh
cargo test --locked -p rsdeseq2 --test gtex_real_wald_parity
cargo test --locked -p rsdeseq2 --test gtex_real_wald_parity full_gtex_wald \
  -- --ignored --nocapture
```

The fixture README records the source data, independent R generation version,
selection policy and regeneration commands. The independent Pasilla holdout also
has committed counts/design and Wald/LRT references; reproduce its checks with
`scripts/validation/check_pasilla_holdout.py` and its [fixture guide](../crates/rsdeseq2/tests/data/pasilla_holdout/README.md).

## Independent stage references

The DESeq2 entry script and its modules live together under
`scripts/references/` and `scripts/references/deseq2/`. Use fresh explicit output
paths to avoid replacing committed expectations:

```sh
Rscript --vanilla scripts/references/generate_deseq2_references.R \
  --output-dir results/reproduced-deseq2-reference \
  --parity-output-dir results/reproduced-deseq2-parity
```

The panels cover normalization/base metadata, supplied-dispersion weighted and
unweighted GLMs, beta priors, dispersion means/trends/MAP, Wald/LRT results,
fitted means/hats and Cook's replacement. Supplied-dispersion references isolate
GLM behavior; they do not establish accuracy for dispersion estimation by themselves.

Committed references are required: absent, empty or malformed files fail tests.
Reference generation is a reviewed maintenance task. Preserve original R/package
versions and platform metadata; never derive expected values from the candidate.
Additional formula, likelihood, posterior, mathematical truth and optimizer
reference generators are grouped in `scripts/references/`.

The bounded optimizer stress generator requires its recorded R 4.6.1 environment:

```sh
OPENBLAS_NUM_THREADS=1 \
  Rscript scripts/references/generate_lbfgsb_synthetic_stress_fixtures.R
```

A fixed hard-row scorer is also available at
`scripts/validation/score_frozen_worst_genes.py`. Its required
[gene/column set](data/wald_frozen_worst100_r461.tsv) prevents selection changes
from hiding error tails. Supply diagnostics for the complete study: missing
fixed genes fail; the scorer does not infer them from a truncated output.

## Platform references and records

[Platform arithmetic](linux_numerical_parity.md) specifies Linux extended
reductions, the pinned BLAS kernel and independent Windows/Intel-macOS CDFs.
The [testing guide](testing-parity.md) lists actual interface and safety checks.

For every new observation, record source commit, executable/package hashes,
R/DESeq2/ashr/apeglm/compiler versions, platform/BLAS, thread counts, input hashes,
settings, raw timings/RSS and result hashes. Keep normalization, gene/trend/MAP
dispersions, coefficients/SEs, Cook's and filtering intermediates so drift can be
localized before examining final p-values.
