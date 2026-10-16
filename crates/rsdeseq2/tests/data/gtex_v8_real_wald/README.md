# GTEx v8 real-data Wald fixture

The source fixture retains all 56,200 gene-count rows and 61 Kidney Cortex
samples from the public GTEx v8 gene-read matrix. The numeric design is
`~ perm_block + condition`, with `conditionB` tested. Original sample IDs,
condition labels and permutation blocks are retained.

The full DESeq2 references were generated with R 4.6.1 and DESeq2 1.52.0,
with parametric dispersions, no beta prior, no count replacement, Cook's
cutoff disabled and independent filtering disabled. `source.json` records
versions, source URLs, the upstream matrix hash and dimensions.

## Routine checks

```bash
cargo test -p rsdeseq2 --test gtex_real_wald_parity
```

The routine Wald fit uses the **512-gene × 61-sample** panel in `compact/`.
Selection covers tails and quantiles of base mean, dispersion, absolute effect,
standard error, test statistic, p-value, maximum count and nonzero sample count;
every reference missingness pattern is represented. Additional background genes
are selected deterministically. `compact/selection.json` records every gene,
selection reasons, original row number and hashes of the four full inputs.

Subsetting changes estimated normalization and the multiple-testing family.
The panel therefore **supplies the original full-cohort size factors and final
dispersions**, and freezes new subset adjusted p-values from R's
`stats::p.adjust(method="BH")`. It does not estimate these from the small panel.
A separate routine check still estimates size factors from all 56,200 genes.

Both replay sizes check all six result fields, exact gene/sample alignment,
finite values and matching missingness. Absolute error bounds are unchanged:
baseMean/lfcSE/pvalue/padj < 1e-6 and log2FoldChange/stat < 1e-5;
size-factor relative error < 1e-10. The original full reference has no `padj`
column, so its BH oracle is calculated independently from the R p-values.

## Full replay

```bash
cargo test -p rsdeseq2 --test gtex_real_wald_parity full_gtex_wald \
  -- --ignored --nocapture
```

This explicit slow test fits every original gene and estimates normalization.
The manual Parity workflow runs it. Full data are also retained for the eleven
biological/random contrast benchmarks; their generation is documented in
`docs/gtex-native-wald-benchmark.md`. The compact fixed-dispersion test does not
replace native dispersion, optimizer, weighted, local-fit or full-contrast gates.

## Reproduce the panel and optimizer oracles

From the repository root, using the pinned reference environment:

```bash
python3 scripts/generate_gtex_compact_fixture.py --rscript /path/to/Rscript
/path/to/Rscript --vanilla scripts/generate_gtex_optim_reference.R "$PWD"
```

The optimizer generator freezes six independent `stats::optim` cases: three
GTEx dispersion/effect/expression extremes in both objective-only compatibility
and analytic-gradient modes. Counts, design, normalization, starting values,
bounds, controls and R answers are checked in under
`crates/rsdeseq2/src/glm/irls/fixtures/`, so no neighboring repository is required.

To redownload upstream data and recreate the full references:

```bash
python3 scripts/generate_gtex_real_wald_fixture.py --rscript /path/to/Rscript
```

This caches the 875 MiB compressed upstream matrix under `target/`. Reference
generation requires DESeq2 1.52.0, jsonlite and BiocManager; the Python generators
use only the standard library. Rebuild the compact panel after regenerating the
full references. Never regenerate references to conceal a candidate regression.
