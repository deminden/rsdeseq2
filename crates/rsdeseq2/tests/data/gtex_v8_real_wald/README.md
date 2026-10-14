# GTEx v8 real-data Wald fixture

This fixture contains all 56,200 gene-count rows for 61 Kidney Cortex samples
from the GTEx v8 gene-read matrix. The sample IDs are the exact IDs shared
with the neighboring `rcompat-lbfgsb` real-data fixture. Their original
condition and permutation-block labels are retained; the numeric design is
`~ perm_block + condition`, with `conditionB` tested. The count data and sample
annotations both come from GTEx v8.

`reference_results.tsv` and `reference_size_factors.tsv` were generated with R
4.6.1 and DESeq2 1.52.0. DESeq2's parametric dispersion fit supplies one
dispersion per gene; the Rust test then runs the fixed-dispersion Wald path
with Cook's cutoff and independent filtering disabled, matching the reference
settings. Rows without finite DESeq2 Wald results are retained in the tables
and excluded from finite-value error summaries.

The checked-in fixture is directly testable without R:

```bash
cargo test -p rsdeseq2 --test gtex_real_wald_parity -- --nocapture
```

To redownload the source matrix and sample annotations, extract the same
counts, and recreate the DESeq2 references, install DESeq2 1.52.0 with R 4.6.1
and run:

```bash
python3 scripts/generate_gtex_real_wald_fixture.py \
  --rscript /path/to/Rscript
```

The generator caches the 875 MiB compressed GTEx matrix under `target/`, so it
is not committed. `source.json` records source URLs, the matrix SHA-256, sample
and gene counts, design, and reference-generation versions. The Python
generator uses only the standard library; reference generation also requires
`jsonlite` and `BiocManager`.
