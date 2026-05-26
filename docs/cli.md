# CLI

The CLI exposes validated file-based entry points for implemented Rust stages.
It accepts primitive TSV inputs and does not parse formulas.

```bash
rsdeseq2 size-factors \
  --counts counts.tsv \
  --method ratio \
  --geometric-means geometric_means.tsv \
  --control-genes 0,2,4 \
  --output size_factors.tsv
rsdeseq2 base-mean \
  --counts counts.tsv \
  --size-factors size_factors.tsv \
  --output base_mean.tsv
rsdeseq2 normalized-counts \
  --counts counts.tsv \
  --size-factors size_factors.tsv \
  --output normalized_counts.tsv
rsdeseq2 vst \
  --counts counts.tsv \
  --design design.tsv \
  --blind=false \
  --fit-type mean \
  --output vst.tsv
rsdeseq2 wald \
  --counts counts.tsv \
  --design design.tsv \
  --size-factors size_factors.tsv \
  --observation-weights observation_weights.tsv \
  --fit-type parametric \
  --contrast-name condition_B_vs_A \
  --lfc-threshold 0.5 \
  --alternative greater \
  --use-t \
  --t-degrees-of-freedom-file wald_t_df.tsv \
  --cooks-cutoff 10 \
  --independent-filtering-alpha 0.05 \
  --independent-filtering-theta 0,0.5,1 \
  --output results.tsv
rsdeseq2 lrt \
  --counts counts.tsv \
  --design design.tsv \
  --reduced-design reduced_design.tsv \
  --normalization-factors normalization_factors.tsv \
  --fit-type parametric \
  --coefficient 1 \
  --disable-cooks-cutoff \
  --disable-independent-filtering \
  --output lrt_results.tsv
```

`counts.tsv` has a leading `gene` column followed by samples. `design.tsv` has
a leading `sample` column followed by numeric design-matrix columns, in the same
sample order as the count matrix:

```text
sample	Intercept	condition_B_vs_A
s1	1	0
s2	1	0
s3	1	1
```

`normalization_factors.tsv` is optional for `base-mean`, `normalized-counts`,
`vst`, `wald`, and `lrt`. When supplied, it preempts estimated size factors and
must have the same gene x sample shape as the count matrix:

```text
gene	s1	s2	s3
gene1	1.0	0.9	1.1
gene2	1.2	1.0	0.8
```

`size_factors.tsv` is also optional for `base-mean`, `normalized-counts`, `vst`,
`wald`, and `lrt`. It is a sample-level table with one positive finite factor per
sample:

```text
sample	size_factor
s1	1.0
s2	0.8
s3	1.2
```

Supply either `--size-factors` or `--normalization-factors`, not both.

`--control-genes` is optional for `size-factors`, `base-mean`,
`normalized-counts`, `vst`, `wald`, and `lrt`. It accepts comma-delimited
zero-based row indices and restricts size-factor estimation to those rows.
It has no effect when `--size-factors` or `--normalization-factors` supplies
the normalization directly.

`geometric_means.tsv` is optional for `size-factors`, `base-mean`,
`normalized-counts`, `vst`, `wald`, and `lrt`. It is a two-column gene/value
table used for frozen size-factor estimation:

```text
gene	geo_mean
gene1	3.2
gene2	0
gene3	8.4
```

The values must be finite and non-negative. As with `--control-genes`,
`--geometric-means` only affects estimated size factors; directly supplied
sample size factors or gene/sample normalization factors preempt estimation.

`observation_weights.tsv` is optional for `base-mean`, `vst`, `wald`, and `lrt`. It
uses the same gene x sample shape as the count matrix and accepts non-negative
finite weights:

```text
gene	s1	s2	s3
gene1	1.0	0.8	1.0
gene2	0.5	1.0	1.0
```

The `wald` and `lrt` commands use the implemented GLM-mu native dispersion,
MAP, Cook's cutoff, replacement/refit, independent-filtering, and result-table
assembly path. `wald` can report a coefficient through `--coefficient`, a
design-column name through `--contrast-name`, or a primitive numeric contrast
through `--contrast 0,1,...` in design-column order. It also accepts
thresholded p-values through
`--lfc-threshold` and `--alternative`, with alternatives `greater-abs`,
`greater-abs-upshot`, `greater-abs2014`, `less-abs`, `greater`, and `less`.
It supports Student t p-values with `--use-t` for residual degrees of freedom
or `--t-degrees-of-freedom` for one scalar value recycled over genes.
`--t-degrees-of-freedom-file` accepts a two-column gene/value TSV for per-gene
degrees of freedom.
Both commands accept `--cooks-cutoff`, `--disable-cooks-cutoff`,
`--disable-independent-filtering`, `--independent-filtering-alpha`, and
`--independent-filtering-theta` for DESeq2-style result filtering control.
Formula construction, wrapper metadata preservation, and unsupported fit types
remain outside the CLI for now. The `lrt` command compares the full design
against the supplied reduced numeric design matrix.
