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
rsdeseq2 rlog \
  --counts counts.tsv \
  --design design.tsv \
  --blind=false \
  --fit-type mean \
  --output rlog.tsv
rsdeseq2 wald \
  --counts counts.tsv \
  --design design.tsv \
  --size-factors size_factors.tsv \
  --observation-weights observation_weights.tsv \
  --fit-type parametric \
  --contrast-positive condition_B_vs_A \
  --contrast-negative Intercept \
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
a leading `sample` column followed by numeric design-matrix columns. For
commands that take both files, design rows are aligned by sample label against
the count-matrix columns:

```text
sample	Intercept	condition_B_vs_A
s1	1	0
s2	1	0
s3	1	1
```

`normalization_factors.tsv` is optional for `base-mean`, `normalized-counts`,
`vst`, `rlog`, `wald`, and `lrt`. When supplied, it preempts estimated size
factors and must have the same gene x sample shape as the count matrix. Rows
and columns are aligned by gene and sample labels:

```text
gene	s1	s2	s3
gene1	1.0	0.9	1.1
gene2	1.2	1.0	0.8
```

`size_factors.tsv` is also optional for `base-mean`, `normalized-counts`, `vst`,
`rlog`, `wald`, and `lrt`. It is a sample-level table with one positive finite
factor per sample. Rows are aligned by sample label against the count-matrix
columns:

```text
sample	size_factor
s1	1.0
s2	0.8
s3	1.2
```

Supply either `--size-factors` or `--normalization-factors`, not both.

`--control-genes` is optional for `size-factors`, `base-mean`,
`normalized-counts`, `vst`, `rlog`, `wald`, and `lrt`. It accepts comma-delimited
zero-based row indices and restricts size-factor estimation to those rows.
It has no effect when `--size-factors` or `--normalization-factors` supplies
the normalization directly.

`geometric_means.tsv` is optional for `size-factors`, `base-mean`,
`normalized-counts`, `vst`, `rlog`, `wald`, and `lrt`. It is a two-column gene/value
table used for frozen size-factor estimation. Rows are aligned by gene label
against the count-matrix rows:

```text
gene	geo_mean
gene1	3.2
gene2	0
gene3	8.4
```

The values must be finite and non-negative. As with `--control-genes`,
`--geometric-means` only affects estimated size factors; directly supplied
sample size factors or gene/sample normalization factors preempt estimation.

`observation_weights.tsv` is optional for `base-mean`, `vst`, `rlog`, `wald`,
and `lrt`. It uses the same gene x sample shape as the count matrix and accepts
non-negative finite weights. Rows and columns are aligned by gene and sample
labels:

```text
gene	s1	s2	s3
gene1	1.0	0.8	1.0
gene2	0.5	1.0	1.0
```

The `wald` and `lrt` commands use the implemented GLM-mu native dispersion,
MAP, Cook's cutoff, replacement/refit, independent-filtering, and result-table
assembly path. `wald` and `lrt` can report a design coefficient by zero-based
`--coefficient` or by `--coefficient-name`. They can also report a
design-column contrast through `--contrast-name`, a positive/negative coefficient
list through `--contrast-positive` and `--contrast-negative`, a factor-level
contrast through `--contrast-factor`, `--contrast-numerator`, and
`--contrast-denominator`, or a primitive numeric contrast through
`--contrast 0,1,...` in design-column order. For LRT, contrast flags only
change the displayed effect-size columns; the statistic and p-values remain
the full-vs-reduced likelihood-ratio test. List contrasts can use
`--contrast-positive-weight` and `--contrast-negative-weight` to override the
default `1` and `-1` weights; matching DESeq2 `listValues`, the positive
weight must be greater than zero and the negative weight must be less than
zero. Factor-level contrasts resolve against existing
design coefficient names, with optional `--contrast-reference`; common
non-reference comparisons can also infer a shared reference from coefficient
names such as `condition_B_vs_A` and `condition_C_vs_A`. Supplying
`--contrast-sample-levels` as a two-column sample/level TSV additionally
enables DESeq2-style factor-level all-zero contrast handling and must be paired
with a factor-level contrast request; sample rows are aligned by label against
the count-matrix columns. For LRT result tables, this cleanup zeroes only the
displayed log2 fold change and keeps the
full-vs-reduced statistic and p-values. The CLI still does not parse formulas.
It also accepts thresholded p-values through
`--lfc-threshold` and `--alternative`, with alternatives `greater-abs`,
`greater-abs-upshot`, `greater-abs2014`, `less-abs`, `greater`, and `less`.
It supports Student t p-values with `--use-t` for residual degrees of freedom
or `--t-degrees-of-freedom` for one scalar value recycled over genes.
`--t-degrees-of-freedom-file` accepts a two-column gene/value TSV for per-gene
degrees of freedom. Rows are aligned by gene label against the count-matrix
rows.
Both commands accept `--cooks-cutoff`, `--disable-cooks-cutoff`,
`--disable-independent-filtering`, `--independent-filtering-alpha`, and
`--independent-filtering-theta` for DESeq2-style result filtering control.
Result metadata sidecars are available with `--result-column-metadata-output`
and `--result-table-metadata-output`. Independent-filtering sidecars are
available with `--independent-filter-metadata-output`,
`--independent-filter-num-rej-output`, and
`--independent-filter-lowess-output`; these require independent filtering to
be enabled.
They can also write optional Cook's sidecar tables: `--cooks-distance-output`
for the Cook's distance matrix, `--cooks-replacement-metadata-output` for
replacement/refit scalar metadata, `--cooks-replacement-row-metadata-output`
for row-level replacement/refit metadata, `--cooks-replaced-counts-output` for
replacement counts, `--cooks-candidate-replacement-counts-output` for
candidate replacement counts, and `--cooks-outlier-cells-output` for the
logical outlier-cell assay. Replacement sidecars require Cook's cutoff to be
enabled so the replacement/refit branch runs.
Formula construction, wrapper metadata preservation, and unsupported fit types
remain outside the CLI for now. The `lrt` command compares the full design
against the supplied reduced numeric design matrix.

The `rlog` command fits the implemented GLM-mu dispersion/MAP stages, estimates
the rlog sample-effect prior from normalized counts, `baseMean`, and `dispFit`,
then writes the transformed gene x sample matrix. It defaults to `--blind=true`
with an intercept-only design; use `--blind=false --design design.tsv` for a
design-aware dispersion workflow. Supplying `--frozen-intercept` with a
two-column gene/value TSV and `--rlog-prior-variance` runs the frozen-intercept
rlog transform after fitting the dispersion state, aligning intercept rows by
gene label.
