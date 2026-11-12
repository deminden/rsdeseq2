# Weighted replacement summaries

Replacement planning previously recomputed unweighted `baseMean` and `baseVar`
even when fitting used observation weights. Result assembly then overwrote
`baseMean` with those planning values. This affected weighted replacement
workflows, including enabled replacement with no refitted rows. Expanded
beta-prior refits also consumed these incorrect planning means.

Original DESeq2's `refitWithoutOutliers` calls `getBaseMeansAndVariances` after
replacement. That routine multiplies normalized counts by the original assay
weights before taking ordinary row means and sample variances. The weights are
not divided by their sum or normalized to each row's maximum for these summaries.
Trimmed replacement counts use normalized counts independently of those weights.

The new standalone `prepare_cooks_replacement_refit_with_weights` reuses the
existing weighted summary kernels. Native, fixed-dispersion, coefficient and
contrast workflows pass original builder weights through one planning helper.
All four expanded-prior replacement routes pass their original input weights.
The existing unweighted public function retains its API and arithmetic. No
optimizer, dispersion-fitting formula, threshold, tolerance or dependency changes.

## Independent checks

`generate_weighted_refit_metadata_reference.R` loads only original DESeq2 1.52.0
under R 4.6.2. Controlled Cook's values isolate replacement from optimizer
variation. Four profiles cover unweighted and raw weighted summaries, sample
size factors, matrix normalization factors and no replacement. Raw weights vary
both within and between rows, catching accidental use of GLM-normalized weights.
The fixture contains 160 rows, original replacement counts, means, variances,
zero flags and replacement counts; candidate outputs do not enter generation.

The Rust check requires exact replacement counts, zero flags and replacement
totals, with mean/variance error at most `1e-10 + 1e-12 * abs(reference)`.
Existing frozen oracles and numerical gates are unchanged. The expanded-prior
test also checks the independent four-sample replacement `41.5 -> 41` and raw
weighted mean `12.8625`, closing its earlier circular comparison against a refit
that reused the same plan.

```sh
Rscript scripts/references/generate_weighted_refit_metadata_reference.R "$FRESH_REFERENCE_JSON"
cargo test -p rsdeseq2 --test weighted_refit_metadata --test results
```

The [local report](../data/weighted-replacement-metadata-local.json) binds the
isolated numerical correction and focused checks. The subsequent fresh-fit
authority report covers combined standalone/wrapped validation separately.
This corrects a missing-weight input in main-code planning; it is not an
L-BFGS-B or locfit precision correction. Final combined release checks and
publication remain separate.
