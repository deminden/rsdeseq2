# Automatic dispersion-mean regression

This compact fixture keeps the eight GTEx v8 Kidney Cortex genes with the
largest pre-correction gene-wise dispersion differences in the male/female
contrast, plus one all-zero gene. It uses the 61 samples from the full
`gtex_v8_real_wald` fixture. The supplied size factors were estimated over the
full 56,200-gene matrix; they must not be re-estimated from these nine rows.

The numeric design is `~ sex` with female as the reference and `sexmale` as the
tested coefficient. DESeq2 automatically uses linear fitted means because the
two distinct design rows span the two coefficients. The references also record
an explicit `linearMu=FALSE` fit and a fit with supplied unit observation
weights, which makes DESeq2 select GLM means even with this saturated design.

`reference_*_mu.tsv` contains the fitted means used by gene-wise dispersion
estimation, before final Wald fitting. These references distinguish the two
branches and cover all-zero expansion, dispersion estimates, and normalization.

Regenerate with R 4.6.1 and DESeq2 1.52.0:

```sh
Rscript scripts/generate_gtex_auto_dispersion_reference.R
cargo test -p rsdeseq2 --test auto_dispersion
```

The full data provenance and GTEx source hash are recorded in the neighboring
[`gtex_v8_real_wald`](../gtex_v8_real_wald/README.md) fixture. `sessionInfo.txt`
records the environment used for the compact reference.
