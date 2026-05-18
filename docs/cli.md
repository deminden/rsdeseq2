# CLI

The CLI currently supports only implemented normalization stages.

```bash
rsdeseq2 size-factors --counts counts.tsv --method ratio --output size_factors.tsv
rsdeseq2 base-mean --counts counts.tsv --method poscounts --output base_mean.tsv
```

The CLI should not grow ahead of the validated Rust and R APIs. Full
differential-expression commands should wait until dispersion and GLM parity
exist.
