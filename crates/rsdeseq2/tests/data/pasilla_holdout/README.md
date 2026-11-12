# Independent Pasilla validation panel

This contains 2,048 genes from all seven samples in `pasilla` 1.40.0, provided by
Wolfgang Huber and Alejandro Reyes. The original experiment is Brooks et al.,
"Conservation of an RNA regulatory map between Drosophila and mammals",
Genome Research 21(2):193–202 (2011), PMID 20921232.

The source package declares `License: LGPL`. Its data are retained here as a
separate validation fixture with source attribution; this notice does not change
the Rust code's license. Original data, license metadata and processing context:
https://bioconductor.org/packages/3.23/data/experiment/src/contrib/pasilla_1.40.0.tar.gz

The archive SHA-256, expression strata and fixed selection rule are recorded in
`provenance.json`. Selection uses expression and gene IDs, never fit errors.
The model adjusts for single-read/paired-end library type and tests treatment.
The LRT removes treatment while retaining library type.

Reproduce from the repository root with the recorded source archive and the
`deseq2-r462-native` environment:

```sh
python scripts/references/generate_pasilla_validation.py --source-archive SOURCE.tar.gz \
  --output-dir crates/rsdeseq2/tests/data/pasilla_holdout
Rscript scripts/references/generate_pasilla_reference.R crates/rsdeseq2/tests/data/pasilla_holdout
```

Keep the independent references unchanged when evaluating a candidate. Use the
holdout at acceptance; do not tune individual optimisation candidates against it.
