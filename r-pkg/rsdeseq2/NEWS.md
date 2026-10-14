# rsdeseq2 0.3.1

* Speed up native Wald quantification with shared Cox-Reid design products,
  matrix-order specialization, reused IRLS vectors and design decompositions,
  cached likelihood terms, and deferred error allocations.
* Preserve the validated 0.3.0 numerical behavior and automatic dispersion
  mean-fitting correction from `9acb44b`. The twenty accepted optimisation
  rounds retain byte-identical results and diagnostics on all eleven frozen
  GTEx contrasts.
* Synchronize the R backend with the optimized Rust core and align package,
  crate, backend, and vendored-core versions at 0.3.1.
* Extend the reproducible GTEx process benchmark with an optional previous
  Rust release, alternating tool order and direct time/peak-RSS ratios.

# rsdeseq2 0.3.0

* Native dispersion workflows select DESeq2's automatic mean-fitting branch.
  Unweighted designs with as many distinct sample rows as coefficients use
  linear means; other designs and supplied observation weights use GLM means.
* The validated correction from `9acb44b` reduces all nine relevant global
  maximum errors by at least 596.7-fold across 11 frozen GTEx Kidney Cortex
  contrasts, with preserved gene coverage and finite/NA patterns.
* Fitted `DESeqDataSet` objects receive Rust `mu` and `cooks` assays, fitted
  row metadata, and per-coefficient MLE and standard-error columns.
* Cook's outlier replacement and refitting run for eligible replicate samples.
  `minReplicatesForReplace` defaults to seven; use an integer of at least three
  or `Inf` to disable replacement. Original counts and Cook's distances are
  retained, with replacement assays and flags attached when applicable.
* Dataset result tables carry Bioconductor column descriptions and result
  metadata, including independent-filtering and replacement/refit information.
* Seeded integration tests cover automatic dispersion means, fitted means,
  Cook's distances, outlier replacement/refitting, and result metadata.

# rsdeseq2 0.2.6

* Added dataset VST/rlog transformations, frozen trend/intercept reuse, and
  selected normal, `ashr`, and `apeglm` shrinkage workflows.
