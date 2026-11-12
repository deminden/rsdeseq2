# Shared fresh-fit state

Cook replacement outputs now expose
`CooksReplacementTestOutput::into_merged_fit()`. The standalone engine consumes
the original/refitted models and returns a row-aligned `DeseqFit`, final
`DeseqResults` and `CooksRefitPlan`. The R adapter calls this same method; its
private row-merging implementation has been removed.

Selected refit rows carry their coefficient covariance, dispersion, full/reduced
likelihood, fitted mean, hat and optimizer diagnostics. Original Cook distances
and maxima, normalization and global trend context remain attached. Newly
all-zero rows receive zero base summaries without inventing another fitted model.
Matrices are consumed and updated without cloning the whole fit. Missing stages
on both models are allowed; inconsistent dimensions/designs and invalid zero-row
indices return errors.

Keep the returned plan and final results with the fit. A bare `ResultRequest`
does not encode replacement context; this API does not by itself establish
replacement-aware parity for every later result query. That broader contract
remains part of C19/C21.

## Corrections found by complete fits

The independent [weighted replacement correction](weighted-replacement-metadata-correction.md)
preserves raw observation weights in base means and variances, including enabled
replacement with no refitted rows. It is committed separately from this refactor.

Unnamed matrices exposed two wrapper defects. Assigning a null first covariance
dimname through `[[1]]` removed a dimension entry; list-preserving assignment now
retains the matrix shape. Native result data frames now use compact automatic R
row names when gene names are absent. Named inputs retain their existing names.

## Local evidence

`check_fresh_fit_identity.R` runs actual public `DESeq()` calls with R numerical
functions blocked and compares their complete exported numerical state with an
independent R-free process. Floating-point values, including missing values,
are compared by their raw 64-bit representations. JSON object keys and container
types are canonicalized; numerical fields are not rounded or discarded.

The 24 cases cover mean/local/parametric fits, Wald and one/two-df LRT,
normalization factors and size factors, raw weights and failed weighted rows,
serial/four-worker execution, actual replacement, enabled-but-no-refit,
unreplaceable samples and newly all-zero rows. Three invalid inputs must reject
on both routes. A small parametric boundary case may succeed identically or
reject on both routes for a dispersion error; the larger parametric case must
succeed. Locally 20 cases succeed and four reject. The driver applies the public
R argument precedence before constructing the lower builder: normalization
factors suppress unused supplied size-factor metadata. Every exported field is
still compared. Unnamed matrix dimensions and automatic data-frame row names
are checked explicitly. Both R CI lanes include this checker.

```sh
cargo build --release --locked -p rsdeseq2 --example evaluate_r_interface
Rscript --vanilla scripts/validation/check_fresh_fit_identity.R \
  target/release/examples/evaluate_r_interface "$R_LIBRARY" "$FRESH_OUTPUT_DIR"
```

The connected Rust checks pass 196 tests. An existing replacement test also
checks merged covariance, reduced likelihood/means and optimizer objectives
against the actual selected source rows, preserves original Cook bits and
rejects a malformed plan. Focused R fit, normalization-state and native-boundary
checks pass; 25 earlier accepted standalone payload sets remain exact. The 210
pre-existing core reference files remain unchanged; the new independent weighted
fixture is additional coverage. Neither numerical gates nor optimizer controls
were relaxed.

## One paired workflow pass

All 20 complete standard workflow payloads remain byte-identical to local
`3d36759`. Summed process time is 88.2470 to 88.1812 seconds; maximum process RSS
is 1133.4063 to 1133.2344 MiB, below the 2 GiB gate. This is essentially unchanged
performance and does not count as an optimization round.

The measured library precedes the final four expanded-prior weighted planning
call-site corrections. Those paths are outside these unweighted standard cases,
but the measured and accepted source manifests remain separate: this pass does
not certify the final combined release source. No selective reruns or relaxed
payload comparator were used.

The [source/artifact report](../data/native-fresh-fit-local.json) binds the accepted
626 runtime paths, driver/library identities and measured-source limitation.
Full warning/option/state coverage, later replacement-aware queries, combined
GTEx/held-out precision, performance/RSS, sanitizer/MSRV/platform/package checks
and a final source bundle remain open. Work and commits stay local; no publication
is authorized.
