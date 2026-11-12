# Numerical evidence and active validation inputs

This directory contains the records used by the current guides and validation tools:

- `ci_0.5.0.json`: passing implementation run/job identities and downloaded log hashes.
- `native_sqp_warning_native_*` and `gtex_comparison_provenance.json`: the raw
  eleven-contrast original-DESeq2 precision/time/RSS comparison and its provenance.
- `review_fixes_2026-11-12_R.tsv` and `r_workflow_provenance.json`: complete R payload
  checks, measured workflow costs and their original source/measurement scope.
- `gtex_dispersion_reference_sha256.tsv`, `gtex_dispersion_improvement.tsv` and
  `release_0.4.0_gtex_*`: unchanged input/reference, precision and export acceptance gates.
- `wald_frozen_worst100_r461.tsv`: fixed genes/columns for the hard-row scorer.

Names, generation versions and hashes identify the recorded observations and
accepted inputs. They do not claim a fresh timing run for a different binary.
These records are distinct from operational environment locks in
[`config/r-reference/`](../../config/r-reference/) and test fixtures under the crates/R package.
