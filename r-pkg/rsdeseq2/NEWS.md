# rsdeseq2 0.5.0 (prepared)

* Run supported fitting, formulas, result queries, transforms and normal/ashr/
  apeglm calculations in shared standalone Rust; R translates and constructs objects.
* Support restricted likelihood and saved-dispersion expressions without R execution.
  Reject custom formula functions and unsupported captured values before execution.
* Preserve formula term/factor order, polynomial grouping and offsets. Bound
  interaction expansion and nesting while simplifying repeated variables.
* Share replacement merging, later queries, warning sequences, optimizer controls,
  selected posterior outputs and portable random-state handling across interfaces.
* Retain coefficient identity for posterior summaries and validate result shapes.
  Construct Rust posterior objects through fitting APIs; their identity is private.
* Add configurable grid budgets, cooperative cancellation and safe R interruption/reuse.
* Borrow retained CLI Cook's matrices for diagnostic export.
* Group tooling by purpose, keep current documentation and reproducibility evidence,
  and exclude OS metadata from Rust/R source packages.
* Pass both R lanes, Rust/Python, five platforms, full GTEx and Linux sanitizer checks.
  R package checks retain the known Rust compiled-symbol warning.

See the supported contract and current R guide for method/option limits.
No 0.5.0 release has been published.
