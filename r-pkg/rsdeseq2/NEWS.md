# rsdeseq2 0.6.2

* Bundle version 0.1.2 of the Rust numerical, formula, and shrinkage crates.
* Use Rust binary128 arithmetic for the Linux arm64 formula profile.

# rsdeseq2 0.6.1

* Update the bundled Rust numerical, formula, and shrinkage crates to 0.1.1.

# rsdeseq2 0.6.0

* With `fitType = "parametric"`, a failed parametric dispersion trend is now
  replaced by the local trend inside the Rust core for `DESeq()` as well as the
  transforms, with DESeq2's message. The fit and transform metadata record the
  trend actually used. Previously `DESeq()` stopped with "dispersion fit did
  not converge".
* Remove the `parityMode` argument of `DESeq()` and `fitRust()`. Only strict
  parity mode was ever implemented; `"fast"` and `"diagnostic"` had no effect.

# rsdeseq2 0.5.0

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
