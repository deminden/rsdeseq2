# Frozen local VST migration: source audit

The 0.4.0 candidate still evaluates `dispersionFunction(object)` in R for a
frozen local VST. It also constructs normalized counts and the sinh grid in R.
This path does not satisfy the shared Rust numerical-engine requirement.

In the pinned DESeq2 1.52.0 environment, `localDispersionFit` creates a closure
whose environment holds a fitted `locfit` object. Its function evaluates that
fit at log means and exponentiates the prediction. DESeq2's local VST evaluates
that closure on a 999-point sinh grid, then integrates, interpolates and anchors
the transform. The Rust core already has local trend fitting and
`vst_local_from_dispersion_grid`, but reconstructing a fit from training data
is not proven equivalent to using the exact saved R fit.

The wrapper also accepts arbitrary closures tagged `fitType = "local"`. The
existing integration tests intentionally exercise a custom analytic closure,
including a normalization-factor matrix. A tag alone therefore cannot identify
a genuine DESeq2 locfit closure. Calling it to discover its values would keep
the R calculation and could execute arbitrary code.

The next implementation needs to:

1. Generate independent fixtures from actual DESeq2 local fits, including the
   all-minimum fallback, serialization/reload, both normalization forms, and
   out-of-range grid predictions. Record the original fit representation and
   grid/output values without invoking rsdeseq2.
2. Determine whether the saved locfit evaluation grid and coefficients can be
   validated and imported into a Rust prediction representation without refitting.
   Check interpolation, extrapolation and error behavior against those fixtures.
3. Move maximum normalized count and grid construction into shared Rust, then
   pass only validated, owned fit primitives across the R boundary. Test direct
   standalone/wrapped bit identity independently from DESeq2 precision.
4. Reject unsupported or unrecognized closures without executing them. Preserve
   support only for representations that can be translated safely and verified;
   document any compatibility change explicitly. The earlier user decision to
   reject custom *formula* functions does not by itself establish what to do
   with custom frozen dispersion functions.
5. Re-run the full source-package, GTEx, cross-interface and paired workflow gates
   before making performance or release claims. Keep the 0.4.0 evidence intact.

This audit establishes the concrete representation problem; it does not complete
the frozen-local VST port, the broader Rust-only migration, or the release audit.
