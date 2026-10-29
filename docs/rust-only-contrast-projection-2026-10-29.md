# Shared standalone contrast projection

This checkpoint moves production R Wald/LRT contrast effects and standard errors into
the standalone Rust engine. R no longer computes these with matrix products,
per-gene covariance loops or square roots. Primitive Wald queries call the same
projection and hypothesis implementation as native fitted contrast queries.
LRT display effects use the raw projection without evaluating a Wald test.
The broader [Rust-only migration](rust-only-wrapper-migration.md) remains active;
package versions remain 0.3.1.

## Implementation and numerical policy

`MatrixView` borrows shape-checked row-major or column-major storage.
`CovarianceView` describes the corresponding per-gene covariance tensor.
`InferenceEstimates` validates aligned array shapes once and can be reused for
multiple requests. Neither view copies nor transposes full beta, SE or covariance
arrays. Layout changes indexing, never the arithmetic implementation or logical
summation order. A projection reuses one coefficient-squared scratch buffer
across genes and allocates its gene-length output vectors.

The shared projection retains the existing native scaled-sum and ordered-product
arithmetic. A single active coefficient uses supplied marginal SEs, as the R
wrapper previously did; covariance remains necessary for general contrasts or
when marginal SEs are absent. This also makes SE-only primitive inputs available
through the public standalone API. Missing variance context gives an explicit
error rather than a fabricated diagonal covariance.

Wald contrasts retain the existing standalone joint invalid-row policy: invalid
estimates, SEs or statistics mask both displayed inference fields. Missing t df
alone retain the statistic, effect and SE while masking the p-value. The wrapper
now follows that same rule. Existing native overflow/cancellation assertions pass
unchanged; the initial independent-field implementation failed two of them and
was corrected without changing their expectations.

Raw projection keeps independently available effects and SEs, including zero
SEs. Native LRT builders and fitted-object contrast queries now use this path,
matching primitive R LRT projection. A missing SE must not discard a valid effect
or alter the stored full-versus-reduced LRT statistic/p-value. The public legacy
result builder accepting an already-computed Wald output remains supported.

The C entry points validate types, complete matrix/tensor shapes, storage lengths,
checked products, df and scalar options before borrowing pointers. Optional
buffers and Wald/Normal/t modes have explicit flags. No pointer is retained and
no R numerical callback or alternate numerical mode is introduced. Inference
argument parsing is shared by the vector and contrast bridge modules.

## Verification

- Rust: **1,040 passed**, zero failures, one existing opt-in full-GTEx test ignored
  in the ordinary suite; the separate full replay passes below. Backend unit
  test: one passed.
- Checked R source package: **Status OK**, **1,426 assertions**, no warnings,
  failures or skips. All **242** packaged source files match current sources,
  allowing only generated DESCRIPTION fields and DCF whitespace.
- **64 direct standalone/wrapped bit-equality cases** pass against the checked
  package. Eighteen new contrast cases cover raw/Normal/t outputs, fractional and
  signed weights, marginal-SE precedence, covariance-only/SE-only inputs, zero
  and missing SEs, and missing df. R storage is column-major and standalone test
  input is row-major. Contrast checks additionally trap R matrix products and
  square roots, alongside the existing statistical-function/import guards.
- Independent small-matrix R references verify effects, covariance SEs and
  directional t inference with the existing strict precision bounds. Asymmetric
  covariance entries expose tensor-index mistakes. Rust tests cover shape
  aliases/overflow, unavailable context and preservation of LRT inference when
  displayed SEs cannot form a Wald test. The **1,169** pinned scalar probability
  references remain unchanged and pass.
- The eleven full GTEx contrasts pass **121 per-contrast and nine global gates**.
  Minimum target/global error reductions remain **664×/597×**. All **22** result
  and diagnostic exports are byte-identical to the preceding shared-probability
  checkpoint. Earlier probability corrections are not compared as if they were
  byte-identical to the original approximation.
- Root/backend Clippy and formatting, Rust 1.97.1, strict Rust documentation and
  vendor synchronization and crate packaging checks pass. C ASan/UBSan and Rust ASan boundary tests
  pass for the final bridge. R, installed dependencies and Rust's standard library
  are not instrumented; leak detection is disabled.

The [validation record](data/rust_only_contrasts_2026-10-29_validation.json)
binds these checks to immutable executables and the checked source archive.
See the [GTEx gates](data/rust_only_contrasts_2026-10-29_gtex_precision.tsv),
[global gates](data/rust_only_contrasts_2026-10-29_gtex_global_precision.tsv) and
[export identity table](data/rust_only_contrasts_2026-10-29_native_output_identity.tsv).
Detailed logs and artifacts remain under
`results/codebase_improvements/rust-only-contrasts-*`.

## Remaining work

Contrast all-zero decisions and remaining Cook's/offset/formula calculations
still need shared Rust paths. The multi-request bridge and `resultsList` have
not been restored; reusable borrowed views alone do not complete C20. Optional
ashr/apeglm/posterior migration and the final C01–C23 audit remain outstanding.
This checkpoint makes no speed/RSS claim and counts no optimization round.
The combined result/batch change still requires the prescribed paired timing
and memory assessment.
