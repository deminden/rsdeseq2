# Version 0.4.0 release candidate

This candidate gathers the changes since 0.3.1 into a minor release. Formula
coding, validation and fitted-state behavior include intentional compatibility
corrections, so this is not a patch-only performance release. Versions agree in
the crate, R package, backend, vendored core and both lockfiles.

## Included changes

- One shared Rust implementation for supported formula expressions, frame row
  selection, categorical preparation, offsets and matrix products. Custom
  formula functions are rejected before execution. R retains symbolic metadata,
  argument translation, identifier alignment and object construction.
- Correct ordered-factor polynomial coding, factor conversion and unused-level
  semantics, higher-order interaction coding and multiplication order. Factor
  labels may merge levels or provide a scalar prefix. Releveling an ordered
  factor is rejected. Ordered coefficients use named or numeric selectors;
  treatment-style level selectors are rejected with an explanatory error.
- Expanded-factor recoding for normal shrinkage uses the standalone Rust engine.
  Native expanded matrices retain ordered polynomial coordinates and use full
  indicators for unordered factors. This representation does not certify the
  unsupported source-ordered formula or failed expanded-prior routes in DESeq2.
- Shared Rust inference, contrasts, Cook's/zero policy, independent filtering and
  BH; bounded result batches and validated reuse of mutable fitted datasets.
- Rank/reduced-model, argument, normalization-offset and worker fixes; fewer
  temporary fit copies and deferred replacement-count allocation.
- Refactored large modules and consolidated repetitive tests while retaining
  independent numerical references. Audited stable dependency updates include
  locfit 0.2.0 and L-BFGS-B 0.3.0.

## Compatibility evidence

The final candidate passes 1,072 Rust tests, three backend tests and an R source
package check with 2,615 assertions, zero failures, zero warnings and zero skips.
Release-profile numerical/CLI checks also pass. Both dependency graphs pass
all-target Clippy with warnings denied and Rust 1.97.1 checks. Vendor/version
integrity and the Python harness suite pass (46 tests, one platform skip).

Direct standalone/wrapped comparisons pass 189 primitive cases, 188 expression
cases, 54 existing native-formula cases and the new 51 factor matrix/frame plus
nine expanded-factor cases. Numeric arrays and primitive metadata are exact;
native and R-style formula labels/order use explicit, checked column mappings.
This does not imply that distinct formula APIs expose identical labels.

The factor fixture generator requires R 4.6.2 / DESeq2 1.52.0 and never loads
rsdeseq2. Fresh regeneration is byte-identical and existing paths are rejected.
An unchanged production-source control from `fada553`, with its old test-only
driver, disagrees with original R in 40 of the 51 matrix cases. Current matrices,
frame classes/levels/references and invalid-input decisions match the fixtures.
Existing frozen fixtures, tolerance bounds and earlier release hash manifests
are unchanged.

All 121 per-contrast and nine global GTEx precision gates pass. All 22 native
result/diagnostic exports are byte-identical to `fada553`. The
[current GTEx comparison](gtex-native-wald-benchmark.md) records the final
candidate's single paired speed/RSS pass and original-DESeq2 comparison.

Fresh native-boundary controls pass four formula and four replacement
unwind/recovery cases under C ASan/UBSan and Rust ASan. The C/Rust sources in
that disposable build match the final candidate. R, preinstalled dependencies
and Rust's standard library are not instrumented; leak detection is disabled.
The ordinary checked library has no test hooks or sanitizer instrumentation.

## Artifacts and paired observations

The checked R source archive is SHA-256
`56ae1c2e617e0ea45a3f6cadd34ece1e654d8e60491cf4830d39f9131c5d0add`.
Its 322 source members match the workspace byte for byte; DESCRIPTION's original
DCF fields match after whitespace unfolding, with only R's generated Packaged
and NeedsCompilation fields added. The tested CLI is
`c41f9357e61ef612ec53bc6f45ae930f05e4782afb30d4710e3b2f2294376dfa`.
The [validation record](data/release_0.4.0_validation.json) binds the checked
source, numerical gates, executable and both measurement panels.

| One paired pass | Previous elapsed | Candidate elapsed | Previous maximum RSS | Candidate maximum RSS |
| --- | ---: | ---: | ---: | ---: |
| Native, 11 contrasts | 38.340 s | 40.670 s | 303.031 MiB | 307.672 MiB |
| R, 20 workflows | 112.764 s | 111.325 s | 1,140.063 MiB | 1,137.766 MiB |

All twenty complete R payloads match byte for byte, with no metadata exception.
Every R case remains below 2 GiB. The native candidate is 6.078% slower in this
pass; the R panel takes 1.276% less elapsed time. Background CPU activity was
observed during measurement. All observations are retained, with no unchanged
case retimed or discarded. These numbers do not certify unchanged speed or an
isolated optimization effect; no optimization round is counted. Publication
readiness still includes the combined performance audit and platform gates.

## Remaining scope and publication

The [Rust-only wrapper migration](rust-only-wrapper-migration.md) is unfinished.
Frozen local VST reuse and optional ashr/apeglm posterior, scaling, tail and
s-value calculations still involve R. General argument/method/warning semantics,
the redundant initial BH pass, C21's authoritative result path and C23's combined
performance/platform audit remain open. This release does not claim a fully
Rust-only R wrapper or complete Bioconductor compatibility.

The local preparation includes the checked R source archive, verified crate,
tested macOS ARM64 executable, checksums, source provenance and release notes.
Publishing still requires the coordinated workflow's reference-R and five
platform-target gates. No remote tag, push or publication is part of this
preparation. See [release procedures](release.md).

## Reproduce the factor checks

```sh
cargo build --release --locked --example evaluate_r_interface
Rscript scripts/generate_native_factor_reference.R /fresh/path/factors.json
Rscript scripts/check_native_factor_identity.R \
  target/release/examples/evaluate_r_interface R_LIBRARY FRESH_OUTPUT_DIR
cargo test --locked --test native_factor --test native_formula
```

Use `deseq2-r462-native` for the original-R generator and a separately installed
candidate source package for the direct comparison. The R CI workflow also runs
all four direct comparison scripts.
