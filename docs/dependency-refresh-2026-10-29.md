# Stable Rust dependency refresh

The 2026-10-29 refresh updates both the standalone and R-wrapper Rust dependency
graphs. The baseline is `a6e376a`. Package versions remain 0.3.1, the declared
minimum Rust version remains 1.97.1, and the larger Rust-only migration continues.

## Direct dependencies

The crates.io registry was checked for all 15 direct dependencies, including
development dependencies. Every direct dependency now resolves to its latest
non-yanked stable release at the time of the audit.

| Dependency | Before | After |
| --- | --- | --- |
| rcompat-locfit | 0.1.5 | 0.2.0 |
| rcompat-lbfgsb | 0.2.2 | 0.3.0 |
| statrs | 0.18.0 | 0.19.1 |
| sha2 | 0.10.9 | 0.11.1 |

The other 11 direct dependencies already resolved to their latest stable
versions. The `same-file` declaration now records its resolved 1.0.6 minimum
explicitly. The [direct dependency audit](data/dependency_refresh_2026-10-29_direct.tsv)
records resolved before/after versions, release dates, checksums, declared minimum
Rust versions and registry source URLs.

Both lockfiles were refreshed and the vendored manifest regenerated. Every
production package/version in the R backend lock is also present in the root
lock. The refresh removes the older duplicate nalgebra 0.33.3 dependency and its
associated older numerical dependencies. Statrs now shares nalgebra 0.35.0 with
the main engine.

The locfit public API remains unchanged. Existing local-trend calls use its new
implementation. L-BFGS-B calls use 0.3.0 with the existing default arithmetic
profile; optimizer controls and scientific options are unchanged. The only
application source adjustment is SHA-256 text encoding: sha2's new digest array
does not implement the previous hexadecimal formatter. The completion marker
still records the same 64 lowercase hexadecimal characters. Its existing CLI
test independently formats each digest byte to verify the published hashes.

## Transitive dependencies

The audit covers **120 registry crate names**, with **124 resolved root graph
entries** including optional and platform-specific packages. Cargo selected the
latest compatible stable versions for the existing upstream requirements.
Seven families remain on older major/minor lines:

| Family | Upstream constraint |
| --- | --- |
| anes | criterion requires `^0.1.4` |
| glam 0.30–0.33 | nalgebra exposes optional conversions for those specific versions |
| itertools 0.13 | criterion and criterion-plot require `^0.13` |
| num-bigint 0.4 | num-rational's optional dependency requires `^0.4.0` |
| r-efi 6 | getrandom requires `^6` for its platform support |
| syn 2 | nalgebra-macros and zerocopy-derive require syn 2 |
| windows-link 0.2 | windows-sys requires `^0.2.1` |

These constraints come from current upstream releases. Forcing their newer major
versions would require upstream API migrations. The
[complete registry audit](data/dependency_refresh_2026-10-29_registry_audit.tsv)
records each resolved version, latest stable version, backend presence and
constraint. Optional/platform entries in this inventory do not imply they are
compiled for the default macOS target.

## Validation

The full Rust suite and checked R source package pass with the original
references and tolerances. Root and backend checks pass on Rust 1.97.1, and both
Clippy checks pass without warnings. The full eleven-contrast GTEx replay passes
121 per-contrast and nine global precision gates, with the same 664× minimum
target and 597× minimum global error reductions. Every one of the 22 result and
diagnostic exports is byte-identical to the pre-refresh checkpoint; all recorded
maximum reference errors are unchanged.

The [native precision panel](data/dependency_refresh_2026-10-29_gtex_precision.tsv),
[global gates](data/dependency_refresh_2026-10-29_gtex_global_precision.tsv) and
[reviewed output hashes](data/dependency_refresh_2026-10-29_native_output_sha256.tsv)
retain the independent reference and exact-output evidence.

## Speed and memory

One accepted paired pass covered all eleven complete 56,200-gene GTEx contrasts,
with the previous binary, updated binary and original DESeq2. Each tool used one
thread. Full result and diagnostic exports were verified separately before
measurement, every timed native result matched its reviewed hash, and the 44
frozen input/reference files remained unchanged. Compilation and other fits did
not overlap the measurement.

| Measurement | Before | Updated | Ratio before/updated |
| --- | ---: | ---: | ---: |
| Native total elapsed time, 11 contrasts | 36.064 s | 36.469 s | 0.9889× |
| Native maximum peak RSS | 302.266 MiB | 302.328 MiB | 0.9998× |
| R workflow total elapsed time, 14 cases | 59.448 s | 60.193 s | 0.9876× |
| R workflow maximum peak RSS | 1129.266 MiB | 1129.781 MiB | 0.9995× |

Native elapsed time increased **1.12%** and R workflow elapsed time **1.25%** in
these single passes. These measurements do not establish a speed gain or its
statistical significance; no optimization round is counted. The individual
outliers are retained and the timing pass was not repeated. Original DESeq2's
aggregate time was **192.971 s**, giving **5.291×** DESeq2/updated-Rust speedup.

All **14 serialized R workflow payloads are byte-identical**, including fitted
state and metadata. The panel covers Wald/LRT results, ten-result extraction at
1,200 and 56,200 genes, all three shrinkage methods, VST, replacement with and
without refits, and four rlog sizes. Every case remains below the registered
2 GiB peak-RSS budget. Maximum RSS above is the maximum across each panel, rather
than the memory of a single common case.

See the [native timings](data/dependency_refresh_2026-10-29_native_timings.tsv),
[per-contrast comparisons](data/dependency_refresh_2026-10-29_native_summary.tsv),
[R timings](data/dependency_refresh_2026-10-29_r_timings.tsv) and
[R payload identities](data/dependency_refresh_2026-10-29_r_output_identity.tsv).
Detailed logs and immutable artifacts are retained under
`results/codebase_improvements/dependency-refresh-2026-10-29/`.

## Artifact and interface checks

The final package passes **1,040 Rust tests**, one backend unit test and
**1,426 R assertions**, with one pre-existing opt-in Rust test excluded from the
ordinary suite. The full GTEx replay above covers its separate acceptance scope.
The checked R source archive matches all **242** packaged source files, allowing
only R build-generated DESCRIPTION fields and DCF whitespace. The crate archive
also passes packaging and compilation verification.

All **64 standalone/wrapped bit-equality cases** pass against the checked package,
including guards against R statistical callbacks and targeted R matrix-product
and square-root calculations. Native-boundary tests pass with C ASan/UBSan and
Rust ASan using the updated dependencies. The R interpreter, installed R packages
and Rust standard library are not instrumented, and leak detection is disabled.
Formatting, strict Rust documentation, vendor synchronization and the Python
support-tool suite pass as well.

The [validation record](data/dependency_refresh_2026-10-29_validation.json) binds
these checks to both lockfiles, immutable executables and source archives. The
R reference environment remains `deseq2-r462-native`; existing numerical fixtures,
release manifests and tolerances are preserved. Shared all-zero rules, bounded
result batches and the remaining Rust-only migration still require completion.
