# Latest-stable dependency recheck

The crates.io API was rechecked at **2026-10-30 06:35:50 UTC**, following the
[numerical dependency refresh](dependency-refresh-2026-10-29.md). The baseline
for this follow-up is `af6e9b8`. All **15 direct dependencies**, including
development dependencies, already resolve to their latest non-yanked stable
release. Both standalone Rust and the R wrapper use **rcompat-locfit 0.2.0** and
**rcompat-lbfgsb 0.3.0**. Their controls and arithmetic profile are unchanged.

Six related transitive packages had been published earlier that morning:

| Package | Before | Updated |
| --- | --- | --- |
| js-sys | 0.3.106 | 0.3.107 |
| web-sys | 0.3.106 | 0.3.107 |
| wasm-bindgen | 0.2.129 | 0.2.130 |
| wasm-bindgen-macro | 0.2.129 | 0.2.130 |
| wasm-bindgen-macro-support | 0.2.129 | 0.2.130 |
| wasm-bindgen-shared | 0.2.129 | 0.2.130 |

The root lockfile now includes these releases. Their route is Criterion's
Plotters WebAssembly support. They are absent from the R backend and from the
native macOS build, including native development/build dependencies with all
features enabled. `once_cell` 1.21.4 was removed because the updated wasm-bindgen
no longer requires it. The R backend lockfile is byte-identical to the baseline;
every backend registry package/version remains present in the root lockfile.

The audit examined all 120 previously resolved registry names. The final graph
contains **119 names / 123 registry entries** after the unused package removal.
The seven upstream-constrained families documented in the preceding refresh
remain: anes, older optional glam conversions, itertools, num-bigint, r-efi,
syn 2 and windows-link. No incompatible major versions are forced into their
parents. The [complete registry inventory](data/dependency_recheck_2026-10-30_registry.tsv)
records resolved/latest versions, constraints, checksums and registry sources.
The [validation record](data/dependency_recheck_2026-10-30_validation.json)
records publication times, lockfile identities and native dependency-tree hashes.

## Verification and numerical scope

Locked standalone checks for all targets, the locked R-backend check, vendor
integrity and crate packaging/compilation verification pass. Complete native
dependency trees before/after are byte-identical for both graphs. A locked
WebAssembly dependency-tree query resolves the new versions; this is resolution
evidence, not a claim that this project was compiled or tested on WebAssembly.
No numerical source, direct numerical dependency, R package source or reference
was changed by this follow-up. Package versions remain 0.3.1; Rust minimum
version remains 1.97.1.

The preceding [preparation checkpoint](rust-only-preparation-2026-10-30.md)
retains its 1,060 Rust tests, 1,668 R assertions, 131 direct interface identity
cases, 121+9 GTEx precision gates and 22 identical native exports. Its measured
eleven-contrast native throughput gain is 1.61%, with all twenty paired R
payloads byte-identical and individual memory costs explicitly retained. This
inactive dependency update does not introduce a new timing claim or optimization
round. The comparison environment remains `deseq2-r462-native`; its independently
pinned R reference packages are unchanged. The broader Rust-only migration and
final combined performance/memory audit remain open.
