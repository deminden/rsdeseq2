# Validated development checkpoint — 2026-11-12

All **23 checkpoints** in the [original improvement plan](codebase-improvement-plan-2026-10-22.md)
are accepted for the [documented supported contract](../native-supported-contract.md).
The implementation is pushed to `dev` at `d1ffc9813565352044bbd214fcc436efcd3b05b6`.
This report and the final documentation update change no executable, test,
workflow, dependency or package sources after that CI-tested revision.

The development source has not been tagged, merged to the main branch, released
on GitHub or published to a package registry. Package version 0.4.0 remains a
development version after the already published 0.4.0; a future release needs
version coordination and a new publication instruction.

## Current checks

| Gate | Result and source |
| --- | --- |
| Linux Rust, formatting and Clippy | [Passed](https://github.com/deminden/rsdeseq2/actions/runs/43215474416): 1,155 Rust tests; the existing full-data opt-in test is checked separately below. All 49 Python tooling tests run; two requiring R skip in this Rust-only lane. |
| Pinned R 4.6.2 / original packages | [Passed](https://github.com/deminden/rsdeseq2/actions/runs/43215474305): source package check, independently regenerated Linux formula/SQP/uniform references, actual standalone/wrapped identity with R calculations blocked, and both real-R workflow-comparison gate tests. |
| Current R lane | Passed in the same R run: package and actual-interface checks. This lane remains advisory by policy; its result here is successful. |
| Five target panels | [Passed](https://github.com/deminden/rsdeseq2/actions/runs/43215474471): 101 release-profile numerical/CLI tests on each Linux x86_64/ARM64, macOS Intel/ARM64 and Windows x86_64 target; one existing opt-in test excluded per panel. |
| Full GTEx replay | [Passed](https://github.com/deminden/rsdeseq2/actions/runs/43215558030): all 56,200 genes and six Wald fields. The duplicate ordinary-oracle job was deliberately skipped by `full_only`; the complete Rust lane runs those tests. |
| Minimum compiler | Local Rust 1.97.1 workspace/all-targets/locked check repeated successfully after the Linux correction. |
| Vendor and source identity | Core/R vendor check passes. The exact successful pinned-CI R archive is retained; all 430 included source files match, allowing only standard DESCRIPTION build metadata/whitespace. |

Both Linux R package checks retain **one compiled-code warning** for
`_exit`/`abort`/`exit` symbols in linked Rust code. This is visible in the logs;
no zero-warning Linux package claim is made. Numerical tests and actual interface
checks pass. The preceding macOS package check separately recorded 2,867 passing
assertions with Status OK and no warnings/skips; that report keeps its own scope.

The [machine-readable CI proof](../data/dev-validation-2026-11-12.json) records every
run, successful job, log hash, reference source and measurement identity.

## Failures found and corrected by actual Linux execution

The complete formula-fitting tests initially selected the macOS snapshot on
Linux. Nine of 32 original-R profiles differ by at most `1.11e-16` across these
platforms. The new Linux snapshot preserves original R's numerical fields
exactly; both fresh and committed-reference comparisons retain exact bits.
The portable elementary-function panel also incorrectly selected a nonexistent
Linux filename. Its checker now selects the common reference explicitly and
reports missing required files directly.

Fresh original ashr then exposed a half-uniform posterior-SD error exceeding the
existing `2e-6` bound. Original Linux R rounds a density square and addition
separately, while the macOS reference uses a fused operation. Rust now preserves
the original platform policy in the shared helper. This is a core arithmetic
correction, with no R calculation or optimizer change. Both fresh independent
accuracy and exact standalone/wrapped identity pass on actual Linux.
See the [Linux arithmetic notes](../linux_numerical_parity.md),
[formula proof](../data/dev-ci-formula-2026-11-12.json) and
[posterior proof](../data/dev-ci-uniform-2026-11-12.json). No existing reference
value, precision limit or macOS arithmetic was changed to close either gate.

## Precision, speed and memory scope

The rebuilt macOS executable is byte-identical to the measured executable:
`447ff5b2df1d2f55b94d4940c2c6507cb5063c859e4d6e3c667bae3c8982a780`.
The [existing single paired pass](../reports/gtex-native-wald-benchmark.md) therefore retains
its exact executable identity: **43.328 s Rust / 242.121 s original DESeq2
(5.59× time ratio)**, maximum process RSS **307.469 / 1,768.922 MiB**.
All eleven contrasts, 121 per-contrast plus nine global precision gates and 22
accepted native exports remain intact. The independent held-out Wald/LRT exports
remain exact on the accepted native path.

The earlier twenty macOS R workflows measured **89.638 s**, maximum RSS
**1,133.391 MiB**, and exact payload preservation under the 2 GiB budget. These
remain observations of their recorded package. The only runtime-source change
since then is the Linux x86_64 density branch; the macOS expression and default
quantification path are unchanged. No Linux speed/memory result is inferred.
External load was present, all observations were retained, and no new timing
pass, causal speed improvement or optimization round is claimed.

The [local acceptance](final-plan-local-acceptance.md) and
[warning follow-up](native-sqp-warning-platform-gates.md) retain their scoped
sanitizer, authority, boundary, API and workflow evidence. Current cross-platform
execution supplements those checks; it does not broaden their instrumentation.

## Final local bundle

`dist/rsdeseq2-0.4.0-local-<final-revision>.zip` contains the verified crate,
byte-for-byte CI-checked R source, macOS ARM64 executable, tracked repository
source archive, provenance and SHA-256 checksums. Its verification record binds
the final documentation revision to the CI-tested code revision and checks that
only the main README and docs changed between them. It also verifies every
runtime source path, package source file, repository archive member and bundled
artifact. Personal files, build caches and the untracked optimization journal
are excluded. The generated verification record and raw logs stay under
`results/dev-ci/`; no generated bundle is committed.

This closes C23's documentation and development-artifact acceptance. A future
publication still runs the coordinated tag gates and uses a separately agreed
version. Branch CI and this local bundle do not authorize that publication.
