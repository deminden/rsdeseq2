# SQP warning parity and fresh platform references

The final local audit exposed a gap in the next Linux validation lane: it
regenerated older formula references but did not regenerate the new 32-case
complete-formula matrix panel or nine complete/nine primitive sparse SQP
profiles. Those original-only generators are now required in the pinned R lane.
Fresh outputs are saved as CI artifacts, tested against Rust with the unchanged
exact formula/`2e-6` SQP bounds, and used for actual standalone/wrapped checks.
The R lane also records the exact built archive name for checking/upload rather
than relying on a wildcard that could select an unrelated tarball.
The ordinary Rust suite still uses its committed fixtures; CI additionally
supplies fresh files through explicit test-only environment variables.

Strengthening this gate found a real diagnostic omission. The earlier
[sparse SQP report](native-sparse-sqp-development.md) claimed original warning
parity, but its numerical checks did not assert the complete original warning
sequence. On `short_sqp` and `sparse_one_sqp`, Rust emitted ashr's nonconvergence
warning while omitting the preceding mixsqp warning. Exact native/wrapped
agreement could not detect an omission shared by both interfaces.

`location_fit.rs` now emits the mixsqp warning before the ashr warning for a
nonconverged SQP fit. The standalone engine owns the message and original default
line wrapping; R does not format or compute it. Other optimizers retain their
own warning policy. The independent original-reference test now compares the
full warning list as well as posterior moments, s-values, likelihoods and
convergence. No precision tolerance or optimizer calculation changes.

Both generators record platform, R and BLAS/LAPACK provenance. Deliberately
corrupting one matrix value or posterior target makes the fresh-reference tests
fail, establishing that they consume the fresh files rather than silently
falling back to committed data. Saved fixtures remain unchanged.

The [local evidence report](data/native-sqp-warning-local.json) records the
current package, numerical preservation and actual-interface checks. The
[preceding combined measurements](final-plan-local-acceptance.md) retain their
binary/package identities and background-load limitation; they are not relabeled
as measurements of an untested artifact.

The corrected source passes **1,154 Rust tests**, **2,867 R assertions** with
`R CMD check` Status OK, and **49 Python tests** with the real installed-package
oracle gate enabled. The 101-test macOS release-profile panel, full 56,200-gene
test, MSRV 1.97.1 and current scoped native sanitizers pass. All 22 frozen GTEx
exports and eight held-out exports are unchanged. For the nine complete/nine
primitive sparse SQP profiles, every non-warning output field is exact against
the preceding executable; only the expected warning-list addition changes.

One final paired pass measures **43.328 s** current native versus **46.476 s**
before this warning correction and **242.121 s** original DESeq2 (**5.59×**
original/current total-time ratio). Native maximum RSS is **307.469 / 308.000 /
1,768.922 MiB**, respectively. Twenty R workflows take **89.638 versus 89.423 s**
(**+0.240%** elapsed), with maximum RSS **1,133.391 versus 1,133.859 MiB**; all
remain below the predeclared 2 GiB budget. All twenty complete R payloads are
byte-identical under strict comparison, with no correction exception enabled.

The measurement task ran panels sequentially after all other task builds/fits
ended. Other user, compiler, virtualization and system workloads remained active;
load snapshots and every observation are retained. No case was discarded or
retimed, no causal small speed improvement is claimed, and no optimization round
is counted. The [current comparison](gtex-native-wald-benchmark.md) records all
contrasts, stage times and original-package precision; the README compares memory
only against original DESeq2.

## Reproduction

Choose fresh output files/directories and use the pinned native R environment.

```sh
"$R_SCRIPT" --vanilla scripts/generate_formula_fitting_reference.R "$FRESH_FORMULA_JSON"
"$R_SCRIPT" --vanilla scripts/generate_sparse_sqp_reference.R "$FRESH_SQP_JSON"
RSDESEQ2_FORMULA_FITTING_REFERENCE="$FRESH_FORMULA_JSON" \
RSDESEQ2_SPARSE_SQP_REFERENCE="$FRESH_SQP_JSON" \
  cargo test --locked -p rsdeseq2 --test formula_fitting --test sqp_controls
"$R_SCRIPT" --vanilla scripts/check_location_shrinkage_identity.R \
  "$STANDALONE_DRIVER" "$R_LIBRARY" "$FRESH_LOCATION_DIR" "$FRESH_SQP_JSON"
"$R_SCRIPT" --vanilla scripts/check_native_ashr_optimizers.R \
  "$STANDALONE_DRIVER" "$R_LIBRARY" "$FRESH_OPTIMIZER_DIR" - "$FRESH_SQP_JSON"
```

These commands were first tested locally on macOS. The later
[dev validation](dev-validation-2026-11-12.md) records successful current Linux
and five-target execution. The [Linux notes](linux_numerical_parity.md) document
two additional failures found by that execution and their independent fixes.
Tags and publication remain on hold.

## Check platforms without publishing

The `Platform checks` workflow now runs the same eight release-profile numerical
and CLI integration suites on Linux x86-64/ARM64, macOS Intel/ARM64 and Windows
x86-64 for branch pushes, pull requests or manual checks. The existing release
workflow and its provenance writer share `scripts/check_release_target.py`'s
test list. The CI-only workflow has read-only permissions and contains no
packaging, artifact upload, release, tag or registry operation.
The [local platform-path audit](data/platform-validation-local.json) records the
workflow lint, matching matrix, native runner results and environment limits.

For a locally runnable target, use:

```sh
python scripts/check_release_target.py --target aarch64-apple-darwin
```

This ARM Mac cannot execute Intel binaries in its current configuration, lacks a
Linux C cross-compiler and has no local Linux/Windows runtime. The remaining
platform gates required suitable runners. The user subsequently authorized
`dev` branch checks; all five panels now pass as recorded in the final validation.
Tags and publication remain separate actions and are still on hold.
