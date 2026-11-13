# Release preparation

The current source version is **0.5.0**. Its implementation is prepared for `main` and has
passed [source CI](testing-parity.md#ci-status). No 0.5.0 tag, GitHub release or
registry publication has been created. [Release notes](releases/0.5.0.md)
describe the prepared distribution.

## Local source checks

The Rust crate, R package, R backend, vendored core and own lockfile entries must
agree on the version. Check the committed vendor without repairing it:

```sh
scripts/release/sync_r_vendor.sh --check
cargo fmt --all -- --check
cargo clippy --locked --workspace --all-targets -- -D warnings
cargo test --locked --workspace
cargo +1.97.1 check --locked --workspace --all-targets
python3 scripts/validation/check_docs.py
python3 -m unittest discover -s scripts/tests -p 'test_*.py'
cargo package --locked -p rsdeseq2
```

Run the R package check, exact interface checks and full GTEx replay from
[testing](testing-parity.md). Numerical changes also require the eleven-contrast
precision/result/diagnostic gates and matched time/RSS procedure in
[reproducibility](reproducibility.md) and [benchmarks](benchmarks.md).
Normal acceptance replays independent fixtures; generating new expected values
is a separate reviewed maintenance operation.

Review staged files. Exclude caches, generated objects/archives, datasets, local
journals, OS metadata and machine-specific paths. Both package inventories must
exclude `.DS_Store`. Preserve measurement/reference identities rather than
relabelling them as measurements of a new binary.

R package checks retain the documented compiled-symbol warning for Rust runtime
`abort`/`exit` entry points. Passing CI is not a warning-free or CRAN-acceptance claim.

## CI and artifact binding

Branch pushes run Rust/Python, both R lanes and five release-profile platform
panels. The Parity workflow supplies the full GTEx replay. These jobs do not
publish releases. Windows and Intel macOS generate independent R uniform-CDF
references; pinned Linux R regenerates formula/SQP references and runs actual
interruption and C ASan/UBSan plus Rust ASan checks.

When publication is explicitly requested, the matching version tag triggers
`.github/workflows/release.yml`. Through `validate-release.yml` it requires
Rust/Python, pinned R, full GTEx and the minimum compiler check. Current R is
advisory in that workflow; its actual result remains visible. Each binary target
also runs its release-profile numerical/CLI panel.

CI uploads the R source tarball it checked and the verified crate archive with
a manifest binding their hashes to the revision. Publishing verifies that
manifest and tag/version agreement. Crate publication stages version-pinned
README links and requires repackaging to reproduce the checked archive bytes.
Binary provenance records compiler, commit and executable/archive hashes.

Binary archives use `scripts/release/prepare_release_docs.py` to generate a
platform-specific README with executable instructions, local license notices,
and documentation links pinned to the release tag. The same helper renders
GitHub release notes with absolute tag-pinned links; missing notes fail packaging
instead of silently generating a different release description. Repository CI
checks these generated documents without building binaries or contacting GitHub.

To inspect the documentation locally, use a fresh ignored staging directory:

```sh
mkdir -p results/release-docs
cp LICENSE THIRD_PARTY_NOTICES.md results/release-docs/
python3 scripts/release/prepare_release_docs.py --tag v0.5.0 binary \
  --directory results/release-docs --target x86_64-unknown-linux-gnu
python3 scripts/release/prepare_release_docs.py --tag v0.5.0 notes \
  --output results/release-docs/release-notes.md
```

These commands only prepare local files; they do not create a tag or release.

## Publication procedure

Only perform these steps after an explicit publication request:

1. Select committed source, confirm an unused version/tag and align all own
   package versions. Complete the source and numerical gates above.
2. Update the README installation/publication status for the intended release,
   then verify CI for that final commit, current notes, package contents and
   artifact identities.
3. Push the matching `vX.Y.Z` tag. The coordinated workflow validates source
   before crate publication and binary packaging/upload.
4. Confirm checked source artifacts, five target archives, checksums and
   provenance records are attached to the matching release.

Crate publishing uses `publish-crates.yml` and `CARGO_REGISTRY_TOKEN`.
`release-binaries.yml` distributes Linux x86_64/aarch64, macOS x86_64/aarch64
and Windows x86_64 archives.

If only the final GitHub upload fails after validation and all binary builds
pass, run **Release binaries** manually from the default branch with the existing
`release_tag` and original coordinated `source_run_id`. Recovery verifies the
tag, original checks and hashes and uploads already checked archives. It does
not move the tag or republish the crate.
