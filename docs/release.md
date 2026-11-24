# Release preparation

The checkout contains **0.6.2**. See the [version notes](releases/0.6.2.md)
and [local acceptance results](testing-parity.md#ci-status). Release workflows
check the tagged source and bind artifacts to its commit.

## Local source checks

Run the source and package checks before staging a release. The vendor
check verifies version agreement across the Rust crates, R package,
R backend and their lockfiles.

```sh
scripts/release/sync_r_vendor.sh --check
cargo fmt --all -- --check
cargo clippy --locked --workspace --all-targets -- -D warnings
cargo test --locked --workspace
cargo +1.96.0 check --locked --workspace --all-targets
python3 scripts/validation/check_docs.py
python3 -m unittest discover -s scripts/tests -p 'test_*.py'
cargo package --locked --all-features -p rsdeseq2
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

Record package-check warnings with the release results. The local R 4.6.2
macOS arm64 check reports the Rust runtime's `_abort` import.

## CI and artifact binding

Branch pushes run Rust/Python, three R lanes and five release-profile platform
panels. The Parity workflow supplies the full GTEx replay. These jobs do not
publish releases. Windows and Intel macOS generate independent R uniform-CDF
references; pinned Linux R regenerates formula/SQP references and runs actual
interruption and C ASan/UBSan plus Rust ASan checks.

When publication is explicitly requested, the matching version tag triggers
`.github/workflows/release.yml`. Through `validate-release.yml` it requires
Rust/Python, pinned R, full GTEx and the minimum compiler check. Current R and the previous Bioconductor release are advisory in that
workflow. The previous release records its DESeq2 version and package
check results as artifacts. Each binary target
also runs its release-profile numerical/CLI panel.

CI uploads the checked R source tarball and packaged Rust crate with a
manifest binding their hashes to the revision. Library tests run from the
unpacked Rust crate against its locked, published dependencies with the host R
numeric profile.
Publishing verifies the manifest and tag/version agreement.
Crate publication stages version-pinned
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
python3 scripts/release/prepare_release_docs.py --tag v0.6.2 binary \
  --directory results/release-docs --target aarch64-unknown-linux-gnu
python3 scripts/release/prepare_release_docs.py --tag v0.6.2 notes \
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

If an asset build fails because its release harness needs correction, keep the
published crate and immutable tag intact. Correct the workflow on the default
branch, then dispatch **Release binaries** with the original `release_tag`,
`source_run_id`, and `rebuild_binaries=true`. This mode checks out the tag, verifies
the original checked-source manifest, reruns all five target panels and requires
all rebuilt jobs to pass before upload. Original source validation and crate
publication must already have succeeded. Windows and Intel macOS regenerate
the same independent R 4.6.2 CDF references used by normal platform CI, at the
existing bounds. Only the asset workflow is repaired; tagged numerical source
and registry packages remain immutable.

## Publish the rcompat crates

The `rcompat-vX.Y.Z` tag runs the [rcompat release workflow](../.github/workflows/release-rcompat.yml).
It checks the three package versions, tests their packaged archives, and
publishes `rcompat-numeric`, `rcompat-formula`, then `rcompat-shrinkage` to
crates.io. Give each package the tag's version and update the workspace,
fuzz and R backend lockfiles before pushing the tag.
