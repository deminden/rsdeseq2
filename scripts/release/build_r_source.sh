#!/usr/bin/env sh
# Build the self-contained R source package (hardening 8.3).
#
# The Rust core is not committed under r-pkg/rsdeseq2/src/rust/vendor. This
# script generates it from crates/, vendors the third-party crates with
# `cargo vendor` so the package builds with `cargo build --offline` (as CRAN
# and Bioconductor require), runs `R CMD build`, and removes the generated
# trees again. Usage, from anywhere:
#   scripts/release/build_r_source.sh [output-directory]
set -eu
repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
out_dir=${1:-$(pwd)}
mkdir -p "$out_dir"
out_dir=$(CDPATH= cd -- "$out_dir" && pwd)
rust_dir="$repo_root/r-pkg/rsdeseq2/src/rust"

cleanup() {
  rm -rf "$rust_dir/vendor-crates" "$rust_dir/cargo"
}
trap cleanup EXIT

python3 "$repo_root/scripts/release/sync_r_vendor.py"
# cargo vendor prints the source replacement; configure passes the file to
# `cargo build --offline --config cargo/config.toml` when it exists. Cargo
# resolves its relative `directory` against the parent of `cargo/`, as for
# a `.cargo/config.toml`, without shipping a hidden directory.
mkdir -p "$rust_dir/cargo"
(cd "$rust_dir" && cargo vendor --locked --versioned-dirs vendor-crates) > "$rust_dir/cargo/config.toml"
(cd "$out_dir" && R CMD build "$repo_root/r-pkg/rsdeseq2")
