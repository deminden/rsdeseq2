#!/usr/bin/env sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
vendor="$repo_root/r-pkg/rsdeseq2/src/rust/vendor/rsdeseq2"
rm -rf "$vendor/src"
mkdir -p "$vendor"
cp -R "$repo_root/crates/rsdeseq2/src" "$vendor/src"
rm -f "$vendor/src/main.rs"
cat > "$vendor/README.md" <<'EOF'
# Vendored Rust Core

This package contains the Rust core used to build the R backend. Regenerate it
from `crates/rsdeseq2` with `scripts/sync_r_vendor.sh` before packaging a
release. Do not edit this copy directly.
EOF
cat > "$vendor/Cargo.toml" <<'EOF'
[package]
name = "rsdeseq2"
version = "0.3.0"
edition = "2024"
license = "MIT"
rust-version = "1.97.1"

[lib]
name = "rsdeseq2"
path = "src/lib.rs"

[dependencies]
clap = { version = "4.6.7", features = ["derive"] }
csv = "1.4.0"
nalgebra = "0.35.0"
rayon = "1.12.0"
rcompat-lbfgsb = "0.2.2"
rcompat-locfit = "0.1.5"
serde = { version = "1.0.229", features = ["derive"] }
serde_json = "1.0.154"
statrs = "0.18.0"
thiserror = "2.0.21"
EOF
printf 'Synchronized Rust core into %s\n' "$vendor"
