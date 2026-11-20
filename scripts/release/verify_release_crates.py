#!/usr/bin/env python3
"""Build and test checked .crate archives with their locked dependencies."""
import hashlib
import json
from pathlib import Path
import subprocess
import tarfile
import tempfile

ROOT = Path(__file__).resolve().parents[2]


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    packages = sorted((ROOT / "target/package").glob("*.crate"))
    if not packages:
        raise ValueError("No Cargo archives to verify in target/package")

    with tempfile.TemporaryDirectory(prefix="rsdeseq2-release-crates-") as temporary:
        work = Path(temporary)
        vendor = work / "vendor"
        config = work / "vendor.toml"
        config.write_text(subprocess.check_output(
            ["cargo", "vendor", "--locked", str(vendor)], cwd=ROOT, text=True,
        ))

        for archive in packages:
            with tarfile.open(archive) as package:
                package.extractall(vendor, filter="data")
            source = vendor / archive.stem
            # Dependents' packaged lockfiles pin the sibling archive checksum.
            files = {p.relative_to(source).as_posix(): sha256(p)
                     for p in source.rglob("*") if p.is_file()}
            (source / ".cargo-checksum.json").write_text(
                json.dumps({"package": sha256(archive), "files": files}),
            )

        for archive in packages:
            options = [
                "--locked", "--offline", "--all-features",
                "--manifest-path", str(vendor / archive.stem / "Cargo.toml"),
                "--config", str(config), "--target-dir", str(work / "target"),
            ]
            subprocess.run(["cargo", "build", *options], cwd=work, check=True)
            subprocess.run(["cargo", "test", "--lib", *options], cwd=work, check=True)


if __name__ == "__main__":
    main()
