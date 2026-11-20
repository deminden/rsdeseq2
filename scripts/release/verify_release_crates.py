#!/usr/bin/env python3
"""Build and test checked .crate archives with their locked dependencies."""
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import shutil
import subprocess
import tarfile
import tempfile

ROOT = Path(__file__).resolve().parents[2]


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def extract_crate(archive, vendor):
    stem = archive.stem
    source = vendor / stem
    source.mkdir()
    with tarfile.open(archive) as package:
        for member in package:
            parts = PurePosixPath(member.name).parts
            if (
                not parts
                or parts[0] != stem
                or any(part.endswith((" ", ".")) for part in parts)
                or "\\" in member.name
                or ":" in member.name
            ):
                raise ValueError(f"Invalid path in {archive.name}: {member.name!r}")
            target = source.joinpath(*parts[1:])
            if member.isdir():
                target.mkdir(parents=True, exist_ok=True)
            elif member.isfile() and len(parts) > 1:
                target.parent.mkdir(parents=True, exist_ok=True)
                data = package.extractfile(member)
                if data is None:
                    raise ValueError(f"Missing file data in {archive.name}: {member.name!r}")
                with data, target.open("xb") as output:
                    shutil.copyfileobj(data, output)
                target.chmod(0o755 if member.mode & 0o100 else 0o644)
            else:
                raise ValueError(f"Unsupported entry in {archive.name}: {member.name!r}")
    return source


def main():
    packages = sorted((ROOT / "target/package").glob("*.crate"))
    if not packages:
        raise ValueError("No Cargo archives to verify in target/package")

    # The workspace Cargo config selects this profile for tests that use R fixtures.
    test_env = {**os.environ, "RSDESEQ2_NUMERIC_PROFILE": "r-host"}
    with tempfile.TemporaryDirectory(prefix="rsdeseq2-release-crates-") as temporary:
        work = Path(temporary)
        vendor = work / "vendor"
        config = work / "vendor.toml"
        config.write_text(subprocess.check_output(
            ["cargo", "vendor", "--locked", str(vendor)], cwd=ROOT, text=True,
        ))

        for archive in packages:
            source = extract_crate(archive, vendor)
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
            subprocess.run(
                ["cargo", "test", "--lib", *options], cwd=work, env=test_env, check=True,
            )


if __name__ == "__main__":
    main()
