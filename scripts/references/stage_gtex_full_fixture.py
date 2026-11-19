#!/usr/bin/env python3
"""Restore the two large GTEx reference tables from a pinned Git revision."""

import argparse
from hashlib import sha256 as new_sha256
import json
from pathlib import Path
import subprocess
from tempfile import NamedTemporaryFile

ROOT = Path(__file__).resolve().parents[2]
RELATIVE_FIXTURE = Path("crates/rsdeseq2/tests/data/gtex_v8_real_wald")
FIXTURE = ROOT / RELATIVE_FIXTURE
SOURCE_REVISION = "2d86916c089f55bcddb33aee1af9b079f2b9384e"
TABLES = ("counts.tsv", "reference_results.tsv")


def sha256(path):
    digest = new_sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def stage(output_dir, fetch):
    manifest = json.loads((FIXTURE / "compact/selection.json").read_text())
    output_dir.mkdir(parents=True, exist_ok=True)
    for name in TABLES:
        expected = manifest["source_sha256"][name]
        output = output_dir / name
        if output.is_file() and sha256(output) == expected:
            continue
        if fetch:
            known = subprocess.run(
                ["git", "cat-file", "-e", f"{SOURCE_REVISION}^{{commit}}"],
                cwd=ROOT, check=False, stderr=subprocess.DEVNULL,
            )
            if known.returncode:
                subprocess.run(
                    ["git", "fetch", "--no-tags", "--depth=1", "origin", SOURCE_REVISION],
                    cwd=ROOT, check=True,
                )
            fetch = False
        with NamedTemporaryFile(dir=output_dir, delete=False) as handle:
            temporary = Path(handle.name)
            try:
                subprocess.run(
                    ["git", "show", f"{SOURCE_REVISION}:{(RELATIVE_FIXTURE / name).as_posix()}"],
                    cwd=ROOT, check=True, stdout=handle,
                )
                handle.close()
                if sha256(temporary) != expected:
                    raise ValueError(f"SHA-256 mismatch for {name} in {SOURCE_REVISION}")
                temporary.replace(output)
            finally:
                temporary.unlink(missing_ok=True)
        print(f"Verified {name} from {SOURCE_REVISION}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=FIXTURE)
    parser.add_argument("--fetch", action="store_true", help="fetch the pinned commit if absent")
    args = parser.parse_args()
    stage(args.output_dir, args.fetch)


if __name__ == "__main__":
    main()
