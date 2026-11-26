#!/usr/bin/env python3
"""Write the small committed samples of R's collation orders and the oracle
summary.

Inputs (default /tmp/rcollate): root_single.tsv, root_strings.tsv,
tailoring_<locale>.tsv, <locale>_single.tsv and expansion_<locale>.tsv,
written by the probe scripts in this directory. The oracles are not
committed (about 450 MB); their SHA-256, sizes and the verification counts go
to
docs/provenance/observations/collation_r462_macos/oracle_summary.tsv, and a
fixed-seed sample of each tailoring order goes to
crates/rcompat-collate/tests/data/sorted_tailorings.tsv for the crate's tests.
A sample of a sorted list is still in R's order: two sampled strings tie
exactly when every string between them ties.

Usage: sample_collation_orders.py [--input DIR] [--per-locale N]
"""

from __future__ import annotations

import argparse
import hashlib
import random
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SUMMARY = ROOT / "docs/provenance/observations/collation_r462_macos/oracle_summary.tsv"
SAMPLE = ROOT / "crates/rcompat-collate/tests/data/sorted_tailorings.tsv"
LOCALES = [
    "sv", "fi", "nb", "nn", "es", "gl", "pl", "cs", "sk", "hu", "hr", "bs", "sl", "ro", "lt", "lv",
    "et", "tr", "az", "ru", "uk", "be", "bg", "sr", "mk", "de_collation_phonebook",
    "da", "is", "el", "fr_CA",
]


def rows(path: Path):
    out = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("# completely ignorable"):
            break
        if line.startswith("#") or not line:
            continue
        f = line.split("\t")
        out.append((f[0], f[1] if len(f) > 1 else ""))
    return out


def header_locale(path: Path) -> str:
    for line in path.read_text(encoding="utf-8").splitlines():
        if "requested " in line:
            return line.split("requested ")[1].split(";")[0]
    return "root"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--input", type=Path, default=Path("/tmp/rcollate"))
    parser.add_argument("--per-locale", type=int, default=150)
    args = parser.parse_args()
    files = ["root_single.tsv", "root_strings.tsv", "expansion_root.tsv", "tailoring_root.tsv"]
    for loc in LOCALES:
        files += [f"tailoring_{loc}.tsv", f"{loc}_single.tsv", f"expansion_{loc}.tsv"]
    summary = ["# oracle file\tstrings\tadjacent pairs\tties\tSHA-256 (regenerate with the probe scripts)"]
    for name in files:
        path = args.input / name
        r = rows(path)
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        ties = sum(1 for _, rel in r if rel == "=")
        summary.append(f"{name}\t{len(r)}\t{len(r) - 1}\t{ties}\t{digest}")
    SUMMARY.write_text("\n".join(summary) + "\n", encoding="utf-8")
    rng = random.Random(20261127)
    out = ["# Samples of R 4.6.2 (ICU 78.3) orders under tailored collations, written by",
           "# scripts/provenance/sample_collation_orders.py from probe_collation_tailorings.R.",
           "# A line \"@ locale\" starts a list; then code points and the relation to the next."]
    for loc in LOCALES:
        path = args.input / f"tailoring_{loc}.tsv"
        r = rows(path)
        group, g = [], 0
        for _, rel in r:
            group.append(g)
            if rel != "=":
                g += 1
        pick = sorted(rng.sample(range(len(r)), args.per_locale))
        out.append(f"@ {header_locale(path)}")
        for a, b in zip(pick, pick[1:] + [None]):
            rel = "" if b is None else ("=" if group[a] == group[b] else "<")
            out.append(f"{r[a][0]}\t{rel}")
    SAMPLE.write_text("\n".join(out) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
