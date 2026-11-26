#!/usr/bin/env python3
"""Write the strings of the expansion oracle of R's tailored collations.

A character that sorts as several primary weights (a ligature, a squared or
parenthesised form, a stored expansion) is compared at the tertiary level with
the plain spelling of its primaries, element by element. Under a collation
that sorts uppercase first (Danish) the case of each element decides those
comparisons, which the single-code-point oracles cannot show. For every
member of a primary group that R sorts as several atoms
(`root_expansions.tsv`), this script writes the member, its spelling in
lowercase atoms, the spelling with every cased atom in uppercase and with only
the first in uppercase, and the compatibility decomposition of the member.
The uppercase letters are those whose lowercase form (from Rust, in
`decompositions.tsv`) is the atom. `probe_collation_oracle.R sort` then
records R's order of the strings under each locale.

Inputs (default /tmp/rcollate): root_single.tsv, root_expansions.tsv and
decompositions.tsv. Output: expansion_input.tsv.

Usage: collation_expansion_strings.py [--input DIR]
"""

from __future__ import annotations

import argparse
from pathlib import Path


def parse(h: str) -> list[int]:
    return [int(x, 16) for x in h.split()]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--input", type=Path, default=Path("/tmp/rcollate"))
    args = parser.parse_args()
    spelling: dict[int, list[int]] = {}
    for line in (args.input / "root_expansions.tsv").read_text(encoding="utf-8").splitlines():
        if line.startswith("#") or not line:
            continue
        head, atoms = line.split("\t")[:2]
        spelling[int(head, 16)] = parse(atoms)
    nfkd: dict[int, list[int]] = {}
    upper: dict[int, int] = {}
    for line in (args.input / "decompositions.tsv").read_text(encoding="utf-8").splitlines():
        f = line.split("\t")
        c = int(f[0], 16)
        nfkd[c] = parse(f[2])
        low = parse(f[3])
        if len(low) == 1 and low[0] != c:
            upper.setdefault(low[0], c)
    # Primary groups of R's single-code-point order.
    groups: list[list[int]] = [[]]
    for line in (args.input / "root_single.tsv").read_text(encoding="utf-8").splitlines():
        if line.startswith("# completely ignorable"):
            break
        if line.startswith("#"):
            continue
        cp, rel = line.split("\t")[:2]
        groups[-1].append(int(cp, 16))
        if rel == "<":
            groups.append([])
    strings: set[tuple[int, ...]] = set()
    for members in groups:
        head = next((c for c in members if c in spelling), None)
        if head is None:
            continue
        atoms = spelling[head]
        cased = [i for i, a in enumerate(atoms) if a in upper]
        variants = [atoms, [upper.get(a, a) for a in atoms]]
        if cased:
            first = list(atoms)
            first[cased[0]] = upper[atoms[cased[0]]]
            variants.append(first)
        for v in variants:
            strings.add(tuple(v))
        for c in members:
            strings.add((c,))
            if c in nfkd and nfkd[c] != [c]:
                strings.add(tuple(nfkd[c]))
    # Tab and line feed are the level-test separators of the oracle.
    lines = sorted(" ".join(f"{c:04X}" for c in s) for s in strings if not {9, 10} & set(s))
    (args.input / "expansion_input.tsv").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"strings {len(lines)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
