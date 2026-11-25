#!/usr/bin/env python3
"""Check that provenance statements in tracked sources are covered by the ledger.

Every tracked source file must match a row of docs/provenance/ledger.md, and
every provenance statement (phrases such as "adapted from", names of upstream
projects, licence file names) must sit in a file whose row is resolved
(`clean` or `notice kept`). The check makes statements traceable; it is not a
reason to delete or soften them.

Usage: scripts/validation/check_provenance.py [--list] [--report]
  --list    print every provenance statement with its ledger row
  --report  print unresolved rows and exit 0
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
LEDGER = ROOT / "docs/provenance/ledger.md"
RESOLVED = {"clean", "notice kept"}
STATUSES = RESOLVED | {"to replace", "to remove", "owner to confirm"}

# Tracked files whose text is code or code documentation.
SOURCE_GLOBS = [
    "crates/**/*.rs",
    "crates/**/*-LICENSE",
    "crates/*/THIRD_PARTY_NOTICES.md",
    "THIRD_PARTY_NOTICES.md",
    "r-pkg/rsdeseq2/R/*.R",
    "r-pkg/rsdeseq2/src/rust/src/**/*.rs",
    "r-pkg/rsdeseq2/tests/testthat/*.R",
    "r-pkg/rsdeseq2/inst/COPYRIGHTS",
]

STATEMENT = re.compile(
    r"adapted from|adaptation of|\bports? of\b|\bported from|\bportions? of\b"
    r"|translat(?:ed|ion) (?:of|from)|copied (?:from|verbatim)"
    r"|follows the (?:fortran|source|c code|algorithm)|\bfollows R's\b"
    r"|openblas|\blinpack\b|\blapack\b|\bdqrdc2?\b|\bglibc\b|\bboost\b|\bstatrs\b"
    r"|\blbfgspp\b|\bmixsqp\b|\bslatec\b|\bnetlib\b|r-source|\bnmath\b|src/appl"
    r"|\bqf\.c\b|external/|copyright \(c\)|-LICENSE\b|\bGPL\b|\bLGPL\b|\bBSD\b",
    re.IGNORECASE,
)


def glob_to_regex(pattern: str) -> re.Pattern[str]:
    out = []
    i = 0
    while i < len(pattern):
        if pattern.startswith("**/", i):
            out.append("(?:.*/)?")
            i += 3
        elif pattern.startswith("**", i):
            out.append(".*")
            i += 2
        elif pattern[i] == "*":
            out.append("[^/]*")
            i += 1
        else:
            out.append(re.escape(pattern[i]))
            i += 1
    return re.compile("".join(out) + r"\Z")


def ledger_rows() -> list[dict[str, str]]:
    rows = []
    in_rows = False
    for line in LEDGER.read_text().splitlines():
        if line.startswith("## "):
            in_rows = line.strip() == "## Rows"
            continue
        if not in_rows or not line.startswith("|") or line.startswith("|---"):
            continue
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if len(cells) != 5 or cells[0] == "Path":
            continue
        path, computes, source, evidence, status = cells
        if status not in STATUSES:
            raise SystemExit(f"ledger: unknown status {status!r} for {path}")
        rows.append(
            {
                "path": re.sub(r"`([^`]*)`", r"\1", path),
                "external": not path.startswith("`") or not path.endswith("`"),
                "status": status,
                "source": source,
            }
        )
    return rows


def tracked_sources() -> list[str]:
    out = subprocess.run(
        ["git", "ls-files"],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.split("\n")
    patterns = [glob_to_regex(g) for g in SOURCE_GLOBS]
    return sorted(p for p in out if p and any(r.match(p) for r in patterns))


def row_for(path: str, rows: list[dict[str, str]]) -> dict[str, str] | None:
    best = None
    for row in rows:
        if row["external"]:
            continue
        if glob_to_regex(row["path"]).match(path):
            if best is None or len(row["path"]) > len(best["path"]):
                best = row
    return best


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--list", action="store_true")
    parser.add_argument("--report", action="store_true")
    args = parser.parse_args()

    rows = ledger_rows()
    problems = []
    for path in tracked_sources():
        row = row_for(path, rows)
        if row is None:
            problems.append(f"{path}: no ledger row")
            continue
        text = (ROOT / path).read_text(errors="replace")
        hits = [
            (n, m.group(0))
            for n, line in enumerate(text.splitlines(), 1)
            for m in [STATEMENT.search(line)]
            if m
        ]
        if path.endswith("-LICENSE") and not hits:
            hits = [(1, "licence file")]
        for n, phrase in hits:
            if args.list:
                print(f"{path}:{n}: {phrase!r} -> {row['path']} [{row['status']}]")
            if row["status"] not in RESOLVED:
                problems.append(f"{path}:{n}: {phrase!r} under unresolved row {row['path']} ({row['status']})")
    unresolved = [r for r in rows if r["status"] not in RESOLVED]
    for row in unresolved:
        problems.append(f"ledger row unresolved: {row['path']} ({row['status']})")
    confirm = sum(1 for r in rows if r["status"] == "owner to confirm")

    for p in problems:
        print(p, file=sys.stderr)
    print(
        f"{len(problems)} problem(s); {len(unresolved)} unresolved row(s), "
        f"{confirm} marked 'owner to confirm'",
        file=sys.stderr,
    )
    if args.report:
        return 0
    return 1 if problems else 0


if __name__ == "__main__":
    raise SystemExit(main())
