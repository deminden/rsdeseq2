#!/usr/bin/env python3
"""Fail when two top-level modules of one workspace crate reference each other.

Edges are `crate::<module>` paths in production sources (comments, strings and
`#[cfg(test)]` items excluded), as counted by `hardening_metrics.py`.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location('hardening_metrics', ROOT / 'scripts/validation/hardening_metrics.py')
METRICS = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(METRICS)

def check(pairs: set[tuple[str, str]]) -> list[str]:
    return [f'mutual dependency: {a} <-> {b}' for a, b in sorted(pairs)]


def main() -> int:
    pairs = set(METRICS.mutual_pairs(METRICS.production_sources()))
    problems = check(pairs)
    for problem in problems:
        print(problem)
    if not problems:
        print('module dependency check passed (zero mutual pairs)')
    return 1 if problems else 0


if __name__ == '__main__':
    raise SystemExit(main())
