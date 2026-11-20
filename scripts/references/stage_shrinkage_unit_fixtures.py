#!/usr/bin/env python3
"""Extract independent leaf references for the standalone shrinkage crate."""

import argparse
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
REFERENCE_DATA = ROOT / "crates/rsdeseq2/tests/data"
UNIT_DATA = ROOT / "crates/rcompat-shrinkage/tests/data"


def fixture_contents():
    fixtures = (
        ("beta_binomial_posterior_r462.json", "analytic_leaves", "beta_binomial_analytic_leaves_r462.json"),
        ("posterior_analysis_r462.json", "grid_leaves", "posterior_grid_leaves_r462.json"),
    )
    for source, key, output in fixtures:
        reference = json.loads((REFERENCE_DATA / source).read_text())
        leaves = reference[key]
        if not isinstance(leaves, list) or not leaves:
            raise ValueError(f"{source} has no {key}")
        fixture = {"source_fixture": source, "source": reference["source"], key: leaves}
        yield output, json.dumps(fixture, indent=2) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--check", action="store_true", help="check the tracked unit fixtures")
    action.add_argument("--output-dir", type=Path, help="write a new copy of both unit fixtures")
    args = parser.parse_args()
    if args.check:
        for name, contents in fixture_contents():
            if (UNIT_DATA / name).read_text() != contents:
                parser.error(f"{name} differs from the independent source fixture")
    else:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        for name, contents in fixture_contents():
            with (args.output_dir / name).open("x") as output:
                output.write(contents)


if __name__ == "__main__":
    main()
