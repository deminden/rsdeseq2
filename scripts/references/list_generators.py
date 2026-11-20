#!/usr/bin/env python3
"""List repository generators, their outputs and the checks that consume them."""

import argparse
import csv
from pathlib import Path
import sys

DIRECTORY = Path(__file__).resolve().parent
MANIFEST = DIRECTORY / "MANIFEST.tsv"
OTHER_OUTPUTS = DIRECTORY / "OTHER_OUTPUTS.tsv"


def read_rows(path):
    with path.open(newline="") as source:
        return list(csv.DictReader(
            (line for line in source if not line.startswith("#")), delimiter="\t",
        ))


def generator_rows():
    return sorted(
        (row for path in (MANIFEST, OTHER_OUTPUTS) for row in read_rows(path)
         if row["generator"] != "-"),
        key=lambda row: (row["generator"], row["pattern"]),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("filter", nargs="?", default="", help="show rows containing this text")
    args = parser.parse_args()
    search = args.filter.casefold()
    writer = csv.DictWriter(
        sys.stdout, fieldnames=("generator", "pattern", "check", "note"),
        delimiter="\t", lineterminator="\n",
    )
    writer.writeheader()
    for row in generator_rows():
        if any(search in field.casefold() for field in row.values()):
            writer.writerow(row)


if __name__ == "__main__":
    main()
