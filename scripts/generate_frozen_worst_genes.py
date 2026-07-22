#!/usr/bin/env python3
"""Generate a frozen worst-gene fixture from parity diagnostics."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


RESULT_COLUMNS = ["baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]
OUTPUT_COLUMNS = ["rank", "gene", "result_column", "baseline_abs_error"]


def positive_int(text: str) -> int:
    value = int(text)
    if value <= 0:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return value


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--diagnostics", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--rows", type=positive_int, default=100)
    return parser.parse_args()


def parse_finite_error(row: dict[str, str], column: str, row_number: int) -> float | None:
    text = row[column]
    try:
        value = float(text)
    except ValueError as error:
        raise ValueError(f"row {row_number}: {column} is not numeric: {text!r}") from error
    return value if math.isfinite(value) else None


def load_ranked_rows(path: Path) -> list[tuple[float, str, str]]:
    metric_columns = [f"{column}_abs" for column in RESULT_COLUMNS]
    required = {"gene", "max_abs", *metric_columns}
    ranked: list[tuple[float, str, str]] = []

    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = set(reader.fieldnames or [])
        missing = sorted(required - fields)
        if missing:
            raise ValueError(f"diagnostics table is missing required columns: {', '.join(missing)}")

        for row_number, row in enumerate(reader, start=2):
            gene = row["gene"]
            if not gene:
                raise ValueError(f"row {row_number}: gene is empty")
            max_abs = parse_finite_error(row, "max_abs", row_number)
            if max_abs is None or max_abs < 0.0:
                raise ValueError(f"row {row_number}: max_abs must be finite and non-negative")

            finite_metrics: list[tuple[float, str]] = []
            for result_column, metric_column in zip(RESULT_COLUMNS, metric_columns):
                value = parse_finite_error(row, metric_column, row_number)
                if value is not None:
                    if value < 0.0:
                        raise ValueError(
                            f"row {row_number}: {metric_column} must be non-negative"
                        )
                    finite_metrics.append((value, result_column))
            if not finite_metrics:
                raise ValueError(f"row {row_number}: no finite result-column errors")

            baseline_abs_error, result_column = max(
                finite_metrics,
                key=lambda item: (item[0], -RESULT_COLUMNS.index(item[1])),
            )
            if not math.isclose(max_abs, baseline_abs_error, rel_tol=1e-12, abs_tol=0.0):
                raise ValueError(
                    f"row {row_number}: max_abs {max_abs:.17g} does not match greatest "
                    f"finite result-column error {baseline_abs_error:.17g}"
                )
            ranked.append((max_abs, gene, result_column))

    return sorted(ranked, key=lambda item: (-item[0], item[1], item[2]))


def main() -> int:
    args = parse_args()
    ranked = load_ranked_rows(args.diagnostics)
    if len(ranked) < args.rows:
        raise ValueError(
            f"diagnostics table has {len(ranked)} rows; cannot write {args.rows}-row fixture"
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(OUTPUT_COLUMNS)
        for rank, (baseline_abs_error, gene, result_column) in enumerate(
            ranked[: args.rows], start=1
        ):
            writer.writerow([rank, gene, result_column, f"{baseline_abs_error:.17g}"])
    print(f"wrote {args.rows} rows to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
