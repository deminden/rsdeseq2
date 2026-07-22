#!/usr/bin/env python3
"""Score a parity diagnostics table against a frozen worst-gene baseline."""

from __future__ import annotations

import argparse
import csv
import math
import statistics
from collections import Counter
from pathlib import Path


FIXTURE_COLUMNS = {"rank", "gene", "result_column", "baseline_abs_error"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fixture", type=Path, required=True)
    parser.add_argument("--diagnostics", type=Path, required=True)
    parser.add_argument("--min-median-improvement", type=float, default=10.0)
    parser.add_argument(
        "--max-error-ratio",
        type=float,
        default=1.0,
        help="maximum candidate/baseline worst-error ratio (default: 1.0)",
    )
    parser.add_argument(
        "--report-only",
        action="store_true",
        help="print measurements without applying pass/fail comparisons",
    )
    return parser.parse_args()


def require_columns(
    fieldnames: list[str] | None,
    required: set[str],
    table_name: str,
) -> None:
    columns = fieldnames or []
    duplicates = sorted(column for column, count in Counter(columns).items() if count > 1)
    if duplicates:
        raise ValueError(
            f"{table_name} table has duplicate columns: {', '.join(duplicates)}"
        )
    missing = sorted(required - set(columns))
    if missing:
        raise ValueError(
            f"{table_name} table is missing required columns: {', '.join(missing)}"
        )


def index_unique_genes(
    rows: list[dict[str, str]],
    table_name: str,
) -> dict[str, dict[str, str]]:
    indexed: dict[str, dict[str, str]] = {}
    first_rows: dict[str, int] = {}
    for row_number, row in enumerate(rows, start=2):
        gene = row["gene"]
        if not gene:
            raise ValueError(f"{table_name} row {row_number}: gene is empty")
        if gene in indexed:
            raise ValueError(
                f"{table_name} table has duplicate gene {gene!r} at rows "
                f"{first_rows[gene]} and {row_number}"
            )
        indexed[gene] = row
        first_rows[gene] = row_number
    return indexed


def parse_non_negative_error(text: str | None, description: str) -> float:
    if text is None:
        raise ValueError(f"{description} is not numeric: None")
    try:
        value = float(text)
    except ValueError as error:
        raise ValueError(f"{description} is not numeric: {text!r}") from error
    if not math.isfinite(value) or value < 0.0:
        raise ValueError(f"{description} must be finite and non-negative")
    return value


def format_missing_genes(genes: list[str]) -> str:
    shown = genes[:10]
    suffix = (
        f", ... (+{len(genes) - len(shown)} more)"
        if len(genes) > len(shown)
        else ""
    )
    return ", ".join(shown) + suffix


def main() -> int:
    args = parse_args()
    with args.fixture.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require_columns(reader.fieldnames, FIXTURE_COLUMNS, "fixture")
        baseline = list(reader)

    if len(baseline) != 100:
        raise ValueError(f"expected 100 frozen rows, found {len(baseline)}")
    baseline_by_gene = index_unique_genes(baseline, "fixture")

    result_columns: set[str] = set()
    for row_number, frozen in enumerate(baseline, start=2):
        result_column = frozen["result_column"]
        if not result_column:
            raise ValueError(f"fixture row {row_number}: result_column is empty")
        result_columns.add(f"{result_column}_abs")

    with args.diagnostics.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require_columns(reader.fieldnames, {"gene", *result_columns}, "diagnostics")
        diagnostics = list(reader)

    if not diagnostics:
        raise ValueError("diagnostics table is empty")
    candidate_by_gene = index_unique_genes(diagnostics, "diagnostics")

    missing_genes = sorted(set(baseline_by_gene) - set(candidate_by_gene))
    if missing_genes:
        noun = "gene" if len(missing_genes) == 1 else "genes"
        raise ValueError(
            f"diagnostics table is missing {len(missing_genes)} frozen {noun}: "
            f"{format_missing_genes(missing_genes)}"
        )

    old_errors: list[float] = []
    new_errors: list[float] = []
    pair_improvements: list[float] = []

    for row_number, frozen in enumerate(baseline, start=2):
        gene = frozen["gene"]
        old = parse_non_negative_error(
            frozen["baseline_abs_error"],
            f"fixture row {row_number}: baseline_abs_error for frozen gene {gene!r}",
        )
        column = f"{frozen['result_column']}_abs"
        new = parse_non_negative_error(
            candidate_by_gene[gene][column],
            f"diagnostics {column} for frozen gene {gene!r}",
        )
        old_errors.append(old)
        new_errors.append(new)
        pair_improvements.append(old / max(new, 1e-300))

    old_median = statistics.median(old_errors)
    new_median = statistics.median(new_errors)
    median_improvement = old_median / max(new_median, 1e-300)
    mean_improvement = statistics.mean(old_errors) / max(statistics.mean(new_errors), 1e-300)
    improved = sum(new < old for old, new in zip(old_errors, new_errors))
    at_least_ten = sum(ratio >= 10.0 for ratio in pair_improvements)

    print(f"frozen_rows\t{len(baseline)}")
    print(f"candidate_rows_present\t{len(baseline)}")
    print(f"baseline_median_abs_error\t{old_median:.17g}")
    print(f"candidate_median_abs_error\t{new_median:.17g}")
    print(f"median_improvement\t{median_improvement:.9g}")
    print(f"mean_improvement\t{mean_improvement:.9g}")
    print(f"rows_improved\t{improved}")
    print(f"rows_improved_at_least_10x\t{at_least_ten}")
    print(f"baseline_max_abs_error\t{max(old_errors):.17g}")
    print(f"candidate_max_abs_error\t{max(new_errors):.17g}")
    max_error_ratio = max(new_errors) / max(max(old_errors), 1e-300)
    print(f"max_error_ratio\t{max_error_ratio:.9g}")

    if args.report_only:
        return 0

    failed = False
    if median_improvement < args.min_median_improvement:
        print(
            f"FAIL: median improvement {median_improvement:.6g}x is below "
            f"{args.min_median_improvement:.6g}x"
        )
        failed = True
    if max_error_ratio > args.max_error_ratio:
        print(
            f"FAIL: maximum error ratio {max_error_ratio:.6g} exceeds "
            f"{args.max_error_ratio:.6g}"
        )
        failed = True
    if failed:
        return 1
    print(
        f"PASS: median improvement is at least {args.min_median_improvement:g}x "
        f"and maximum error ratio is at most {args.max_error_ratio:g}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
