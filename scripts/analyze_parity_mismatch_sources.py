#!/usr/bin/env python3
"""Attribute real-data parity tail rows to upstream numeric stages.

This script consumes the diagnostics TSV produced by `scripts/real_data_parity.py`.
It does not call R; optional DESeq2-side optimizer fixture summaries are read
from generated offline fixture directories when available.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
from collections import Counter
from pathlib import Path
from typing import Iterable


DEFAULT_DIAGNOSTICS = Path(
    "results/benchmarks/real_data_parity_non_lbfgsb_start_probe_diagnostics.tsv"
)
DEFAULT_LBFGSB_ROOT = Path("results/fixtures/lbfgsb_hard_real_2026-06-01")
DEFAULT_SE_COVARIANCE_ROOT = Path("results/fixtures/se_covariance_hard_real_2026-06-05")
METRICS = [
    "baseMean_abs",
    "log2FoldChange_abs",
    "lfcSE_abs",
    "stat_abs",
    "pvalue_abs",
    "padj_abs",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--diagnostics", type=Path, default=DEFAULT_DIAGNOSTICS)
    parser.add_argument("--lbfgsb-root", type=Path, default=DEFAULT_LBFGSB_ROOT)
    parser.add_argument("--se-covariance-root", type=Path, default=DEFAULT_SE_COVARIANCE_ROOT)
    parser.add_argument("--top", type=int, default=12)
    return parser.parse_args()


def parse_float(text: str | None) -> float:
    if text is None or text == "" or text == "NA":
        return math.nan
    return float(text)


def finite(value: float) -> bool:
    return math.isfinite(value)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_optional_gene_set(path: Path) -> set[str]:
    if not path.exists():
        return set()
    with path.open(newline="") as fh:
        return {row["gene"] for row in csv.DictReader(fh, delimiter="\t")}


def read_optional_fixture_summary(root: Path, contrast: str) -> dict[str, dict[str, str]]:
    path = root / "contrasts" / contrast / "lbfgsb" / "all_gene_optimizer_summary.tsv.gz"
    if not path.exists():
        return {}
    with gzip.open(path, "rt", newline="") as fh:
        return {row["gene"]: row for row in csv.DictReader(fh, delimiter="\t")}


def optimizer_ran(row: dict[str, str], prefix: str = "") -> bool:
    return finite(parse_float(row.get(f"{prefix}rustBetaOptimIter")))


def bool_text(row: dict[str, str], column: str) -> bool:
    return row.get(column, "").lower() == "true"


def stage_label(row: dict[str, str], fixed_dispersion_se_genes: set[str]) -> str:
    any_optim = optimizer_ran(row) or optimizer_ran(row, "refit")
    replaced = bool_text(row, "replace") or bool_text(row, "refitReplace")
    beta_abs = parse_float(row.get("log2FoldChange_abs"))
    se_abs = parse_float(row.get("lfcSE_abs"))
    stat_abs = parse_float(row.get("stat_abs"))
    pvalue_abs = parse_float(row.get("pvalue_abs"))
    padj_abs = parse_float(row.get("padj_abs"))
    if replaced and any_optim:
        return "replacement optimizer/refit"
    if replaced:
        return "replacement/refit"
    if any_optim and finite(beta_abs) and beta_abs >= max(se_abs, 1e-12):
        return "optimizer beta target"
    if any_optim:
        return "optimizer-routed row"
    if finite(se_abs) and se_abs >= 1e-8 and stat_abs > beta_abs:
        if row["gene"] in fixed_dispersion_se_genes:
            return "MAP-dispersion propagated SE/stat tail"
        return "non-optimizer SE/stat tail"
    if finite(padj_abs) and padj_abs > 1e-8 and pvalue_abs < 1e-8:
        return "BH propagation"
    return "non-optimizer residual"


def top_rows(rows: Iterable[dict[str, str]], column: str, n: int) -> list[dict[str, str]]:
    return sorted(
        rows,
        key=lambda row: parse_float(row.get(column)) if finite(parse_float(row.get(column))) else -1.0,
        reverse=True,
    )[:n]


def markdown_table(headers: list[str], rows: list[list[object]]) -> None:
    print("| " + " | ".join(headers) + " |")
    print("| " + " | ".join("---" for _ in headers) + " |")
    for row in rows:
        print("| " + " | ".join(str(cell) for cell in row) + " |")


def fmt(value: float) -> str:
    if not finite(value):
        return "NA"
    return f"{value:.3e}"


def main() -> None:
    args = parse_args()
    diagnostics = read_tsv(args.diagnostics)
    fixture_cache: dict[str, dict[str, dict[str, str]]] = {}
    fixed_dispersion_se_genes = read_optional_gene_set(
        args.se_covariance_root / "selected_diagnostics.tsv"
    )

    print(f"# Parity Mismatch Source Report\n")
    print(f"Diagnostics: `{args.diagnostics}`")
    print(f"Rows inspected: {len(diagnostics)}\n")
    if fixed_dispersion_se_genes:
        print(
            "Focused fixed-dispersion SE/covariance fixture: "
            f"`{args.se_covariance_root}` ({len(fixed_dispersion_se_genes)} genes)\n"
        )

    print("## Top-Row Stage Counts\n")
    count_rows: list[list[object]] = []
    for column in ["max_abs", "log2FoldChange_abs", "lfcSE_abs", "stat_abs", "pvalue_abs", "padj_abs"]:
        selected = top_rows(diagnostics, column, min(args.top, len(diagnostics)))
        labels = Counter(stage_label(row, fixed_dispersion_se_genes) for row in selected)
        count_rows.append([column, len(selected), ", ".join(f"{key}: {value}" for key, value in labels.items())])
    markdown_table(["ranked by", "top rows", "stage labels"], count_rows)

    print("\n## Largest Rows\n")
    table_rows: list[list[object]] = []
    for row in top_rows(diagnostics, "max_abs", min(args.top, len(diagnostics))):
        contrast = row["contrast"]
        fixture_cache.setdefault(
            contrast, read_optional_fixture_summary(args.lbfgsb_root, contrast)
        )
        fixture = fixture_cache[contrast].get(row["gene"], {})
        de_optim_beta = parse_float(fixture.get("conditionBetaUseOptim"))
        replaced = bool_text(row, "replace") or bool_text(row, "refitReplace")
        rust_beta = parse_float(row.get("fitBeta") or row.get("refitBeta"))
        ref_beta = parse_float(row.get("ref_log2FoldChange"))
        rust_vs_fixture = "NA"
        ref_vs_fixture = "NA"
        if not replaced and finite(de_optim_beta) and finite(rust_beta):
            rust_vs_fixture = fmt(abs(rust_beta - de_optim_beta))
        if not replaced and finite(de_optim_beta) and finite(ref_beta):
            ref_vs_fixture = fmt(abs(ref_beta - de_optim_beta))
        table_rows.append(
            [
                row["gene"],
                row["worst_column"],
                stage_label(row, fixed_dispersion_se_genes),
                fmt(parse_float(row.get("max_abs"))),
                fmt(parse_float(row.get("log2FoldChange_abs"))),
                fmt(parse_float(row.get("lfcSE_abs"))),
                fmt(parse_float(row.get("stat_abs"))),
                "yes" if optimizer_ran(row) or optimizer_ran(row, "refit") else "no",
                rust_vs_fixture,
                ref_vs_fixture,
            ]
        )
    markdown_table(
        [
            "gene",
            "worst",
            "stage label",
            "max abs",
            "LFC abs",
            "SE abs",
            "stat abs",
            "optim",
            "Rust beta vs fixture",
            "ref beta vs fixture",
        ],
        table_rows,
    )

    print("\n## Interpretation\n")
    print(
        "- If `Rust beta vs fixture` matches the LFC error while `ref beta vs fixture` is near zero, "
        "the first visible divergence is the Rust optimizer target for that row."
    )
    print(
        "- Rows labeled `MAP-dispersion propagated SE/stat tail` have no beta-optimizer diagnostics "
        "and are present in the focused fixed-dispersion fixture. In that fixture, injecting DESeq2's "
        "final dispersions makes betaSE/covariance machine-tight, so the visible statistic gap is "
        "propagated from upstream MAP dispersion line-search sensitivity."
    )
    print(
        "- Rows labeled `non-optimizer SE/stat tail` have no beta-optimizer diagnostics but are not "
        "covered by the focused fixed-dispersion fixture; generate a targeted fixture before assigning "
        "the source to covariance arithmetic."
    )
    print(
        "- Rows labeled `BH propagation` generally inherit adjusted-p-value movement from other "
        "p-value tail rows rather than having a large local beta/statistic mismatch."
    )


if __name__ == "__main__":
    main()
