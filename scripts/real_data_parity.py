#!/usr/bin/env python3
"""Run harsh real-data parity checks against offline DESeq2 outputs.

The script expects a study directory with:

- 01_inputs/{tissue}_raw_counts.tsv.gz
- 02_deseq2_outputs/{tissue}_norm_counts.tsv.gz
- 02_deseq2_outputs/{tissue}_{contrast}_repNN_deseq2_results.tsv.gz
- 02_null_splits/{tissue}_{contrast}_repNN_groups.tsv

It never calls R. Reference outputs are treated as offline fixtures.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import os
import shutil
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np


DEFAULT_STUDY_ROOT = Path(
    os.environ.get(
        "RSDESEQ2_REAL_DATA_ROOT",
        "/home/den/bio/decor_method_study",
    )
)


@dataclass
class CommandStats:
    elapsed_s: float
    max_rss_kb: int | None
    swaps: int | None


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", newline="")
    return path.open("r", newline="")


def run_timed(cmd: list[str]) -> CommandStats:
    timed = ["/usr/bin/time", "-v", *cmd]
    start = time.monotonic()
    proc = subprocess.run(timed, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    elapsed = time.monotonic() - start
    if proc.returncode != 0:
        sys.stderr.write(proc.stdout)
        sys.stderr.write(proc.stderr)
        raise SystemExit(proc.returncode)
    max_rss = None
    swaps = None
    for line in proc.stderr.splitlines():
        line = line.strip()
        if line.startswith("Maximum resident set size"):
            max_rss = int(line.rsplit(":", 1)[1].strip())
        elif line.startswith("Swaps:"):
            swaps = int(line.rsplit(":", 1)[1].strip())
    return CommandStats(elapsed_s=elapsed, max_rss_kb=max_rss, swaps=swaps)


def read_header(path: Path) -> list[str]:
    with open_text(path) as fh:
        return next(csv.reader(fh, delimiter="\t"))


def parse_float(text: str) -> float:
    if text == "NA" or text == "":
        return math.nan
    return float(text)


def finite_abs_diff(a: float, b: float) -> float:
    if math.isnan(a) and math.isnan(b):
        return 0.0
    if math.isnan(a) or math.isnan(b):
        return math.inf
    return abs(a - b)


def infer_size_factors(raw_counts: Path, normalized_counts: Path) -> dict[str, float]:
    raw_fh = open_text(raw_counts)
    norm_fh = open_text(normalized_counts)
    try:
        raw_reader = csv.reader(raw_fh, delimiter="\t")
        norm_reader = csv.reader(norm_fh, delimiter="\t")
        raw_header = next(raw_reader)
        norm_header = next(norm_reader)
        if raw_header != norm_header:
            raise ValueError("raw and normalized-count headers differ")
        samples = raw_header[1:]
        factors: list[float | None] = [None] * len(samples)
        raw_iter = iter(raw_reader)
        raw_row = next(raw_iter, None)
        for norm_row in norm_reader:
            while raw_row is not None and raw_row[0] != norm_row[0]:
                raw_row = next(raw_iter, None)
            if raw_row is None:
                raise ValueError(f"raw counts missing normalized-count gene {norm_row[0]}")
            for i, (raw_text, norm_text) in enumerate(zip(raw_row[1:], norm_row[1:])):
                if factors[i] is not None:
                    continue
                raw = parse_float(raw_text)
                norm = parse_float(norm_text)
                if raw > 0.0 and norm > 0.0 and math.isfinite(raw) and math.isfinite(norm):
                    factors[i] = raw / norm
            if all(x is not None for x in factors):
                break
        missing = [sample for sample, factor in zip(samples, factors) if factor is None]
        if missing:
            raise ValueError(f"could not infer size factors for {len(missing)} samples")
        return {sample: float(factor) for sample, factor in zip(samples, factors)}
    finally:
        raw_fh.close()
        norm_fh.close()


def write_size_factors(path: Path, size_factors: dict[str, float], samples: Iterable[str] | None = None) -> None:
    selected = list(samples) if samples is not None else list(size_factors)
    with path.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample", "size_factor"])
        for sample in selected:
            writer.writerow([sample, f"{size_factors[sample]:.17g}"])


def compare_size_factors(got: Path, expected: dict[str, float]) -> dict[str, float | int]:
    max_abs = 0.0
    max_rel = 0.0
    n = 0
    with got.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            n += 1
            sample = row["sample"]
            size_field = "sizeFactor" if "sizeFactor" in row else "size_factor"
            diff = abs(float(row[size_field]) - expected[sample])
            max_abs = max(max_abs, diff)
            denom = max(abs(expected[sample]), 1e-300)
            max_rel = max(max_rel, diff / denom)
    return {"n": n, "max_abs": max_abs, "max_rel": max_rel}


def compare_matrix(got: Path, expected: Path) -> dict[str, float | int]:
    max_abs = 0.0
    max_rel = 0.0
    mismatched = 0
    n_values = 0
    with got.open(newline="") as got_fh, open_text(expected) as exp_fh:
        got_reader = csv.reader(got_fh, delimiter="\t")
        exp_reader = csv.reader(exp_fh, delimiter="\t")
        got_header = next(got_reader)
        exp_header = next(exp_reader)
        if got_header != exp_header:
            raise ValueError("matrix headers differ")
        got_iter = iter(got_reader)
        got_row = next(got_iter, None)
        n_rows = 0
        for exp_row in exp_reader:
            while got_row is not None and got_row[0] != exp_row[0]:
                got_row = next(got_iter, None)
            if got_row is None:
                raise ValueError(f"output missing expected gene {exp_row[0]}")
            n_rows += 1
            for got_text, exp_text in zip(got_row[1:], exp_row[1:]):
                got_value = parse_float(got_text)
                exp_value = parse_float(exp_text)
                diff = finite_abs_diff(got_value, exp_value)
                n_values += 1
                if not math.isfinite(diff):
                    mismatched += 1
                    continue
                max_abs = max(max_abs, diff)
                max_rel = max(max_rel, diff / max(abs(exp_value), 1e-300))
    return {
        "n_rows": n_rows,
        "n_values": n_values,
        "mismatched": mismatched,
        "max_abs": max_abs,
        "max_rel": max_rel,
    }


def compare_base_mean(got: Path, expected_norm: Path) -> dict[str, float | int]:
    expected: dict[str, float] = {}
    with open_text(expected_norm) as fh:
        reader = csv.reader(fh, delimiter="\t")
        next(reader)
        for row in reader:
            values = [parse_float(x) for x in row[1:]]
            expected[row[0]] = sum(values) / len(values)
    max_abs = 0.0
    max_rel = 0.0
    n = 0
    with got.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if row["gene"] not in expected:
                continue
            n += 1
            exp = expected[row["gene"]]
            base_field = "baseMean" if "baseMean" in row else "base_mean"
            diff = finite_abs_diff(float(row[base_field]), exp)
            max_abs = max(max_abs, diff)
            max_rel = max(max_rel, diff / max(abs(exp), 1e-300))
    return {"n": n, "max_abs": max_abs, "max_rel": max_rel}


def load_groups(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def retained_samples(groups: list[dict[str, str]]) -> list[str]:
    return [
        row["sample_id"]
        for row in groups
        if row.get("retained") == "TRUE" and row.get("condition") in {"A", "B"}
    ]


def design_rows(groups: list[dict[str, str]], samples: list[str], blocked: bool) -> tuple[list[str], list[list[float]], str]:
    by_sample = {row["sample_id"]: row for row in groups}
    rows = [by_sample[s] for s in samples]
    if blocked:
        blocks = sorted({row["perm_block"] for row in rows if row.get("perm_block") and row["perm_block"] != "NA"})
        columns = ["Intercept", *[f"perm_block_{block}" for block in blocks[1:]], "condition_B_vs_A"]
        matrix = []
        for row in rows:
            values = [1.0]
            values.extend(1.0 if row["perm_block"] == block else 0.0 for block in blocks[1:])
            values.append(1.0 if row["condition"] == "B" else 0.0)
            matrix.append(values)
        rank = np.linalg.matrix_rank(np.asarray(matrix, dtype=float))
        valid_blocks = all(
            {row["condition"] for row in rows if row["perm_block"] == block} == {"A", "B"}
            for block in blocks
        )
        if valid_blocks and rank == len(columns):
            return columns, matrix, "perm_block + condition"
    columns = ["Intercept", "condition_B_vs_A"]
    matrix = [[1.0, 1.0 if row["condition"] == "B" else 0.0] for row in rows]
    return columns, matrix, "condition"


def write_design(path: Path, samples: list[str], columns: list[str], matrix: list[list[float]]) -> None:
    with path.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample", *columns])
        for sample, values in zip(samples, matrix):
            writer.writerow([sample, *[f"{x:.17g}" for x in values]])


def write_retained_counts(raw_counts: Path, path: Path, samples: list[str]) -> None:
    with open_text(raw_counts) as in_fh, path.open("w", newline="") as out_fh:
        reader = csv.reader(in_fh, delimiter="\t")
        writer = csv.writer(out_fh, delimiter="\t", lineterminator="\n")
        header = next(reader)
        index = {sample: i for i, sample in enumerate(header)}
        keep = [index[sample] for sample in samples]
        writer.writerow(["gene", *samples])
        for row in reader:
            writer.writerow([row[0], *[row[i] for i in keep]])


def read_result_genes(path: Path) -> list[str]:
    with open_text(path) as fh:
        return [row["gene"] for row in csv.DictReader(fh, delimiter="\t")]


def write_retained_counts_for_genes(
    raw_counts: Path, path: Path, samples: list[str], genes: list[str]
) -> None:
    gene_set = set(genes)
    seen: set[str] = set()
    with open_text(raw_counts) as in_fh, path.open("w", newline="") as out_fh:
        reader = csv.reader(in_fh, delimiter="\t")
        writer = csv.writer(out_fh, delimiter="\t", lineterminator="\n")
        header = next(reader)
        index = {sample: i for i, sample in enumerate(header)}
        keep = [index[sample] for sample in samples]
        writer.writerow(["gene", *samples])
        for row in reader:
            gene = row[0]
            if gene in gene_set:
                writer.writerow([gene, *[row[i] for i in keep]])
                seen.add(gene)
    missing = [gene for gene in genes if gene not in seen]
    if missing:
        raise ValueError(f"raw counts missing {len(missing)} expected result genes")


def write_uncompressed_counts(raw_counts: Path, path: Path) -> None:
    with open_text(raw_counts) as in_fh, path.open("w", newline="") as out_fh:
        shutil.copyfileobj(in_fh, out_fh)


def compare_results(got: Path, expected: Path) -> dict[str, float | int]:
    columns = ["baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]
    got_by_gene: dict[str, dict[str, str]] = {}
    with got.open(newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            got_by_gene[row["gene"]] = row
    max_by_column = {col: 0.0 for col in columns}
    mismatch_by_column = {col: 0 for col in columns}
    diffs_by_column: dict[str, list[float]] = {col: [] for col in columns}
    n = 0
    with open_text(expected) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            gene = row["gene"]
            if gene not in got_by_gene:
                raise ValueError(f"missing result gene {gene}")
            n += 1
            got_row = got_by_gene[gene]
            for col in columns:
                got_value = parse_float(got_row[col])
                exp_value = parse_float(row[col])
                diff = finite_abs_diff(got_value, exp_value)
                if not math.isfinite(diff):
                    mismatch_by_column[col] += 1
                else:
                    max_by_column[col] = max(max_by_column[col], diff)
                    diffs_by_column[col].append(diff)
    out: dict[str, float | int] = {"n": n}
    for col in columns:
        diffs = sorted(diffs_by_column[col])
        out[f"{col}_max_abs"] = max_by_column[col]
        out[f"{col}_mean_abs"] = sum(diffs) / len(diffs) if diffs else math.nan
        out[f"{col}_median_abs"] = diffs[len(diffs) // 2] if diffs else math.nan
        out[f"{col}_p99_abs"] = diffs[min(len(diffs) - 1, int(0.99 * len(diffs)))] if diffs else math.nan
        out[f"{col}_mismatched"] = mismatch_by_column[col]
    return out


def append_result_diagnostics(
    diagnostics: list[dict[str, object]],
    contrast: str,
    got: Path,
    expected: Path,
    replacement_rows: Path | None,
    fit_diagnostics: Path | None,
    refit_diagnostics: Path | None,
    fit_beta: Path | None,
    fit_beta_se: Path | None,
    fit_beta_optim_start: Path | None,
    refit_beta: Path | None,
    refit_beta_se: Path | None,
    refit_beta_optim_start: Path | None,
    coefficient_index: int,
    limit: int,
) -> None:
    if limit <= 0:
        return
    columns = ["baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]
    expected_by_gene: dict[str, dict[str, str]] = {}
    with open_text(expected) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            expected_by_gene[row["gene"]] = row
    replacement_by_gene: dict[str, dict[str, str]] = {}
    if replacement_rows is not None and replacement_rows.exists():
        with replacement_rows.open(newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                replacement_by_gene[row["gene"]] = row
    fit_by_gene = read_optional_gene_rows(fit_diagnostics)
    refit_by_gene = read_optional_gene_rows(refit_diagnostics)
    fit_beta_by_gene = read_optional_gene_matrix_column(fit_beta, coefficient_index)
    fit_beta_se_by_gene = read_optional_gene_matrix_column(fit_beta_se, coefficient_index)
    fit_beta_optim_start_by_gene = read_optional_gene_matrix_column(fit_beta_optim_start, coefficient_index)
    refit_beta_by_gene = read_optional_gene_matrix_column(refit_beta, coefficient_index)
    refit_beta_se_by_gene = read_optional_gene_matrix_column(refit_beta_se, coefficient_index)
    refit_beta_optim_start_by_gene = read_optional_gene_matrix_column(
        refit_beta_optim_start, coefficient_index
    )
    rows: list[tuple[float, dict[str, object], dict[str, float]]] = []
    with got.open(newline="") as fh:
        for got_row in csv.DictReader(fh, delimiter="\t"):
            gene = got_row["gene"]
            exp_row = expected_by_gene[gene]
            diffs = {
                col: finite_abs_diff(parse_float(got_row[col]), parse_float(exp_row[col]))
                for col in columns
            }
            finite_diffs = [value for value in diffs.values() if math.isfinite(value)]
            max_diff = max(finite_diffs) if finite_diffs else math.inf
            worst_column = max(columns, key=lambda col: diffs[col])
            replacement = replacement_by_gene.get(gene, {})
            fit = fit_by_gene.get(gene, {})
            refit = refit_by_gene.get(gene, {})
            rows.append((
                max_diff,
                {
                    "contrast": contrast,
                    "gene": gene,
                    "worst_column": worst_column,
                    "max_abs": max_diff,
                    **{f"{col}_abs": diffs[col] for col in columns},
                    **{f"rust_{col}": got_row[col] for col in columns},
                    **{f"ref_{col}": exp_row[col] for col in columns},
                    "replace": replacement.get("replace", ""),
                    "refitReplace": replacement.get("refitReplace", ""),
                    "newAllZero": replacement.get("newAllZero", ""),
                    "postRefitMaxCooks": replacement.get("postRefitMaxCooks", ""),
                    "dispGeneEst": fit.get("dispGeneEst", ""),
                    "dispGeneIter": fit.get("dispGeneIter", ""),
                    "dispFit": fit.get("dispFit", ""),
                    "dispMAP": fit.get("dispMAP", ""),
                    "dispersion": fit.get("dispersion", ""),
                    "dispIter": fit.get("dispIter", ""),
                    "dispOutlier": fit.get("dispOutlier", ""),
                    "betaConv": fit.get("betaConv", ""),
                    "fullBetaConv": fit.get("fullBetaConv", ""),
                    "reducedBetaConv": fit.get("reducedBetaConv", ""),
                    "betaIter": fit.get("betaIter", ""),
                    "rustBetaOptimIter": fit.get("rustBetaOptimIter", ""),
                    "rustBetaOptimStartObjective": fit.get("rustBetaOptimStartObjective", ""),
                    "rustBetaOptimObjective": fit.get("rustBetaOptimObjective", ""),
                    "rustBetaOptimGradientNorm": fit.get("rustBetaOptimGradientNorm", ""),
                    "reducedBetaIter": fit.get("reducedBetaIter", ""),
                    "deviance": fit.get("deviance", ""),
                    "maxCooks": fit.get("maxCooks", ""),
                    "refitDispGeneEst": refit.get("dispGeneEst", ""),
                    "refitDispGeneIter": refit.get("dispGeneIter", ""),
                    "refitDispFit": refit.get("dispFit", ""),
                    "refitDispMAP": refit.get("dispMAP", ""),
                    "refitDispersion": refit.get("dispersion", ""),
                    "refitDispIter": refit.get("dispIter", ""),
                    "refitDispOutlier": refit.get("dispOutlier", ""),
                    "refitBetaConv": refit.get("betaConv", ""),
                    "refitFullBetaConv": refit.get("fullBetaConv", ""),
                    "refitReducedBetaConv": refit.get("reducedBetaConv", ""),
                    "refitBetaIter": refit.get("betaIter", ""),
                    "refitRustBetaOptimIter": refit.get("rustBetaOptimIter", ""),
                    "refitRustBetaOptimStartObjective": refit.get("rustBetaOptimStartObjective", ""),
                    "refitRustBetaOptimObjective": refit.get("rustBetaOptimObjective", ""),
                    "refitRustBetaOptimGradientNorm": refit.get("rustBetaOptimGradientNorm", ""),
                    "refitReducedBetaIter": refit.get("reducedBetaIter", ""),
                    "refitDeviance": refit.get("deviance", ""),
                    "refitMaxCooks": refit.get("maxCooks", ""),
                    "fitBeta": fit_beta_by_gene.get(gene, ""),
                    "fitBetaSE": fit_beta_se_by_gene.get(gene, ""),
                    "fitBetaOptimStart": fit_beta_optim_start_by_gene.get(gene, ""),
                    "refitBeta": refit_beta_by_gene.get(gene, ""),
                    "refitBetaSE": refit_beta_se_by_gene.get(gene, ""),
                    "refitBetaOptimStart": refit_beta_optim_start_by_gene.get(gene, ""),
                },
                diffs,
            ))
    selected: dict[tuple[str, str], dict[str, object]] = {}
    for _, row, _ in sorted(rows, key=lambda item: item[0], reverse=True)[:limit]:
        selected[(str(row["contrast"]), str(row["gene"]))] = row
    for column in columns:
        column_rows = sorted(
            rows,
            key=lambda item: item[2][column] if math.isfinite(item[2][column]) else math.inf,
            reverse=True,
        )
        for _, row, _ in column_rows[:limit]:
            selected[(str(row["contrast"]), str(row["gene"]))] = row
    diagnostics.extend(selected.values())


def read_optional_gene_rows(path: Path | None) -> dict[str, dict[str, str]]:
    rows: dict[str, dict[str, str]] = {}
    if path is None or not path.exists():
        return rows
    with path.open(newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            rows[row["gene"]] = row
    return rows


def read_optional_gene_matrix_column(path: Path | None, column_index: int) -> dict[str, str]:
    rows: dict[str, str] = {}
    if path is None or not path.exists():
        return rows
    value_index = column_index + 1
    with path.open(newline="") as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if not row or row[0] == "gene":
                continue
            if value_index < len(row):
                rows[row[0]] = row[value_index]
    return rows


def write_summary(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        return
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def primitive_run(args, binary: Path, tissue: str, tmp: Path, rows: list[dict[str, object]]) -> None:
    raw = args.study_root / "01_inputs" / f"{tissue}_raw_counts.tsv.gz"
    norm = args.study_root / "02_deseq2_outputs" / f"{tissue}_norm_counts.tsv.gz"
    size_factors = infer_size_factors(raw, norm)
    raw_tsv = tmp / f"{tissue}_raw_counts.tsv"
    write_uncompressed_counts(raw, raw_tsv)
    size_path = tmp / f"{tissue}_size_factors.tsv"
    write_size_factors(size_path, size_factors)

    got_sf = tmp / f"{tissue}_sf_out.tsv"
    stats = run_timed([str(binary), "size-factors", "--counts", str(raw_tsv), "--output", str(got_sf)])
    cmp = compare_size_factors(got_sf, size_factors)
    rows.append({
        "tissue": tissue,
        "output": "size_factors",
        "elapsed_s": f"{stats.elapsed_s:.6g}",
        "max_rss_kb": stats.max_rss_kb,
        "swaps": stats.swaps,
        **cmp,
    })

    got_norm = tmp / f"{tissue}_normalized.tsv"
    stats = run_timed([
        str(binary),
        "normalized-counts",
        "--counts",
        str(raw_tsv),
        "--size-factors",
        str(size_path),
        "--output",
        str(got_norm),
    ])
    cmp = compare_matrix(got_norm, norm)
    rows.append({
        "tissue": tissue,
        "output": "normalized_counts",
        "elapsed_s": f"{stats.elapsed_s:.6g}",
        "max_rss_kb": stats.max_rss_kb,
        "swaps": stats.swaps,
        **cmp,
    })

    got_base = tmp / f"{tissue}_base.tsv"
    stats = run_timed([
        str(binary),
        "base-mean",
        "--counts",
        str(raw_tsv),
        "--size-factors",
        str(size_path),
        "--output",
        str(got_base),
    ])
    cmp = compare_base_mean(got_base, norm)
    rows.append({
        "tissue": tissue,
        "output": "base_mean",
        "elapsed_s": f"{stats.elapsed_s:.6g}",
        "max_rss_kb": stats.max_rss_kb,
        "swaps": stats.swaps,
        **cmp,
    })


def contrast_run(
    args,
    binary: Path,
    spec: str,
    tmp: Path,
    rows: list[dict[str, object]],
    diagnostics: list[dict[str, object]],
) -> None:
    parts = spec.split(":")
    if len(parts) != 3:
        raise ValueError("--contrast must be tissue:null_type:rep")
    tissue, null_type, rep_text = parts
    rep = int(rep_text)
    stem = f"{tissue}_{null_type}_rep{rep:02d}"
    raw = args.study_root / "01_inputs" / f"{tissue}_raw_counts.tsv.gz"
    norm = args.study_root / "02_deseq2_outputs" / f"{tissue}_norm_counts.tsv.gz"
    expected = args.study_root / "02_deseq2_outputs" / f"{stem}_deseq2_results.tsv.gz"
    groups_path = args.study_root / "02_null_splits" / f"{stem}_groups.tsv"
    groups = load_groups(groups_path)
    samples = retained_samples(groups)
    columns, matrix, design_kind = design_rows(groups, samples, "blocked_permutation" in null_type)
    counts_path = tmp / f"{stem}_counts.tsv"
    design_path = tmp / f"{stem}_design.tsv"
    got = tmp / f"{stem}_wald.tsv"
    replacement_rows = tmp / f"{stem}_replacement_rows.tsv"
    replacement_meta = tmp / f"{stem}_replacement_meta.tsv"
    fit_diagnostics = tmp / f"{stem}_fit_diagnostics.tsv"
    refit_diagnostics = tmp / f"{stem}_refit_diagnostics.tsv"
    fit_beta = tmp / f"{stem}_fit_beta.tsv"
    fit_beta_se = tmp / f"{stem}_fit_beta_se.tsv"
    fit_beta_optim_start = tmp / f"{stem}_fit_beta_optim_start.tsv"
    refit_beta = tmp / f"{stem}_refit_beta.tsv"
    refit_beta_se = tmp / f"{stem}_refit_beta_se.tsv"
    refit_beta_optim_start = tmp / f"{stem}_refit_beta_optim_start.tsv"
    expected_genes = read_result_genes(expected)
    write_retained_counts_for_genes(raw, counts_path, samples, expected_genes)
    write_design(design_path, samples, columns, matrix)
    cmd = [
        str(binary),
        "wald",
        "--counts",
        str(counts_path),
        "--design",
        str(design_path),
        "--fit-type",
        args.fit_type,
        "--coefficient",
        str(len(columns) - 1),
        "--output",
        str(got),
    ]
    if args.diagnostics_output is not None:
        cmd.extend([
            "--cooks-replacement-row-metadata-output",
            str(replacement_rows),
            "--cooks-replacement-metadata-output",
            str(replacement_meta),
            "--fit-diagnostics-output",
            str(fit_diagnostics),
            "--refit-diagnostics-output",
            str(refit_diagnostics),
            "--fit-beta-output",
            str(fit_beta),
            "--fit-beta-se-output",
            str(fit_beta_se),
            "--fit-beta-optim-start-output",
            str(fit_beta_optim_start),
            "--refit-beta-output",
            str(refit_beta),
            "--refit-beta-se-output",
            str(refit_beta_se),
            "--refit-beta-optim-start-output",
            str(refit_beta_optim_start),
        ])
    if args.contrast_size_factors == "full":
        size_path = tmp / f"{stem}_size_factors.tsv"
        write_size_factors(size_path, infer_size_factors(raw, norm), samples)
        cmd.extend(["--size-factors", str(size_path)])
    elif args.contrast_size_factors != "estimate":
        raise ValueError(f"unknown contrast size-factor mode: {args.contrast_size_factors}")
    stats = run_timed(cmd)
    cmp = compare_results(got, expected)
    append_result_diagnostics(
        diagnostics,
        stem,
        got,
        expected,
        replacement_rows if args.diagnostics_output is not None else None,
        fit_diagnostics if args.diagnostics_output is not None else None,
        refit_diagnostics if args.diagnostics_output is not None else None,
        fit_beta if args.diagnostics_output is not None else None,
        fit_beta_se if args.diagnostics_output is not None else None,
        fit_beta_optim_start if args.diagnostics_output is not None else None,
        refit_beta if args.diagnostics_output is not None else None,
        refit_beta_se if args.diagnostics_output is not None else None,
        refit_beta_optim_start if args.diagnostics_output is not None else None,
        len(columns) - 1,
        args.diagnostics_limit,
    )
    rows.append({
        "contrast": stem,
        "output": "wald_results",
        "design": design_kind,
        "size_factor_mode": args.contrast_size_factors,
        "n_samples": len(samples),
        "elapsed_s": f"{stats.elapsed_s:.6g}",
        "max_rss_kb": stats.max_rss_kb,
        "swaps": stats.swaps,
        **cmp,
    })


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--study-root", type=Path, default=DEFAULT_STUDY_ROOT)
    parser.add_argument("--binary", type=Path, default=Path("target/release/rsdeseq2"))
    parser.add_argument("--tissue", action="append", default=[])
    parser.add_argument("--contrast", action="append", default=[])
    parser.add_argument("--contrast-size-factors", choices=["estimate", "full"], default="estimate")
    parser.add_argument("--fit-type", default="parametric")
    parser.add_argument("--output", type=Path, default=Path("results/benchmarks/real_data_parity.tsv"))
    parser.add_argument("--diagnostics-output", type=Path)
    parser.add_argument("--diagnostics-limit", type=int, default=100)
    parser.add_argument("--keep-workdir", action="store_true")
    args = parser.parse_args()

    if not args.binary.exists():
        resolved = shutil.which(str(args.binary))
        if resolved:
            args.binary = Path(resolved)
        else:
            raise SystemExit(f"missing binary: {args.binary}")
    rows: list[dict[str, object]] = []
    diagnostics: list[dict[str, object]] = []
    tmp = Path(tempfile.mkdtemp(prefix="rsdeseq2-real-parity-"))
    try:
        for tissue in args.tissue:
            primitive_run(args, args.binary, tissue, tmp, rows)
        for spec in args.contrast:
            contrast_run(args, args.binary, spec, tmp, rows, diagnostics)
        write_summary(args.output, rows)
        print(f"wrote {args.output}")
        if args.diagnostics_output is not None:
            write_summary(args.diagnostics_output, diagnostics)
            print(f"wrote {args.diagnostics_output}")
    finally:
        if args.keep_workdir:
            print(f"kept workdir {tmp}")
        else:
            shutil.rmtree(tmp)


if __name__ == "__main__":
    main()
