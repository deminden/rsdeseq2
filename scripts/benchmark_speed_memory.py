#!/usr/bin/env python3
"""Benchmark current rsdeseq2 primitives against DESeq2 for time and max RSS."""

from __future__ import annotations

import argparse
import csv
import math
import random
import statistics
import subprocess
import tempfile
import time
from pathlib import Path


def parse_int_list(value: str) -> list[int]:
    return [int(part.strip()) for part in value.split(",") if part.strip()]


def parse_operation_list(value: str) -> list[str]:
    operations = [part.strip() for part in value.split(",") if part.strip()]
    allowed = {"size-factors", "base-mean"}
    unknown = sorted(set(operations) - allowed)
    if unknown:
        raise argparse.ArgumentTypeError(f"unsupported operations: {', '.join(unknown)}")
    return operations


def write_counts(path: Path, genes: int, samples: int, seed: int) -> None:
    rng = random.Random(seed + genes * 1009 + samples * 917)
    sample_scale = [rng.lognormvariate(0.0, 0.25) for _ in range(samples)]
    scale_center = statistics.geometric_mean(sample_scale)
    sample_scale = [value / scale_center for value in sample_scale]

    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["gene", *[f"s{idx + 1}" for idx in range(samples)]])
        for gene in range(genes):
            baseline = rng.lognormvariate(3.2, 1.0)
            row = []
            for sample in range(samples):
                mean = baseline * sample_scale[sample]
                noise = rng.gauss(0.0, math.sqrt(mean + 1.0) + mean * 0.12)
                value = max(0, int(round(mean + noise)))
                if rng.random() < 0.01:
                    value = 0
                row.append(value)
            if not any(row):
                row[rng.randrange(samples)] = 1
            writer.writerow([f"g{gene + 1}", *row])


def count_matrix_shape(path: Path) -> tuple[int, int]:
    with path.open(newline="") as handle:
        header = handle.readline()
        if not header:
            raise ValueError(f"empty count table: {path}")
        columns = len(header.rstrip("\n\r").split("\t"))
        if columns < 2:
            raise ValueError(f"count table must have at least two columns: {path}")
        rows = sum(1 for _ in handle)
    return rows, columns - 1


def parse_elapsed_seconds(value: str) -> float:
    parts = value.strip().split(":")
    if len(parts) == 3:
        hours, minutes, seconds = parts
        return int(hours) * 3600.0 + int(minutes) * 60.0 + float(seconds)
    if len(parts) == 2:
        minutes, seconds = parts
        return int(minutes) * 60.0 + float(seconds)
    return float(value)


def parse_time_metrics(path: Path) -> tuple[float | None, int | None]:
    elapsed = None
    max_rss = None
    for line in path.read_text().splitlines():
        line = line.strip()
        if line.startswith("Elapsed (wall clock) time"):
            value = line.split("):", 1)[1].strip() if "):" in line else line.rsplit(":", 1)[1].strip()
            elapsed = parse_elapsed_seconds(value)
        elif line.startswith("Maximum resident set size"):
            max_rss = int(line.rsplit(":", 1)[1].strip())
    return elapsed, max_rss


def run_timed(time_bin: str, command: list[str], cwd: Path, metrics_path: Path) -> dict[str, object]:
    started = time.perf_counter()
    proc = subprocess.run(
        [time_bin, "-v", "-o", str(metrics_path), *command],
        cwd=cwd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    wall_elapsed = time.perf_counter() - started
    elapsed, max_rss = (None, None)
    if metrics_path.exists():
        elapsed, max_rss = parse_time_metrics(metrics_path)
    if elapsed is None or elapsed == 0.0:
        elapsed = wall_elapsed
    return {
        "status": "ok" if proc.returncode == 0 else "failed",
        "returncode": proc.returncode,
        "elapsed_s": elapsed,
        "max_rss_kb": max_rss,
        "stderr": proc.stderr.strip().replace("\n", " ")[:300],
        "command": " ".join(command),
    }


def read_numeric_output(path: Path) -> list[float]:
    values = []
    with path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        next(reader)
        for row in reader:
            values.append(float(row[1]))
    return values


def max_abs_diff(lhs: Path, rhs: Path) -> float | None:
    if not lhs.exists() or not rhs.exists():
        return None
    left = read_numeric_output(lhs)
    right = read_numeric_output(rhs)
    if len(left) != len(right):
        return None
    return max((abs(a - b) for a, b in zip(left, right)), default=0.0)


def median_or_blank(values: list[float | int | None]) -> str:
    present = [float(value) for value in values if value not in (None, "")]
    if not present:
        return ""
    return f"{statistics.median(present):.6g}"


def min_or_blank(values: list[float | int | None]) -> str:
    present = [float(value) for value in values if value not in (None, "")]
    if not present:
        return ""
    return f"{min(present):.6g}"


def max_or_blank(values: list[float | int | None]) -> str:
    present = [float(value) for value in values if value not in (None, "")]
    if not present:
        return ""
    return f"{max(present):.6g}"


def mad_or_blank(values: list[float | int | None]) -> str:
    present = [float(value) for value in values if value not in (None, "")]
    if not present:
        return ""
    median = statistics.median(present)
    return f"{statistics.median(abs(value - median) for value in present):.6g}"


def write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "benchmark_id",
        "counts_path",
        "tool",
        "operation",
        "method",
        "genes",
        "samples",
        "repeat",
        "elapsed_s",
        "max_rss_kb",
        "max_abs_diff_vs_deseq2",
        "status",
        "returncode",
        "stderr",
        "command",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_summary(path: Path, rows: list[dict[str, object]]) -> None:
    groups: dict[tuple[object, ...], list[dict[str, object]]] = {}
    for row in rows:
        key = (
            row["tool"],
            row["operation"],
            row["method"],
            row["genes"],
            row["samples"],
        )
        groups.setdefault(key, []).append(row)

    fieldnames = [
        "tool",
        "operation",
        "method",
        "genes",
        "samples",
        "runs",
        "ok_runs",
        "median_elapsed_s",
        "min_elapsed_s",
        "max_elapsed_s",
        "mad_elapsed_s",
        "median_max_rss_kb",
        "min_max_rss_kb",
        "max_max_rss_kb",
        "mad_max_rss_kb",
        "max_abs_diff_vs_deseq2",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for key, group_rows in sorted(groups.items()):
            diff_values = [
                row["max_abs_diff_vs_deseq2"]
                for row in group_rows
                if row["max_abs_diff_vs_deseq2"] != ""
            ]
            elapsed_values = [row["elapsed_s"] for row in group_rows]
            rss_values = [row["max_rss_kb"] for row in group_rows]
            writer.writerow(
                {
                    "tool": key[0],
                    "operation": key[1],
                    "method": key[2],
                    "genes": key[3],
                    "samples": key[4],
                    "runs": len(group_rows),
                    "ok_runs": sum(row["status"] == "ok" for row in group_rows),
                    "median_elapsed_s": median_or_blank(elapsed_values),
                    "min_elapsed_s": min_or_blank(elapsed_values),
                    "max_elapsed_s": max_or_blank(elapsed_values),
                    "mad_elapsed_s": mad_or_blank(elapsed_values),
                    "median_max_rss_kb": median_or_blank(rss_values),
                    "min_max_rss_kb": min_or_blank(rss_values),
                    "max_max_rss_kb": max_or_blank(rss_values),
                    "mad_max_rss_kb": mad_or_blank(rss_values),
                    "max_abs_diff_vs_deseq2": median_or_blank(diff_values),
                }
            )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genes", default="1000,10000", help="comma-separated gene counts")
    parser.add_argument("--samples", default="8,16", help="comma-separated sample counts")
    parser.add_argument(
        "--counts-file",
        type=Path,
        help="existing count TSV with gene IDs in the first column; skips synthetic generation",
    )
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--method", choices=["ratio", "poscounts"], default="ratio")
    parser.add_argument(
        "--operations",
        type=parse_operation_list,
        default=parse_operation_list("size-factors,base-mean"),
    )
    parser.add_argument("--seed", type=int, default=20260523)
    parser.add_argument("--rscript", default="Rscript")
    parser.add_argument("--time-bin", default="/usr/bin/time")
    parser.add_argument("--no-build", action="store_true")
    parser.add_argument(
        "--output",
        default="results/benchmarks/speed_memory.tsv",
        type=Path,
    )
    args = parser.parse_args()

    repo = Path(__file__).resolve().parents[1]
    if not args.no_build:
        subprocess.run(["cargo", "build", "-p", "rsdeseq2", "--release"], cwd=repo, check=True)

    binary = repo / "target" / "release" / "rsdeseq2"
    deseq2_script = repo / "scripts" / "benchmark_deseq2.R"
    rows: list[dict[str, object]] = []

    with tempfile.TemporaryDirectory(prefix="rsdeseq2-bench-") as tmp:
        tmpdir = Path(tmp)
        if args.counts_file is not None:
            counts = args.counts_file.resolve()
            genes, samples = count_matrix_shape(counts)
            datasets = [(counts.stem, counts, genes, samples)]
        else:
            datasets = []
            for genes in parse_int_list(args.genes):
                for samples in parse_int_list(args.samples):
                    counts = tmpdir / f"counts_g{genes}_s{samples}.tsv"
                    write_counts(counts, genes, samples, args.seed)
                    datasets.append((f"g{genes}_s{samples}", counts, genes, samples))

        for dataset_name, counts, genes, samples in datasets:
            safe_dataset_name = "".join(
                char if char.isalnum() or char in "._-" else "_"
                for char in dataset_name
            )
            for operation in args.operations:
                for repeat in range(1, args.repeats + 1):
                    benchmark_id = f"{safe_dataset_name}_{operation}_r{repeat}"
                    rust_out = tmpdir / f"{benchmark_id}_rust.tsv"
                    deseq2_out = tmpdir / f"{benchmark_id}_deseq2.tsv"

                    rust_command = [
                        str(binary),
                        operation,
                        "--counts",
                        str(counts),
                        "--method",
                        args.method,
                        "--output",
                        str(rust_out),
                    ]
                    rust_metrics = tmpdir / f"{benchmark_id}_rust.time"
                    rust = run_timed(args.time_bin, rust_command, repo, rust_metrics)

                    deseq2_command = [
                        args.rscript,
                        str(deseq2_script),
                        "--operation",
                        operation,
                        "--counts",
                        str(counts),
                        "--method",
                        args.method,
                        "--output",
                        str(deseq2_out),
                    ]
                    deseq2_metrics = tmpdir / f"{benchmark_id}_deseq2.time"
                    deseq2 = run_timed(args.time_bin, deseq2_command, repo, deseq2_metrics)

                    diff = max_abs_diff(rust_out, deseq2_out)
                    for tool, result in [("rsdeseq2", rust), ("DESeq2", deseq2)]:
                        rows.append(
                            {
                                "benchmark_id": benchmark_id,
                                "counts_path": str(counts),
                                "tool": tool,
                                "operation": operation,
                                "method": args.method,
                                "genes": genes,
                                "samples": samples,
                                "repeat": repeat,
                                "elapsed_s": (
                                    result["elapsed_s"]
                                    if result["elapsed_s"] is not None
                                    else ""
                                ),
                                "max_rss_kb": (
                                    result["max_rss_kb"]
                                    if result["max_rss_kb"] is not None
                                    else ""
                                ),
                                "max_abs_diff_vs_deseq2": diff if diff is not None else "",
                                "status": result["status"],
                                "returncode": result["returncode"],
                                "stderr": result["stderr"],
                                "command": result["command"],
                            }
                        )

    write_rows(args.output, rows)
    summary = args.output.with_name(args.output.stem + "_summary.tsv")
    write_summary(summary, rows)
    print(f"wrote {args.output}")
    print(f"wrote {summary}")


if __name__ == "__main__":
    main()
