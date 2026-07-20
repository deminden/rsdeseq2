#!/usr/bin/env python3
"""Run publication-scope real-data parity jobs with bounded concurrency."""

from __future__ import annotations

import argparse
import csv
import os
import re
import signal
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path


CONTRAST_RE = re.compile(r"cameraPR_(.+)_go_bp[.]manifest[.]json$")


@dataclass(frozen=True)
class Task:
    kind: str
    name: str
    args: list[str]


@dataclass
class TaskResult:
    task: Task
    output: Path
    diagnostics: Path | None
    elapsed_s: float
    returncode: int
    stderr_tail: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--study-root",
        type=Path,
        default=Path("/home/den/bio/rsfgsea/results/decor_method_study"),
    )
    parser.add_argument(
        "--manifest-dir",
        type=Path,
        default=Path(
            "/home/den/bio/rsfgsea/results/decor_method_study/"
            "decor_publication_dossier/15_external_method_comparison/manifests/cameraPR"
        ),
    )
    parser.add_argument("--binary", type=Path, default=Path("target/release/rsdeseq2"))
    parser.add_argument("--driver", type=Path, default=Path("scripts/real_data_parity.py"))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--diagnostics-output", type=Path, required=True)
    parser.add_argument("--failures-output", type=Path, required=True)
    parser.add_argument(
        "--failure-list",
        type=Path,
        help="TSV produced by this script; rerun the failed task names from it.",
    )
    parser.add_argument("--workdir", type=Path)
    parser.add_argument("--rsdeseq2-workers", type=int, default=9)
    parser.add_argument(
        "--rayon-num-threads",
        type=int,
        help="Set RAYON_NUM_THREADS for each rsdeseq2 subprocess.",
    )
    parser.add_argument("--task-timeout-s", type=float, default=7200.0)
    parser.add_argument("--diagnostics-limit", type=int, default=20)
    parser.add_argument(
        "--max-tasks",
        type=int,
        help="Run only the first N selected tasks, useful for smoke batches.",
    )
    parser.add_argument(
        "--task-offset",
        type=int,
        default=0,
        help="Skip the first N selected tasks before applying --max-tasks.",
    )
    parser.add_argument("--tissue-name", action="append", default=[])
    parser.add_argument("--contrast-name", action="append", default=[])
    parser.add_argument(
        "--contrast-size-factors",
        choices=["estimate", "full"],
        default="estimate",
    )
    parser.add_argument("--skip-tissues", action="store_true")
    parser.add_argument("--skip-contrasts", action="store_true")
    return parser.parse_args()


def discover_tissues(study_root: Path) -> list[str]:
    input_dir = study_root / "01_inputs"
    tissues = []
    for path in input_dir.glob("*_raw_counts.tsv.gz"):
        tissues.append(path.name.removesuffix("_raw_counts.tsv.gz"))
    return sorted(tissues)


def discover_contrasts(manifest_dir: Path) -> list[str]:
    contrasts: set[str] = set()
    for path in manifest_dir.glob("*.json"):
        match = CONTRAST_RE.match(path.name)
        if match:
            contrasts.add(match.group(1))
    return sorted(contrasts)


def contrast_spec(name: str) -> str:
    match = re.match(r"^(.+)_full_blocked_permutation_rep(\d+)$", name)
    if not match:
        raise ValueError(f"unsupported publication contrast name: {name}")
    tissue, rep = match.groups()
    return f"{tissue}:full_blocked_permutation:{int(rep)}"


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def tasks_from_failure_list(path: Path) -> list[Task]:
    tasks = []
    for row in read_tsv(path):
        kind = row.get("kind", "")
        name = row.get("name", "")
        if kind == "tissue":
            tasks.append(Task("tissue", name, ["--tissue", name]))
        elif kind == "contrast":
            tasks.append(Task("contrast", name, ["--contrast", contrast_spec(name)]))
    return tasks


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="") as fh:
        if not fields:
            return
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def run_task(task: Task, args: argparse.Namespace, workdir: Path) -> TaskResult:
    output = workdir / f"{task.kind}_{task.name}.tsv"
    diagnostics = workdir / f"{task.kind}_{task.name}_diagnostics.tsv"
    cmd = [
        sys.executable,
        str(args.driver),
        "--study-root",
        str(args.study_root),
        "--binary",
        str(args.binary),
        "--output",
        str(output),
        *task.args,
    ]
    if task.kind == "contrast":
        cmd.extend(
            [
                "--contrast-size-factors",
                args.contrast_size_factors,
                "--diagnostics-output",
                str(diagnostics),
                "--diagnostics-limit",
                str(args.diagnostics_limit),
            ]
        )
    else:
        diagnostics = None
    start = time.monotonic()
    env = os.environ.copy()
    if args.rayon_num_threads is not None:
        env["RAYON_NUM_THREADS"] = str(args.rayon_num_threads)
    proc = subprocess.Popen(
        cmd,
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        start_new_session=True,
    )
    try:
        _, stderr = proc.communicate(timeout=args.task_timeout_s)
        returncode = proc.returncode
    except subprocess.TimeoutExpired as error:
        try:
            os.killpg(proc.pid, signal.SIGTERM)
        except ProcessLookupError:
            pass
        _, timeout_stderr = proc.communicate()
        returncode = 124
        stderr = (error.stderr or "") + (timeout_stderr or "")
        stderr += f"\npublication parity task timed out after {args.task_timeout_s:.1f}s"
    elapsed_s = time.monotonic() - start
    stderr_tail = "\n".join(stderr.splitlines()[-40:])
    return TaskResult(task, output, diagnostics, elapsed_s, returncode, stderr_tail)


def main() -> None:
    args = parse_args()
    workdir = args.workdir or Path(tempfile.mkdtemp(prefix="rsdeseq2-publication-parity-"))
    workdir.mkdir(parents=True, exist_ok=True)

    tasks: list[Task] = []
    if args.failure_list is not None:
        tasks.extend(tasks_from_failure_list(args.failure_list))
    if not args.skip_tissues and args.failure_list is None:
        tissues = sorted(set(args.tissue_name)) if args.tissue_name else discover_tissues(args.study_root)
        for tissue in tissues:
            tasks.append(Task("tissue", tissue, ["--tissue", tissue]))
    if not args.skip_contrasts and args.failure_list is None:
        contrasts = (
            sorted(set(args.contrast_name))
            if args.contrast_name
            else discover_contrasts(args.manifest_dir)
        )
        for contrast in contrasts:
            tasks.append(Task("contrast", contrast, ["--contrast", contrast_spec(contrast)]))
    if args.task_offset:
        tasks = tasks[args.task_offset :]
    if args.max_tasks is not None:
        tasks = tasks[: args.max_tasks]

    print(
        f"running {len(tasks)} publication parity tasks "
        f"with {args.rsdeseq2_workers} rsdeseq2 workers; workdir={workdir}",
        flush=True,
    )

    results: list[TaskResult] = []
    with ThreadPoolExecutor(max_workers=args.rsdeseq2_workers) as executor:
        future_map = {executor.submit(run_task, task, args, workdir): task for task in tasks}
        for future in as_completed(future_map):
            result = future.result()
            results.append(result)
            status = "ok" if result.returncode == 0 else f"failed:{result.returncode}"
            print(
                f"{status}\t{result.task.kind}\t{result.task.name}\t{result.elapsed_s:.1f}s",
                flush=True,
            )

    rows: list[dict[str, object]] = []
    diagnostics: list[dict[str, object]] = []
    failures: list[dict[str, object]] = []
    for result in sorted(results, key=lambda item: (item.task.kind, item.task.name)):
        if result.returncode == 0:
            rows.extend(read_tsv(result.output))
            if result.diagnostics is not None:
                diagnostics.extend(read_tsv(result.diagnostics))
        else:
            failures.append(
                {
                    "kind": result.task.kind,
                    "name": result.task.name,
                    "elapsed_s": f"{result.elapsed_s:.6g}",
                    "returncode": result.returncode,
                    "stderr_tail": result.stderr_tail,
                }
            )

    write_tsv(args.output, rows)
    write_tsv(args.diagnostics_output, diagnostics)
    write_tsv(args.failures_output, failures)
    print(f"wrote {len(rows)} rows to {args.output}")
    print(f"wrote {len(diagnostics)} diagnostics rows to {args.diagnostics_output}")
    print(f"wrote {len(failures)} failure rows to {args.failures_output}")


if __name__ == "__main__":
    main()
