#!/usr/bin/env python3

import csv
import importlib.util
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest import mock


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "benchmark_speed_memory.py"
SPEC = importlib.util.spec_from_file_location("benchmark_speed_memory", SCRIPT)
BENCHMARK = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(BENCHMARK)


class RProvenanceTests(unittest.TestCase):
    def test_queries_requested_rscript(self):
        completed = subprocess.CompletedProcess(
            args=[],
            returncode=0,
            stdout=(
                "r_version\tR version 4.6.1 (2026-06-24)\n"
                "deseq2_version\t1.52.0\n"
                "r_platform\tx86_64-conda-linux-gnu\n"
                "blas\t/path/to/libblas.so\n"
                "lapack\t/path/to/liblapack.so\n"
            ),
            stderr="",
        )
        with mock.patch.object(BENCHMARK.subprocess, "run", return_value=completed) as run:
            provenance = BENCHMARK.query_r_provenance("/chosen/Rscript")

        self.assertEqual(run.call_args.args[0][0], "/chosen/Rscript")
        self.assertEqual(provenance["r_version"], "R version 4.6.1 (2026-06-24)")
        self.assertEqual(provenance["deseq2_version"], "1.52.0")
        self.assertEqual(provenance["r_platform"], "x86_64-conda-linux-gnu")
        self.assertEqual(provenance["blas"], "/path/to/libblas.so")
        self.assertEqual(provenance["lapack"], "/path/to/liblapack.so")

    def test_probe_failure_preserves_benchmark_compatibility(self):
        completed = subprocess.CompletedProcess(
            args=[], returncode=1, stdout="", stderr="R unavailable"
        )
        with mock.patch.object(BENCHMARK.subprocess, "run", return_value=completed):
            provenance = BENCHMARK.query_r_provenance("Rscript")

        self.assertEqual(
            provenance,
            {field: "" for field in BENCHMARK.R_PROVENANCE_FIELDS},
        )

    def test_legacy_rows_get_blank_additive_provenance_columns(self):
        row = {
            "benchmark_id": "g10_s2_size-factors_r1",
            "counts_path": "/tmp/counts.tsv",
            "tool": "DESeq2",
            "operation": "size-factors",
            "method": "ratio",
            "genes": 10,
            "samples": 2,
            "repeat": 1,
            "elapsed_s": 1.0,
            "max_rss_kb": 100,
            "max_abs_diff_vs_deseq2": 0.0,
            "status": "ok",
            "returncode": 0,
            "stderr": "",
            "command": "Rscript benchmark_deseq2.R",
        }
        with tempfile.TemporaryDirectory() as directory:
            raw = Path(directory) / "raw.tsv"
            summary = Path(directory) / "summary.tsv"
            BENCHMARK.write_rows(raw, [row])
            BENCHMARK.write_summary(summary, [row])

            with raw.open(newline="") as handle:
                raw_row = next(csv.DictReader(handle, delimiter="\t"))
            with summary.open(newline="") as handle:
                summary_row = next(csv.DictReader(handle, delimiter="\t"))

        for field in BENCHMARK.R_PROVENANCE_FIELDS:
            self.assertEqual(raw_row[field], "")
            self.assertEqual(summary_row[field], "")


if __name__ == "__main__":
    unittest.main()
