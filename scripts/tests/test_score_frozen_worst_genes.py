#!/usr/bin/env python3

import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCORER = ROOT / "scripts" / "score_frozen_worst_genes.py"


class FrozenWorstGeneScorerTests(unittest.TestCase):
    def run_case(
        self,
        candidate_error: float,
        minimum: float,
        worst_error: float | None = None,
        report_only: bool = False,
        fixture_rows: list[list[object]] | None = None,
        diagnostic_rows: list[list[object]] | None = None,
        fixture_header: list[str] | None = None,
        diagnostics_header: list[str] | None = None,
    ) -> subprocess.CompletedProcess[str]:
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            fixture = directory / "fixture.tsv"
            diagnostics = directory / "diagnostics.tsv"
            with fixture.open("w", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(
                    fixture_header
                    or ["rank", "gene", "result_column", "baseline_abs_error"]
                )
                writer.writerows(
                    fixture_rows
                    or [[rank, f"g{rank}", "stat", "1"] for rank in range(1, 101)]
                )
            with diagnostics.open("w", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(diagnostics_header or ["gene", "max_abs", "stat_abs"])
                writer.writerows(
                    diagnostic_rows
                    or [
                        [
                            f"g{rank}",
                            worst_error
                            if rank == 1 and worst_error is not None
                            else candidate_error,
                            worst_error
                            if rank == 1 and worst_error is not None
                            else candidate_error,
                        ]
                        for rank in range(1, 101)
                    ]
                )
            command = [
                sys.executable,
                str(SCORER),
                "--fixture",
                str(fixture),
                "--diagnostics",
                str(diagnostics),
                "--min-median-improvement",
                str(minimum),
            ]
            if report_only:
                command.append("--report-only")
            return subprocess.run(
                command,
                text=True,
                capture_output=True,
                check=False,
            )

    def test_passes_at_threshold(self):
        result = self.run_case(0.1, 10.0)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn(
            "PASS: median improvement is at least 10x and maximum error ratio is at most 1",
            result.stdout,
        )

    def test_fails_below_threshold(self):
        result = self.run_case(0.2, 10.0)
        self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
        self.assertIn("FAIL: median improvement 5x", result.stdout)

    def test_fails_when_worst_error_regresses(self):
        result = self.run_case(0.01, 10.0, worst_error=1.01)
        self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
        self.assertIn("FAIL: maximum error ratio 1.01 exceeds 1", result.stdout)

    def test_report_only_returns_measurements_without_comparison(self):
        result = self.run_case(0.2, 10.0, worst_error=1.01, report_only=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("candidate_median_abs_error\t0.2", result.stdout)
        self.assertNotIn("PASS:", result.stdout)
        self.assertNotIn("FAIL:", result.stdout)

    def test_rejects_truncated_diagnostics_instead_of_estimating_cutoff(self):
        rows = [[f"g{rank}", 0.01, 0.01] for rank in range(1, 100)]
        # Put a large final value in an intentionally unsorted table. The old
        # scorer silently charged this value as a cutoff for absent g100.
        rows[-1] = ["g99", 50.0, 50.0]
        result = self.run_case(0.01, 10.0, diagnostic_rows=rows)

        self.assertNotEqual(result.returncode, 0)
        self.assertIn("diagnostics table is missing 1 frozen gene: g100", result.stderr)

    def test_rejects_duplicate_diagnostics_genes(self):
        rows = [[f"g{rank}", 0.1, 0.1] for rank in range(1, 101)]
        rows[-1][0] = "g1"
        result = self.run_case(0.1, 10.0, diagnostic_rows=rows)

        self.assertNotEqual(result.returncode, 0)
        self.assertIn("diagnostics table has duplicate gene 'g1'", result.stderr)

    def test_rejects_duplicate_fixture_genes(self):
        rows = [[rank, f"g{rank}", "stat", "1"] for rank in range(1, 101)]
        rows[-1][1] = "g1"
        result = self.run_case(0.1, 10.0, fixture_rows=rows)

        self.assertNotEqual(result.returncode, 0)
        self.assertIn("fixture table has duplicate gene 'g1'", result.stderr)

    def test_rejects_missing_required_columns(self):
        result = self.run_case(
            0.1,
            10.0,
            diagnostics_header=["gene", "max_abs"],
            diagnostic_rows=[[f"g{rank}", 0.1] for rank in range(1, 101)],
        )

        self.assertNotEqual(result.returncode, 0)
        self.assertIn(
            "diagnostics table is missing required columns: stat_abs", result.stderr
        )

    def test_rejects_duplicate_columns(self):
        result = self.run_case(
            0.1,
            10.0,
            diagnostics_header=["gene", "stat_abs", "stat_abs"],
        )

        self.assertNotEqual(result.returncode, 0)
        self.assertIn("diagnostics table has duplicate columns: stat_abs", result.stderr)


if __name__ == "__main__":
    unittest.main()
