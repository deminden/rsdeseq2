#!/usr/bin/env python3

import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
GENERATOR = ROOT / "scripts" / "generate_frozen_worst_genes.py"
RESULT_COLUMNS = ["baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]


class FrozenWorstGeneGeneratorTests(unittest.TestCase):
    def run_generator(
        self,
        rows: list[dict[str, object]],
        count: int = 100,
        fieldnames: list[str] | None = None,
    ) -> tuple[subprocess.CompletedProcess[str], list[dict[str, str]]]:
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            diagnostics = directory / "diagnostics.tsv"
            output = directory / "nested" / "fixture.tsv"
            fields = fieldnames or [
                "gene",
                "max_abs",
                *[f"{column}_abs" for column in RESULT_COLUMNS],
            ]
            with diagnostics.open("w", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    delimiter="\t",
                    fieldnames=fields,
                    lineterminator="\n",
                )
                writer.writeheader()
                writer.writerows(rows)
            result = subprocess.run(
                [
                    sys.executable,
                    str(GENERATOR),
                    "--diagnostics",
                    str(diagnostics),
                    "--output",
                    str(output),
                    "--rows",
                    str(count),
                ],
                text=True,
                capture_output=True,
                check=False,
            )
            generated: list[dict[str, str]] = []
            if output.exists():
                with output.open(newline="") as handle:
                    generated = list(csv.DictReader(handle, delimiter="\t"))
            return result, generated

    @staticmethod
    def diagnostic(gene: str, errors: list[object]) -> dict[str, object]:
        finite = [float(value) for value in errors if str(value).lower() not in {"nan", "inf"}]
        row: dict[str, object] = {"gene": gene, "max_abs": max(finite)}
        row.update(
            {f"{column}_abs": value for column, value in zip(RESULT_COLUMNS, errors)}
        )
        return row

    def test_sorts_rows_and_selects_greatest_finite_result_error(self):
        rows = [
            self.diagnostic("gene_z", [0.1, 0.2, 0.3, 0.4, "nan", "inf"]),
            self.diagnostic("gene_b", [0.5, 0.1, 0.1, 0.1, 0.1, 0.1]),
            self.diagnostic("gene_a", [0.1, 0.5, 0.1, 0.1, 0.1, 0.1]),
        ]
        result, generated = self.run_generator(rows, count=3)

        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(
            generated,
            [
                {
                    "rank": "1",
                    "gene": "gene_a",
                    "result_column": "log2FoldChange",
                    "baseline_abs_error": "0.5",
                },
                {
                    "rank": "2",
                    "gene": "gene_b",
                    "result_column": "baseMean",
                    "baseline_abs_error": "0.5",
                },
                {
                    "rank": "3",
                    "gene": "gene_z",
                    "result_column": "stat",
                    "baseline_abs_error": "0.40000000000000002",
                },
            ],
        )

    def test_rejects_missing_gene_or_max_abs_column(self):
        for missing in ("gene", "max_abs"):
            with self.subTest(missing=missing):
                fields = [
                    "gene",
                    "max_abs",
                    *[f"{column}_abs" for column in RESULT_COLUMNS],
                ]
                fields.remove(missing)
                result, generated = self.run_generator([], count=1, fieldnames=fields)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn(f"missing required columns: {missing}", result.stderr)
                self.assertEqual(generated, [])

    def test_rejects_max_abs_inconsistent_with_result_errors(self):
        row = self.diagnostic("gene_a", [0.1, 0.2, 0.3, 0.4, 0.5, 0.6])
        row["max_abs"] = 0.7
        result, generated = self.run_generator([row], count=1)

        self.assertNotEqual(result.returncode, 0)
        self.assertIn("does not match greatest finite result-column error", result.stderr)
        self.assertEqual(generated, [])

    def test_rejects_request_larger_than_diagnostics(self):
        row = self.diagnostic("gene_a", [0.1, 0.2, 0.3, 0.4, 0.5, 0.6])
        result, generated = self.run_generator([row], count=2)

        self.assertNotEqual(result.returncode, 0)
        self.assertIn("cannot write 2-row fixture", result.stderr)
        self.assertEqual(generated, [])


if __name__ == "__main__":
    unittest.main()
