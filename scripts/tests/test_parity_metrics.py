#!/usr/bin/env python3
"""Tests for significance-aware parity metrics."""

import importlib.util
import math
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location('parity_metrics', ROOT / 'scripts/validation/parity_metrics.py')
METRICS = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(METRICS)


def table(rows):
    return {gene: {'pvalue': str(p), 'padj': str(q), 'stat': str(s)} for gene, p, q, s in rows}


class ParityMetricTests(unittest.TestCase):
    def test_identical_tables_have_no_flips_and_full_overlap(self):
        rows = table([(f'g{i}', 10 ** -(i % 7), 10 ** -(i % 5), i - 50) for i in range(120)])
        metrics = METRICS.significance_metrics(rows, rows)
        self.assertEqual(metrics['flips_padj_0.05'], 0)
        self.assertEqual(metrics['top100_overlap'], 1.0)
        self.assertAlmostEqual(metrics['stat_spearman'], 1.0)
        self.assertEqual(metrics['pvalue_small_max_log10'], 0.0)
        METRICS.require_no_flips(metrics, 'identical')

    def test_flips_and_log10_differences_are_counted(self):
        expected = table([('a', 1e-6, 0.04, 5), ('b', 0.5, 0.6, 0.1), ('c', 1e-4, 'NA', 3)])
        got = table([('a', 1e-5, 0.06, 4), ('b', 0.5, 0.6, 0.2), ('c', 1e-4, 0.01, 3)])
        metrics = METRICS.significance_metrics(got, expected)
        self.assertEqual(metrics['flips_padj_0.05'], 2)  # a loses, c gains
        self.assertEqual(metrics['flips_padj_0.1'], 1)
        self.assertAlmostEqual(metrics['pvalue_small_max_log10'], 1.0)
        with self.assertRaises(SystemExit):
            METRICS.require_no_flips(metrics, 'test')

    def test_missing_reference_gene_is_an_error(self):
        with self.assertRaises(ValueError):
            METRICS.significance_metrics(table([('a', 1, 1, 0)]), table([('b', 1, 1, 0)]))

    def test_spearman_handles_ties_and_nan(self):
        self.assertTrue(math.isnan(METRICS.spearman([1.0], [2.0])))
        self.assertAlmostEqual(METRICS.spearman([1, 2, 2, 3], [10, 20, 20, 30]), 1.0)


if __name__ == '__main__':
    unittest.main()
