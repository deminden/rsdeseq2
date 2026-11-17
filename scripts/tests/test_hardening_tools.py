#!/usr/bin/env python3
"""Tests for the hardening metrics and paired timing helpers."""

import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def load(name, relative):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


METRICS = load('hardening_metrics', 'scripts/validation/hardening_metrics.py')
TIMING = load('paired_gtex_rust_timing', 'scripts/benchmarks/paired_gtex_rust_timing.py')


class StripTests(unittest.TestCase):
    def test_comments_and_strings_do_not_count(self):
        text = 'let a = x.unwrap(); // y.unwrap()\n/* z.unwrap() */ let s = "q.unwrap()";\n'
        stripped = METRICS.strip_comments_and_strings(text)
        self.assertEqual(stripped.count('.unwrap()'), 1)
        self.assertEqual(stripped.count('\n'), text.count('\n'))

    def test_raw_strings_and_chars_are_blanked(self):
        text = 'let r = r#"a.unwrap() "quoted""#; let c = \'{\'; f.expect("x");\n'
        stripped = METRICS.strip_comments_and_strings(text)
        self.assertNotIn('a.unwrap()', stripped)
        self.assertEqual(stripped.count('.expect('), 1)
        self.assertNotIn("'{'", stripped)

    def test_cfg_test_items_are_removed(self):
        text = 'fn prod() { a.unwrap(); }\n#[cfg(test)]\nmod tests {\n fn t() { b.unwrap(); { } }\n}\nfn after() {}\n'
        stripped = METRICS.strip_cfg_test(text)
        self.assertEqual(stripped.count('.unwrap()'), 1)
        self.assertIn('fn after', stripped)
        self.assertEqual(stripped.count('\n'), text.count('\n'))

    def test_cfg_test_use_declaration_is_removed(self):
        stripped = METRICS.strip_cfg_test('#[cfg(test)]\nuse std::fmt;\nfn keep() {}\n')
        self.assertNotIn('std::fmt', stripped)
        self.assertIn('fn keep', stripped)


class ModuleCycleTests(unittest.TestCase):
    def test_mutual_references_are_pairs(self):
        src = METRICS.CRATE / 'src'
        sources = {
            src / 'a.rs': 'use crate::b::X;',
            src / 'b' / 'mod.rs': 'fn f() { crate::a::g(); }',
            src / 'c.rs': 'use crate::a::Y;',
        }
        self.assertEqual(METRICS.mutual_pairs(sources), [('a', 'b')])


class TimingTests(unittest.TestCase):
    def test_contrast_layout_requires_eleven_directories(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for name in TIMING.CONTRASTS:
                (root / name).mkdir()
            for index in range(1, 8):
                (root / 'random_splits' / f'blocked_random_rep{index:02d}').mkdir(parents=True)
            with self.assertRaises(SystemExit):
                TIMING.contrast_dirs(root)
            (root / 'random_splits' / 'unblocked_random_rep01').mkdir()
            self.assertEqual(len(TIMING.contrast_dirs(root)), 11)

    def test_measured_command_reports_failure(self):
        with tempfile.TemporaryDirectory() as tmp:
            log = Path(tmp) / 'log.txt'
            elapsed, rss = TIMING.run_measured([sys.executable, '-c', 'pass'], {}, log)
            self.assertGreaterEqual(elapsed, 0.0)
            self.assertGreater(rss, 0)
            with self.assertRaises(SystemExit):
                TIMING.run_measured([sys.executable, '-c', 'raise SystemExit(3)'], {}, log)


if __name__ == '__main__':
    unittest.main()
