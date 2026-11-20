"""Every test-data file is listed in scripts/references/MANIFEST.tsv."""

import fnmatch
import hashlib
import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from scripts.references.list_generators import DIRECTORY, MANIFEST, generator_rows, read_rows
from scripts.references.stage_shrinkage_unit_fixtures import UNIT_DATA, fixture_contents

DATA = ROOT / "crates/rsdeseq2/tests/data"


def data_files():
    return sorted(
        p.relative_to(ROOT).as_posix()
        for directory in (DATA, UNIT_DATA)
        for p in directory.rglob("*")
        if p.is_file()
    )


class ReferenceManifestTest(unittest.TestCase):
    def test_every_data_file_has_a_row(self):
        rows = read_rows(MANIFEST)
        unmatched = [f for f in data_files() if not any(fnmatch.fnmatchcase(f, r["pattern"]) for r in rows)]
        self.assertEqual(unmatched, [], "add these to scripts/references/MANIFEST.tsv")

    def test_every_row_matches_a_file(self):
        files = data_files()
        stale = [r["pattern"] for r in read_rows(MANIFEST) if not any(fnmatch.fnmatchcase(f, r["pattern"]) for f in files)]
        self.assertEqual(stale, [])

    def test_files_without_generators_have_an_explanation(self):
        for row in read_rows(MANIFEST):
            with self.subTest(pattern=row["pattern"]):
                if row["generator"] == "-":
                    self.assertTrue(row["note"], "a missing generator needs a note")

    def test_catalog_lists_each_reference_generator_with_output_and_check(self):
        entries = generator_rows()
        registered = {row["generator"] for row in entries}
        expected = {
            path.relative_to(ROOT).as_posix()
            for pattern in ("generate_*", "stage_*", "freeze_*")
            for path in DIRECTORY.glob(pattern)
            if path.is_file() and path.suffix in (".R", ".py")
        }
        self.assertEqual(expected - registered, set())
        for row in entries:
            with self.subTest(generator=row["generator"], output=row["pattern"]):
                self.assertTrue((ROOT / row["generator"]).is_file())
                self.assertTrue(row["pattern"])
                self.assertTrue(row["check"])
                self.assertNotEqual(row["check"], "-")

    def test_shrinkage_leaves_match_the_independent_fixtures(self):
        for name, contents in fixture_contents():
            with self.subTest(fixture=name):
                self.assertEqual((UNIT_DATA / name).read_text(), contents)

    def test_student_posterior_reference_uses_the_current_inputs(self):
        source = DATA / "student_uniform_posterior_r462.json"
        fixture = json.loads((DATA / "student_uniform_mathematical.json").read_text())
        self.assertEqual(fixture["source"]["inputs"], source.name)
        self.assertEqual(
            fixture["source"]["input_sha256"], hashlib.sha256(source.read_bytes()).hexdigest(),
        )


if __name__ == "__main__":
    unittest.main()
