"""Every test-data file is listed in scripts/references/MANIFEST.tsv."""

import csv
import fnmatch
import unittest
from pathlib import Path

from scripts.references.stage_shrinkage_unit_fixtures import UNIT_DATA, fixture_contents

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "crates/rsdeseq2/tests/data"
MANIFEST = ROOT / "scripts/references/MANIFEST.tsv"


def load_rows():
    lines = [line for line in MANIFEST.read_text().splitlines() if not line.startswith("#")]
    return list(csv.DictReader(lines, delimiter="\t"))


def data_files():
    return sorted(
        p.relative_to(ROOT).as_posix()
        for directory in (DATA, UNIT_DATA)
        for p in directory.rglob("*")
        if p.is_file()
    )


class ReferenceManifestTest(unittest.TestCase):
    def test_every_data_file_has_a_row(self):
        rows = load_rows()
        unmatched = [f for f in data_files() if not any(fnmatch.fnmatchcase(f, r["pattern"]) for r in rows)]
        self.assertEqual(unmatched, [], "add these to scripts/references/MANIFEST.tsv")

    def test_every_row_matches_a_file(self):
        files = data_files()
        stale = [r["pattern"] for r in load_rows() if not any(fnmatch.fnmatchcase(f, r["pattern"]) for f in files)]
        self.assertEqual(stale, [])

    def test_generators_exist(self):
        for row in load_rows():
            with self.subTest(pattern=row["pattern"]):
                if row["generator"] == "-":
                    self.assertTrue(row["note"], "a missing generator needs a note")
                else:
                    self.assertTrue((ROOT / row["generator"]).is_file(), row["generator"])

    def test_shrinkage_leaves_match_the_independent_fixtures(self):
        for name, contents in fixture_contents():
            with self.subTest(fixture=name):
                self.assertEqual((UNIT_DATA / name).read_text(), contents)


if __name__ == "__main__":
    unittest.main()
