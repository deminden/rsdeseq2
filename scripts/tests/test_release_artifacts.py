"""Publication must reject mutated packages and validation from another commit."""
from io import BytesIO
import json
import os
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import tarfile
import tempfile
import unittest
from unittest import mock

import scripts.release.release_artifact_manifest as release
from scripts.release.verify_release_crates import extract_crate


class ReleaseArtifacts(unittest.TestCase):
    def test_package_hash_and_filename_checks(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            package = root/'rsdeseq2.crate'; package.write_bytes(b'checked source')
            hashes = {package.name: release.sha256(package)}
            release.check_files(root, hashes)
            package.write_bytes(b'changed source')
            with self.assertRaisesRegex(ValueError, 'identity differs'):
                release.check_files(root, hashes)
            with self.assertRaisesRegex(ValueError, 'identity differs'):
                release.check_files(root, {'../rsdeseq2.crate':'ignored'})

    def test_validation_from_another_revision_is_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root/'validation.json').write_text(json.dumps(dict(source_commit='old',trigger_sha='new')))
            with mock.patch.object(sys,'argv',['manifest','verify',str(root)]), mock.patch.object(release,'require_frozen_tree'), mock.patch.object(release,'source_identity',return_value='new'), mock.patch.dict(os.environ, {'GITHUB_SHA':'new'}):
                with self.assertRaisesRegex(ValueError, 'different source revision'):
                    release.main()


class CrateExtraction(unittest.TestCase):
    def test_extracts_files_beneath_the_crate_name(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            archive = root / "example-1.2.3.crate"
            payload = b"#!/bin/sh\n"
            with tarfile.open(archive, "w:gz") as package:
                script = tarfile.TarInfo("example-1.2.3/tools/check.sh")
                script.size = len(payload)
                script.mode = 0o755
                package.addfile(script, BytesIO(payload))
            source = extract_crate(archive, root)
            script = source / "tools/check.sh"
            self.assertEqual(script.read_bytes(), payload)
            if os.name != "nt":
                self.assertTrue(script.stat().st_mode & 0o100)

    def test_rejects_members_outside_the_crate_and_links(self):
        for name, kind in [
            ("example-1.2.3/../outside", tarfile.REGTYPE),
            ("example-1.2.3/.. /outside", tarfile.REGTYPE),
            ("wrong-root/file", tarfile.REGTYPE),
            ("/example-1.2.3/file", tarfile.REGTYPE),
            (r"example-1.2.3/dir\file", tarfile.REGTYPE),
            ("example-1.2.3/link", tarfile.SYMTYPE),
            ("example-1.2.3/link", tarfile.LNKTYPE),
        ]:
            with self.subTest(name=name, kind=kind), tempfile.TemporaryDirectory() as temp:
                root = Path(temp)
                archive = root / "example-1.2.3.crate"
                with tarfile.open(archive, "w:gz") as package:
                    member = tarfile.TarInfo(name)
                    member.type = kind
                    member.linkname = "../outside"
                    package.addfile(member)
                with self.assertRaises(ValueError):
                    extract_crate(archive, root)
                self.assertFalse((root / "outside").exists())
