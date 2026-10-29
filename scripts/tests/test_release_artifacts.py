"""Publication must reject mutated packages and validation from another commit."""
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import release_artifact_manifest as release


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
