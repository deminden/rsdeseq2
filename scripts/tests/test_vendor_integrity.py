"""Vendor drift must fail without fixing the checkout being checked."""
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import sync_r_vendor as vendor


class VendorIntegrity(unittest.TestCase):
    def test_installed_notices_are_complete_and_read_only_check_rejects_drift(self):
        with tempfile.TemporaryDirectory() as temp:
            source, destination = Path(temp)/'NOTICE', Path(temp)/'COPYRIGHTS'
            source.write_text('Required third-party notice\n')
            for content in (None, 'incomplete notice\n'):
                if content is not None: destination.write_text(content)
                before = destination.read_bytes() if destination.exists() else None
                with self.assertRaisesRegex(ValueError, 'third-party notices differ'):
                    vendor.sync_notice(source, destination, check=True)
                self.assertEqual(before, destination.read_bytes() if destination.exists() else None)
            vendor.sync_notice(source, destination)
            vendor.sync_notice(source, destination, check=True)
            self.assertEqual(source.read_bytes(), destination.read_bytes())

    def test_check_rejects_missing_extra_and_changed_sources_without_writing(self):
        for mutation in ('missing', 'extra', 'changed'):
            with self.subTest(mutation=mutation), tempfile.TemporaryDirectory() as temp:
                expected, actual = Path(temp)/'expected', Path(temp)/'actual'
                expected.mkdir(); actual.mkdir()
                (expected/'lib.rs').write_text('reference')
                (actual/'lib.rs').write_text('reference')
                vendor.check_tree(expected, actual)
                if mutation == 'missing': (actual/'lib.rs').unlink()
                elif mutation == 'extra': (actual/'extra.rs').write_text('obsolete')
                else: (actual/'lib.rs').write_text('stale')
                before = vendor.tree_files(actual)
                with self.assertRaisesRegex(ValueError, 'R vendor differs'):
                    vendor.check_tree(expected, actual)
                self.assertEqual(before, vendor.tree_files(actual))

    def test_mismatched_package_versions_fail_before_generation(self):
        for r_version, backend_version in [('0.3.2','0.3.1'), ('0.3.1','0.3.2')]:
            with self.subTest(r_version=r_version, backend_version=backend_version), tempfile.TemporaryDirectory() as temp:
                root = Path(temp)
                description = root/'r-pkg/rsdeseq2/DESCRIPTION'
                description.parent.mkdir(parents=True)
                description.write_text('Version: '+r_version+'\n')
                packages = [[dict(name='rsdeseq2',version='0.3.1')], [dict(name='rsdeseq2-r',version=backend_version)]]
                with mock.patch.object(vendor, 'ROOT', root), mock.patch.object(vendor, 'metadata', side_effect=packages), mock.patch.object(sys, 'argv', ['sync', '--check']), mock.patch.object(vendor.shutil, 'copytree') as copy:
                    with self.assertRaisesRegex(ValueError, 'versions must agree'):
                        vendor.main()
                    copy.assert_not_called()
