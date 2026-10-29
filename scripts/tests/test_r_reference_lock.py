"""Pinned R environments reject incomplete graphs and changed source archives."""
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import install_r_reference as reference


class RReferenceLock(unittest.TestCase):
    def test_install_order_respects_transitive_dependencies(self):
        rows = [dict(package='fit', dependencies='stats;matrix'),
                dict(package='matrix', dependencies='stats'), dict(package='stats', dependencies='')]
        self.assertEqual([r['package'] for r in reference.dependency_order(rows)], ['stats','matrix','fit'])
        for broken in ([dict(package='a', dependencies='missing')],
                       [dict(package='a', dependencies='b'),dict(package='b', dependencies='a')]):
            with self.assertRaises(ValueError):
                reference.dependency_order(broken)

    def test_changed_cached_archive_is_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            archive = Path(temp)/'package.tar.gz'
            archive.write_bytes(b'original archive')
            expected = reference.hashlib.sha256(archive.read_bytes()).hexdigest()
            reference.require_source_hash(archive, expected)
            archive.write_bytes(b'changed archive')
            with self.assertRaisesRegex(ValueError, 'hash differs'):
                reference.require_source_hash(archive, expected)

    def test_cran_archive_fallback_preserves_source_identity(self):
        row = dict(package='ashr', version='2.2-63', url='https://cran.r-project.org/src/contrib/ashr_2.2-63.tar.gz')
        self.assertEqual(reference.source_urls(row)[-1], 'https://cran.r-project.org/src/contrib/Archive/ashr/ashr_2.2-63.tar.gz')
