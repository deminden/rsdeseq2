"""Pinned R environments reject incomplete graphs and changed source archives."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import json
import os
import tempfile
import unittest
from unittest import mock

import scripts.environment.install_r_reference as reference


class RReferenceLock(unittest.TestCase):
    def test_check_only_uses_requested_library_without_dropping_existing_paths(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            manifest = root/'sources.tsv'
            manifest.write_text('package\tversion\turl\tsha256\tdependencies\n' +
                'brio\t1.1.5\thttps://example.org/brio.tar.gz\t'+'0'*64+'\t\n')
            environment = root/'environment.json'
            environment.write_text(json.dumps(dict(r_version='4.6.2')))
            library = root/'requested library'
            args = ['install_r_reference','--check-installed','--rscript','Rscript',
                    '--manifest',str(manifest),'--environment',str(environment),
                    '--library',str(library),'--cache',str(root/'cache')]
            with mock.patch.object(sys,'argv',args), mock.patch.dict(os.environ,{'R_LIBS_USER':'existing-library'}), mock.patch.object(reference.subprocess,'run') as run:
                reference.main()
            self.assertEqual(run.call_args.kwargs['env']['R_LIBS_USER'],
                             str(library.resolve())+os.pathsep+'existing-library')
            self.assertEqual(run.call_args.args[0][-1], 'check')
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

    def test_mirror_fallback_still_requires_the_locked_archive_hash(self):
        expected = b'locked package source'
        row = dict(package='pkgbuild', version='1.4.8',
                   url='https://cran.r-project.org/src/contrib/pkgbuild_1.4.8.tar.gz',
                   sha256=reference.hashlib.sha256(expected).hexdigest())
        for payload in [expected, b'changed package source']:
            with self.subTest(payload=payload), tempfile.TemporaryDirectory() as temp:
                destination = Path(temp)/'source.tar.gz'
                def download(command, **kwargs):
                    if command[-1] == row['url']:
                        raise reference.subprocess.CalledProcessError(22, command)
                    Path(command[command.index('--output')+1]).write_bytes(payload)
                with mock.patch.object(reference.subprocess, 'run', side_effect=download) as run:
                    if payload == expected:
                        reference.download_source(row, destination)
                        self.assertEqual(destination.read_bytes(), expected)
                    else:
                        with self.assertRaisesRegex(ValueError, 'hash differs'):
                            reference.download_source(row, destination)
                        self.assertFalse(destination.exists())
                self.assertEqual(run.call_count, 2)
                self.assertIn('cloud.r-project.org', run.call_args.args[0][-1])
