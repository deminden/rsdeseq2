"""Reject broken navigation and documentation unsuitable for release downloads."""
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts.release.prepare_release_docs import write_binary_readme, write_release_notes
from scripts.validation.check_docs import ROOT, check_markdown, source_version


class DocumentationChecks(unittest.TestCase):
    def test_local_paths_references_and_anchors_ignore_code(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root / 'other page.md').write_text(
                '# `Rust` Usage\n# Repeat\n# Repeat\n<a id="manual"></a>\n', encoding='utf-8')
            page = root / 'README.md'
            page.write_text('''# Home
[usage](<other page.md#rust-usage>)
[repeat][topic]
[topic]: other%20page.md#repeat-1
[manual](other%20page.md#manual)
[home](#home)
`[not a link](absent.md)`
```md
[not a reference][absent]
[not a link](absent.md)
```
~~~md
# Not a heading
~~~
''', encoding='utf-8')
            self.assertEqual([], check_markdown(page, root))

    def test_missing_paths_anchors_references_and_escape_are_rejected(self):
        cases = (
            ('[x](absent.md)', 'missing target'),
            ('[x](#absent)', 'missing anchor'),
            ('[x][absent]', 'undefined reference'),
            ('[x][]', 'undefined reference'),
            ('[x](../outside.md)', 'escapes documentation tree'),
        )
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            page = root / 'README.md'
            for content, diagnostic in cases:
                with self.subTest(content=content):
                    page.write_text(content, encoding='utf-8')
                    self.assertTrue(any(diagnostic in error for error in check_markdown(page, root)))

    def test_package_links_must_be_contained_or_pinned_to_matching_tag(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            package = root / 'package'
            package.mkdir()
            (root / 'guide.md').write_text('# Guide\n', encoding='utf-8')
            page = package / 'README.md'
            for target, diagnostic in (
                ('../guide.md', 'escapes documentation tree'),
                ('https://github.com/deminden/rsdeseq2/blob/main/guide.md', 'not pinned'),
                ('https://github.com/deminden/rsdeseq2/blob/v0.4.0/guide.md', 'not pinned'),
                ('https://github.com/deminden/rsdeseq2/blob/v0.5.0/absent.md', 'missing target'),
            ):
                with self.subTest(target=target):
                    page.write_text(f'[guide]({target})', encoding='utf-8')
                    errors = check_markdown(page, root, package_root=package, tag='v0.5.0')
                    self.assertTrue(any(diagnostic in error for error in errors))
            page.write_text('[guide](https://github.com/deminden/rsdeseq2/blob/v0.5.0/guide.md#guide)', encoding='utf-8')
            self.assertEqual([], check_markdown(page, root, package_root=package, tag='v0.5.0'))

    def test_generated_download_docs_and_release_notes_work_outside_checkout(self):
        tag = 'v' + source_version(ROOT)
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            for notice in ('LICENSE', 'THIRD_PARTY_NOTICES.md'):
                (directory / notice).write_bytes((ROOT / notice).read_bytes())
            for target, command in (
                ('x86_64-unknown-linux-gnu', './rsdeseq2 --version'),
                ('aarch64-apple-darwin', './rsdeseq2 --version'),
                ('x86_64-pc-windows-msvc', r'.\rsdeseq2.exe --version'),
            ):
                with self.subTest(target=target):
                    write_binary_readme(ROOT, directory, tag, target)
                    text = (directory / 'README.md').read_text(encoding='utf-8')
                    self.assertIn(command, text)
                    self.assertNotIn('cargo install', text)
                    self.assertIn(f'/blob/{tag}/docs/cli.md', text)
            output = directory / 'release-notes.md'
            write_release_notes(ROOT, output, tag)
            text = output.read_text(encoding='utf-8')
            self.assertNotIn('](../', text)
            self.assertIn(f'/blob/{tag}/docs/testing-parity.md#ci-status', text)
            self.assertEqual([], check_markdown(output, ROOT, package_root=directory, tag=tag))

    def test_generation_rejects_wrong_version_and_missing_packaged_notice(self):
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            with self.assertRaisesRegex(ValueError, 'match the source version'):
                write_binary_readme(ROOT, directory, 'v0.0.0', 'x86_64-unknown-linux-gnu')
            self.assertFalse((directory / 'README.md').exists())
            with self.assertRaisesRegex(ValueError, 'missing target.*LICENSE'):
                write_binary_readme(ROOT, directory, 'v' + source_version(ROOT), 'x86_64-unknown-linux-gnu')


if __name__ == '__main__':
    unittest.main()
