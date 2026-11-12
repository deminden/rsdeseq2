#!/usr/bin/env python3
"""Check repository Markdown and generated release documentation without network access.

Supports the repository's inline/reference links, ATX headings and HTML anchors.
External sites and historical revisions are not fetched. Current main/dev/tag
links into this repository are resolved against the checkout instead.
"""
import argparse
from pathlib import Path
import re
import subprocess
import sys
import tempfile
from urllib.parse import unquote, urlsplit

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

ROOT = Path(__file__).resolve().parents[2]
REPOSITORY = 'https://github.com/deminden/rsdeseq2'
LINK = re.compile(r'\]\((<[^>\n]+>|[^\s)]+)(?:\s+["\'][^\n]*?["\'])?\)')
REFERENCE = re.compile(r'^ {0,3}\[([^]\n]+)\]:\s*(<[^>\n]+>|\S+)', re.M)


def prose(text, *, strip_inline=True):
    """Exclude fenced and inline code from link parsing, preserving line breaks."""
    lines = []
    fence = None
    for line in text.splitlines():
        match = re.match(r'^ {0,3}(`{3,}|~{3,})', line)
        if fence:
            if re.fullmatch(r' {0,3}' + re.escape(fence[0]) + '{' + str(len(fence)) + r',}\s*', line):
                fence = None
            lines.append('')
        elif match:
            fence = match[1]
            lines.append('')
        else:
            lines.append(line)
    text = '\n'.join(lines)
    return re.sub(r'(`+).*?\1', '', text) if strip_inline else text


def anchors(text):
    # Heading code spans contribute their text to the GitHub anchor.
    text = prose(text, strip_inline=False)
    result = set(re.findall(r'(?:id|name)=["\']([^"\']+)["\']', text))
    counts = {}
    for match in re.finditer(r'^ {0,3}#{1,6}\s+(.+?)\s*#*$', text, re.M):
        label = re.sub(r'\[([^]]+)\]\([^)]*\)', r'\1', match[1])
        label = re.sub(r'<[^>]+>', '', label)
        slug = re.sub(r'[^\w\- ]', '', label.lower()).replace(' ', '-')
        suffix = counts.get(slug, 0)
        counts[slug] = suffix + 1
        result.add(slug + (f'-{suffix}' if suffix else ''))
    return result


def source_version(root):
    manifest = (root / 'crates/rsdeseq2/Cargo.toml').read_text(encoding='utf-8')
    return re.search(r'^version\s*=\s*"([^"]+)"', manifest, re.M)[1]


def check_markdown(path, root, *, package_root=None, tag=None):
    """Return diagnostics; packaged relative links must exist inside the package."""
    path, root = path.resolve(), root.resolve()
    text = prose(path.read_text(encoding='utf-8'))
    refs = {key.casefold(): value for key, value in REFERENCE.findall(text)}
    errors = []
    for label, reference in re.findall(r'\[([^]\n]+)\]\[([^]\n]*)\]', text):
        if (reference or label).casefold() not in refs:
            errors.append(f'{path.name}: undefined reference [{reference or label}]')
    targets = LINK.findall(text) + list(refs.values())
    for raw in targets:
        url = urlsplit(raw.strip('<>'))
        if url.scheme or url.netloc:
            prefix = '/deminden/rsdeseq2/blob/'
            if url.netloc != 'github.com' or not url.path.startswith(prefix):
                continue
            revision, separator, relative = url.path[len(prefix):].partition('/')
            if not separator:
                errors.append(f'{path.name}: invalid repository link {raw}')
                continue
            if package_root and revision != tag:
                errors.append(f'{path.name}: release link is not pinned to {tag}: {raw}')
                continue
            if revision not in {'main', 'dev', tag}:
                continue
            target = root / unquote(relative)
            boundary = root
        else:
            target = path.parent / unquote(url.path) if url.path else path
            boundary = package_root.resolve() if package_root else root
        target = target.resolve()
        if not target.is_relative_to(boundary):
            errors.append(f'{path.name}: link escapes documentation tree: {raw}')
        elif not target.exists():
            errors.append(f'{path.name}: missing target {raw}')
        elif url.fragment and target.suffix.lower() == '.md':
            if unquote(url.fragment) not in anchors(target.read_text(encoding='utf-8')):
                errors.append(f'{path.name}: missing anchor {raw}')
    return errors


def check_repository(root):
    tracked = subprocess.check_output(['git', 'ls-files', '-z'], cwd=root).decode().split('\0')
    paths = [root / name for name in tracked if name.lower().endswith('.md')]
    tag = 'v' + source_version(root)
    errors = []
    for path in paths:
        if not path.is_file():
            errors.append(f'{path.relative_to(root)}: tracked documentation is missing')
        else:
            errors.extend(f'{path.relative_to(root)}: {error}' for error in check_markdown(path, root, tag=tag))
    # Exercise the same generation and validation that binary packaging uses.
    from scripts.release.prepare_release_docs import write_binary_readme, write_release_notes
    with tempfile.TemporaryDirectory() as temp:
        output = Path(temp)
        for target in ('x86_64-unknown-linux-gnu', 'x86_64-pc-windows-msvc'):
            (output / 'LICENSE').write_bytes((root / 'LICENSE').read_bytes())
            (output / 'THIRD_PARTY_NOTICES.md').write_bytes((root / 'THIRD_PARTY_NOTICES.md').read_bytes())
            write_binary_readme(root, output, tag, target)
        write_release_notes(root, output / 'release-notes.md', tag)
    return paths, errors


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.parse_args()
    paths, errors = check_repository(ROOT)
    if errors:
        raise SystemExit('\n'.join(errors))
    print(f'Documentation checks passed: {len(paths)} Markdown files and generated release documentation.')


if __name__ == '__main__':
    main()
