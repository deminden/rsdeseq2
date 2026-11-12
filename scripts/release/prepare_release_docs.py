#!/usr/bin/env python3
"""Generate and validate tag-pinned documentation for binary downloads and releases."""
import argparse
from pathlib import Path
import sys
from urllib.parse import quote, unquote, urlsplit

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts.validation.check_docs import LINK, REFERENCE, REPOSITORY, check_markdown, source_version

ROOT = Path(__file__).resolve().parents[2]


def require_tag(root, tag):
    if tag != 'v' + source_version(root):
        raise ValueError('Documentation tag must match the source version')


def binary_readme(tag, target):
    windows = target.endswith('-msvc')
    command = r'.\rsdeseq2.exe' if windows else './rsdeseq2'
    language = 'powershell' if windows else 'sh'
    base = f'{REPOSITORY}/blob/{tag}'
    return f'''# rsdeseq2 {tag.removeprefix('v')}

This archive contains the compiled CLI for `{target}`. R and a Rust compiler
are not required to run it.

## Run

Extract the archive, open a terminal in the extracted directory, then run:

```{language}
{command} --version
{command} --help
{command} wald --help
```

Use your own count and numeric design tables:

```{language}
{command} wald --counts counts.tsv --design design.tsv --fit-type mean --coefficient 1 --output results.tsv
```

`--fit-type mean` explicitly selects a constant mean dispersion trend. Choose
the trend appropriate to your analysis. Coefficient indices are zero-based.
See the [CLI guide]({base}/docs/cli.md) for input schemas and options.
You can place the executable in a directory on your PATH for use elsewhere.

## Documentation

- [Supported workflows and limits]({base}/docs/compatibility.md)
- [Precision, speed and memory]({base}/README.md#validation-and-performance)
- [R interface]({base}/docs/r-wrapper.md)
- [Release notes]({base}/docs/releases/{tag.removeprefix('v')}.md)

For bugs, report the version, target and command through
[GitHub Issues]({REPOSITORY}/issues).

## License

Project code is MIT-licensed. Adapted numerical components additionally carry
BSD-3-Clause and BSD-2-Clause-Views terms. Keep [LICENSE](LICENSE) and
[third-party notices](THIRD_PARTY_NOTICES.md) with redistributed binaries.
'''


def validate(path, root, directory, tag):
    errors = check_markdown(path, root, package_root=directory, tag=tag)
    if errors:
        raise ValueError('\n'.join(errors))


def write_binary_readme(root, directory, tag, target):
    require_tag(root, tag)
    path = directory / 'README.md'
    path.write_text(binary_readme(tag, target), encoding='utf-8')
    validate(path, root, directory, tag)


def write_release_notes(root, output, tag):
    require_tag(root, tag)
    source = root / 'docs/releases' / f'{tag.removeprefix("v")}.md'
    text = source.read_text(encoding='utf-8')
    errors = check_markdown(source, root, tag=tag)
    if errors:
        raise ValueError('\n'.join(errors))

    def pinned(raw):
        url = urlsplit(raw.strip('<>'))
        if url.scheme or url.netloc:
            return raw
        target = (source.parent / unquote(url.path)).resolve() if url.path else source.resolve()
        relative = target.relative_to(root.resolve()).as_posix()
        return (f'{REPOSITORY}/blob/{tag}/{quote(relative)}'
                + (f'?{url.query}' if url.query else '')
                + (f'#{url.fragment}' if url.fragment else ''))

    # Release notes use ordinary inline/reference links, not links inside code.
    text = LINK.sub(lambda m: m[0].replace(m[1], pinned(m[1]), 1), text)
    text = REFERENCE.sub(lambda m: m[0].replace(m[2], pinned(m[2]), 1), text)
    output.write_text(text, encoding='utf-8')
    validate(output, root, output.parent, tag)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tag', required=True)
    sub = parser.add_subparsers(dest='mode', required=True)
    binary = sub.add_parser('binary')
    binary.add_argument('--directory', required=True, type=Path)
    binary.add_argument('--target', required=True)
    notes = sub.add_parser('notes')
    notes.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    if args.mode == 'binary':
        write_binary_readme(ROOT, args.directory, args.tag, args.target)
    else:
        write_release_notes(ROOT, args.output, args.tag)


if __name__ == '__main__':
    main()
