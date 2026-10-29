#!/usr/bin/env python3
"""Generate or check the R package's vendored core without changing check input."""
import argparse
import json
from pathlib import Path
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
VENDOR = ROOT / 'r-pkg/rsdeseq2/src/rust/vendor/rsdeseq2'


def metadata(manifest):
    return json.loads(subprocess.check_output([
        'cargo', 'metadata', '--manifest-path', str(manifest), '--locked',
        '--no-deps', '--format-version', '1'], text=True))['packages']


def vendor_manifest(package):
    lines = ['[package]']
    for key, value in [('name', package['name']), ('version', package['version']),
                       ('edition', package['edition']), ('license', package['license']),
                       ('rust-version', package['rust_version'])]:
        lines.append(f'{key} = {json.dumps(value)}')
    lines += ['', '[lib]', 'name = "rsdeseq2"', 'path = "src/lib.rs"']
    sections = {}
    for dep in package['dependencies']:
        if dep['kind'] == 'dev':
            continue
        if dep.get('path'):
            raise ValueError('Vendoring a path dependency requires an explicit packaging rule')
        kind = 'build-dependencies' if dep['kind'] == 'build' else 'dependencies'
        section = kind if dep['target'] is None else f'target.{json.dumps(dep["target"])}.{kind}'
        fields = [f'version = {json.dumps(dep["req"])}']
        if dep['rename']:
            fields.append(f'package = {json.dumps(dep["name"])}')
        if dep['features']:
            fields.append(f'features = {json.dumps(dep["features"])}')
        if dep['optional']:
            fields.append('optional = true')
        if not dep['uses_default_features']:
            fields.append('default-features = false')
        if dep.get('registry'):
            raise ValueError('Custom dependency registries require an explicit packaging rule')
        sections.setdefault(section, []).append(f'{dep["rename"] or dep["name"]} = {{ '+', '.join(fields)+' }')
    for name in sorted(sections):
        lines += ['', f'[{name}]', *sorted(sections[name])]
    if package.get('features'):
        lines += ['', '[features]']
        lines += [f'{name} = {json.dumps(values)}' for name, values in sorted(package['features'].items())]
    return '\n'.join(lines)+'\n'


def tree_files(root):
    return {str(p.relative_to(root)): p.read_bytes() for p in root.rglob('*') if p.is_file()}


def check_tree(expected, actual):
    wanted, got = tree_files(expected), tree_files(actual)
    changed = sorted(key for key in wanted.keys() | got.keys() if wanted.get(key) != got.get(key))
    if changed:
        raise ValueError('R vendor differs; run scripts/sync_r_vendor.sh:\n'+'\n'.join(changed))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true', help='read-only source, manifest and version check')
    args = parser.parse_args()
    package = next(p for p in metadata(ROOT/'Cargo.toml') if p['name'] == 'rsdeseq2')
    version = package['version']
    description = (ROOT/'r-pkg/rsdeseq2/DESCRIPTION').read_text().splitlines()
    r_version = next(line.split(':', 1)[1].strip() for line in description if line.startswith('Version:'))
    backend = next(p for p in metadata(ROOT/'r-pkg/rsdeseq2/src/rust/Cargo.toml') if p['name'] == 'rsdeseq2-r')
    if r_version != version or backend['version'] != version:
        raise ValueError('Rust crate, R package and R backend versions must agree')
    with tempfile.TemporaryDirectory(prefix='rsdeseq2-vendor-') as temp:
        expected = Path(temp)/'rsdeseq2'
        shutil.copytree(ROOT/'crates/rsdeseq2/src', expected/'src')
        (expected/'src/main.rs').unlink(missing_ok=True)
        (expected/'Cargo.toml').write_text(vendor_manifest(package))
        (expected/'README.md').write_text(
            '# Vendored Rust Core\n\n'
            'This package contains the Rust core used to build the R backend. Regenerate it\n'
            'from `crates/rsdeseq2` with `scripts/sync_r_vendor.sh` before packaging a\n'
            'release. Do not edit this copy directly.\n')
        if args.check:
            check_tree(expected, VENDOR)
            print('R vendor, manifest and package versions match.')
        else:
            shutil.copytree(expected/'src', VENDOR/'src', dirs_exist_ok=True)
            # Remove stale source paths as well as copying new ones.
            for path in list((VENDOR/'src').rglob('*')):
                if path.is_file() and not (expected/'src'/path.relative_to(VENDOR/'src')).exists():
                    path.unlink()
            for name in ('Cargo.toml', 'README.md'):
                shutil.copyfile(expected/name, VENDOR/name)
            print(f'Synchronized Rust core into {VENDOR}')


if __name__ == '__main__':
    main()
