#!/usr/bin/env python3
"""Generate the R package's copy of the Rust core, or check its inputs.

The copy (r-pkg/rsdeseq2/src/rust/vendor) is generated, not committed: R's
configure runs this script in a git checkout, and
scripts/release/build_r_source.sh runs it before `R CMD build`. `--check`
verifies, without writing, that the copy can be generated and that versions
and notices agree.
"""
import argparse
import json
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import shutil
import subprocess
import re
import tempfile

ROOT = Path(__file__).resolve().parents[2]
VENDOR_ROOT = ROOT / 'r-pkg/rsdeseq2/src/rust/vendor'
VENDOR = VENDOR_ROOT / 'rsdeseq2'


def metadata(manifest):
    return json.loads(subprocess.check_output([
        'cargo', 'metadata', '--manifest-path', str(manifest), '--locked',
        '--no-deps', '--format-version', '1'], text=True))['packages']


# Features of workspace crates the R backend never enables. They are emptied,
# and the optional dependencies only they use are left out of the vendored
# manifests,
# so an offline build against the vendored crates (`cargo vendor` vendors only
# what Cargo.lock resolves) does not look for them.
UNUSED_FEATURES = {'rsdeseq2': {'cli'}}


def vendor_manifest(package):
    dropped = UNUSED_FEATURES.get(package['name'], set())
    # A dropped feature stays declared, empty, so `cfg(feature = ...)` checks pass.
    features = {name: ([] if name in dropped else [v for v in values if v not in dropped])
                for name, values in package.get('features', {}).items()}
    used = {v[4:] for values in features.values() for v in values if v.startswith('dep:')}
    dropped_deps = {v[4:] for name in dropped for v in package.get('features', {}).get(name, [])
                    if v.startswith('dep:')} - used
    package = dict(package, features=features,
                   dependencies=[d for d in package['dependencies']
                                 if (d['rename'] or d['name']) not in dropped_deps])
    lines = ['[package]']
    for key, value in [('name', package['name']), ('version', package['version']),
                       ('edition', package['edition']), ('license', package['license']),
                       ('rust-version', package['rust_version'])]:
        lines.append(f'{key} = {json.dumps(value)}')
    lib_name = package['name'].replace('-', '_')
    lines += ['exclude = [".DS_Store", "**/.DS_Store"]', '', '[lib]', f'name = "{lib_name}"', 'path = "src/lib.rs"']
    sections = {}
    for dep in package['dependencies']:
        if dep['kind'] == 'dev':
            continue
        kind = 'build-dependencies' if dep['kind'] == 'build' else 'dependencies'
        section = kind if dep['target'] is None else f'target.{json.dumps(dep["target"])}.{kind}'
        fields = [f'version = {json.dumps(dep["req"])}']
        if dep.get('path'):
            # Workspace crates are vendored side by side.
            fields.insert(0, f'path = {json.dumps("../" + dep["name"])}')
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


def vendored_crates(workspace):
    """`rsdeseq2` and the workspace crates it depends on through paths."""
    names, pending = [], ['rsdeseq2']
    while pending:
        name = pending.pop()
        if name in names:
            continue
        names.append(name)
        pending += [dep['name'] for dep in workspace[name]['dependencies']
                    if dep.get('path') and dep['kind'] != 'dev']
    return sorted(names)


def tree_files(root):
    return {str(p.relative_to(root)): p.read_bytes() for p in root.rglob('*') if p.is_file() and p.name != '.DS_Store'}


def check_tree(expected, actual):
    wanted, got = tree_files(expected), tree_files(actual)
    changed = sorted(key for key in wanted.keys() | got.keys() if wanted.get(key) != got.get(key))
    if changed:
        raise ValueError('R vendor differs; run scripts/release/sync_r_vendor.sh:\n'+'\n'.join(changed))


def sync_notice(source, destination, check=False):
    """Ship identical license/notice files in packages; check is read-only."""
    wanted = source.read_bytes()
    if check:
        if not destination.is_file() or destination.read_bytes() != wanted:
            raise ValueError(f'Package third-party notices differ: {destination.name}; run scripts/release/sync_r_vendor.sh')
    else:
        destination.write_bytes(wanted)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true', help='read-only source, manifest and version check')
    args = parser.parse_args()
    workspace = {p['name']: p for p in metadata(ROOT/'Cargo.toml')}
    package = workspace['rsdeseq2']
    version = package['version']
    description = (ROOT/'r-pkg/rsdeseq2/DESCRIPTION').read_text().splitlines()
    r_version = next(line.split(':', 1)[1].strip() for line in description if line.startswith('Version:'))
    # Read, not resolved: the backend's path dependencies are the generated copy.
    backend_manifest = (ROOT/'r-pkg/rsdeseq2/src/rust/Cargo.toml').read_text()
    backend_version = re.search(r'^version\s*=\s*"([^"]+)"', backend_manifest, re.M).group(1)
    if r_version != version or backend_version != version:
        raise ValueError('Rust crate, R package and R backend versions must agree')
    crates = vendored_crates(workspace)
    if any(workspace[name]['version'] != version for name in crates):
        raise ValueError('Workspace crate versions must agree')
    with tempfile.TemporaryDirectory(prefix='rsdeseq2-vendor-') as temp:
        for name in crates:
            source = Path(workspace[name]['manifest_path']).parent
            expected = Path(temp)/name
            shutil.copytree(source/'src', expected/'src', ignore=shutil.ignore_patterns('.DS_Store'))
            (expected/'src/main.rs').unlink(missing_ok=True)
            (expected/'Cargo.toml').write_text(vendor_manifest(workspace[name]))
            (expected/'README.md').write_text(
                '# Vendored Rust Core\n\n'
                f'This package contains the Rust core used to build the R backend. Regenerate it\n'
                f'from `crates/{name}` with `scripts/release/sync_r_vendor.sh` before packaging a\n'
                'release. Do not edit this copy directly.\n')
        expected_root = Path(temp)
        if not args.check:
            VENDOR_ROOT.mkdir(parents=True, exist_ok=True)
            for path in VENDOR_ROOT.iterdir():
                if path.is_dir() and path.name not in crates:
                    shutil.rmtree(path)
            for name in crates:
                target = VENDOR_ROOT/name
                shutil.copytree(expected_root/name/'src', target/'src', dirs_exist_ok=True)
                # Remove stale source paths as well as copying new ones.
                for path in list((target/'src').rglob('*')):
                    if path.is_file() and path.name != '.DS_Store' and not (expected_root/name/'src'/path.relative_to(target/'src')).exists():
                        path.unlink()
                for file in ('Cargo.toml', 'README.md'):
                    shutil.copyfile(expected_root/name/file, target/file)
            print(f'Synchronized Rust core ({", ".join(crates)}) into {VENDOR_ROOT}')
        for name in ('LICENSE', 'THIRD_PARTY_NOTICES.md'):
            # Every published workspace crate ships the notices.
            for crate in sorted(workspace):
                sync_notice(ROOT/name, Path(workspace[crate]['manifest_path']).parent/name, args.check)
        sync_notice(ROOT/'THIRD_PARTY_NOTICES.md', ROOT/'r-pkg/rsdeseq2/inst/COPYRIGHTS', args.check)
        if args.check:
            print('R vendor manifests, notices and package versions match.')


if __name__ == '__main__':
    main()
