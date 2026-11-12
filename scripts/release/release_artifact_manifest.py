#!/usr/bin/env python3
"""Bind checked release artifacts to their source commit and verify downloads."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import subprocess
from scripts.release.stage_release_crate import require_frozen_tree

ROOT = Path(__file__).resolve().parents[2]


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_file_sha256(relative):
    # Git blob bytes avoid checkout CRLF differences on Windows release runners.
    content = subprocess.check_output(['git','show','HEAD:'+relative], cwd=ROOT)
    return hashlib.sha256(content).hexdigest()


def source_identity():
    return subprocess.check_output(['git','rev-parse','HEAD'], cwd=ROOT, text=True).strip()


def check_files(root, files):
    for name, expected in files.items():
        if Path(name).name != name or not (root/name).is_file() or sha256(root/name) != expected:
            raise ValueError(f'Release artifact identity differs: {name}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=['write','verify'])
    parser.add_argument('directory', type=Path)
    args = parser.parse_args()
    require_frozen_tree()
    root = args.directory
    path = root/'validation.json'
    if args.mode == 'write':
        files = {p.name: sha256(p) for p in root.iterdir() if p.is_file() and p.name != path.name}
        if len([p for p in files if p.endswith('.crate')]) != 1 or len([p for p in files if p.endswith('.tar.gz')]) != 1:
            raise ValueError('Exactly one checked Rust crate and R source package are required')
        check_files(root, files)
        data = dict(source_commit=source_identity(), trigger_sha=os.environ.get('GITHUB_SHA', source_identity()),
                    release_tag=os.environ.get('GITHUB_REF_NAME'), files=files,
                    rustc=subprocess.check_output(['rustc','--version'], text=True).strip(),
                    r_reference_lock_sha256=source_file_sha256('config/r-reference/r-reference-sources-2026-10-29.tsv'),
                    r_reference_environment_sha256=source_file_sha256('config/r-reference/r-reference-environment-2026-10-29.json'),
                    required_checks=['rust','r-reference','frozen-parity','msrv'])
        path.write_text(json.dumps(data, indent=2, sort_keys=True)+'\n')
    else:
        data = json.loads(path.read_text())
        if data['source_commit'] != source_identity() or data['trigger_sha'] != os.environ.get('GITHUB_SHA', source_identity()):
            raise ValueError('Validated artifacts belong to a different source revision')
        if data['release_tag'] != os.environ.get('GITHUB_REF_NAME'):
            raise ValueError('Validated artifacts belong to a different release tag')
        if data['r_reference_lock_sha256'] != source_file_sha256('config/r-reference/r-reference-sources-2026-10-29.tsv'):
            raise ValueError('R reference source lock changed')
        if data['r_reference_environment_sha256'] != source_file_sha256('config/r-reference/r-reference-environment-2026-10-29.json'):
            raise ValueError('R reference environment changed')
        check_files(root, data['files'])
    print('Release artifact identities verified.')


if __name__ == '__main__':
    main()
