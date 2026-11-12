#!/usr/bin/env python3
"""Prepare the sole permitted tracked staging edit before crate packaging."""
import json
import os
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import subprocess

ROOT = Path(__file__).resolve().parents[2]


def require_frozen_tree():
    changed = subprocess.check_output(['git','diff','HEAD','--name-only'], cwd=ROOT, text=True).splitlines()
    untracked = subprocess.check_output(['git','ls-files','--others','--exclude-standard','--',
        'crates','r-pkg','scripts','config','.cargo'], cwd=ROOT, text=True).splitlines()
    generated = {'r-pkg/rsdeseq2/src/Makevars', 'r-pkg/rsdeseq2/src/Makevars.win'}
    if set(untracked) - generated:
        raise ValueError('Release staging contains untracked source: '+', '.join(untracked))
    if set(changed) - {'crates/rsdeseq2/README.md'}:
        raise ValueError('Release staging contains unvalidated tracked changes: '+', '.join(changed))


def main():
    packages = json.loads(subprocess.check_output(['cargo','metadata','--locked','--no-deps','--format-version','1'], cwd=ROOT, text=True))['packages']
    version = next(p['version'] for p in packages if p['name'] == 'rsdeseq2')
    tag = os.environ.get('RELEASE_TAG', os.environ.get('GITHUB_REF_NAME', ''))
    if tag != 'v'+version:
        raise ValueError(f'Release tag {tag!r} does not match crate version {version}')
    require_frozen_tree()
    path = ROOT/'crates/rsdeseq2/README.md'
    path.write_text(path.read_text().replace('/blob/main/', f'/blob/{tag}/'))
    require_frozen_tree()


if __name__ == '__main__':
    main()
