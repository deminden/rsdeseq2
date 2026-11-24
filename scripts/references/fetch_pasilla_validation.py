#!/usr/bin/env python3
"""Fetch the independently licensed Pasilla fixture for local validation."""

import hashlib
import json
from pathlib import Path
import tempfile
from urllib.request import urlopen


ROOT = Path(__file__).resolve().parents[2]
FIXTURE = Path('crates/rsdeseq2/tests/data/pasilla_holdout')
ARCHIVE = 'https://raw.githubusercontent.com/deminden/rsdeseq2/f1c9e49f08b293b62dd0c8133ca9ad7e6655bbb6'
DATA_FILES = ('counts.tsv', 'design.tsv', 'reference_lrt.tsv', 'reference_wald.tsv')


def main():
    output = ROOT / FIXTURE
    expected = json.loads((output / 'sha256.json').read_text())
    for name in DATA_FILES:
        target = output / name
        if target.is_file():
            if hashlib.sha256(target.read_bytes()).hexdigest() != expected[name]:
                raise ValueError(f'{target} differs from the frozen fixture')
            continue
        with urlopen(f'{ARCHIVE}/{FIXTURE.as_posix()}/{name}', timeout=30) as response:
            data = response.read()
        if hashlib.sha256(data).hexdigest() != expected[name]:
            raise ValueError(f'Downloaded {name} differs from the frozen fixture')
        with tempfile.TemporaryDirectory(dir=output) as temporary:
            downloaded = Path(temporary) / name
            downloaded.write_bytes(data)
            downloaded.replace(target)
    print(f'Pasilla validation files verified in {output}')


if __name__ == '__main__':
    main()
