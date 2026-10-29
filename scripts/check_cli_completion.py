#!/usr/bin/env python3
"""Verify every file in an rsdeseq2 CLI completion manifest."""
import argparse
import hashlib
import json
from pathlib import Path


def verify(path):
    path = Path(path)
    original = path.read_bytes()
    manifest = json.loads(original)
    if manifest.get('schema_version') != 1 or not manifest.get('generation'):
        raise ValueError('Unknown or incomplete generation metadata')
    files = manifest.get('files')
    if not isinstance(files, list) or not files:
        raise ValueError('Completion manifest contains no files')
    seen = set()
    for item in files:
        output = Path(item['path'])
        identity = output.resolve(strict=True)
        if not output.is_absolute() or identity in seen:
            raise ValueError('Completion paths must be absolute and unique')
        seen.add(identity)
        digest = hashlib.sha256()
        length = 0
        with output.open('rb') as handle:
            for block in iter(lambda: handle.read(1024*1024), b''):
                digest.update(block)
                length += len(block)
        if length != item['bytes'] or digest.hexdigest() != item['sha256']:
            raise ValueError(f'Incomplete or changed output: {output}')
    if path.read_bytes() != original:
        raise ValueError('Generation changed during verification')
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('manifest', type=Path)
    args = parser.parse_args()
    manifest = verify(args.manifest)
    print(f"Verified generation {manifest['generation']}: {len(manifest['files'])} files")


if __name__ == '__main__':
    main()
