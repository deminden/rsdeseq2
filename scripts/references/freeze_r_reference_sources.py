#!/usr/bin/env python3
"""Freeze exact installed R package source URLs and hashes from a package TSV."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
import hashlib
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import tarfile
import subprocess


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--packages', type=Path, required=True, help='package/version/repository/dependencies TSV')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--cache', type=Path, default=Path('results/r-reference-sources'))
    args = parser.parse_args()
    args.cache.mkdir(parents=True, exist_ok=True)
    with args.packages.open() as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))

    def freeze(row):
        name, version = row['package'], row['version']
        filename = f'{name}_{version}.tar.gz'
        path = args.cache/filename
        repository = row.get('repository') or ('Bioconductor '+row['url'].split('/packages/')[1].split('/')[0] if 'bioconductor.org/packages/' in row.get('url','') else 'CRAN')
        if repository.startswith('Bioconductor'):
            bioc = repository.split()[-1]
            urls = [f'https://bioconductor.org/packages/{bioc}/{kind}/src/contrib/{filename}'
                    for kind in ('bioc', 'data/annotation', 'data/experiment', 'workflows')]
        else:
            urls = [f'https://cran.r-project.org/src/contrib/{filename}',
                    f'https://cran.r-project.org/src/contrib/Archive/{name}/{filename}']
        errors = []
        for url in urls:
            try:
                subprocess.run(["curl", "--fail", "--location", "--silent", "--show-error", "--max-time", "60", "--output", str(path), url], check=True, capture_output=True)
                with tarfile.open(path) as archive:
                    description = archive.extractfile(f'{name}/DESCRIPTION').read().decode()
                if f'Version: {version}' not in description.splitlines():
                    raise ValueError('archive version differs')
                print(f'{name} {version}: frozen', flush=True)
                return dict(package=name, version=version, url=url, sha256=digest(path), dependencies=row['dependencies'])
            except (OSError, ValueError, tarfile.TarError, subprocess.CalledProcessError) as error:
                errors.append(f'{url}: {error}')
        raise ValueError('\n'.join(errors))

    with ThreadPoolExecutor(max_workers=4) as pool:
        frozen = list(pool.map(freeze, rows))
    with args.output.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=['package','version','url','sha256','dependencies'], delimiter='\t')
        writer.writeheader(); writer.writerows(frozen)


if __name__ == '__main__':
    main()
