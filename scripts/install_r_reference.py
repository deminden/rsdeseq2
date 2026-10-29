#!/usr/bin/env python3
"""Restore or verify the pinned numerical R reference environment."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
from urllib.parse import urlparse

ROOT = Path(__file__).resolve().parents[1]


def load_manifest(path):
    with path.open() as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    names = [row['package'] for row in rows]
    if not rows or len(set(names)) != len(names):
        raise ValueError('Empty or duplicate R package lock')
    for row in rows:
        if not re.fullmatch(r'[A-Za-z][A-Za-z0-9.]*', row['package']) or not re.fullmatch(r'[0-9][0-9A-Za-z.-]*', row['version']):
            raise ValueError('Invalid locked package identity')
        if not re.fullmatch(r'[0-9a-f]{64}', row['sha256']) or urlparse(row['url']).scheme != 'https':
            raise ValueError('Invalid locked source hash or URL')
    return rows


def dependency_order(rows):
    remaining = {row['package']: row for row in rows}
    known = set(remaining)
    for row in rows:
        if set(filter(None, row['dependencies'].split(';'))) - known:
            raise ValueError('Dependency is absent from the R source lock')
    ordered = []
    while remaining:
        ready = sorted(name for name, row in remaining.items()
                       if not set(filter(None, row['dependencies'].split(';'))) & remaining.keys())
        if not ready:
            raise ValueError('Cyclic R package dependency lock')
        ordered.extend(remaining.pop(name) for name in ready)
    return ordered


def require_source_hash(path, expected):
    if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
        raise ValueError(f'Locked R source hash differs: {path}')


def source_urls(row):
    urls = [row['url']]
    if urlparse(row['url']).netloc == 'cran.r-project.org':
        archive = f'https://cran.r-project.org/src/contrib/Archive/{row["package"]}/{row["package"]}_{row["version"]}.tar.gz'
        if archive not in urls:
            urls.append(archive)
    return urls


def download_source(row, destination):
    temp = destination.with_suffix('.download')
    for url in source_urls(row):
        try:
            subprocess.run(['curl','--fail','--location','--silent','--show-error','--max-time','120','--output',str(temp),url], check=True)
        except subprocess.CalledProcessError:
            continue
        require_source_hash(temp, row['sha256'])
        temp.replace(destination)
        return
    raise ValueError(f'Locked R source could not be downloaded: {row["package"]}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, default=ROOT/'docs/data/r-reference-sources-2026-10-29.tsv')
    parser.add_argument('--cache', type=Path, default=ROOT/'results/r-reference-sources')
    parser.add_argument('--library', type=Path, default=ROOT/'results/r-reference-library')
    parser.add_argument('--rscript', default=shutil.which('Rscript'))
    parser.add_argument('--environment', type=Path, default=ROOT/'docs/data/r-reference-environment-2026-10-29.json')
    parser.add_argument('--check-installed', action='store_true')
    args = parser.parse_args()
    if not args.rscript:
        parser.error('Rscript is required')
    environment = json.loads(args.environment.read_text())
    rows = dependency_order(load_manifest(args.manifest))
    args.cache.mkdir(parents=True, exist_ok=True)
    for row in rows:
        source = args.cache/f'{row["package"]}_{row["version"]}.tar.gz'
        if not args.check_installed:
            if not source.exists():
                download_source(row, source)
            require_source_hash(source, row['sha256'])
        row['archive'] = str(source.resolve())
    ordered = args.cache/'installation-order.tsv'
    with ordered.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader(); writer.writerows(rows)
    subprocess.run([args.rscript, '--vanilla', str(ROOT/'scripts/install_locked_r_packages.R'),
                    str(ordered.resolve()), str(args.library.resolve()), environment['r_version'],
                    'check' if args.check_installed else 'install'], check=True)


if __name__ == '__main__':
    main()
