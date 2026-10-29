#!/usr/bin/env python3
"""Select a fixed independent Pasilla validation panel without inspecting fit results."""
import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import tarfile

SOURCE_SHA256 = '75122d45eb2c415d3e94d1b9bc9e4f6978c4504ef7e70f6bfeccb20547cb3db5'
SOURCE_URL = 'https://bioconductor.org/packages/3.23/data/experiment/src/contrib/pasilla_1.40.0.tar.gz'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-archive', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    if hashlib.sha256(args.source_archive.read_bytes()).hexdigest() != SOURCE_SHA256:
        raise ValueError('Pasilla source archive differs from the frozen source')
    with tarfile.open(args.source_archive) as archive:
        counts = list(csv.DictReader(io.StringIO(archive.extractfile('pasilla/inst/extdata/pasilla_gene_counts.tsv').read().decode()), delimiter='\t'))
        annotation = list(csv.DictReader(io.StringIO(archive.extractfile('pasilla/inst/extdata/pasilla_sample_annotation.csv').read().decode())))
    samples = list(counts[0])[1:]
    by_sample = {row['file'].removesuffix('fb'): row for row in annotation}
    bins = [[] for _ in range(5)]
    for row in counts:
        mean = sum(int(row[s]) for s in samples)/len(samples)
        index = 0 if mean == 0 else 1 if mean < 10 else 2 if mean < 100 else 3 if mean < 1000 else 4
        bins[index].append(row['gene_id'])
    selected = set()
    for group, quota in zip(bins, [64,496,496,496,496]):
        if len(group) < quota: raise ValueError('Source expression stratum is too small')
        selected.update(sorted(group, key=lambda gene: hashlib.sha256(('20261029:'+gene).encode()).digest())[:quota])
    args.output_dir.mkdir(parents=True, exist_ok=True)
    with (args.output_dir/'counts.tsv').open('w',newline='') as handle:
        writer=csv.writer(handle,delimiter='\t',lineterminator='\n')
        writer.writerow(['gene',*samples])
        writer.writerows([row['gene_id'],*[row[s] for s in samples]] for row in counts if row['gene_id'] in selected)
    with (args.output_dir/'design.tsv').open('w',newline='') as handle:
        writer=csv.writer(handle,delimiter='\t',lineterminator='\n')
        writer.writerow(['sample','Intercept','typesingle_read','conditiontreated'])
        writer.writerows([s,1,int(by_sample[s]['type']=='single-read'),int(by_sample[s]['condition']=='treated')] for s in samples)
    (args.output_dir/'provenance.json').write_text(json.dumps(dict(source=SOURCE_URL,source_sha256=SOURCE_SHA256,
        source_package='pasilla 1.40.0',source_license='LGPL',genes=len(selected),samples=samples,
        selection='Expression strata 0/(0,10)/[10,100)/[100,1000)/[1000,infinity); SHA256(20261029:gene) within each stratum; quotas 64/496/496/496/496',
        purpose='Independent validation; not used for optimization tuning'),indent=2)+'\n')


if __name__ == '__main__':
    main()
