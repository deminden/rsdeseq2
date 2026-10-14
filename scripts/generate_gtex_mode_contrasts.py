#!/usr/bin/env python3
"""Prepare biological and random GTEx contrasts and optional offline DESeq2 references."""
from __future__ import annotations
import argparse
import csv
import os
import random
import shutil
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SEED = 20261014

def write_design(folder, names, metadata, values):
    folder.mkdir(parents=True, exist_ok=True)
    with (folder/'design_matrix.tsv').open('w',newline='') as handle:
        writer=csv.writer(handle,delimiter='\t',lineterminator='\n')
        writer.writerow(['sample',*names])
        for row in metadata:
            writer.writerow([row['sample'],*values(row)])

def full_counts(source, destination):
    if not destination.exists():
        try:
            os.link(source,destination)
        except OSError:
            shutil.copyfile(source,destination)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fixture-dir',type=Path,default=ROOT/'crates/rsdeseq2/tests/data/gtex_v8_real_wald')
    parser.add_argument('--output-dir',type=Path,default=ROOT/'results/benchmarks/gtex_mode_contrasts')
    parser.add_argument('--rscript',default=shutil.which('Rscript'))
    parser.add_argument('--counts-only',action='store_true')
    args=parser.parse_args()
    with (args.fixture_dir/'sample_metadata.tsv').open() as handle:
        metadata=list(csv.DictReader(handle,delimiter='\t'))
    count_source=args.fixture_dir/'counts.tsv'
    contrasts=[]
    block_names=['Intercept','perm_blockmale_3','perm_blockmale_4','conditionB']
    block=lambda row:[1,int(row['perm_block']=='male_3'),int(row['perm_block']=='male_4')]
    d=args.output_dir/'condition_blocked'
    write_design(d,block_names,metadata,lambda row:[*block(row),int(row['condition']=='B')])
    full_counts(count_source,d/'counts.tsv');contrasts.append(d)
    d=args.output_dir/'sex_female_vs_male'
    write_design(d,['Intercept','sexmale'],metadata,lambda row:[1,int(row['sex']=='male')])
    full_counts(count_source,d/'counts.tsv');contrasts.append(d)
    d=args.output_dir/'male_age4_vs_age3'
    males=[r for r in metadata if r['sex']=='male']
    write_design(d,['Intercept','male_age4'],males,lambda row:[1,int(row['age_bin']=='4')])
    with count_source.open() as source,(d/'counts.tsv').open('w',newline='') as output:
        reader=csv.reader(source,delimiter='\t');header=next(reader)
        samples=[r['sample'] for r in males];columns=[header.index(s) for s in samples]
        writer=csv.writer(output,delimiter='\t',lineterminator='\n')
        writer.writerow(['gene',*samples])
        for row in reader:writer.writerow([row[0],*[row[i] for i in columns]])
    contrasts.append(d)
    rng=random.Random(SEED)
    strata={name:[r for r in metadata if r['perm_block']==name] for name in sorted({r['perm_block'] for r in metadata})}
    for rep in range(1,6):
        selected=set()
        for rows in strata.values():
            selected.update(rng.sample([r['sample'] for r in rows],sum(r['condition']=='B' for r in rows)))
        d=args.output_dir/'random_splits'/f'blocked_random_rep{rep:02d}'
        write_design(d,[*block_names[:-1],'randomB'],metadata,lambda row:[*block(row),int(row['sample'] in selected)])
        full_counts(count_source,d/'counts.tsv');contrasts.append(d)
    for rep in range(1,4):
        selected=set(rng.sample([r['sample'] for r in metadata],sum(r['condition']=='B' for r in metadata)))
        d=args.output_dir/'random_splits'/f'unblocked_random_rep{rep:02d}'
        write_design(d,['Intercept','randomB'],metadata,lambda row:[1,int(row['sample'] in selected)])
        full_counts(count_source,d/'counts.tsv');contrasts.append(d)
    for d in contrasts:
        print(d,flush=True)
        if not args.counts_only:
            if not args.rscript:raise RuntimeError('pass --rscript for the pinned offline DESeq2 environment')
            subprocess.run([args.rscript,'--vanilla',str(ROOT/'scripts/generate_gtex_mode_contrast_references.R'),str(d)],check=True)
    return 0

if __name__=='__main__':
    raise SystemExit(main())
