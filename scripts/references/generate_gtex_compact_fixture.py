#!/usr/bin/env python3
"""Select a deterministic GTEx Wald panel; retain full-cohort normalization.

No Rust output is used for selection or reference generation. The original
56,200-gene fixture remains available for full replay and contrast generation.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import subprocess

ROOT = Path(__file__).resolve().parents[2]

def read(path):
    with path.open(newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))

def select_genes(counts, reference, size=512):
    if len(counts) != len(reference) or [r['gene'] for r in counts] != [r['gene'] for r in reference]:
        raise ValueError('count/reference gene alignment mismatch')
    if len({r['gene'] for r in counts}) != len(counts):
        raise ValueError('duplicate genes')
    reasons = {}
    def add(i, reason):
        reasons.setdefault(i, set()).add(reason)
    # Cover each missingness pattern, including all-zero rows.
    masks = {}
    for i, r in enumerate(reference):
        mask = tuple(r[c] in ('NA', 'NaN', '') for c in ('dispersion', 'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue'))
        masks.setdefault(mask, []).append(i)
    for mask, indices in sorted(masks.items()):
        for i in indices[:8]: add(i, 'missingness_' + ''.join(str(int(v)) for v in mask))
    metrics = {c: [(i, abs(float(r[c]))) for i, r in enumerate(reference)
                   if r[c] not in ('NA', 'NaN', '') and math.isfinite(float(r[c]))]
               for c in ('baseMean', 'dispersion', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue')}
    metrics['nonzero_samples'] = [(i, sum(int(v)>0 for k,v in r.items() if k!='gene')) for i,r in enumerate(counts)]
    metrics['maximum_count'] = [(i, max(int(v) for k,v in r.items() if k!='gene')) for i,r in enumerate(counts)]
    for name, values in metrics.items():
        ranked = sorted(values, key=lambda v: (v[1], v[0]))
        for i,_ in ranked[:8] + ranked[-8:]: add(i, name+'_tail')
        for q in range(33): add(ranked[(len(ranked)-1)*q//32][0], name+'_quantile')
    if len(reasons) > size:
        raise ValueError(f'panel size {size} cannot retain {len(reasons)} mandatory coverage genes')
    # Uniform positions spread extra genes across the original matrix; deterministic
    # and independent of optimizer/platform differences.
    for i in range(len(counts)):
        candidate = i * 7919 % len(counts)
        if len(reasons) == min(size, len(counts)): break
        add(candidate, 'background')
    return [(i, sorted(why)) for i,why in sorted(reasons.items())]

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fixture-dir', type=Path, default=ROOT/'crates/rsdeseq2/tests/data/gtex_v8_real_wald')
    parser.add_argument('--rscript', required=True)
    args = parser.parse_args()
    counts, reference = (read(args.fixture_dir/name) for name in ('counts.tsv','reference_results.tsv'))
    selected = select_genes(counts, reference)
    out = args.fixture_dir/'compact'; out.mkdir(exist_ok=True)
    for name, rows in (('counts.tsv', counts), ('reference_results.tsv', reference)):
        with (out/name).open('w', newline='') as h:
            w=csv.DictWriter(h,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n')
            w.writeheader(); w.writerows(rows[i] for i,_ in selected)
    subprocess.run([args.rscript,'--vanilla',str(ROOT/'scripts/references/generate_gtex_compact_reference.R'),str(out)],check=True)
    manifest = dict(selection_version=1, source_genes=len(counts), selected_genes=len(selected), samples=len(counts[0])-1,
        normalization='full-cohort reference_size_factors.tsv; no subset re-estimation',
        padj='R stats::p.adjust(method="BH") on selected reference p-values',
        source_sha256={name:hashlib.sha256((args.fixture_dir/name).read_bytes()).hexdigest()
                       for name in ('counts.tsv','reference_results.tsv','reference_size_factors.tsv','design_matrix.tsv')},
        genes=[dict(gene=counts[i]['gene'],source_row=i,reasons=why) for i,why in selected])
    (out/'selection.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(f'{len(counts)} -> {len(selected)} genes, all {len(counts[0])-1} samples')

if __name__=='__main__': main()
