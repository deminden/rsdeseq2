#!/usr/bin/env python3
"""Significance-aware parity of the eleven GTEx contrasts against DESeq2.

Runs the result-only native Wald command (Cook's cutoff and independent
filtering disabled, as in the frozen references) for each contrast and
compares with `reference_results.tsv` using `parity_metrics`: log10
differences of small p-values, Spearman correlation of the statistic, top-k
overlap and significance flips. Any flip fails (exit status 1) unless
`--report-only` is given.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from scripts.benchmarks.paired_gtex_rust_timing import contrast_dirs, verify_hashes  # noqa: E402
from scripts.validation.parity_metrics import rows_by_gene, significance_metrics  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--binary', type=Path, default=ROOT / 'target/release/rsdeseq2')
    parser.add_argument('--fixture-root', type=Path,
                        default=ROOT / 'results/benchmarks/gtex_mode_contrasts')
    parser.add_argument('--reference-hashes', type=Path,
                        default=ROOT / 'docs/data/gtex_dispersion_reference_sha256.tsv')
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--threads', type=int, default=2)
    parser.add_argument('--report-only', action='store_true')
    args = parser.parse_args()
    dirs = contrast_dirs(args.fixture_root)
    verify_hashes(dirs, args.reference_hashes)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, RAYON_NUM_THREADS=str(args.threads))
    table = []
    for d in dirs:
        result = args.output_dir / f'{d.name}.tsv'
        subprocess.run([str(args.binary.resolve()), 'wald', '--counts', str(d / 'counts.tsv'),
                        '--design', str(d / 'design_matrix.tsv'), '--output', str(result),
                        '--disable-cooks-cutoff', '--disable-independent-filtering'],
                       check=True, env=env, stdout=subprocess.DEVNULL)
        with result.open(newline='') as got_handle, (d / 'reference_results.tsv').open(newline='') as ref_handle:
            metrics = significance_metrics(rows_by_gene(csv.DictReader(got_handle, delimiter='\t')),
                                           rows_by_gene(csv.DictReader(ref_handle, delimiter='\t')))
        table.append({'contrast': d.name, **metrics})
        print(d.name, json.dumps(metrics), flush=True)
    with (args.output_dir / 'significance.tsv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(table[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(table)
    flips = sum(row['flips_padj_0.05'] + row['flips_padj_0.1'] for row in table)
    print(f'total significance flips: {flips}')
    return 0 if flips == 0 or args.report_only else 1


if __name__ == '__main__':
    raise SystemExit(main())
