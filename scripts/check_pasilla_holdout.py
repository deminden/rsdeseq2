#!/usr/bin/env python3
"""Milestone check of the frozen independent panel; never regenerate references.

The absolute bounds below reuse the existing GTEx parity test bounds. Dispersion
and adjusted-p errors are reported; exact baseline identity also guards them.
"""
import argparse
import json
import os
from pathlib import Path
import subprocess

from benchmark_gtex_native_wald import sha256
from check_gtex_dispersion_improvement import compare

ROOT = Path(__file__).resolve().parents[1]
BOUNDS = dict(baseMean=1e-6, log2FoldChange=1e-5, lfcSE=1e-6, stat=1e-5, pvalue=1e-6)
COLUMNS = (*BOUNDS, 'dispersion', 'padj')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--baseline-binary', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--fixture', type=Path, default=ROOT/'crates/rsdeseq2/tests/data/pasilla_holdout')
    args = parser.parse_args()
    frozen = json.loads((args.fixture/'sha256.json').read_text())
    def verify():
        if any(sha256(args.fixture/name) != value for name, value in frozen.items()):
            raise ValueError('Independent panel differs from its frozen identity')
    verify()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    reduced = args.output_dir/'reduced.tsv'
    reduced.write_text('\n'.join('\t'.join(line.split('\t')[:3]) for line in
                      (args.fixture/'design.tsv').read_text().splitlines())+'\n')
    binaries = dict(baseline=args.baseline_binary.resolve(), candidate=args.binary.resolve())
    identities = {label:sha256(path) for label,path in binaries.items()}
    env = dict(os.environ, RAYON_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
    report = dict(status='in_progress', absolute_bounds=BOUNDS, fixture_sha256=frozen,
                  binaries={label:dict(path=str(path),sha256=identities[label]) for label,path in binaries.items()}, tests={})
    report_path = args.output_dir/'acceptance.json'
    report_path.write_text(json.dumps(report,indent=2)+'\n')
    for test in ['wald', 'lrt']:
        outputs = {}
        for label,binary in binaries.items():
            prefix = args.output_dir/f'{test}_{label}'
            outputs[label] = [prefix.with_suffix('.tsv'),prefix.with_suffix('.diagnostics.tsv'),
                              prefix.with_suffix('.beta.tsv'),prefix.with_suffix('.se.tsv')]
            result,diagnostics,beta,se = outputs[label]
            command = [str(binary), test, '--counts', str(args.fixture/'counts.tsv'),
                       '--design', str(args.fixture/'design.tsv'), '--disable-cooks-cutoff',
                       '--disable-independent-filtering', '--output', str(result),
                       '--fit-diagnostics-output', str(diagnostics), '--fit-beta-output', str(beta),
                       '--fit-beta-se-output', str(se)]
            if test == 'lrt': command += ['--reduced-design',str(reduced)]
            subprocess.run(command,check=True,env=env)
        exact = all(sha256(a) == sha256(b) for a,b in zip(outputs['baseline'],outputs['candidate']))
        metrics = compare(outputs['candidate'][0],args.fixture/f'reference_{test}.tsv',COLUMNS)
        oracle = all(row['n_finite'] > 0 and row['finite_na_mismatches'] == 0 and
                     (row['column'] not in BOUNDS or row['max_abs'] < BOUNDS[row['column']]) for row in metrics)
        report['tests'][test] = dict(exact_baseline=exact,oracle_passed=oracle,metrics=metrics)
        report_path.write_text(json.dumps(report,indent=2)+'\n')
        if not exact or not oracle: raise ValueError(f'{test} independent milestone check failed; see {report_path}')
    verify()
    if any(sha256(path) != identities[label] for label,path in binaries.items()):
        raise ValueError('Executable changed during independent validation')
    report['status'] = 'complete'
    report_path.write_text(json.dumps(report,indent=2)+'\n')
    print('Independent Wald/LRT panel passed; result, diagnostic, beta and SE exports match baseline exactly.')


if __name__ == '__main__':
    main()
