#!/usr/bin/env python3
"""Regenerate exact Linux formula oracles with original base R, independently of rsdeseq2."""
import argparse
import json
import os
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[1]
PANELS = {
    'formula_expression': 'formula_expression_r462',
    'extended_formula': 'formula_expression_extended_r462',
    'model_matrix': 'model_matrix_r462',
    'native_formula': 'native_formula_r462',
    'numeric_transform': 'numeric_transform_r462',
    'native_factor': 'native_factor_r462',
    'formula_fitting': 'formula_fitting_r462',
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rscript', default='Rscript')
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=False)
    env = dict(os.environ, LC_ALL='C.UTF-8', OPENBLAS_CORETYPE='Haswell', OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
    for generator, stem in PANELS.items():
        output = args.output_dir / (stem + '.json')
        subprocess.run([args.rscript, '--vanilla', str(ROOT / 'scripts' / ('generate_' + generator + '_reference.R')), str(output.resolve())], env=env, check=True)
        actual = json.loads(output.read_text())
        expected = json.loads((ROOT / 'crates/rsdeseq2/tests/data' / (stem + '_linux_x86_64.json')).read_text())
        metadata_field = 'source' if generator == 'formula_fitting' else 'metadata'
        if actual[metadata_field]['platform'] != 'x86_64-pc-linux-gnu':
            raise ValueError('Linux x86_64 R is required to regenerate these platform oracles')
        # Session/provenance fields differ from the original independently generated
        # snapshots. Every input, expected bit, attribute, and rejection stays exact.
        actual.pop(metadata_field)
        expected.pop(metadata_field)
        if actual != expected:
            raise ValueError('Independent R oracle differs: ' + stem)
        print(stem + ': independently regenerated inputs and expected values match exactly')


if __name__ == '__main__':
    main()
