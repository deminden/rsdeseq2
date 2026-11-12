#!/usr/bin/env python3
"""Run the same target/profile numerical panel before and during a release."""
import argparse
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[1]
RELEASE_TESTS = (
    'model_validity', 'wald_reference', 'lrt_reference', 'dispersion_reference',
    'gtex_real_wald_parity', 'thread_invariance', 'cli', 'frozen_local_vst',
    # Small independent panels exercise target-sensitive posterior arithmetic,
    # sparse optimizer state, uncertainty propagation and formula column order.
    'bounded_posterior', 'sqp_controls', 'heavy_tail', 'posterior_analysis',
    'formula_fitting',
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--target', required=True)
    args = parser.parse_args()
    command = ['cargo', 'test', '--release', '--locked', '-p', 'rsdeseq2',
               '--target', args.target]
    for name in RELEASE_TESTS:
        command.extend(['--test', name])
    subprocess.run(command, cwd=ROOT, check=True)


if __name__ == '__main__':
    main()
