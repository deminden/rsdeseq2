#!/usr/bin/env python3
"""Run the same target/profile numerical panel before and during a release."""
import argparse
import os
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import subprocess

ROOT = Path(__file__).resolve().parents[2]
RELEASE_TESTS = (
    'model_validity', 'wald_reference', 'lrt_reference', 'dispersion_reference',
    'gtex_real_wald_parity', 'thread_invariance', 'cli', 'frozen_local_vst',
    # Small independent panels exercise target-sensitive posterior arithmetic,
    # sparse optimizer state, uncertainty propagation and formula column order.
    'bounded_posterior', 'sqp_controls', 'heavy_tail', 'posterior_analysis',
    'formula_fitting',
)
LINUX_ARM64_FORMULA_TESTS = (
    'formula_expression', 'model_matrix', 'native_factor', 'native_formula',
    'numeric_transform',
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--target', required=True)
    args = parser.parse_args()
    build = ['--release', '--locked', '-p', 'rsdeseq2', '--target', args.target]
    linux_arm64 = args.target == 'aarch64-unknown-linux-gnu'
    tests = RELEASE_TESTS + LINUX_ARM64_FORMULA_TESTS if linux_arm64 else RELEASE_TESTS
    host_env = {**os.environ, 'RSDESEQ2_NUMERIC_PROFILE': 'r-host'}
    portable_env = {**os.environ, 'RSDESEQ2_NUMERIC_PROFILE': 'portable'}
    subprocess.run(
        ['cargo', 'test', *build, '--test', 'integration', '--',
         *(f'{name}::' for name in tests)],
        cwd=ROOT, check=True, env=host_env,
    )
    subprocess.run(
        ['cargo', 'test', *build, '--test', 'integration', '--',
         'bounded_posterior::'],
        cwd=ROOT, check=True, env=portable_env,
    )
    subprocess.run(
        ['cargo', 'run', *build, '--example', 'golden', '--', '--check'],
        cwd=ROOT, check=True, env=portable_env,
    )
    if linux_arm64:
        subprocess.run(
            ['cargo', 'run', *build, '--example', 'golden', '--', '--check'],
            cwd=ROOT, check=True, env=host_env,
        )


if __name__ == '__main__':
    main()
