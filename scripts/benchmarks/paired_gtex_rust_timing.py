#!/usr/bin/env python3
"""Paired Rust-only timing over the eleven full GTEx Wald contrasts.

Runs a candidate and an optional baseline `rsdeseq2` binary on every contrast
with the result-only command used by `benchmark_gtex_native_wald.py`, alternating
which binary goes first per contrast and repetition. Each process is measured
for wall time and peak RSS (`wait4` rusage, no external `time` binary needed).

Acceptance follows docs/development.md: a round counts only when the summed
candidate throughput is at least 2% higher than the baseline's over the same
paired pass, and the candidate's result files are byte-identical to the
baseline's (IDENTICAL tasks) unless `--allow-output-change` is given.

The fixture root must hold the reviewed contrast bundle; its input/reference
hashes are verified against docs/data/gtex_dispersion_reference_sha256.tsv.
This runner does not compare against DESeq2 or produce accepted R timings.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import platform
import statistics
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CONTRASTS = ('condition_blocked', 'sex_female_vs_male', 'male_age4_vs_age3')
REQUIRED_GAIN = 1.02


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def contrast_dirs(root: Path) -> list[Path]:
    dirs = [root / name for name in CONTRASTS]
    dirs += sorted((root / 'random_splits').glob('*_random_rep*'))
    if len(dirs) != 11 or not all(d.is_dir() for d in dirs):
        raise SystemExit(f'expected 11 contrast directories under {root}')
    return dirs


def verify_hashes(dirs: list[Path], manifest: Path) -> None:
    with manifest.open(newline='') as handle:
        expected = {(r['contrast'], r['file']): r['sha256']
                    for r in csv.DictReader(handle, delimiter='\t')}
    for d in dirs:
        for name in ('counts.tsv', 'design_matrix.tsv'):
            if sha256(d / name) != expected.get((d.name, name)):
                raise SystemExit(f'input hash differs from the frozen manifest: {d / name}')


def run_measured(cmd: list[str], env: dict[str, str], log: Path) -> tuple[float, int]:
    """Return (wall seconds, peak RSS bytes) for one child process."""
    with log.open('w') as handle:
        start = time.perf_counter()
        process = subprocess.Popen(cmd, env=env, stdout=handle, stderr=subprocess.STDOUT)
        _, status, usage = os.wait4(process.pid, 0)
        elapsed = time.perf_counter() - start
    process.returncode = os.waitstatus_to_exitcode(status)
    if process.returncode != 0:
        raise SystemExit(f'command failed ({process.returncode}): {" ".join(cmd)}; see {log}')
    scale = 1 if sys.platform == 'darwin' else 1024  # ru_maxrss: bytes on macOS, KiB on Linux
    return elapsed, usage.ru_maxrss * scale


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--binary', type=Path, required=True, help='candidate rsdeseq2 binary')
    parser.add_argument('--baseline-binary', type=Path, help='baseline rsdeseq2 binary')
    parser.add_argument('--fixture-root', type=Path,
                        default=ROOT / 'results/benchmarks/gtex_mode_contrasts')
    parser.add_argument('--reference-hashes', type=Path,
                        default=ROOT / 'docs/data/gtex_dispersion_reference_sha256.tsv')
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--repeats', type=int, default=1)
    parser.add_argument('--threads', type=int, default=1, help='RAYON_NUM_THREADS for each run')
    parser.add_argument('--allow-output-change', action='store_true',
                        help='NUMERIC tasks: do not require byte-identical result files')
    args = parser.parse_args()
    if args.repeats < 1 or args.threads < 1:
        parser.error('repeats and threads must be positive')

    dirs = contrast_dirs(args.fixture_root)
    verify_hashes(dirs, args.reference_hashes)
    tools = {'candidate': args.binary.resolve()}
    if args.baseline_binary:
        tools['baseline'] = args.baseline_binary.resolve()
    binary_hashes = {name: sha256(path) for name, path in tools.items()}
    if binary_hashes.get('baseline') == binary_hashes['candidate']:
        parser.error('candidate and baseline binaries are byte-identical')
    args.output_dir.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env['RAYON_NUM_THREADS'] = str(args.threads)

    rows = []
    for index, d in enumerate(dirs):
        out = args.output_dir / d.name
        out.mkdir(exist_ok=True)
        for repeat in range(1, args.repeats + 1):
            names = list(tools)
            if (index + repeat) % 2 and len(names) == 2:
                names.reverse()
            for name in names:
                result = out / f'{name}_rep{repeat:02d}.tsv'
                cmd = [str(tools[name]), 'wald', '--counts', str(d / 'counts.tsv'),
                       '--design', str(d / 'design_matrix.tsv'), '--output', str(result),
                       '--disable-cooks-cutoff', '--disable-independent-filtering']
                elapsed, rss = run_measured(cmd, env, out / f'{name}_rep{repeat:02d}.log')
                rows.append(dict(contrast=d.name, tool=name, repetition=repeat, order=names.index(name),
                                 elapsed_s=f'{elapsed:.4f}', peak_rss_mib=f'{rss / 2**20:.2f}',
                                 result_sha256=sha256(result)))
                print(f'{d.name} {name} rep{repeat}: {elapsed:.2f}s {rss / 2**20:.1f} MiB', flush=True)

    with (args.output_dir / 'timings.tsv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)

    summary = {}
    for name in tools:
        mine = [r for r in rows if r['tool'] == name]
        per_contrast = {d.name: statistics.median(float(r['elapsed_s']) for r in mine if r['contrast'] == d.name)
                        for d in dirs}
        summary[name] = dict(
            binary_sha256=binary_hashes[name],
            summed_median_s=round(sum(per_contrast.values()), 4),
            max_peak_rss_mib=max(float(r['peak_rss_mib']) for r in mine),
            median_peak_rss_mib=statistics.median(float(r['peak_rss_mib']) for r in mine),
        )
    verdict = None
    if 'baseline' in tools:
        identical = all(
            len({r['result_sha256'] for r in rows if r['contrast'] == d.name}) == 1 for d in dirs)
        gain = summary['baseline']['summed_median_s'] / summary['candidate']['summed_median_s']
        verdict = dict(throughput_ratio=round(gain, 4), required=REQUIRED_GAIN,
                       outputs_identical=identical,
                       accepted=gain >= REQUIRED_GAIN and (identical or args.allow_output_change))
    report = dict(platform=platform.platform(), machine=platform.machine(),
                  cpu_count=os.cpu_count(), threads=args.threads, repeats=args.repeats,
                  timing_scope='process startup, input parsing, full native Wald fit and result TSV output',
                  summary=summary, verdict=verdict)
    (args.output_dir / 'summary.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))
    if verdict and not verdict['outputs_identical'] and not args.allow_output_change:
        return 1
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
