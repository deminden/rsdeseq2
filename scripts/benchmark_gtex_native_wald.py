#!/usr/bin/env python3
"""Sequential process time and peak-RSS measurements for eleven GTEx Wald contrasts.

Includes process startup, count/design loading, fitting and result-table output
for both tools. Frozen inputs and DESeq2 references are read-only. Runs alternate
tool order across repetitions; no fits run concurrently.
"""
from __future__ import annotations
import argparse
import csv
import json
import os
import platform
import shutil
import statistics
import subprocess
import time
from pathlib import Path
from check_gtex_dispersion_improvement import compare, RESULT_COLUMNS, write_rows, sha256

ROOT = Path(__file__).resolve().parents[1]

def peak_rss_bytes(path):
    """Normalize native macOS bytes and GNU time KiB to bytes."""
    for line in path.read_text().splitlines():
        line = line.strip()
        if line.endswith('maximum resident set size'):
            return int(line.split()[0])
        if line.startswith('Maximum resident set size (kbytes):'):
            return int(line.rsplit(':',1)[1].strip()) * 1024
    raise ValueError('Peak resident memory was not reported in '+str(path))

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fixture-root',type=Path,default=ROOT/'results/benchmarks/gtex_mode_contrasts')
    parser.add_argument('--output-dir',type=Path,default=ROOT/'results/benchmarks/gtex_native_wald_speed')
    parser.add_argument('--binary',type=Path,default=ROOT/'target/release/rsdeseq2')
    parser.add_argument('--rscript',default=shutil.which('Rscript'))
    parser.add_argument('--baseline-binary',type=Path,help='Optional previous Rust release for a direct release comparison')
    parser.add_argument('--baseline-label',default='previous release')
    parser.add_argument('--validated-results-root',type=Path,help='Frozen native results directory; require exact hashes for the current Rust binary')
    parser.add_argument('--repeats',type=int,default=3)
    parser.add_argument('--rust-threads',type=int,default=1)
    parser.add_argument('--blas-threads',type=int,default=1)
    parser.add_argument('--time-bin',default='/usr/bin/time',
                        help='Native macOS time or GNU time executable for per-process peak RSS')
    parser.add_argument('--contrasts',nargs='+',help='Optional named contrast subset')
    args = parser.parse_args()
    if not args.rscript or min(args.repeats,args.rust_threads,args.blas_threads) < 1:
        parser.error('Rscript and positive repeat/thread counts are required')
    contrasts = [args.fixture_root/n for n in ('condition_blocked','sex_female_vs_male','male_age4_vs_age3')]
    contrasts += sorted((args.fixture_root/'random_splits').glob('*_random_rep*'))
    if len(contrasts) != 11:
        parser.error('Generate all eleven contrasts with generate_gtex_mode_contrasts.py first')
    if args.contrasts:
        selected = set(args.contrasts)
        if not selected <= {d.name for d in contrasts}:
            parser.error('Unknown contrast name')
        contrasts = [d for d in contrasts if d.name in selected]
    args.output_dir.mkdir(parents=True,exist_ok=True)
    time_args = ['-l'] if platform.system()=='Darwin' and Path(args.time_bin).name=='time' else ['-v']
    env = os.environ.copy()
    env.update(RAYON_NUM_THREADS=str(args.rust_threads),OPENBLAS_NUM_THREADS=str(args.blas_threads),
               OMP_NUM_THREADS=str(args.blas_threads),MKL_NUM_THREADS=str(args.blas_threads))
    frozen = []
    for d in contrasts:
        for name in ('counts.tsv','design_matrix.tsv','reference_results.tsv','reference_fit_diagnostics.tsv'):
            frozen.append(dict(contrast=d.name,file=name,sha256=sha256(d/name)))
    write_rows(args.output_dir/'input_reference_sha256.tsv',frozen)
    binary_hash = sha256(args.binary)
    baseline_hash = sha256(args.baseline_binary) if args.baseline_binary else None
    subprocess.run([args.rscript,'--vanilla',str(ROOT/'scripts/benchmark_gtex_native_wald.R'),
                    '--provenance'],check=True,env=env,
                   stdout=(args.output_dir/'R-session.txt').open('w'))
    (args.output_dir/'provenance.json').write_text(json.dumps(dict(
        platform=platform.platform(),machine=platform.machine(),python=platform.python_version(),
        rust_threads=args.rust_threads,blas_threads=args.blas_threads,repeats=args.repeats,
        binary_sha256=binary_hash,rscript=args.rscript,
        baseline_binary_sha256=baseline_hash,baseline_label=args.baseline_label if args.baseline_binary else None,
        tool_order='alternates by contrast/repetition; rotates all three tools when a baseline is supplied',
        timing_scope='process startup, input parsing, full native Wald fit and result TSV output',
        memory_metric='per-process maximum resident set size, normalized to MiB (2^20 bytes)',
        memory_scope='entire fitting process, including interpreter/libraries, input data and native allocations',
        time_command=[args.time_bin,*time_args],
        cpu_model=subprocess.check_output(['sysctl','-n','machdep.cpu.brand_string'],text=True).strip()
                  if platform.system()=='Darwin' else platform.processor()),indent=2)+'\n')
    timings=[]; summary=[]; parity=[]
    for contrast_index,d in enumerate(contrasts):
        out=args.output_dir/d.name;out.mkdir(exist_ok=True)
        for repeat in range(1,args.repeats+1):
            tools=('rsdeseq2','DESeq2','rsdeseq2_baseline') if args.baseline_binary else ('rsdeseq2','DESeq2')
            offset=(contrast_index+repeat-1)%len(tools)
            order=tools[offset:]+tools[:offset]
            for tool in order:
                output=out/f'{tool}_rep{repeat:02d}.tsv'
                if tool in ('rsdeseq2','rsdeseq2_baseline'):
                    binary=args.binary if tool=='rsdeseq2' else args.baseline_binary
                    cmd=[str(binary.resolve()),'wald','--counts',str(d/'counts.tsv'),
                         '--design',str(d/'design_matrix.tsv'),'--output',str(output),
                         '--disable-cooks-cutoff','--disable-independent-filtering']
                else:
                    cmd=[args.rscript,'--vanilla',str(ROOT/'scripts/benchmark_gtex_native_wald.R'),
                         str(d/'counts.tsv'),str(d/'design_matrix.tsv'),str(output)]
                with (out/f'{tool}_rep{repeat:02d}.log').open('w') as log:
                    resources=out/f'{tool}_rep{repeat:02d}.resources.txt'
                    start=time.perf_counter()
                    subprocess.run([args.time_bin,*time_args,'-o',str(resources),*cmd],
                                   env=env,check=True,stdout=log,stderr=subprocess.STDOUT)
                    elapsed=time.perf_counter()-start
                if tool=='rsdeseq2' and args.validated_results_root:
                    expected=args.validated_results_root/d.name/'results.tsv'
                    if sha256(output)!=sha256(expected):
                        raise ValueError('Numerical regression versus frozen native output: '+str(output))
                rss=peak_rss_bytes(resources)
                timings.append(dict(contrast=d.name,tool=tool,repetition=repeat,elapsed_s=elapsed,
                                    peak_rss_bytes=rss,peak_rss_mib=rss/(1024**2),
                                    result_sha256=sha256(output)))
                write_rows(args.output_dir/'timings.tsv',timings)
                for metric in compare(output,d/'reference_results.tsv',RESULT_COLUMNS):
                    if metric['finite_na_mismatches']:
                        raise ValueError('finite/NA mismatch in '+str(output))
                    parity.append(dict(contrast=d.name,tool=tool,repetition=repeat,**metric))
                print(f'{d.name}: {tool} repeat {repeat}/{args.repeats}: '
                      f'{elapsed:.2f}s, peak RSS {rss/(1024**2):.2f} MiB',flush=True)
        rust=[r['elapsed_s'] for r in timings if r['contrast']==d.name and r['tool']=='rsdeseq2']
        reference=[r['elapsed_s'] for r in timings if r['contrast']==d.name and r['tool']=='DESeq2']
        rust_memory=[r['peak_rss_mib'] for r in timings if r['contrast']==d.name and r['tool']=='rsdeseq2']
        reference_memory=[r['peak_rss_mib'] for r in timings if r['contrast']==d.name and r['tool']=='DESeq2']
        with (d/'counts.tsv').open() as h:
            samples=len(next(csv.reader(h,delimiter='\t')))-1
        genes=next(r['n_genes'] for r in parity if r['contrast']==d.name)
        summary.append(dict(contrast=d.name,n_samples=samples,n_genes=genes,repeats=args.repeats,
            rust_threads=args.rust_threads,blas_threads=args.blas_threads,
            rust_median_s=statistics.median(rust),deseq2_median_s=statistics.median(reference),
            speedup=statistics.median(reference)/statistics.median(rust),
            rust_peak_rss_mib=statistics.median(rust_memory),
            deseq2_peak_rss_mib=statistics.median(reference_memory),
            memory_ratio=statistics.median(reference_memory)/statistics.median(rust_memory),
            rust_min_peak_rss_mib=min(rust_memory),rust_max_peak_rss_mib=max(rust_memory),
            deseq2_min_peak_rss_mib=min(reference_memory),deseq2_max_peak_rss_mib=max(reference_memory),
            rust_min_s=min(rust),rust_max_s=max(rust),deseq2_min_s=min(reference),deseq2_max_s=max(reference)))
        if args.baseline_binary:
            baseline=[r for r in timings if r['contrast']==d.name and r['tool']=='rsdeseq2_baseline']
            baseline_s=statistics.median(r['elapsed_s'] for r in baseline)
            baseline_rss=statistics.median(r['peak_rss_mib'] for r in baseline)
            summary[-1].update(rust_baseline_median_s=baseline_s,
                rust_release_speedup=baseline_s/statistics.median(rust),
                rust_baseline_peak_rss_mib=baseline_rss,
                rust_release_memory_ratio=baseline_rss/statistics.median(rust_memory))
        write_rows(args.output_dir/'summary.tsv',summary)
        write_rows(args.output_dir/'parity.tsv',parity)
    for row in frozen:
        d=next(d for d in contrasts if d.name==row['contrast'])
        if sha256(d/row['file'])!=row['sha256']:
            raise ValueError('Frozen input/reference was modified')
    if sha256(args.binary)!=binary_hash:
        raise ValueError('Binary changed during timing')
    if args.baseline_binary and sha256(args.baseline_binary)!=baseline_hash:
        raise ValueError('Baseline binary changed during timing')
    print('All frozen inputs and references unchanged.',flush=True)

if __name__=='__main__':
    main()
