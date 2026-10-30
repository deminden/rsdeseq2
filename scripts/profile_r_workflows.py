#!/usr/bin/env python3
"""One paired R workflow/RSS pass with exact output checks and a native R runtime."""
import argparse
import csv
import json
import os
import signal
from pathlib import Path
import subprocess
import time
from benchmark_gtex_native_wald import peak_rss_bytes, sha256

ROOT = Path(__file__).resolve().parents[1]
CASES = [
    ('wald',1200,16), ('multi_results',1200,16), ('lrt',1200,16),
    ('multi_results',56200,16),
    ('batch_results',1200,16), ('batch_results',56200,16),
    ('two_group_results',56200,16),
    ('offsets_size',56200,16), ('offsets_matrix',56200,16),
    ('shrink_normal',1200,16), ('shrink_ashr',1200,16), ('shrink_apeglm',1200,16),
    ('vst',1200,16), ('replacement',40,8), ('replacement_none',1200,16),
    ('replacement_none',56200,16),
    ('rlog',1024,16), ('rlog',2048,32), ('rlog',2048,64), ('rlog',4096,96),
]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rscript', required=True)
    parser.add_argument('--baseline-library', type=Path, required=True)
    parser.add_argument('--candidate-library', type=Path)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--cases', nargs='+', help='WORKFLOW_GENESxSAMPLES case names')
    parser.add_argument('--allow-design-origin-migration', action='store_true',
                        help='Validate formula origin and normalize only the schema 1→2 metadata addition')
    parser.add_argument('--time-bin', default='/usr/bin/time')
    parser.add_argument('--max-rss-mib', type=float, default=2048)
    args = parser.parse_args()
    if os.name != 'posix': parser.error('Workflow RSS profiling requires macOS or GNU time on Linux')
    args.output_dir.mkdir(parents=True, exist_ok=True)
    cases = [(f'{w}_{g}x{s}',w,g,s) for w,g,s in CASES]
    if args.cases:
        selected = set(args.cases)
        if selected - {case[0] for case in cases}:
            parser.error('Unknown workflow case')
        cases = [case for case in cases if case[0] in selected]
    libraries = dict(baseline=args.baseline_library.resolve())
    if args.candidate_library: libraries['candidate'] = args.candidate_library.resolve()
    env = dict(os.environ, RAYON_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1')
    provenance = dict(rscript=args.rscript, libraries={k:str(v) for k,v in libraries.items()}, repeats=1,
                      max_rss_mib=args.max_rss_mib, status='in_progress', paired=bool(args.candidate_library),
                      allow_design_origin_migration=args.allow_design_origin_migration,
                      comparison_script_sha256=sha256(ROOT/'scripts/compare_r_workflow_payloads.R'),
                      workflow_script_sha256=sha256(ROOT/'scripts/profile_r_workflow.R'),
                      package_files={label:{str(p.relative_to(lib)):sha256(p) for p in (lib/'rsdeseq2').rglob('*') if p.is_file()} for label,lib in libraries.items()})
    (args.output_dir/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    rows = []
    for index,(name,workflow,genes,samples) in enumerate(cases):
        order = list(libraries)
        if index % 2: order.reverse()
        outputs = {}
        for label in order:
            lib = libraries[label]; out = args.output_dir/f'{name}_{label}'
            resources = out.with_suffix('.resources.txt'); output = out.with_suffix('.rds')
            command = [args.rscript,'--vanilla',str(ROOT/'scripts/profile_r_workflow.R'),str(lib),workflow,str(genes),str(samples),str(output),str(out.with_suffix('.stages.tsv'))]
            start = time.perf_counter()
            with out.with_suffix('.log').open('w') as log:
                timed_command = [args.time_bin,'-l' if os.uname().sysname == 'Darwin' else '-v','-o',str(resources),*command]
                with subprocess.Popen(timed_command,env=env,stdout=log,stderr=subprocess.STDOUT,start_new_session=True) as proc:
                    try:
                        status = proc.wait(timeout=300)
                    except subprocess.TimeoutExpired:
                        os.killpg(proc.pid,signal.SIGKILL)
                        proc.wait()
                        raise
                    if status: raise subprocess.CalledProcessError(status,timed_command)
            elapsed = time.perf_counter()-start
            rss = peak_rss_bytes(resources)/(1024**2)
            if rss > args.max_rss_mib: raise ValueError('Workflow exceeded its pre-registered RSS budget')
            rows.append(dict(case=name,tool=label,elapsed_s=elapsed,peak_rss_mib=rss,output_sha256=sha256(output)))
            outputs[label] = output
            print(f'{name}: {label} {elapsed:.3f}s {rss:.2f} MiB',flush=True)
        if len(outputs) == 2:
            subprocess.run([args.rscript,'--vanilla',str(ROOT/'scripts/compare_r_workflow_payloads.R'),
                            str(outputs['baseline']),str(outputs['candidate']),workflow,
                            'ORIGIN_V2' if args.allow_design_origin_migration else 'STRICT',
                            str(args.output_dir/f'{name}_comparison.tsv')],check=True,env=env)
        with (args.output_dir/'runs.tsv').open('w',newline='') as handle:
            writer=csv.DictWriter(handle,fieldnames=list(rows[0]),delimiter='\t');writer.writeheader();writer.writerows(rows)
    provenance['status']='complete'
    (args.output_dir/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')


if __name__ == '__main__':
    main()
