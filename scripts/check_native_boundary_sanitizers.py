#!/usr/bin/env python3
"""Build a separate ASan/UBSan R bridge and run isolated malformed-call tests.

Requires the selected R's Clang compiler and its dynamic ASan runtime. Rust
address instrumentation uses RUSTC_BOOTSTRAP only in this disposable build;
normal packages and the supported stable/MSRV build are unchanged. The R
interpreter, preinstalled R dependencies and Rust standard library are not
rebuilt with instrumentation. Leak detection is disabled for this R process.
"""
import argparse
import json
import os
import platform
from pathlib import Path
import shlex
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rscript', type=Path, required=True)
    parser.add_argument('--r-library', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    report = dict(status='in_progress', c='address,undefined', rust='address',
                  rust_stdlib_instrumented=False, r_interpreter_instrumented=False,
                  preinstalled_r_dependencies_instrumented=False, leak_detection=False, replacement_unwind_hooks=True, formula_unwind_hooks=True)
    report_path = out/'sanitizers.json'
    def record():
        report_path.write_text(json.dumps(report, indent=2)+'\n')
    record()
    rscript = args.rscript.resolve()
    r = rscript.with_name('R')
    compiler = shlex.split(subprocess.check_output([str(r), 'CMD', 'config', 'CC'], text=True))
    linux = sys.platform.startswith('linux')
    if linux and not any('clang' in Path(part).name for part in compiler):
        clang = shutil.which('clang')
        if clang is None:
            raise RuntimeError('Linux boundary instrumentation requires clang')
        compiler = [clang]
    try:
        resources = Path(subprocess.check_output([*compiler, '--print-resource-dir'], text=True).strip())
    except subprocess.CalledProcessError:
        report.update(status='unavailable', reason='Selected R compiler does not expose Clang sanitizer resources')
        record()
        raise
    machine = {'AMD64': 'x86_64', 'arm64': 'aarch64'}.get(platform.machine(), platform.machine())
    pattern = f'libclang_rt.asan-{machine}.so' if linux else 'libclang_rt.asan*_dynamic.dylib'
    runtimes = sorted(resources.rglob(pattern))
    if not runtimes:
        report.update(status='unavailable', reason='Selected Clang has no supported dynamic ASan runtime')
        record()
        raise RuntimeError(report['reason'])
    report['asan_runtime'] = str(runtimes[0])
    # An explicit target prevents Cargo from instrumenting host proc macros:
    # loading their ASan dylibs into the ordinary rustc process is unsupported.
    target = subprocess.check_output(['rustc', '-vV'], text=True)
    target = next(line.split(': ', 1)[1] for line in target.splitlines() if line.startswith('host: '))
    report['target'] = target
    source = out/'source'
    if source.exists():
        raise ValueError('Use a new output directory; instrumented source is immutable per run')
    shutil.copytree(ROOT/'r-pkg/rsdeseq2', source,
                    ignore=shutil.ignore_patterns('target', 'Makevars', 'Makevars.win', '*.o', '*.so', '*.dll', '*.a'))
    template = source/'src/Makevars.in'
    template.write_text(template.read_text().replace('rust/target/release/', f'rust/target/{target}/release/'))
    library = out/'library'
    library.mkdir()
    makevars = out/'Makevars.sanitizers'
    shared_runtime = ' -shared-libasan' if linux else ''
    makevars.write_text(('CC = '+shlex.join(compiler)+'\n' if linux else '')+'CFLAGS += -fsanitize=address,undefined -fno-omit-frame-pointer -DRSDESEQ2_TEST_REPLACEMENT_UNWIND -DRSDESEQ2_TEST_FORMULA_UNWIND\n'
                       'CXXFLAGS += -fsanitize=address,undefined -fno-omit-frame-pointer\n'
                       'LDFLAGS += -fsanitize=address,undefined'+shared_runtime+'\n')
    env = dict(os.environ, R_MAKEVARS_USER=str(makevars), RUSTC_BOOTSTRAP='1',
               R_LIBS_USER=os.pathsep.join([str(library), str(args.r_library.resolve())]),
               RUSTFLAGS='-Zsanitizer=address -Cforce-frame-pointers=yes', CARGO_BUILD_JOBS='2',
               CARGO_BUILD_TARGET=target,
               OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', RAYON_NUM_THREADS='1')
    report['rustflags'] = env['RUSTFLAGS']
    record()
    try:
        with (out/'install.log').open('w') as log:
            subprocess.run([str(r), 'CMD', 'INSTALL', '--preclean', '--no-test-load',
                            '--library='+str(library), str(source)], env=env, check=True,
                           stdout=log, stderr=subprocess.STDOUT, timeout=600)
        fingerprints = list((source/f'src/rust/target/{target}/release').rglob('lib-rsdeseq2_r.json'))
        # Cargo 1.100 stores these below build/CRATE/HASH/fingerprint;
        # earlier Cargo uses .fingerprint/CRATE-HASH. The disposable source
        # is fresh, so both layouts must describe this instrumented build.
        report['rust_fingerprint_evidence'] = [
            dict(path=str(p.relative_to(source)), rustflags=json.loads(p.read_text()).get('rustflags', []))
            for p in fingerprints]
        report['rust_address_instrumentation_verified'] = any(
            '-Zsanitizer=address' in json.loads(p.read_text()).get('rustflags', []) for p in fingerprints)
        dll = library/'rsdeseq2/libs/rsdeseq2.so'
        symbols = subprocess.check_output(['nm', '-u', str(dll)], text=True)
        report['asan_symbols_verified'] = '__asan_' in symbols
        report['ubsan_symbols_verified'] = '__ubsan_' in symbols
        if not all(report[key] for key in ['rust_address_instrumentation_verified', 'asan_symbols_verified', 'ubsan_symbols_verified']):
            raise ValueError('Instrumentation was not present in the built artifact')
        record()
        r_home = subprocess.check_output([str(rscript), '--vanilla', '-e', 'cat(R.home())'], text=True).strip()
        r_exec = Path(r_home)/'bin/exec/R'
        library_paths = os.pathsep.join([str(Path(r_home)/'lib'), str(rscript.parent.parent/'lib')])
        # Rscript enters R through its shell wrapper, where macOS strips DYLD_*.
        # Enter the interpreter directly and use the same path for child tests.
        env.update(R_HOME=r_home,
                   DYLD_FALLBACK_LIBRARY_PATH=library_paths,
                   RSDESEQ2_SANITIZER_LIBRARY_PATH=library_paths,
                   RSDESEQ2_NATIVE_BOUNDARY_R_EXEC=str(r_exec),
                   RSDESEQ2_SANITIZER_RUNTIME=str(runtimes[0]),
                   RSDESEQ2_REQUIRE_REPLACEMENT_HOOKS='1',
                   RSDESEQ2_REQUIRE_FORMULA_HOOKS='1',
                   ASAN_OPTIONS='detect_leaks=0:abort_on_error=1:detect_stack_use_after_return=1',
                   UBSAN_OPTIONS='halt_on_error=1:print_stacktrace=1')
        preload_var = 'LD_PRELOAD' if linux else 'DYLD_INSERT_LIBRARIES'
        library_var = 'LD_LIBRARY_PATH' if linux else 'DYLD_FALLBACK_LIBRARY_PATH'
        env[preload_var] = str(runtimes[0])
        env[library_var] = os.pathsep.join([str(runtimes[0].parent), library_paths])
        env['RSDESEQ2_SANITIZER_PRELOAD_VAR'] = preload_var
        env['RSDESEQ2_SANITIZER_LIBRARY_VAR'] = library_var
        env['RSDESEQ2_SANITIZER_LIBRARY_PATH'] = env[library_var]
        expression = ('library(rsdeseq2); testthat::test_dir(commandArgs(TRUE)[1], '
                      'filter="native-boundary",load_package="installed",package="rsdeseq2",'
                      'reporter="summary",stop_on_failure=TRUE)')
        with (out/'tests.log').open('w') as log:
            subprocess.run([str(r_exec), '--vanilla', '--slave', '-e', expression, '--args',
                            str(ROOT/'r-pkg/rsdeseq2/tests/testthat')], env=env, check=True,
                           stdout=log, stderr=subprocess.STDOUT, timeout=180)
    except Exception as error:
        report.update(status='failed', error=str(error))
        record()
        raise
    tests_log = (out/'tests.log').read_text()
    if 'Replacement unwind: 4 forced errors released their owners' not in tests_log:
        report.update(status='failed', error='Replacement unwind controls were not exercised')
        record()
        raise ValueError(report['error'])
    report['replacement_unwind_error_cases_passed'] = 4
    if 'Formula unwind: 4 forced errors released their owners' not in tests_log:
        report.update(status='failed', error='Formula unwind controls were not exercised')
        record()
        raise ValueError(report['error'])
    report['formula_unwind_error_cases_passed'] = 4
    report.update(status='passed')
    record()
    print('Instrumented native boundary tests passed; scope recorded in '+str(report_path))


if __name__ == '__main__':
    main()
