#!/usr/bin/env python3
"""Check actual R SIGINT containment and reuse after native posterior cancellation."""
import argparse
import os
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import signal
import subprocess
import tempfile
import time

SCRIPT = r'''
args <- commandArgs(TRUE)
suppressPackageStartupMessages(library(rsdeseq2, lib.loc=args[[1L]]))
request <- getFromNamespace(".rsdeseq2_formula_request", "rsdeseq2")
y <- matrix(as.double(15:34), nrow=1L)
x <- cbind(1, diag(20L)[,1:7,drop=FALSE])
arguments <- list(Y=y,x=x,param=0.1,coef=8L,method="nbinomCR",
  interval.type="credible",ngrid=50,ngrid.nuis=5,max.grid.evaluations=NULL)
limited <- arguments
limited$max.grid.evaluations <- 1
rejected <- tryCatch(request(list(operation="apeglm_shrinkage",arguments=limited)),error=identity)
stopifnot(inherits(rejected,"error"),grepl("budget",conditionMessage(rejected)))
writeLines("READY")
flush(stdout())
result <- tryCatch(request(list(operation="apeglm_shrinkage",arguments=arguments)),error=identity)
stopifnot(inherits(result,"error"),grepl("cancelled",conditionMessage(result)))
# A second native request must succeed in the same interpreter after cleanup.
small <- arguments
small$interval.type <- "laplace"
again <- request(list(operation="apeglm_shrinkage",arguments=small))
stopifnot(is.list(again),length(again$posterior_mode)==1L)
writeLines("CANCELLED_AND_REUSED")
'''


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rscript', type=Path, required=True)
    parser.add_argument('--r-library', type=Path, required=True)
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix='rsdeseq2-interrupt-') as directory:
        script = Path(directory) / 'interrupt.R'
        script.write_text(SCRIPT)
        env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', RAYON_NUM_THREADS='1')
        process = subprocess.Popen([str(args.rscript), '--vanilla', str(script), str(args.r_library.resolve())],
                                   stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, env=env)
        try:
            # Startup includes Bioconductor loading. A bounded process timer also
            # prevents waiting forever for READY if loading hangs.
            import threading
            timer = threading.Timer(60, process.kill)
            timer.start()
            ready = process.stdout.readline().strip()
            if ready != 'READY':
                raise RuntimeError('R did not reach the native query: ' + ready + process.stderr.read())
            time.sleep(0.5)
            if process.poll() is not None:
                raise RuntimeError('Native grid terminated before the interrupt probe')
            process.send_signal(signal.SIGINT)
            stdout, stderr = process.communicate(timeout=30)
            if process.returncode != 0 or 'CANCELLED_AND_REUSED' not in stdout:
                raise RuntimeError(f'R cancellation/reuse failed ({process.returncode}):\n{stdout}\n{stderr}')
            print('Actual SIGINT cancelled the posterior and the same R process completed a second native fit')
        finally:
            timer.cancel()
            if process.poll() is None:
                process.kill()
            process.communicate()


if __name__ == '__main__':
    main()
