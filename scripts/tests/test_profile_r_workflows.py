"""Precision gates for the explicitly opted-in fitted provenance migration."""
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
RSCRIPT = os.environ.get('RSDESEQ2_TEST_RSCRIPT') or shutil.which('Rscript')


@unittest.skipUnless(RSCRIPT, 'Rscript is needed for workflow payload comparison')
class WorkflowPayloadComparison(unittest.TestCase):
    def test_only_validated_origin_metadata_may_change(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            source = r'''
args <- commandArgs(TRUE); directory <- args[[1]]
cd <- data.frame(condition=factor(c("A","A","B","B")))
state <- structure(list(schemaVersion=1L,design=model.matrix(~condition,cd),
  colData=cd,rawPvalue=c(.1,.2),beta=matrix(1:4,2)),class="rsdeseq2RustFit")
a <- list(counts=matrix(1:8,2),fitted=state,output=list(values=c(.1,.2),metadata=list(rs=state)))
saveRDS(a,file.path(directory,"a.rds"))
b <- a
for (field in c("fitted","nested")) {
  s <- if(field=="fitted")b$fitted else b$output$metadata$rs
  s$schemaVersion <- 2L;s$designOrigin <- list(kind="formula",variables="condition")
  if(field=="fitted")b$fitted <- s else b$output$metadata$rs <- s
}
saveRDS(b,file.path(directory,"valid.rds"))
variants <- list(
  wrong_origin=function(x) {x$fitted$designOrigin$variables <- "batch"; x},
  wrong_schema=function(x) {x$fitted$schemaVersion <- 3L; x},
  changed_pvalue=function(x) {x$fitted$rawPvalue[1] <- .10000000000000002; x},
  changed_output=function(x) {x$output$values[1] <- .10000000000000002; x},
  unknown_field=function(x) {x$fitted$unknown <- TRUE; x},
  corrupt_nested=function(x) {x$output$metadata$rs$beta[1] <- 99L; x},
  wrong_matrix=function(x) {x$fitted$design[1,2] <- 1; x})
for(name in names(variants))saveRDS(variants[[name]](b),file.path(directory,paste0(name,".rds")))
'''
            subprocess.run([RSCRIPT, '--vanilla', '-e', source, str(directory)], check=True,
                           capture_output=True, text=True)
            cases = [('a', 'STRICT', True), ('valid', 'STRICT', False),
                     ('valid', 'ORIGIN_V2', True)]
            cases += [(name, 'ORIGIN_V2', False) for name in
                      ['wrong_origin', 'wrong_schema', 'changed_pvalue', 'changed_output',
                       'unknown_field', 'corrupt_nested', 'wrong_matrix']]
            for name, mode, succeeds in cases:
                with self.subTest(name=name, mode=mode):
                    report = directory/'report.tsv'
                    report.unlink(missing_ok=True)
                    result = subprocess.run([RSCRIPT, '--vanilla',
                        str(ROOT/'scripts/compare_r_workflow_payloads.R'),
                        str(directory/'a.rds'), str(directory/f'{name}.rds'),
                        'two_group_results', mode, str(report)], capture_output=True, text=True)
                    self.assertEqual(result.returncode == 0, succeeds, result.stdout+result.stderr)
                    if succeeds and mode == 'ORIGIN_V2':
                        self.assertEqual(len(report.read_text().splitlines()), 3)
