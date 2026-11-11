"""Precision gates for the explicitly opted-in fitted provenance migration."""
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
RSCRIPT = os.environ.get('RSDESEQ2_TEST_RSCRIPT') or shutil.which('Rscript')
R_LIBRARY = os.environ.get('RSDESEQ2_TEST_RLIB')


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

    @unittest.skipUnless(R_LIBRARY, 'Set RSDESEQ2_TEST_RLIB for the installed-package SQP gate')
    def test_sqp_exception_requires_original_accuracy_and_other_fields_exact(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            source = r'''
args <- commandArgs(TRUE); directory <- args[[1]]
library(rsdeseq2,lib.loc=normalizePath(args[[2]]))
set.seed(7142)
counts <- matrix(rnbinom(60*8,mu=80,size=5),nrow=60,
  dimnames=list(paste0("g",1:60),paste0("s",1:8)))
design <- cbind(Intercept=1,conditionB=rep(c(0,1),each=4))
fit <- rsdeseq2::DESeq(counts,design=design,fitType="mean",
  sizeFactors=rep(1,8),minReplicatesForReplace=Inf,nproc=1)
raw <- rsdeseq2::results(fit,name=tail(fit$resultsNames,1))
oracle <- ashr::ash(raw$log2FoldChange,raw$lfcSE,mixcompdist="normal",method="shrink")$result
output <- data.frame(baseMean=raw$baseMean,log2FoldChange=oracle$PosteriorMean,
  lfcSE=oracle$PosteriorSD,pvalue=raw$pvalue,padj=raw$padj,row.names=rownames(raw))
a <- list(counts=counts,fitted=fit,output=output)
b <- a; a$output$log2FoldChange <- a$output$log2FoldChange+2e-8
b$output$log2FoldChange <- b$output$log2FoldChange+1e-8
saveRDS(a,file.path(directory,"a.rds"));saveRDS(b,file.path(directory,"valid.rds"))
variants <- list(
  wrong_mean=function(x) {x$output$log2FoldChange[1] <- oracle$PosteriorMean[1]+1e-3;x},
  wrong_sd=function(x) {x$output$lfcSE[1] <- oracle$PosteriorSD[1]+1e-3;x},
  missing_mean=function(x) {x$output$log2FoldChange[1] <- NA_real_;x},
  changed_pvalue=function(x) {x$output$pvalue[1] <- x$output$pvalue[1]+1e-10;x},
  changed_fit=function(x) {x$fitted$beta[1] <- x$fitted$beta[1]+1e-10;x},
  changed_counts=function(x) {x$counts[1] <- x$counts[1]+1;x},
  changed_names=function(x) {rownames(x$output) <- rev(rownames(x$output));x},
  extra_column=function(x) {x$output$unknown <- 1;x},
  extra_payload=function(x) {x$unknown <- 1;x})
for(name in names(variants))saveRDS(variants[[name]](b),file.path(directory,paste0(name,".rds")))
'''
            subprocess.run([RSCRIPT, '--vanilla', '-e', source, str(directory), R_LIBRARY],
                           check=True, capture_output=True, text=True)
            cases = [('valid', 'shrink_ashr', 'STRICT', False),
                     ('valid', 'shrink_ashr', 'ASHR_SQP', True),
                     ('valid', 'wald', 'ASHR_SQP', False)]
            cases += [(name, 'shrink_ashr', 'ASHR_SQP', False) for name in
                      ['wrong_mean', 'wrong_sd', 'missing_mean', 'changed_pvalue',
                       'changed_fit', 'changed_counts', 'changed_names', 'extra_column',
                       'extra_payload']]
            for name, workflow, mode, succeeds in cases:
                with self.subTest(name=name, workflow=workflow, mode=mode):
                    report = directory/'report.tsv'
                    report.unlink(missing_ok=True)
                    result = subprocess.run([RSCRIPT, '--vanilla',
                        str(ROOT/'scripts/compare_r_workflow_payloads.R'),
                        str(directory/'a.rds'), str(directory/f'{name}.rds'), workflow,
                        mode, str(report), R_LIBRARY], capture_output=True, text=True)
                    self.assertEqual(result.returncode == 0, succeeds, result.stdout+result.stderr)
                    if succeeds:
                        self.assertEqual(len(report.read_text().splitlines()), 3)
