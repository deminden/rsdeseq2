"""Failure-boundary tests for the current GTEx parity and timing tools."""
import csv
import json
import math
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
import check_gtex_dispersion_improvement as parity
import benchmark_gtex_native_wald as benchmark
import generate_gtex_mode_contrasts as contrasts
import generate_gtex_compact_fixture as compact


def write(path, rows, columns=('gene','value')):
    with path.open('w',newline='') as h:
        writer=csv.writer(h,delimiter='\t',lineterminator='\n')
        writer.writerow(columns);writer.writerows(rows)


class ParityIntegrityTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.actual=Path(self.temp.name)/'actual.tsv'
        self.reference=Path(self.temp.name)/'reference.tsv'
        write(self.reference,[('g1','1'),('g2','NA')])

    def compare(self):
        return parity.compare(self.actual,self.reference,('value',))[0]

    def test_checks_values_and_missingness_in_both_directions(self):
        for rows,error,mismatches in [([('g1','1.125'),('g2','NA')],.125,0),
                                     ([('g1','NA'),('g2','2')],0,2)]:
            with self.subTest(rows=rows):
                write(self.actual,rows)
                result=self.compare()
                self.assertEqual(result['max_abs'],error)
                self.assertEqual(result['finite_na_mismatches'],mismatches)

    def test_rejects_nonfinite_numbers_instead_of_treating_them_as_na(self):
        for value in ('Inf','-Inf','inf','nan','1e999'):
            with self.subTest(value=value):
                write(self.actual,[('g1','1'),('g2',value)])
                with self.assertRaisesRegex(ValueError,'non-finite'):self.compare()

    def test_rejects_duplicate_missing_empty_or_misaligned_gene_tables(self):
        for rows in [[],[('g1','1'),('g1','2'),('g2','NA')],
                     [('g1','1')],[('g1','1'),('wrong','NA')],[('','1'),('g2','NA')]]:
            with self.subTest(rows=rows):
                write(self.actual,rows)
                with self.assertRaises(ValueError):self.compare()

    def test_rejects_duplicate_missing_columns_and_malformed_widths(self):
        for text in ['gene\tvalue\tvalue\ng1\t1\t2\n',
                     'gene\tother\ng1\t1\n','gene\tvalue\ng1\n',
                     'gene\tvalue\ng1\t1\textra\n']:
            with self.subTest(text=text):
                self.actual.write_text(text)
                with self.assertRaises(ValueError):self.compare()

    def test_all_na_has_zero_error_but_no_finite_evidence(self):
        write(self.actual,[('g1','NA')]);write(self.reference,[('g1','NA')])
        result=self.compare()
        self.assertEqual(result['max_abs'],0)
        self.assertEqual(result['n_finite'],0)

    def test_error_reduction_handles_zero_without_inventing_improvement(self):
        self.assertEqual(parity.reduction_factor(0,0),1)
        self.assertEqual(parity.reduction_factor(1,0),math.inf)
        self.assertEqual(parity.reduction_factor(2,1),2)
        for before,after in [(1,-1),(math.nan,1),(1,math.inf),(-1,1)]:
            with self.subTest(before=before,after=after):
                with self.assertRaises(ValueError):parity.reduction_factor(before,after)

    def test_global_gates_require_finite_coverage_and_no_missingness_drift(self):
        rows=[dict(category=category,column=col,before_max_abs=1,after_max_abs=.001,n_finite=10,finite_na_mismatches=0)
              for category,cols in [('dispersion_stage',parity.STAGE_COLUMNS),('result',('log2FoldChange','lfcSE','stat','pvalue','padj'))]
              for col in cols]
        self.assertTrue(all(r['passed'] for r in parity.global_maxima(rows,100)))
        for field,value in [('n_finite',0),('finite_na_mismatches',1)]:
            with self.subTest(field=field):
                changed=[dict(r) for r in rows];changed[0][field]=value
                self.assertFalse(parity.global_maxima(changed,100)[0]['passed'])
        # An invalid non-maximal row must not disappear inside max().
        for bad in (math.nan, math.inf, -1):
            with self.subTest(bad=bad):
                changed=[dict(r) for r in rows]
                changed.append(dict(rows[0], after_max_abs=bad))
                with self.assertRaises(ValueError):parity.global_maxima(changed,100)
        for factor in (0,-1,math.inf,math.nan):
            with self.assertRaises(ValueError):parity.global_maxima(rows,factor)

    def test_hash_manifest_is_required_and_exact(self):
        actual=[dict(contrast='c',file='f',sha256='abc')]
        path=Path(self.temp.name)/'hashes.tsv'
        with self.assertRaises(FileNotFoundError):parity.verify_reference_hashes(path,actual)
        for rows in [[],[('c','f','wrong')],[('c','f','abc'),('c','f','abc')],[('other','f','abc')]]:
            with self.subTest(rows=rows):
                write(path,rows,('contrast','file','sha256'))
                with self.assertRaises(ValueError):parity.verify_reference_hashes(path,actual)
        write(path,[('c','f','abc')],('contrast','file','sha256'))
        parity.verify_reference_hashes(path,actual)


class TimingIntegrityTests(unittest.TestCase):
    def test_parses_native_memory_units_and_rejects_unreported_memory(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'time.txt'
            for text,expected in [('  2048 maximum resident set size\n',2048),
                                  ('Maximum resident set size (kbytes): 2048\n',2097152)]:
                path.write_text(text)
                self.assertEqual(benchmark.peak_rss_bytes(path),expected)
            for text in ('no memory\n','-1 maximum resident set size\n','0 maximum resident set size\n'):
                path.write_text(text)
                with self.assertRaises(ValueError):benchmark.peak_rss_bytes(path)

    def test_binary_and_output_hash_guards_detect_mutations(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'file';path.write_bytes(b'original')
            digest=parity.sha256(path)
            benchmark.require_hash(path,digest,'baseline binary')
            path.write_bytes(b'changed')
            for label in ('baseline binary','binary','validated numerical output'):
                with self.assertRaisesRegex(ValueError,label):benchmark.require_hash(path,digest,label)


class CompactFixtureTests(unittest.TestCase):
    def test_panel_is_reproducible_and_retains_extremes_and_missingness(self):
        directory=ROOT/'crates/rsdeseq2/tests/data/gtex_v8_real_wald'
        counts=compact.read(directory/'counts.tsv');reference=compact.read(directory/'reference_results.tsv')
        manifest=json.loads((directory/'compact/selection.json').read_text())
        selected=compact.select_genes(counts,reference)
        self.assertEqual([dict(gene=counts[i]['gene'],source_row=i,reasons=why) for i,why in selected],manifest['genes'])
        self.assertEqual(len(selected),512)
        self.assertEqual(len(counts[0])-1,61)
        for name,digest in manifest['source_sha256'].items():self.assertEqual(parity.sha256(directory/name),digest)
        self.assertEqual([counts[i] for i,_ in selected],compact.read(directory/'compact/counts.tsv'))
        small=compact.read(directory/'compact/reference_results.tsv')
        for actual,(i,_) in zip(small,selected):
            self.assertEqual({k:v for k,v in actual.items() if k!='padj'},reference[i])
        self.assertTrue(any(r['dispersion']=='NA' for r in small))
        self.assertEqual(min(float(r['pvalue']) for r in small if r['pvalue']!='NA'), min(float(r['pvalue']) for r in reference if r['pvalue']!='NA'))
        self.assertTrue(any(float(r['pvalue'])>.9 for r in small if r['pvalue']!='NA'))
        with self.assertRaises(ValueError):compact.select_genes(counts,reference,size=1)
        with self.assertRaises(ValueError):compact.select_genes(counts,reference[::-1])


class ContrastGenerationTests(unittest.TestCase):
    def test_random_contrasts_are_deterministic_aligned_and_preserve_block_totals(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);fixture=root/'source';fixture.mkdir()
            metadata=[]
            for block,sex,age in [('female_3','female','3'),('male_3','male','3'),('male_4','male','4')]:
                for i in range(4):metadata.append((f'{block}_{i}',sex,age,block,'B' if i<2 else 'A'))
            write(fixture/'sample_metadata.tsv',metadata,('sample','sex','age_bin','perm_block','condition'))
            samples=[r[0] for r in metadata]
            write(fixture/'counts.tsv',[['g1',*range(12)],['g2',*range(20,32)]],['gene',*samples])
            outputs=[]
            for name in ('first','second'):
                output=root/name
                with mock.patch.object(sys,'argv',['generate','--fixture-dir',str(fixture),'--output-dir',str(output),'--counts-only']),mock.patch('builtins.print'):
                    contrasts.main()
                outputs.append(output)
            paths=sorted(outputs[0].rglob('design_matrix.tsv'))
            self.assertEqual(len(paths),11)
            for path in paths:
                relative=path.relative_to(outputs[0])
                self.assertEqual(path.read_bytes(),(outputs[1]/relative).read_bytes())
                design=compact.read(path);counts=compact.read(path.parent/'counts.tsv')
                self.assertEqual([r['sample'] for r in design],list(counts[0])[1:])
                if path.parent.name.startswith('blocked_random'):
                    for block in ('female_3','male_3','male_4'):
                        self.assertEqual(sum(int(r['randomB']) for r in design if r['sample'].startswith(block)),2)
                elif path.parent.name.startswith('unblocked_random'):self.assertEqual(sum(int(r['randomB']) for r in design),6)
            male=compact.read(outputs[0]/'male_age4_vs_age3/counts.tsv')
            self.assertEqual(list(male[0].values())[1:],list(map(str,range(4,12))))

    def test_existing_counts_cannot_silently_remain_stale(self):
        with tempfile.TemporaryDirectory() as directory:
            source=Path(directory)/'source';target=Path(directory)/'target'
            source.write_bytes(b'correct');target.write_bytes(b'stale')
            with self.assertRaisesRegex(ValueError,'differ'):contrasts.full_counts(source,target)
            target.write_bytes(b'correct');contrasts.full_counts(source,target)


class AcceptedBenchmarkContract(unittest.TestCase):
    def test_missing_accuracy_baseline_errors_before_running_tools(self):
        with mock.patch.object(sys, 'argv', ['benchmark', '--rscript', '/mock/Rscript']), mock.patch.object(benchmark.subprocess, 'run') as run:
            with self.assertRaises(SystemExit) as raised:
                benchmark.main()
            self.assertEqual(raised.exception.code, 2)
            run.assert_not_called()

    def test_wrong_finite_output_cannot_match_accepted_baseline(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp)/'result.tsv'
            path.write_text('gene\tpvalue\ng1\t0.01\n')
            expected = benchmark.sha256(path)
            path.write_text('gene\tpvalue\ng1\t0.02\n')
            with self.assertRaisesRegex(ValueError, 'hash differs'):
                benchmark.require_hash(path, expected, 'validated numerical output')

    def test_all_missing_and_mismatched_metrics_are_not_accepted(self):
        for finite, missing in [(0, 0), (1, 1)]:
            with self.subTest(finite=finite, missing=missing), self.assertRaises(ValueError):
                benchmark.validate_parity_metrics([dict(n_finite=finite, finite_na_mismatches=missing)])
        benchmark.validate_parity_metrics([dict(n_finite=10, finite_na_mismatches=0)])

    def test_native_baseline_is_checked_against_its_reviewed_manifest(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            d = root/'contrast'; d.mkdir()
            output = d/'results.tsv'; output.write_text('gene\tpvalue\ng1\t0.01\n')
            manifest = root/'hashes.tsv'
            with manifest.open('w', newline='') as handle:
                writer = csv.DictWriter(handle, fieldnames=['contrast','file','release_sha256'], delimiter='\t')
                writer.writeheader()
                writer.writerow(dict(contrast='contrast', file='results.tsv', release_sha256=benchmark.sha256(output)))
            benchmark.frozen_native_hashes(root, manifest, [d])
            output.write_text('gene\tpvalue\ng1\t0.02\n')
            with self.assertRaisesRegex(ValueError, 'hash differs'):
                benchmark.frozen_native_hashes(root, manifest, [d])

if __name__ == "__main__":
    unittest.main()
