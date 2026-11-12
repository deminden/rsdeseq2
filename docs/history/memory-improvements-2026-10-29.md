# Memory and workflow improvement evidence — 2026-10-29

The memory milestone is implemented; the broader plan remains in progress. New measurements use `deseq2-r462-native` (R 4.6.2,
DESeq2 1.52.0), one Rayon worker and one BLAS/OpenMP worker. Historical numerical
references and published timing records retain their original provenance.

## Baseline and workload selection

The pre-optimization checkpoint is `ff863e2`. Its executable and installed R
package are preserved under ignored `results/codebase_improvements/m2-baseline`.
[Baseline provenance](../data/memory_workflow_baseline_2026-10-29.json) records the
workload script and executable hashes. Initial process profiles were:

| R workflow | Genes × samples | Elapsed seconds | Peak RSS MiB |
| --- | ---: | ---: | ---: |
| Fit and ten selected results | 1,200 × 16 | 2.774 | 856.42 |
| Fit and design-aware rlog | 2,048 × 64 | 5.246 | 910.02 |

These are single baseline profiles, not speedup evidence. The rlog stage itself
used 3.071 seconds; the process measurement includes startup, fitting and output.
One covariance array at that rlog size is 66.02 MiB. Buffer estimates are distinct
from measured process RSS. The registered per-process workflow budget is 2 GiB.

Whole-process sampling of the blocked GTEx fit points to Cox-Reid matrix
preparation, logarithms and GLM means/IRLS; covariance/hat calculations and
allocation also appear. Sampling is diagnostic and is not used as a throughput
ratio. The sampled result table remains byte-identical to the saved release.

Criterion workloads now exercise 512-gene/16-sample IRLS and dispersions, fitted
contrast extraction, a 2,048-gene pipeline and a 256-gene/24-sample rlog transform.
`cargo bench ... -- --test` executes these workloads once as a validity check.
They replace the four-count placeholder and normalization-only pipeline scope.

## Paired workflow checks

`profile_r_workflows.py` runs separate baseline/candidate R processes in
alternating order. It checks complete saved numerical output with zero tolerance,
records internal stage time and whole-process peak RSS, and stops on errors,
changed outputs, timeouts or the RSS budget. The cases cover repeated results,
LRT, three shrinkage methods, VST, actual/no replacement and an rlog size grid.
For example, after activating the recorded R environment:

```sh
python scripts/benchmarks/profile_r_workflows.py --rscript "$CONDA_PREFIX/bin/Rscript" \
  --baseline-library results/codebase_improvements/m2-baseline/R \
  --candidate-library PATH/TO/CANDIDATE/LIBRARY \
  --output-dir results/benchmarks/memory-candidate
```

Use one accepted paired eleven-contrast native benchmark for a combined candidate,
with the reviewed native-output manifest. Count a speed round only at 2% higher
summed throughput. Memory improvements can be reported separately from that bar.
No improvement is claimed by the baseline profiles above.

## Independent validation panel

A Pasilla panel contains 2,048 genes and all seven samples, adjusting for library
type and testing treatment. Counts are selected by fixed expression strata and
hashed gene IDs, without inspecting fit differences. Wald/LRT references were
frozen with native R 4.6.2 before the memory changes. See the fixture's README,
source archive hash and generators for reproduction. Evaluate this panel at the
milestone; keep individual candidate tuning on the existing development cases.

## Implemented retention and ownership changes

C15 moves consumed fit arrays instead of copying them for fallback routing,
reduced-model assembly and replacement merging. Raw Wald/LRT statistics and
refit diagnostics are merged with the corresponding rows. C16 flattens ordered
batches with an 8 MiB temporary-output target and a 2,048-row cap; parallel errors
are reported by lowest gene index. There are no concurrent unsafe writes.

C17 computes each row's covariance and all fallback/finite checks unchanged,
then discards covariance that mean fitting and high-level rlog do not expose.
Public full-fit/rlog-fit functions retain their covariance. Weighted/unweighted
and forced-optimizer tests compare every retained field exactly. C18 passes an
explicit absent replacement buffer when replacement is disabled; C protection
and Rust pointer handling cover that representation. Enabled replacement retains
its existing buffer contract. An opaque cross-language handle for the small
no-refit buffer has not earned its added lifecycle complexity in this panel.

## Paired measurements

All thirteen baseline/candidate saved workflow outputs compare with zero
numerical tolerance and matching attributes. The complete fitted state is saved,
including replacement and raw-statistic fields. Each process stays below the
registered 2 GiB RSS budget. These are single alternating-order pairs, not
confidence intervals or repeated timing estimates.

| rlog genes × samples | Baseline seconds | Candidate seconds | Baseline MiB | Candidate MiB | RSS decrease |
| --- | ---: | ---: | ---: | ---: | ---: |
| 1,024 × 16 | 2.062 | 2.089 | 829.39 | 829.19 | 0.02% |
| 2,048 × 32 | 2.684 | 2.593 | 818.00 | 769.53 | 5.93% |
| 2,048 × 64 | 5.201 | 5.195 | 910.28 | 687.25 | 24.50% |
| 4,096 × 96 | 21.658 | 21.329 | 1,651.55 | 738.09 | 55.31% |

[All workflow measurements](../data/memory_workflows_paired_2026-10-29.tsv) include
fits, repeated results, shrinkage, VST and replacement. Small workflows remain
within 6.6 MiB of baseline RSS; fixed interpreter/library memory dominates them.
The first seven pairs completed before a profiling-only path-with-spaces bug
interrupted replacement-fixture loading. That loader was fixed and the six
remaining pairs ran once in a separate segment. Completed measurements were not
repeated. Raw segment provenance and resource logs are retained locally.

The independent Pasilla Wald/LRT panel passes the existing GTEx absolute bounds,
with identical missingness. Both tests' result tables, diagnostics, coefficients
and SE exports are byte-identical to the preserved native baseline. The
[holdout acceptance record](../data/pasilla_holdout_acceptance_2026-10-29.json)
includes frozen reference and executable hashes. References were not regenerated.
Run `scripts/validation/check_pasilla_holdout.py --help` for milestone reproduction.

The Rust workspace passed 1,022 tests with one existing ignored full-data case;
the R suite passes in the native R 4.6.2 environment. A memory improvement does
not count as a 2% total-quantification speed round. The one paired eleven-contrast pass completed against native DESeq2/R 4.6.2.
Summed process times were 36.216 seconds candidate, 35.768 seconds preserved Rust
baseline and 193.093 seconds DESeq2. The candidate is 5.33× DESeq2 and 1.24%
slower than the Rust baseline in this single pass; this is not a counted speed
round or evidence of a statistically established timing difference. Peak CLI RSS
falls on all eleven contrasts, from baseline/candidate ratios of 1.040–1.213;
candidate RSS is 213.02–301.95 MiB versus DESeq2 1,600.84–1,878.39 MiB.

All eleven candidate result exports match the reviewed v0.3.1 output hashes.
Direct native-R comparison has zero missingness mismatches in every column.
Worst absolute differences across 618,200 gene/contrast rows are: baseMean
5.59e-9, dispersion 1.24e-4, log2FC 7.15e-6, SE 2.75e-5, statistic 1.66e-4,
p-value 3.74e-5 and adjusted p-value 2.87e-5. These are errors against DESeq2,
not changes introduced by the memory candidate; native output bytes are unchanged.

Raw committed evidence: [summary](../data/memory_gtex_summary_2026-10-29.tsv),
[timings](../data/memory_gtex_timings_2026-10-29.tsv),
[direct parity](../data/memory_gtex_direct_parity_2026-10-29.tsv),
[acceptance](../data/memory_gtex_acceptance_2026-10-29.json),
[provenance](../data/memory_gtex_provenance_2026-10-29.json) and
[R session](../data/memory_gtex_R-session_2026-10-29.txt). The full frozen-reference
comparison remains in [parity](../data/memory_gtex_parity_2026-10-29.tsv).
