# Memory and workflow improvement evidence — 2026-10-29

Work is in progress. New measurements use `deseq2-r462-native` (R 4.6.2,
DESeq2 1.52.0), one Rayon worker and one BLAS/OpenMP worker. Historical numerical
references and published timing records retain their original provenance.

## Baseline and workload selection

The pre-optimization checkpoint is `ff863e2`. Its executable and installed R
package are preserved under ignored `results/codebase_improvements/m2-baseline`.
[Baseline provenance](data/memory_workflow_baseline_2026-10-29.json) records the
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
python scripts/profile_r_workflows.py --rscript "$CONDA_PREFIX/bin/Rscript" \
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
