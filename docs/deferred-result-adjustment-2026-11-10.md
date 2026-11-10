# Deferred initial builder adjustment

After the frozen-local VST checkpoint `5fef4f1`, private native builder stages
now construct raw Wald/LRT rows before Cook's policy and final independent
filtering. Their initial BH calculation was discarded: the final filtering
operation recomputes adjustment even when filtering is disabled.

The shared row kernels accept private assembly options. Public primitive
constructors and public fitted convenience queries still return adjusted
p-values immediately. The private unadjusted adapters are used only at the six
builder sites that finalize row policy afterward. This removes the discarded
initial constructor adjustment; it does not assert that every intermediate
finalization throughout all public builder workflows has been eliminated.
Wald probabilities, LRT statistics, dispersion/GLM fitting, masking, filtering
and arithmetic order are unchanged. Existing weak constructor assertions now
check numerical p-values and exact expected BH relationships.

## Validation

The [validation record](data/deferred_adjustment_validation_2026-11-10.json)
binds the numerical source tree, checked R source/library and measured CLI.
Local checks pass: 1,077 Rust tests, three backend tests, 2,634 R assertions
with `R CMD check` status OK, and 47 Python tests with one platform skip.
Both dependency graphs pass Clippy and MSRV 1.97.1; vendored code matches.
All 321 R/C/Rust files in the checked source archive match the working source.

Full GTEx replay retains all 22 result/diagnostic hashes and passes all 121
per-contrast plus nine global gates. The independent Pasilla Wald/LRT panel
passes its existing bounds and retains all eight exports. All twenty paired
R workflow payloads remain byte-identical. Complete frozen VST output and
999-point grids agree across interfaces in eight cases while R calculations
are blocked; the existing 188 expression, 54 native-parser and 60 factor/expanded
cases also pass unchanged. These checks retain independent original-DESeq2
accuracy tests in addition to cross-interface equality.

## One paired performance pass

The candidate combines frozen-local VST migration and the builder change above.
The baseline uses the 0.4.0 numerical source `d64d2f2`. Measurements use the same
local compiler, Apple M5, native R 4.6.2/DESeq2 1.52.0 and one Rust/BLAS worker.
Each case runs once per tool with alternating/rotating order. Timings include
process startup, input and output. RSS covers each entire process. No case was
discarded or repeated; this pass does not estimate timing variation or isolate
the contribution of either individual change.

| Panel | Release baseline time | Combined candidate time | Baseline peak RSS | Candidate peak RSS |
| --- | ---: | ---: | ---: | ---: |
| Native eleven contrasts | 35.911 s | 35.515 s | 301.891 MiB | 301.938 MiB |
| Twenty R workflows | 87.055 s | 87.170 s | 1,131.188 MiB | 1,131.250 MiB |

The native total shows **1.116% higher throughput** (1.103% lower elapsed time).
This is below the user's 2% requirement, so **zero optimization rounds count**.
Original DESeq2 takes 188.803 s with maximum RSS 1,764.391 MiB in the same native
panel: a 5.316× time ratio. All R workflows stay below the 2 GiB budget.
The earlier frozen-local prototype's slower measurements remain in its
[audit](frozen-local-vst-audit-2026-10-31.md); this new observation does not
erase them or certify the full plan's final combined performance gate.

Raw retained tables: [native summary](data/deferred_adjustment_native_2026-11-10.tsv),
[native individual runs](data/deferred_adjustment_native_runs_2026-11-10.tsv),
and [R workflow runs](data/deferred_adjustment_r_workflows_2026-11-10.tsv).
The public README keeps its released comparison against original DESeq2.

## Reproduction

Build both numerical sources with the same compiler. Keep the baseline binary
and installed R library separate from the candidate. Generate the documented
GTEx fixtures and validate exact candidate exports before accepted timing:

```sh
python3 scripts/benchmark_gtex_native_wald.py \
  --binary CANDIDATE_CLI --baseline-binary BASELINE_CLI \
  --baseline-label '0.4.0 numerical source d64d2f2' \
  --validated-results-root VALIDATED_GTEX_OUTPUT \
  --rscript R462SCRIPT --repeats 1 --rust-threads 1 --blas-threads 1 \
  --output-dir /tmp/deferred-adjustment-native
python3 scripts/profile_r_workflows.py \
  --rscript R462SCRIPT --baseline-library BASELINE_LIBRARY \
  --candidate-library CANDIDATE_LIBRARY --output-dir /tmp/deferred-adjustment-r
```

See [the contrast-generation documentation](gtex-native-wald-benchmark.md) for fixture
generation. Optional posterior engines, general method/warning semantics and
the final result-authority/performance/artifact audit remain open; this
checkpoint does not complete C19/C21/C23 or change the published 0.4.0 release.
