# Replacement execution requires eligible samples

Rust previously entered a second dispersion/GLM fit when the raw Cook planner
found outlier candidates, even if no sample was eligible for replacement. The
counts were unchanged. The R adapter hid the candidate flags in that situation,
so ordinary output comparisons did not expose the unnecessary execution.

Original DESeq2's top-level workflow skips replacement/refitting if no design
cell reaches `minReplicatesForReplace`. The shared engine now uses
`CooksRefitPlan::requires_refit()` for execution across all 16 native,
fixed-dispersion and expanded-prior fitting branches. It combines the candidate
refit condition with sample eligibility. Result merging retains the original
table when no samples are eligible; a merged fit no longer requires a nonexistent
refitted model in that case. Later requests use the same eligibility predicate.

The raw planner still reports outlier candidates in `replace`, `n_refit`,
`refit_rows`, `should_refit` and its metadata. Those public fields and primitive
calculations are unchanged. Use `requires_refit()` to decide execution; candidate
metadata does not imply an actual second fit. R's existing public no-replacement
flags and fit schema are unchanged.

This skips constructing an unnecessary second fitted model. It does not change
count replacement arithmetic, fitting objectives, solvers, thresholds, dependencies
or numerical tolerances. Eligibility is resolved outside per-gene merge loops.
No new total-speed or process-memory claim is made before combined measurement.

## Independent execution evidence

`generate_replacement_execution_reference.R` uses original DESeq2 1.52.0 under
R 4.6.2. It traces `refitWithoutOutliers` and gene-wise dispersion estimation for
Wald/LRT, scalar/matrix normalization and weighted/unweighted profiles. With four
replicates per cell and a minimum of seven, all eight profiles execute zero
refits and one gene-estimation pass. Their complete public results are identical
to replacement-disabled fits. Two positive controls with a minimum of three
execute one refit and two gene-estimation passes, confirming the instrumentation.

```sh
Rscript scripts/references/generate_replacement_execution_reference.R "$FRESH_REFERENCE_JSON"
cargo test --locked -p rsdeseq2 --test replacement_execution --test results
```

A single regression matrix applies the eight profiles to native/fixed-dispersion
and coefficient/contrast routes. All 32 cases previously constructed an
unnecessary refit. They now preserve original result rows and fitted matrix bits,
return no refit, and support subsequent replacement-aware requests. Each case
requires real candidate rows, so an empty-plan shortcut cannot make it pass.
The existing expanded-prior no-refit test also covers all four coefficient/
contrast and normalization routes without duplicating its fitting fixture.

The focused Rust checks pass 216 tests; the final loop-invariant cleanup is also
checked by the 54 directly connected tests. Root/backend Clippy, formatting and
vendor checks pass. All 212 previously committed core reference files remain
unchanged; the additional original execution fixture is separate.

The rebuilt standalone and actual wrapped engine agree exactly on all 31 fresh
fit requests and 113 later queries with R numerical calls blocked. Every complete
payload is also unchanged against `1846c48`, including missing-value bits and
filtering metadata. The [local report](../data/replacement-execution-local.json)
binds the 626 runtime paths and accepted execution/identity artifacts.

This closes the ineligible-refit finding. C19/C21 still require the broader
supported-contract and factor/formula query audit; C23 still requires final
combined precision, performance/RSS, platform/package and source-bundle checks.
No optimization round is counted and no publication is authorized.
