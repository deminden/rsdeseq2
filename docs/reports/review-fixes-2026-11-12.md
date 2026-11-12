# Review corrections after the focused refactor

Development source after `f6099e8` addresses the seven-item review. All current
Rust, R, five-platform and full-GTEx CI checks pass at `8de67af1fd3a136c5b205c2927b13fa52e5ac9dd`.
Independent affected-platform R CDF checks and Linux instrumented boundary tests
also pass. Changes are on `dev`; no tag, release or registry publication is part
of this work.

## Posterior query contract

A fitted posterior retains the coefficient identity owning its probability
summaries. `summary_coefficient()` reports it. `coefficient_summary()` rejects
a different coefficient when sign or threshold s-values are present, and checks
related matrix/vector shapes before indexing. Existing mutable result fields
remain available, but the private identity means downstream code must obtain a
posterior from the fitting APIs rather than constructing it with a struct literal.
This is an intentional source-compatibility transition for development after
0.4.0. No fitted arithmetic or original-reference bounds change.

## Formula expansion and CLI ownership

Formula interactions remove repeated variable identities while preserving their
first occurrence. `x*x`, `x:x`, `x/x`, and equivalent grouped formulas simplify
as in independent base R 4.6.2 checks. Repeated star operands are removed before
subset enumeration. Additive groups and dot placeholders keep their cross-variable
interactions: `(x+y)*(x+y)` and `.:.` include the main effects and `x:y`.
Resolved dot powers and Cartesian products use the same preflight bounds.
Generated interaction expansion is limited to 16,384 terms, input to 1 MiB and
parenthesis nesting to 128 levels; excessive requests return an error before
unbounded expansion. Formula parsing now has a real private module with explicit
parent-facing functions. CLI execution, adapters, readers and writers
also have compiler-visible private modules. Public entry points remain stable.

The CLI borrows Cook's distances from the retained original fit. Workflows
without a retained fit continue to own their diagnostic matrix separately.
The three unconditional matrix clones have been removed.

## Posterior work controls

`PosteriorAnalysisOptions::estimated_grid_evaluations()` reports batch-wide grid
likelihood work, including parameter uncertainty. MAP-only, Laplace and unselected
queries report zero. `PosteriorWorkControl` defaults to a budget of 10,000,000 grid
evaluations; callers can set another `max_grid_evaluations`, or `None` for unlimited
work. R's apeglm arguments accept `max.grid.evaluations` through the same Rust
policy; explicit `NULL` removes the limit. Work is checked before fitting.

`PosteriorCancellation` supports shared cooperative cancellation. Fits check
between genes and grid evaluation polls every 256 likelihood calls. A frontend
can attach a polling callback confined to the originating thread. R's callback
uses `R_ToplevelExec` to contain its interrupt jump and return a cancellation
status to Rust; R never jumps across an active Rust frame. R polls are throttled
to 50 ms. Cancellation returns an error rather than partial scientific results.

## Platform coverage

The existing release-target script additionally runs bounded posterior,
sparse SQP, heavy-tail, uncertainty and formula-order suites using their existing
independent fixtures and bounds. The sanitizer driver selects dynamic Clang ASan
on Darwin or Linux and instruments C with ASan/UBSan and Rust with ASan. The pinned
Linux R CI lane runs its malformed-call and forced-unwind controls and retains
logs. Cargo 1.100 moves fingerprints into `build/CRATE/HASH/fingerprint`;
the driver now checks that layout as well as the earlier `.fingerprint` layout
and retains the actual compiler-flag evidence. R, dependencies and Rust's standard library remain uninstrumented and
leak detection remains disabled, as recorded in the driver evidence.

## Independent platform rounding

Expanded platform coverage exposed one ultra-narrow uniform CDF discrepancy on
Windows and Intel macOS. The prior reference was generated on ARM macOS. An
independent R 4.6.2 Windows probe returns a component CDF of exactly 0.5, whereas
ARM R returns 0.50113895216400917 for the same three nearly coincident normal log
probabilities. The original mixture weights agree. The existing 2e-8 CDF bound
cannot establish parity across those different reference arithmetic schedules.

`generate_uniform_cdf_reference.R` evaluates the uniform-prior probability
model independently using only base R, without loading rsdeseq2 or selecting
Rust output. On ARM its 6,500 values match the complete frozen original-ashr
CDF panel within 1.12e-16. Windows and Intel macOS CI generate a fresh CDF panel
with pinned R 4.6.2. Rust checks every entry, its input values and coverage, at
the unchanged 2e-8 bound. The frozen original fixtures remain unchanged, and all
other likelihood, weight, mean and SD comparisons remain active.

For these targets, reproduce the CDF reference before running the panel (use an
absolute output path because Cargo executes integration tests from the crate):

```sh
Rscript --vanilla scripts/references/generate_uniform_cdf_reference.R "$PWD/results/uniform-cdf.tsv"
RSDESEQ2_UNIFORM_CDF_REFERENCE="$PWD/results/uniform-cdf.tsv" \
  cargo test --locked -p rsdeseq2 --test bounded_posterior
```

## Local validation

- Full Rust checkpoint: 1,157 pass, one existing opt-in case ignored. Subsequent
  resolved-dot changes pass all 68 design tests and the formula-order panel.
- Full R package checkpoint: 2,867 assertions, no test failures/warnings/skips;
  the package check retains the known compiled-symbol warning.
- Nine actual standalone/wrapped identity panels pass with R numerical calls
  blocked, including complete matrix/formula fits, later queries, posterior
  outputs, frozen-local transformations, expressions and factors.
- Actual SIGINT cancels an active native grid; the same R process then completes
  another native fit. The same probe checks configured grid-budget rejection.
  Actual-fit tests also verify that a frontend polling callback is never invoked
  from a worker thread, while a shared cancellation token reaches workers.
- Forty-nine Python tooling tests pass with the installed R/original-ashr gate
  enabled and no skips. Rust 1.100 Clippy and Rust 1.97.1 checks pass locally.
- Local Rust 1.100 ASan/UBSan boundary checks pass with the new fingerprint
  layout verified, including four replacement-owner and four formula-owner
  forced-unwind controls. Instrumentation scope remains explicit above.

## One paired performance pass

The existing eleven-contrast native and twenty-workflow R procedures use the
saved `f6099e8` artifacts as baseline. Cases alternate execution order; every
observation is retained. Own fits/builds do not overlap measurements; external
system load remains uncontrolled. No new optimization round is counted.

| Panel | Baseline | Candidate | Observed difference |
| --- | ---: | ---: | ---: |
| Eleven native contrasts | 73.115 s | 71.492 s | -2.22% |
| Maximum native process RSS | 306.609 MiB | 305.891 MiB | -0.719 MiB |
| Twenty R workflows | 180.460 s | 180.779 s | +0.18% |
| Maximum R process RSS | 1,136.703 MiB | 1,159.000 MiB | +22.297 MiB; below 2 GiB budget |

All 22 native result/diagnostic exports per executable and all twenty complete R
payloads are byte-identical. Forty-four native input/reference hashes remain
unchanged. Native times include result and diagnostic serialization and must
not replace the README's original-DESeq2 result-only comparisons. Development
memory deltas remain in this report. The single pass does not establish exact
speed preservation or explain the observed memory variation.

## Final development checks

| Workflow | Result |
| --- | --- |
| [Parity](https://github.com/deminden/rsdeseq2/actions/runs/43311556521) | All 56,200 GTEx genes and six Wald fields passed. The duplicate ordinary-oracle job is intentionally skipped; Rust CI runs those tests. |
| [Rust](https://github.com/deminden/rsdeseq2/actions/runs/43311494102) | Passed formatting, Clippy, 1,159 Rust tests and Python integrity checks. The existing opt-in full GTEx case runs separately. |
| [Platform checks](https://github.com/deminden/rsdeseq2/actions/runs/43311494595) | All five target/profile panels passed, including independent R 4.6.2 uniform CDF references on Windows and Intel macOS. |
| [R](https://github.com/deminden/rsdeseq2/actions/runs/43311494448) | Both pinned R 4.6.2 and current R passed package and native/wrapped identity checks. The pinned lane also passed actual interruption, independent Linux references, ASan/UBSan and the strict original-ashr gates. |

The Linux sanitizer artifact verifies Rust address instrumentation using the
Cargo 1.100 fingerprint layout, C address/undefined-behavior symbols, and all
eight forced-owner cleanup controls. The independent x86 CDF artifacts retain
their platform, R version, every input/CDF value and SHA-256 digest. Production
numerical and wrapper sources remain identical to the paired-measurement source
`f05cbd1`; the subsequent changes affect tests, reference generation, checks and
documentation. No frozen oracle or precision bound was changed.

Linux package checks retain the known compiled-symbol warning from the Rust
static library. Tags and publication remain on hold. The final documentation commit changes only README/docs after the
CI-tested source; earlier artifact identities remain unchanged.

Machine-readable [evidence](../data/review_fixes_2026-11-12.json),
[native observations](../data/review_fixes_2026-11-12_native.tsv) and
[R observations](../data/review_fixes_2026-11-12_R.tsv) retain every case. Raw logs,
resources, reference outputs and complete payloads are retained under the ignored
`results/review-fixes-2026-11-12/` directory. The final public-wrapper shrinkage
group also passes 71 assertions, including explicit unlimited-budget transport
and configured-budget rejection alongside original posterior parity.
