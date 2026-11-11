# Weight-failure state after replacement

Original DESeq2 initially marks a row with insufficient weighted design rank as
`allZero` to skip fitting. If another row triggers replacement, its
`refitWithoutOutliers` routine recomputes base metadata for every row. That resets
`allZero` to count-zero status before any selected rows are refitted. The
unrelated failed-weight row still has missing fitted statistics and convergence
diagnostics; it is not a newly fitted row.

The shared merged-fit implementation previously kept that initial skip flag on
unselected rows. It now refreshes count-zero flags in place when samples are
eligible and replacement occurs, then merges selected refit rows as before.
`weights_fail` stays separate and preserves missing diagnostics in later stored
queries even when `allZero` has been reset. Disabled replacement and enabled
replacement with no marked rows retain their initial flags. This also handles
replacement that creates only new-zero rows and needs no GLM refit.

No count fitting, likelihood, optimizer, dispersion formula, numerical bound or
dependency changes. The zero-flag update uses existing storage. The correction
belongs to main-code fitted-state handling, not locfit or L-BFGS-B.

## Independent original checks

`generate_weighted_replacement_state_reference.R` uses only original DESeq2
1.52.0 under R 4.6.2. Eight actual fits cover Wald and LRT, sample and matrix
normalization factors, disabled/no-refit replacement, ordinary replacement and
new-zero-only replacement. A fixed unrelated row fails the weight-rank check.
The generator asserts actual replacement totals and validates public result
columns and the correct Wald/LRT convergence columns before checking missingness.
It never imports candidate results.

The committed fixture records complete zero/failure flags, replacement totals,
the failed row's raw-weight base mean/variance and missing diagnostics. The
previous engine fails four state cases; the corrected engine passes all eight.
Counts/flags/totals and missingness are exact; base summaries retain the earlier
`1e-10 + 1e-12 * abs(reference)` bound. Both immediate and later standalone
results keep failed statistics, p-values and convergence missing.

```sh
Rscript scripts/generate_weighted_replacement_state_reference.R "$FRESH_REFERENCE_JSON"
cargo test --locked -p rsdeseq2 --test weighted_refit_metadata
```

The existing test module reuses its numeric decoder for this panel. All 211
previously committed core reference files are unchanged; the new fixture is
additional. The earlier exploratory probe established the flag discrepancy;
its generic metadata statistic/p-value lookup was insufficient to prove
missingness. This accepted generator checks the real public columns and rejects
missing or empty columns, avoiding a vacuous missingness assertion.

## Combined local evidence

The fresh-fit checker now includes weighted rank failures in both tests, with
replacement disabled, no refits, ordinary refits and new-zero-only replacement.
It asserts the intended failed-row flag transition and missing diagnostics.
All **31 fresh-fit requests and 113 later queries** agree across the rebuilt
standalone/R interfaces, including raw floating-point bits, logical fields,
filtering metadata and acceptance decisions. Locally 27 fits succeed and four
controlled cases reject. R numerical functions remain blocked.

All 153 connected Rust tests, root/backend Clippy, formatting, vendor integrity
and the production R numerical-boundary guard pass. The
[source/artifact report](data/weighted-replacement-state-local.json) binds the
626 runtime paths and before/after evidence. It adds no speed/RSS claim or
optimization round. C19/C21/C23 remain open for the broader supported-contract
audit and final combined validation/source bundle. No publication is authorized.
