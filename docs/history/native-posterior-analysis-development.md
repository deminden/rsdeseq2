# Local posterior intervals, contrasts and parameter uncertainty

This unreleased checkpoint adds optional apeglm posterior queries to the shared
Rust engine, following the [random-start checkpoint](native-random-start-development.md).
It does not complete C19/C21/C23 or change the released 0.4.0 artifacts.

## Shared engine

`NegativeBinomialPosteriorOptions::analysis` accepts `PosteriorAnalysisOptions`:
optional Laplace, highest-density or equal-tail intervals; a coefficients-by-
contrasts matrix; one parameter standard error per gene; and an optional grid
likelihood override. Results carry natural-unit intervals and contrast modes/SDs,
plus monotonicity, tail-area and evaluation-count diagnostics. Contrasts use the
full row covariance, including off-diagonal terms. A zero contrast has zero mode
and uncertainty for a valid row. Invalid Laplace rows retain missing summaries.
MAP-only methods retain missing intervals/contrasts without evaluating a grid.

Laplace intervals share the native normal quantile already used for random
starts. Grid queries retain apeglm's Simpson normalization, rectangular mass
queries, exponential edge approximation, HPD mesh selection and nuisance-grid
ordering. Correlations are used for contrasts; the source grid itself integrates
a rectangular nuisance mesh without correlation weights. MAP and Laplace SDs
are preserved when a grid changes intervals/sign summaries. A nonmonotone grid
returns missing interval/sign summaries and an explicit diagnostic.

Nuisance combinations are streamed: the engine retains one coefficient-grid
accumulator and one candidate vector, rather than their Cartesian product.
Parameter uncertainty mixes normalized scalar-parameter scenarios using normal
weights. As in the original helper, `ngrid.nuis=1` evaluates nuisance coefficients
at their MAP and ignores the parameter variation; parameter-uncertainty requests
still return missing sign rates/s-values. Parameter SD and a threshold cannot be
combined. Parameter uncertainty currently requires one scalar parameter per gene.
All grid shape/count arithmetic is checked before fitting.

Specialized fitting retains its specialized likelihood. An explicitly supplied
restricted `log.lik` can override the **grid** likelihood without changing MAP
fitting, matching the original method separation. Restricted regression
expressions now support primitive numeric one-dimensional indexing, including
`param[1]` and `param[-1]`. Matrix single indices read column order. Indices must
be finite integers; mixed positive/negative, logical/string and multidimensional
indices are unsupported. This extension does not permit custom formula functions,
class dispatch or execution of R functions/bindings.

The R adapter imports `interval.type`, `interval.level`, `ngrid`, `nsd`,
`ngrid.nuis`, `nsd.nuis`, `contrasts` and `param.sd`. The direct primitive bridge
returns the extra matrices; `lfcShrink()` retains its coefficient result-table
shape and uses native grid sign summaries for its s-values. Extra intervals or
posterior contrasts are not added as result-table columns. With no analysis
request, grids/extra matrices are not allocated and the established default
posterior path is retained.

## Explicit original-package corrections

Original apeglm 1.34.0 has three independently reproduced edge failures:

* Combining a grid interval with a threshold stores `thresh` but reads
  `threshold` during assembly, producing a replacement-length error. Rust returns
  the grid threshold probability and applies the established sign policy.
* A single-coefficient grid with nuisance integration constructs an empty
  Cartesian product and fails. Rust evaluates its one zero-dimensional nuisance
  combination. A closed-form Gaussian posterior checks mode, SD and interval.
* A beta-binomial credible interval whose lower probability falls before the
  first grid cell drops its lower endpoint through R's zero index and recycles
  the upper endpoint. Rust inverts the same fitted exponential tail. The frozen
  original values remain present; a separately recorded R root-finder reference
  verifies the corrected bound.

These corrections are documented differences, not changes to old oracles or
precision bounds. Other original mesh approximations remain observable, including
missing probabilities when a queried threshold lies beyond the right grid edge.

## Independent validation

The original-only generator uses R 4.6.2/apeglm 1.34.0 in
`deseq2-r462-native`. Forty profiles cover 120 gene rows across all three interval
types, correlated designs, weights/offsets, alternate priors, narrow/fine grids,
nuisance-at-MAP, parameter uncertainty, NB/beta-binomial methods and explicit grid
likelihood overrides. Sixteen fixed-mode Gaussian queries isolate probability,
tail and interval arithmetic from optimizer/Hessian variation.

Existing MAP/SD/sign gates remain unchanged. New mesh intervals have a declared
one-cell comparison allowance because HPD membership can jump at nearly tied
heights; fixed-grid probability/tail comparisons use 2e-10. These interval bounds
do not establish exact quantile integration or relax existing numerical gates.

All forty complete analysis outputs match bit for bit between actual R/native
and independent R-free processes with R probability, quantile, matrix, optimizer
and random routines blocked. Four public 128-gene workflows cover HPD,
equal-tail, parameter uncertainty and a custom grid likelihood: maximum original-
apeglm errors are 2.31e-14 log2 coefficient, 2.08e-13 log2 SD and 4.48e-16 s-value.
The existing 1,200-gene default/sign/threshold/normalization checks also pass.

Immutable `e27a10a` driver replay preserves complete default and 164 random-start
posterior batches, plus all thirty RNG profiles, exactly. The extracted normal
quantile changes no accepted draws or serialized states. Focused Rust tests,
forty R likelihood-import assertions, core/backend Clippy and vendor integrity
pass. [Machine-readable evidence](../data/native-posterior-analysis-local.json)
binds source files and local artifacts.

## Reproduction and remaining scope

```sh
Rscript scripts/references/generate_posterior_analysis_reference.R "$FRESH_REFERENCE_JSON"
cargo test -p rsdeseq2 --lib fixed_grid_queries
cargo test -p rsdeseq2 --test posterior_analysis --test regression_posterior --test posterior_rng --test heavy_tail
cargo build --release --example evaluate_r_interface
Rscript scripts/validation/check_native_posterior_analysis.R "$DRIVER" "$R_LIBRARY" "$FRESH_OUTPUT" "$PAIRED_FITS_RDS"
```

Fresh paths protect prior evidence. The optional fourth checker argument supplies
a cached count-fit snapshot for the four public comparisons, avoiding another
count fit. `native_posterior_check_helpers.R` shares the statistical firewall and
primitive encodings with the general-regression checker.

Advanced ashr methods/options, other optimizer/RNG policies, complete warning/
authority reconciliation and final combined precision/speed/memory/sanitizer/
platform/source-bundle gates remain open. This opt-in query checkpoint makes no
quantification-speed or RSS claim; unchanged paths are not retimed and no 2%
optimization round is counted. All work stays local until publication is authorized.
