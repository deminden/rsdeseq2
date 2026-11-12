# R wrapper

Version 0.5.0 runs supported fitting, inference, transformations and shrinkage in
the standalone Rust core. R translates arguments, aligns identifiers and builds
R/Bioconductor objects; C marshals buffers. The [contract](native-supported-contract.md)
maps every export to its numerical owner. [CI](testing-parity.md#ci-status)
checks original-package accuracy and exact standalone/wrapped equality.

## Fitting and stored objects

`DESeq()` and `fitRust()` accept counts with an explicit numeric design, a list
containing `counts`/`design`, or a supported `DESeqDataSet`. They run native
dispersion estimation and Wald or full/reduced LRT fitting. Formula metadata is
translated in R; expressions and design values are computed in Rust.
`estimateDispersionsRust()` returns final dispersions from the same pipeline.

Matrix/list inputs return `rsdeseq2RustFit`, retaining counts, design, coefficients,
covariance, normalization, test and result state. Dataset fits return the dataset
with the Rust fit in metadata, `mu`/original `cooks` assays, and implemented row
metadata: base mean/variance, all-zero flags, trend/final dispersion, convergence,
maximum Cook's distance and coefficient MLE/SE columns.

Use `rsdeseq2::results()` and `rsdeseq2::resultsNames()` for the stored fit. The
package does not register a replacement method for `DESeq2::results()`. Missing
model/count/covariance/provenance cannot be reconstructed from an arbitrary
result table. Changed fitted sample inputs require refitting; supported gene
subsets and reordering map to original fitted rows.

Normalization-factor matrices take precedence over size factors. Supplied weights,
normalization and formula offsets are forwarded to Rust. Offsets compose once
from retained base normalization, including repeated fits. Older objects without
required provenance need the original input or a refit. Formula fits canonicalize
unused sample-factor metadata to ones; explicit matrix fits retain supplied metadata.

## Results and replacement

`results()` supports named coefficients and documented numeric, list and
factor-level contrasts. `contrast` takes precedence over `name`. Wald threshold
alternatives and supplied t degrees of freedom use shared Rust probabilities.
LRT coefficient/contrast selection changes the reported LFC/SE while retaining
the model-comparison statistic/p-value; nondefault Wald-only options fail.

`cooksCutoff=TRUE` selects the automatic cutoff, `FALSE`/`NULL` disables masking,
and a nonnegative finite value supplies the cutoff. Unknown/unused arguments
fail. Independent filtering and BH adjustment are recomputed for each request,
including gene subsets. Raw and weighted all-zero/failure flags remain distinct.

Replacement requires candidate rows and eligible samples. The default threshold
is seven replicates per design cell; `minReplicatesForReplace` accepts an integer
at least three or `Inf` to disable replacement. Original counts/Cook's assays
are preserved. Replaced datasets get `replaceCounts`, `replaceCooks` and a
replacement flag; results retain eligible-sample/refit metadata. No eligible
samples means no unnecessary second fit. The low-count Cook's rule requires
matching formula/factor and design context and uses original first-tied maxima.

Standalone later queries can retain the fit and replacement plan and use
`results_for_request_with_replacement`. Primitive fits reject requests needing
context they do not contain.

`resultsList()` validates shared fitted input once and evaluates independent
queries in chunks of at most 64 requests (default 32):

```r
tables <- rsdeseq2::resultsList(dds, list(
  forward = list(contrast = c("condition", "B", "A")),
  reverse = list(contrast = c("condition", "A", "B"), alpha = 0.05),
  unmasked = list(name = tail(rsdeseq2::resultsNames(dds), 1L), cooksCutoff = FALSE)
), batchSize = 32L)
```

Entries override common arguments, including explicit `NULL`. Names/order,
classes, values and metadata match individual calls; each query keeps its own
Cook's/filtering decisions. Supported primitive Wald/LRT and list inputs are accepted.

## Formulas

Supported native formulas include declared factor/numeric columns, interactions,
subtraction, additive groups, dot expansion, offsets and built-in transforms
such as `log`, `scale`, `poly`, `I`, factor conversion and `relevel`.
Complete fits preserve original R term order, factor coding and polynomial
matrix groups. Repeated variables simplify before bounded expansion. Limits are
1 MiB input, 128 parenthesis levels and 16,384 generated interaction terms.

Custom formula functions and splines are unsupported and rejected without
execution. Arbitrary R classes and the full R formula grammar are outside scope.
Matching reduced-LRT offsets are accepted; conflicting offsets fail. Standalone
callers use `formula_design_for_fitting_from_model_frame()` for the complete-fit
layout; expanded-design compatibility helpers expose their own explicit layout.

## Transformations

Dataset `vst()`, `varianceStabilizingTransformation()` and `rlog()` return
`DESeqTransform`. They support blind/design-aware fitting, normalization factors,
weights and offsets. Fast VST selects its subset using ordinary normalized-count
means even with weights, propagates weights into fitting and errors if fewer
than `nsub` genes qualify. Design-aware rlog can reuse stored `dispFit`; frozen
gene intercepts and a shared sample prior variance are supported.

Full VST reuses saved parametric/mean coefficients or the standard
one-dimensional quadratic Gaussian DESeq2 locfit profile. Rust constructs the
999-point grid, predicts dispersions and performs integration/interpolation.
Restricted custom dispersion functions are translated without executing R:

```r
trend <- local({
  asymptote <- 0.15
  function(mu) asymptote + 0.3 / sqrt(base::pmax(mu, 1e-8))
})
attr(trend, "fitType") <- "local"
```

A function must have one mean argument and one supported expression with stored
primitive constants. Qualified built-ins avoid masking. Custom calls, multiple
statements, active bindings, delayed values and unresolved globals fail before
execution. Saved profiles outside the documented representation are unsupported.

## Shrinkage

`lfcShrink()` consumes a dataset fitted by `rsdeseq2::DESeq()` and returns
`DESeqResults` with method-specific prior metadata. Exact coefficient names and
unambiguous DESeq2-cleaned aliases are supported. `type="normal"` accepts named/
indexed coefficients and supported expanded-design factor, numeric/list contrasts.
Numeric contrasts require supplied results. List contrasts average selected
coefficients and need supplied results when those names are absent from the
unshrunk design. Normal shrinkage rejects formula interactions. Ordered factors
use polynomial coordinates; choose named coefficients or numeric contrasts rather
than treatment-style factor selectors. Transformed-factor selectors remain
unsupported when original DESeq2 cannot resolve them.

Positive normal `lfcThreshold` recomputes posterior Wald tests and filtering/BH.
Thresholded ashr/apeglm summaries use false-sign-or-small s-values.

### ashr

All supported mixture fitting, likelihood evaluation, posterior moments/SDs,
s-values and selected outputs run in Rust. The declared normal-likelihood prior
grids, normal/t observations, canonical normal-mixture/log-F and Poisson/binomial
likelihoods and restricted likelihood expressions are supported within their
explicit prior/likelihood combinations.

`optmethod` accepts `mixSQP`, `mixEM`, `w_mixEM`, `mixVBEM`, `cxxMixSquarem` and
native `mixIP`. EM supports K=1 and higher-order RRE/MPE cycling, including
`control=list(K=2, method="mpe", square=TRUE)` where applicable. Variational
unit/positive fractional priors are supported. Native `mixIP` uses a KKT-checked
log-barrier solver with `rtol`; it does not claim original Mosek execution parity.
SQP exposes active-set, direction, line-search, regularization and optional
`tol.svd` controls, including `maxiter.activeset=NULL`.

Fixed-prior (`fixg=TRUE`) fitting skips unused numerical controls but validates
the optimizer name and prior. Active fitting validates active controls. Numeric
`outputlevel` and exact named subsets select returned fields; omitted fields stay
absent. Positive-threshold requests require `PosteriorMean` and `fitted_g`, even
at level four. `mode="estimate"` profiles the fitted mixture in Rust, preserving
the supported original two-element mode behavior and explicit native intervals.

Canonical likelihood objects import stored primitive parameters without calling
closures. Noise-SD matrices and row-specific degrees of freedom keep observation
alignment. Equal total normal-mixture weight/SD lengths flatten into a global
mixture; otherwise component lists recycle against SD columns, preserving warnings.
Normal-mixture/log-F likelihoods retain their supported uniform-prior combinations
and original missingness policy. Poisson identity/log and binomial identity/logit
links use native mode/grid/FDR defaults and intensity/probability moments.

Custom `lik` lists can supply `name`, `lcdfFUN(x)`, `lpdfFUN(x)` and optional
`etruncFUN(a,b)`/`e2truncFUN(a,b)`. Each function is one supported arithmetic,
comparison, elementary-math, `ifelse`, `pmin` or `pmax` expression with primitive
scalar captures. Rust integrates the density when a first moment is absent;
without a second moment posterior SD stays missing. Reductions, arrays, arbitrary
calls, active bindings and unevaluated promises are rejected. This expression
contract does not enable custom functions in formulas.

### apeglm

Supported methods are `general`, `nbinomR`, `nbinomCR`, `nbinomC`, `nbinomC*`,
`betabinR`, `betabinCR`, `betabinC` and `betabinC*`. Public `lfcShrink` accepts a
named/indexed coefficient; its `contrast=` coefficient selector is unsupported.
Posterior covariance contrasts inside the engine are a separate option.

General `log.lik=function(y,x,beta,param,offset)` is translated from one restricted
expression. It supports the elementary subset, primitive one-dimensional indexing,
`%*%`, `as.numeric`, `as.double`, `lgamma`, `dnorm`, `dnbinom`, `dbinom` and `dt`.
Densities may be qualified with `stats::`; stored captures must be primitive
scalars. Canonical `apeglm::logLikNB` is recognized without execution.

BFGS/L-BFGS-B run in Rust when refinement requires them. MAP-only C methods ignore
optimizer, bounds and `log.link`; CR prefits defer controls until refinement is
needed. General/beta-binomial methods ignore `log.link`. Empty/missing/nonfinite
scalar bounds are unbounded, extra entries drop, and bounds apply only to exact
`optim.method="L-BFGS-B"`. CR with `log.link=FALSE` rejects all-zero rows.
Other refinement families are unsupported.

`interval.type="HPD"`/`"credible"`, `interval.level`, `ngrid`, `nsd`, `ngrid.nuis`
and `nsd.nuis` configure native grid analysis. Scalar `param.sd` requires a grid
interval without a threshold and leaves sign summaries missing. The `contrasts`
option requests native posterior covariance contrasts. Primitive analysis returns
interval/contrast matrices; public `lfcShrink()` retains its selected-coefficient
table shape. Probability summaries retain coefficient ownership and reject
mismatched coefficient/shape queries.

Grid work defaults to ten million likelihood evaluations per batch.
`max.grid.evaluations` changes the limit; explicit `NULL` removes it. Work limits
are checked before fitting and do not reduce precision.

### Random starts and interruption

MAP-only `nbinomC*`/`betabinC*` import `.Random.seed` and write back the advanced
state using native Mersenne-Twister or L'Ecuyer-CMRG with Inversion normals.
Missing state is initialized in Rust; other RNG configurations fail explicitly.
Use `set.seed()` for reproducible calls. R random-generation functions are not
called by the production engine.

Cancellation returns an error, not partial results. Native grid polling supports
shared cancellation; frontend callbacks run only on their originating thread.
R interruption is contained with `R_ToplevelExec`, so no R jump crosses an active
Rust frame. Interpreter reuse after cancellation is tested.

## Primitive results and alignment

Primitive helpers include size factors, normalized counts/base metadata,
Cook's masking, result-table assembly, independent filtering and diagnostic schemas.
Legacy `native` arguments remain accepted; they do not select R calculations.
Metadata/schema helpers assemble objects without reconstructing absent fits.

`waldFitRust()`/`lrtFitRust()` and `nbinomWaldTestRust()`/`nbinomLRTRust()` package
already-computed fields. They do not fit an unfitted dataset. Wald requests use
supplied beta/SE/covariance and required contrast context. LRT requests retain
supplied statistics/probabilities or obtain chi-square tails from Rust with `lrtDf`.
`resultsTableRust()` assembles supplied vectors and computes BH if `padj` is absent;
it does not fit models, replace counts or independently apply Cook's/filtering policy.
`applyIndependentFilteringRust()` returns filtered BH, theta/rejection counts,
selected threshold and lowess metadata from Rust.

Named vectors/matrices align to unique gene/sample/coefficient names. Duplicate,
missing or ambiguous names fail. Lookup uses exact names first, then unambiguous
R-cleaned aliases. Stored `factorReferences`, declared `factorLevels` and compatible
`colData` establish character-contrast references; declared factor order precedes
observed sample order. Named count/design/weight/normalization and covariance
axes preserve their required alignment. See installed R help for supported object
field aliases and complete argument signatures.

`baseMetadataRust()` multiplies normalized counts by raw nonnegative observation
weights before ordinary row means/sample variances; it does not modify `mcols()`.
`applyCooksCutoffRust()` can use caller-supplied low-count context but does not
infer formula semantics from a primitive table. Diagnostics helpers expose
implemented fields rather than every Bioconductor convention.

## Execution and installation

`nproc` controls reusable local worker pools for fitting, replacement and transforms.
`strict`, `fast` and `diagnostic` currently select the same numerical algorithm.
The R package statically links its vendored core and requires Cargo for source
installation. Before packaging core changes, run `scripts/release/sync_r_vendor.sh`
and its read-only `--check` mode from the repository root.
