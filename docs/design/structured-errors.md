# Design note 2.1 — structured error variants

Status: **implemented in steps 2.1b and 2.1c** (2026-11-18), on Denis's request
"do error redesign".

## Problem

`DeseqError` has 13 variants. Five of them (`InvalidOptions`,
`InvalidDispersion`, `InvalidSizeFactors`, `InvalidCounts`,
`UnsupportedFeature`) carry one free-text `reason`. There are about 700 such
messages, and 417 of them are `InvalidOptions`. A caller can only tell
"bad parameter value", "bad formula", "unknown coefficient" and "internal bug"
apart by parsing English. The crate does this itself in one place:
`contrasts::aliases` matches on a message prefix to decide whether a
coefficient search failed.

## Constraints

- **Messages stay byte-identical.** 190 R `expect_error` checks and 158 Rust
  `.contains(...)` checks match message text, and so may users' scripts.
  Every new variant formats exactly the text its sites produced before.
  This makes the change IDENTICAL for numeric output and for messages.
- **`ErrorKind` stays the coarse, stable classification.** Each new variant
  maps to the kind its sites had, except internal invariants, which get their
  own kind. `ErrorKind` itself is unreleased, so adding kinds costs nothing.
- `DeseqError` is `#[non_exhaustive]`, so adding variants is not a breaking
  change. Moving existing errors out of `InvalidOptions` is a behaviour
  change for code that matches that variant; the CHANGELOG says so.

## Variants added

| Variant | Replaces | Display (unchanged) | `kind()` |
| --- | --- | --- | --- |
| `InvalidParameter { domain, parameter, requirement }` | "*name* must be finite / positive / non-negative / finite and positive / finite and non-negative / non-empty / finite and within [0, 1]" in `InvalidOptions`, `InvalidDispersion`, `InvalidSizeFactors`, `InvalidCounts` | `"{domain prefix}: {parameter} must be {requirement}"` | `domain` |
| `InvalidFormula { reason }` | `InvalidOptions` whose message starts with `formula` | `"invalid options: {reason}"` | `Formula` (new) |
| `CoefficientNotFound { reason, candidates }` | the "not present in coefficient names" family | `"invalid options: {reason}"` | `Options` |
| `Internal { what }` | `invariant(..)`, which returned `InvalidOptions` | `"invalid options: internal invariant violated: {what}"` | `Internal` (new) |

`Requirement` is a public `#[non_exhaustive]` enum whose `Display` is the
exact phrase. `domain` is an `ErrorKind` limited to `Options`, `Dispersion`,
`Normalization` and `Counts`; its prefix is the old variant's prefix.

## Step 2.1c — `MissingPrerequisite` (messages reworded, approved by Denis)

`MissingPrerequisite { domain, operation, prerequisite }` formats
`"{domain prefix}: {operation} requires {prerequisite}"` and replaces the 119
"requires" messages. 62 already had that form and keep their text. The other
57 were reworded into it:

| Before | After |
| --- | --- |
| normal/truncated-normal mixtures require a normal observation likelihood | a normal or truncated-normal mixture requires a normal observation likelihood |
| model matrix factors require valid one-based codes and at least two declared levels | each model matrix factor requires valid one-based codes and at least two declared levels |
| scale controls require one finite value per column | scale control requires one finite value per column |
| scale controls require one value per column | scale control requires one value per column |
| polynomials require finite values and positive degree | poly() requires finite values and positive degree |
| orthogonal polynomials require nonmissing finite input and positive degree | an orthogonal polynomial requires nonmissing finite input and positive degree |
| orthogonal polynomial contrasts require at most 95 levels | an orthogonal polynomial contrast requires at most 95 levels |
| parallel extrema require numeric or logical vectors | pmin/pmax requires numeric or logical vectors |
| parallel extrema require at least one input | pmin/pmax requires at least one input |
| poly coefs require names | poly(coefs =) requires names |
| numeric formula operations require numeric covariates | a numeric formula operation requires numeric covariates |
| result requests require exactly one stored Wald or LRT test | a result request requires exactly one stored Wald or LRT test |
| replacement rows require a refitted model | replacement requires a refitted model |
| fitted-object contrast results require either Wald or LRT output, not both | a fitted-object contrast result requires either Wald or LRT output, not both |
| fitted-object contrast results require stored Wald or LRT output | a fitted-object contrast result requires stored Wald or LRT output |
| at least one factor is required | a factor design requires at least one factor |
| at least one factor or numeric covariate is required | an additive design requires at least one factor or numeric covariate |
| coefficient names are required to resolve coefficient names | coefficient name resolution requires coefficient names |
| at least 3 replicates are required for outlier replacement | outlier replacement requires at least 3 replicates |
| coefficient names are required to resolve named contrasts | named contrast resolution requires coefficient names |
| stored design is required for a result request | a result request requires a stored design |
| stored LRT output is required for LRT contrast results | an LRT contrast result requires stored LRT output |
| stored {field} is required for fitted-object results | a fitted-object result requires stored {field} |
| Cook's distances are required before replacement refit | replacement refit requires Cook's distances |
| original dispersion function is required before replacement refit | replacement refit requires the original dispersion function |
| original dispersion prior variance is required before replacement refit | replacement refit requires the original dispersion prior variance |
| original log-dispersion variance is required before replacement refit | replacement refit requires the original log-dispersion variance |
| a fitted dispersion trend is required before VST | VST requires a fitted dispersion trend |
| fitted dispersion trend values are required before rlog | rlog requires fitted dispersion trend values |
| fitted dispersion trend values are required before frozen rlog | frozen rlog requires fitted dispersion trend values |
| Cook's distances are required before factor-level Cook's filtering | factor-level Cook's filtering requires Cook's distances |
| MAP dispersions are required before Wald fitting | Wald fitting requires MAP dispersions |
| MAP dispersions are required before LRT fitting | LRT fitting requires MAP dispersions |
| a fitted fast-VST dispersion trend is required | fast VST requires a fitted fast-VST dispersion trend |
| gene-wise dispersions are required before MAP fitting | MAP fitting requires gene-wise dispersions |
| fitted dispersion trend is required before MAP fitting | MAP fitting requires a fitted dispersion trend |
| fitted means are required before MAP fitting | MAP fitting requires fitted means |
| gene-wise dispersions are required before trend fitting | trend fitting requires gene-wise dispersions |
| --rlog-prior-variance is required with --frozen-intercept | --frozen-intercept requires --rlog-prior-variance |
| count gene/sample names are required to align *X* | aligning *X* requires count gene/sample names |

`kind()` is `domain`, as before.

## Design, contrast and GLM diagnostics

| Variant | Fields | Kind |
| --- | --- | --- |
| `DesignRankDeficient` | Validation context, rank, column count and zero-column labels | `Options` |
| `LrtDesignNotNested` | Full design rank and combined rank | `Options` |
| `AmbiguousCoefficient` | Requested names or levels and matching columns in model-matrix order | `Options` |
| `AmbiguousName` | Name scope, requested or conflicting declared name, matching entries | Original `Formula` or `Options` |
| `NameNotFound` | Name scope and requested name | Original `Formula` or `Options` |
| `AllZeroGene` | Zero-based fit row and original count-matrix row | `Counts` |
| `InvalidGeneValue` | Operation, zero-based gene and optional sample, rejected value | Original domain |
| `InvalidNegativeBinomialValue` | Likelihood stage, optional zero-based gene and sample, rejected value | `Numeric` for means and weights, `Dispersion` for the other stages |
| `CoxReidColumnSumOverflow` | Zero-based gene, sample and design column, finite term and overflowing partial sum | `Numeric` |

These variants preserve the messages of the errors they replace. An intercept
fit reports the sample for a rejected per-sample calculation. A failed
gene-wide reduction can have no sample.

Negative-binomial matrix likelihoods and the cached IRLS likelihood retain
the original count-matrix row. Row likelihoods have no gene index. For a
weighted term or partial sum, the sample identifies the operation that
overflowed. The Cox–Reid overflow's message still displays the finite term
from the original diagnostic. Its `value` field contains the actual sum.

Name lookup preserves exact-name priority. Ambiguous lookups list matching
entries in input order. Declaration checks stop at the first collision and
report the two entries involved. For factor levels, the scope carries the
factor name when known and distinguishes contrast operands, reference levels,
declared levels and the zero-based index of an observed sample.

## Not in this step
- `scoped_invalid_options` sites keep their `"{scope}: {reason}"` text.
- The free-text variants stay; nothing is deprecated.
