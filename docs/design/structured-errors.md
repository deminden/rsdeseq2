# Structured errors

Match `DeseqError` variants for field-level information and use `kind()`
for coarse classification. The enum is non-exhaustive. Coordinates in
diagnostics are zero-based unless a field documents another convention.

| Variant | Fields | Kind |
| --- | --- | --- |
| `InvalidParameter` | Error domain, parameter name and a typed `Requirement` | Supplied domain |
| `MissingPrerequisite` | Error domain, operation and missing prerequisite | Supplied domain |
| `InvalidFormula` | Formula rejection reason | `Formula` |
| `CoefficientNotFound` | Search reason and available candidates | `Options` |
| `Internal` | Violated invariant | `Internal` |

## Design, contrast and GLM diagnostics

| Variant | Fields | Kind |
| --- | --- | --- |
| `DesignRankDeficient` | Validation context, rank, column count and zero-column labels | `Options` |
| `LrtDesignNotNested` | Full design rank and combined rank | `Options` |
| `AmbiguousCoefficient` | Requested names or levels and matching columns in model-matrix order | `Options` |
| `AmbiguousName` | Name scope, requested or conflicting declared name, matching entries | `Formula` or `Options` |
| `NameNotFound` | Name scope and requested name | `Formula` or `Options` |
| `AllZeroGene` | Zero-based fit row and original count-matrix row | `Counts` |
| `InvalidGeneValue` | Operation, zero-based gene and optional sample, rejected value | Operation's error kind |
| `InvalidCalculation` | Operation, rejected value and available gene or coefficient position | Operation's error kind |
| `InvalidDispersionArithmetic` | Calculation context, rejected value and optional zero-based gene, coefficient or reduction element | `Dispersion` |
| `NoUsableDispersionTrendRows` | Minimum dispersion and number of rows checked | `Dispersion` |
| `InvalidNegativeBinomialValue` | Likelihood stage, optional zero-based gene and sample, rejected value | `Numeric` for means and weights, `Dispersion` for the other stages |
| `CoxReidColumnSumOverflow` | Zero-based gene, sample and design column, finite term and overflowing partial sum | `Numeric` |

An intercept fit reports the sample for a rejected per-sample calculation. A failed
gene-wide reduction can have no sample.

Negative-binomial matrix likelihoods and the cached IRLS likelihood retain
the original count-matrix row. Row likelihoods have no gene index. For a
weighted term or partial sum, the sample identifies the operation that
overflowed. The Cox–Reid overflow's message displays the finite term. Its
`value` field contains the actual sum.

Name lookup preserves exact-name priority. Ambiguous lookups list matching
entries in input order. Declaration checks stop at the first collision and
report the two entries involved. For factor levels, the scope carries the
factor name when known and distinguishes contrast operands, reference levels,
declared levels and the zero-based index of an observed sample.

Fitted dispersion checks report the failing input row for batched
predictions. A failed constant fit can be diagnosed before any row
is evaluated. Beta prior quantile failures include the coefficient
when called through the matrix fit. Trigamma keeps the original input
in its message and exposes the failed arithmetic result in `value`.

## Posterior queries

- `PosteriorCancelled` covers both explicit cancellation and an interrupt
  reported by the calling thread.
- `PosteriorGridBudgetExceeded` reports estimated evaluations and the budget.
  `PosteriorGridOverflow` records the requested grid sizes and the
  calculation whose count overflowed.
- `InvalidPosteriorControl` carries the rejected control and value.
  `InvalidGeneValue` carries invalid per-gene parameter uncertainties.
- `PosteriorShapeMismatch` identifies the array, axis and actual and
  expected dimensions. `PosteriorSummaryCoefficientMismatch` identifies
  the requested coefficient and the owner of the probability summaries.
- Regression errors identify missing parameters, incompatible likelihood
  matrix operands, unsupported optimizers and genes without a compiled prefit.

## Mixture optimizers

`InvalidOptimizerCount` reports the control name and rejected value.
`CompiledEmNumericMethodRequired` retains the supplied text values.
`EmPolynomialOrderExceeded` reports the requested order and retained
component count. `MixtureInputLengthMismatch` reports the expected and
actual length of the affected input.

Mixture `InvalidCalculation` errors carry the rejected numeric value.
The location can identify an observation, a component, or both axes of a
likelihood matrix entry. The standalone validation transport retains
its original text when it serializes these errors.

## Formula resources

`FormulaLimitExceeded` reports the resource, limit and observed value.
An expansion preflight can provide a lower bound because it stops counting
as soon as the limit is exceeded. Arithmetic overflow can leave the observed
value unavailable.

`UnsupportedFormulaFunction` retains the function spelling and the required
namespace when the function is recognized. `InvalidFormulaTransportLength`
reports the declared count, per-element minimum and remaining bytes.
`InvalidModelMatrixResponse` distinguishes an out-of-range index from a
factor supplied as a numeric response.

## Analysis and results

`InvalidRun` identifies the incompatible settings through `RunConflict`.
The analysis checks incompatible contrast metadata and linear-mu
replacement before evaluating a supplied formula. A beta-prior coefficient
request is checked against the standard design before constructing its
contrast.

`WaldOptionsOnStoredLrt` identifies the incompatible Wald settings.
`ReplacementRowOutOfBounds` and `EmptySampleLevel` retain the rejected
row. Beta-binomial validation reports the gene and sample for invalid
concentrations, successes or trial counts and the expected parameter shape
for mismatched matrices.

`InvalidPValue` retains the adjustment or filtering context, zero-based
index and rejected value. `UnwrittenOutput` retains the final output path
when an output transaction detects a missing staged write.
`LrtOffsetsMismatch` reports the zero-based sample of the first difference
and each formula's offset there, when available.

`MatrixShapeMismatch` reports the expected and actual `(rows, columns)`
for matrices such as normalization factors, observation weights and
fitted means. Callers can distinguish a transposed input from the expected
shape when both have the same number of elements. The error kind is `Shape`.
