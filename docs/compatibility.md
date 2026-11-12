# Compatibility

This guide describes version 0.5.0. The Rust library, CLI and R package share
one numerical engine. R translates inputs and constructs objects; it does not
calculate statistics or provide a fallback numerical implementation.

## Supported workflows

| Area | Supported behavior |
| --- | --- |
| Normalization | `ratio` and `poscounts`, supplied size factors/geometric means, control genes, normalized counts and weighted base metadata |
| Dispersion | Automatic linear/GLM mean selection, Cox–Reid adjustment, parametric/mean/local trends, prior estimation, MAP shrinkage and outlier handling |
| GLMs | Supplied-dispersion and native-dispersion negative-binomial Wald/LRT fits, observation weights, ridge/beta-prior refits and bounded optimizer fallback |
| Results | Named coefficients, supported numeric/list/factor contrasts, Wald thresholds and supplied t degrees of freedom, Cook's policies, independent filtering and BH adjustment |
| Stored fits | Input/provenance validation, gene subsets/reordering, replacement-aware later queries and bounded result batches |
| Transforms | Mean/parametric/local VST, strict fast-VST subsets, blind/design-aware rlog, frozen intercepts and supported saved dispersion curves |
| Shrinkage | Normal, ashr and apeglm methods within the [supported contract](native-supported-contract.md) |
| R objects | Matrix/list fits, supported `DESeqDataSet` fits, `DESeqResults` and `DESeqTransform` outputs |

The [CLI guide](cli.md) specifies its numeric file interface. The [R guide](r-wrapper.md)
describes formula, fitted-object and shrinkage options. The [contract](native-supported-contract.md)
maps every R export to its Rust owner and accuracy tests.

## Formula and state semantics

Complete formula fits preserve R term order, factor coding, polynomial groups
and supported offsets. Use `formula_design_for_fitting_from_model_frame()` for
standalone designs matching that layout. The `expanded_formula_design*`
compatibility helpers expose their own explicit expanded layout.

Formula functions are restricted to declared built-ins. Custom functions are
rejected without execution. Restricted likelihood and frozen-dispersion
expressions have separate translation contracts; accepting those expressions
does not enable arbitrary functions in formulas.

Normalization-factor matrices take precedence over sample size factors. Formula
workflows canonicalize unused sample-factor metadata to ones; explicit matrix
workflows retain the supplied sample metadata. Formula offsets are composed
once from retained base normalization. Matching reduced-LRT offsets are accepted;
conflicting offsets fail.

Changed fitted sample inputs require a refit. Supported gene subsets preserve
their mapping to fitted rows and recompute subset-specific BH/filtering results.
Primitive result constructors package existing fits; they cannot reconstruct
missing counts, design, covariance or replacement context.

Cook's replacement requires both candidate genes and eligible replicate samples.
Raw all-zero flags and weighted failure flags retain their separate meanings
after replacement. The low-count Cook's rule requires matching formula/factor
and design context, including first-tied maximum handling.

## Numerical expectations

Accuracy uses independent original-R/DESeq2/ashr/apeglm references and
field-specific absolute/relative bounds. Exact equality between standalone
and R-wrapped execution is checked separately. This does not imply bitwise
equality to original DESeq2 for every dataset.

`strict`, `fast` and `diagnostic` currently use the same numerical algorithm.
Worker-count invariance is tested. Platform-specific reference arithmetic is
documented in [numerical parity](linux_numerical_parity.md).

Optimizer trajectories, dispersion line-search boundaries and local smoothing
can amplify small floating-point differences. Exact optimizer replay alone does
not establish whole-workflow parity. Low-residual-df dispersion-prior estimation
does not reproduce R's seeded Monte Carlo/loess realization exactly. Compare
normalization, gene/trend/MAP dispersions, coefficients, SEs, Cook's decisions
and filtering before diagnosing a final p-value difference.

The `poscounts` method sums positive log counts and divides by all samples,
not just the nonzero samples. Rust control-gene indices are zero-based; numeric
R indices are one-based. See [algorithms](algorithms.md) for statistical formulas.

## Unsupported scope

- Complete Bioconductor mutation/metadata semantics, arbitrary original fitted
  objects, plotting and convenience APIs outside the documented exports.
- glmGamPoi fitting and a guarantee for every weighted, unstable or local-fit
  design beyond the independently tested branches. Low-level DESeq2 ridge fits
  on rows marked `weightsFail` differ from supported high-level failure handling.
- Full R formula/class syntax, custom formula functions and splines.
- Saved local-dispersion profiles outside the standard one-dimensional quadratic
  Gaussian locfit representation, or closures outside restricted expressions.
- Normal-shrinkage formula interactions and unsupported transformed/ordered
  treatment-style factor selectors. Numeric normal contrasts require supplied results.
- Public apeglm `contrast=` coefficient selection, refinement optimizers other
  than BFGS/L-BFGS-B, and random-start configurations other than Mersenne-Twister
  or L'Ecuyer-CMRG with Inversion normals.
- Prior/likelihood combinations and active controls rejected by the shared Rust
  contract. Native `mixIP` is KKT checked; it does not reproduce an unavailable
  original Mosek execution.

Current [CI evidence](testing-parity.md#ci-status) and the
[GTEx comparison](reports/gtex-native-wald-benchmark.md) establish the tested
inputs and settings. They do not certify the complete DESeq2 API.
