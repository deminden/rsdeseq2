# Development after 0.4.0

* Run canonical normal-mixture and log-F observation likelihoods entirely in
  shared Rust, including per-observation parameters, beta/gamma tails and
  integration. Import only stored primitive captures; retain source missingness,
  prior restrictions and near-one CDF cancellation behavior.

* Add native SQP active-set, line-search, Hessian-regularization and optional SVD
  controls, including automatic active-set limits. Retain independent original
  and exact-SVD references without changing accepted default outputs.

* Estimate ashr prior modes with shared Rust likelihood profiling. Reuse fitting
  preparation without posterior moments, preserve source search bounds and
  warnings, and expose explicit native intervals and search diagnostics.

* Select ashr posterior fields and fitted priors with numeric output levels or
  exact named subsets in Rust. Preserve omitted columns and require posterior
  means/prior for thresholded s-values; retain original numerical references.

* Add native cycling RRE/MPE and squared polynomial acceleration, plus a separate
  interior-point mixture solver. Retain independent original cycling and
  high-precision KKT references; direct Mosek comparison is not claimed.

* Add native ordinary, weighted, variational and compiled-style EM mixture
  fitting with K=1 acceleration controls and unit/fractional variational priors.
  Correct compiled second-probe convergence and diagnostics; retain independent
  original-package references and unchanged numerical gates.

* Add shared Rust Laplace/grid intervals, scalar parameter uncertainty and
  covariance-based posterior contrasts. Restricted grid likelihood expressions
  support primitive numeric subscripts without executing R functions. Preserve
  default/random-start output bits and document three original apeglm grid bugs.

* Run default ashr and negative-binomial apeglm shrinkage in shared Rust,
  including prior adaptation, posterior summaries, scale conversion and s-values.
  Restricted ashr likelihood expressions are imported without R execution.
  Advanced optional methods remain under development and fail explicitly.
* Retain original apeglm iteration/stopping behavior with a standalone
  MIT-licensed LBFGSpp algorithm adaptation. Ship its attribution with Rust source
  and installed R packages; the read-only vendor gate checks notice integrity.
* Add shared Rust general apeglm and analytic `nbinomR` fitting, with BFGS and
  bounded L-BFGS-B. Translate restricted five-argument likelihood expressions,
  matrix products and supported density calls without invoking R callbacks.
  Custom formula functions remain unsupported.
* Add native beta-binomial logit shrinkage: `betabinR`, `betabinCR` and MAP-only
  `betabinC`, with analytic scores, zero-success rows and concentration/trial
  validation. Preserve native optimizer warnings through the R interface.
* Correct shared BFGS metric resets, convergence probes and iteration counting;
  retain all previous posterior reference values and precision bounds.
* Add MAP-only `nbinomC*` and `betabinC*` with native random starts, portable
  MT19937/L'Ecuyer-CMRG Inversion state import/export and post-error state updates.
  Preserve all integer state-word bits. R random-generation functions are not used.

* Import standard saved DESeq2 local-dispersion coefficients without refitting.
  Construct frozen VST grids and evaluate/interpolate the curves in shared Rust.
* Translate restricted custom dispersion expressions into Rust. Reject custom
  calls, active bindings and delayed expressions before execution.
* Borrow count/normalization matrices when preparing the 999-point frozen grid,
  avoiding a full normalized matrix transfer through the primitive transport.

# rsdeseq2 0.4.0

* Match Linux original-R extended reductions, polynomial BLAS rounding and
  per-argument sum precision in Rust. Independently regenerate six exact Linux
  formula reference panels under a fixed Haswell reference kernel.
* Preserve quantification precision and paired workflow payloads. Reference and
  compatibility R CI, full GTEx replay and all five release targets pass.
* Retain the visible R CMD check warning for Rust runtime termination symbols.

* Share formula arithmetic, model-frame preparation, offsets and model-matrix
  calculations between the standalone Rust core and the R interface. Reject
  custom formula functions before execution; supported built-ins use Rust.
* Correct ordered-factor polynomial coding, factor conversion and unused-level
  handling, higher-order interaction coding and multiplication order. Reject
  treatment-style level selectors for ordered factors; use polynomial coefficient
  names or numeric contrasts. Move expanded-factor recoding into Rust.
* Evaluate Wald hypotheses, LRT probabilities, Cook's cutoffs, independent
  filtering and multiple-testing adjustments in the shared Rust engine. Preserve
  unmasked/refitted state, validate cached inputs, and apply masking before BH.
* Harden input, rank, reduced-model, normalization-offset and worker validation.
  Reduce fit copies, bound temporary results and defer replacement-count output
  allocation until a refit actually produces it.
* Refactor large implementation and test modules, consolidate repetitive tests,
  and add independent original-R fixtures plus direct cross-interface bit checks.
  Upgrade locfit, L-BFGS-B and other dependencies to the audited stable versions.
* This release does not complete the Rust-only R-wrapper migration. Frozen local
  VST reuse and optional ashr/apeglm posterior calculations still involve R.
  See the repository migration plan for the remaining work and supported scope.

# rsdeseq2 0.3.1

* Speed up native Wald quantification with shared Cox-Reid design products,
  matrix-order specialization, reused IRLS vectors and design decompositions,
  cached likelihood terms, and deferred error allocations.
* Preserve the validated 0.3.0 numerical behavior and automatic dispersion
  mean-fitting correction from `9acb44b`. The twenty accepted optimisation
  rounds retain byte-identical results and diagnostics on all eleven frozen
  GTEx contrasts.
* Synchronize the R backend with the optimized Rust core and align package,
  crate, backend, and vendored-core versions at 0.3.1.
* Extend the reproducible GTEx process benchmark with an optional previous
  Rust release, alternating tool order and direct time/peak-RSS ratios.

# rsdeseq2 0.3.0

* Native dispersion workflows select DESeq2's automatic mean-fitting branch.
  Unweighted designs with as many distinct sample rows as coefficients use
  linear means; other designs and supplied observation weights use GLM means.
* The validated correction from `9acb44b` reduces all nine relevant global
  maximum errors by at least 596.7-fold across 11 frozen GTEx Kidney Cortex
  contrasts, with preserved gene coverage and finite/NA patterns.
* Fitted `DESeqDataSet` objects receive Rust `mu` and `cooks` assays, fitted
  row metadata, and per-coefficient MLE and standard-error columns.
* Cook's outlier replacement and refitting run for eligible replicate samples.
  `minReplicatesForReplace` defaults to seven; use an integer of at least three
  or `Inf` to disable replacement. Original counts and Cook's distances are
  retained, with replacement assays and flags attached when applicable.
* Dataset result tables carry Bioconductor column descriptions and result
  metadata, including independent-filtering and replacement/refit information.
* Seeded integration tests cover automatic dispersion means, fitted means,
  Cook's distances, outlier replacement/refitting, and result metadata.

# rsdeseq2 0.2.6

* Added dataset VST/rlog transformations, frozen trend/intercept reuse, and
  selected normal, `ashr`, and `apeglm` shrinkage workflows.
