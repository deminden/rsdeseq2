# Native numerical authority: current audit

The original [C01–C23 plan](codebase-improvement-plan-2026-10-22.md) and the strict
[one-engine requirement](rust-only-wrapper-migration.md) remain the scope. This
is an in-progress audit, not a declaration that C19/C21/C23 are complete.

| Area | Inspected evidence | Disposition / next action |
| --- | --- | --- |
| Production R statistical calls | AST scan of all 39 R files, now eight negative controls; quoted templates excluded; no direct R-evaluation/statistical calls found in the inspected C bridge | Removed unused Cook's arithmetic and the later [duplicate LRT df calculation](native-fit-route-audit.md). The guard now catches statistical dimension subtraction too. Runtime firewalls and complete workflow validation remain separate requirements. |
| Large-count probability/density precision | Original-only panel, independent 100-digit densities and actual native/wrapped identity | [Corrected](native-large-shape-development.md): 299 baseline failing values become zero, with unchanged bounds. No optimizer is involved. |
| Fitted result policy and bulk bridge | Inspected result-request/prepared-contrast paths, 189 primitive cases, 24 actual fresh-fit cases and 86 later queries | [Later-query checkpoint](native-replacement-query-development.md): the standalone request evaluator accepts retained replacement context, shares eligibility and missing-value formatting, and preserves exact numerical bits. Factor/formula query context, full warning/option/state coverage and combined current-source checks remain open. |
| Weight-rank failure combined with replacement | Original DESeq2 fit with supplied factors and matching native input | Confirmed open state defect: an unrelated rank-failed row's `allZero` resets after actual replacement in original DESeq2 but remains true in the merged fit. Both retain missing p-values. Correct and independently cover before declaring fitted-state parity. |
| Canonical normal-mixture captures | Original 2.2-63 constructor inspection, eight independent kernel and 18 complete profiles | [Implemented and checked](native-capture-policy-development.md): matrix captures flatten globally at equal total lengths; otherwise source-ordered component recycling and observation alignment run in Rust. Complete native/wrapped outputs agree exactly. Nonmultiple-recycling warning multiplicity remains in the broader warning audit. |
| Compiled negative-binomial option dispatch | Twelve independent original profiles with traced refinement calls, hard inputs and actual bridge identity | [Implemented and checked](native-capture-policy-development.md): MAP-only methods ignore unused controls; CR defers controls until refinement, uses source BFGS/L-BFGS-B policy and retains zero-row rejection. All 21 earlier payload sets remain exact. Other refinement optimizer families remain explicitly unsupported. This is a main-code contract correction, not a demonstrated dependency defect. |
| General/beta-binomial argument dispatch | 21 original profiles, actual refinement tracing, 76 R assertions and complete interface comparisons | [Implemented and checked](native-regression-method-policy-development.md): unused controls are deferred/ignored, scalar bounds share one Rust normalization policy, and exact-name versus abbreviated optimizer dispatch follows source behavior. All 24 earlier payload sets remain exact. |
| Other advanced method/warning/error contracts | Known explicit unsupported cases, prior per-method reports, portable MT/CMRG Inversion state handling | Complete the parameter/capture/default/warning and supported-state inventory; retain explicit, documented errors for invalid or unsupported contracts. Do not introduce an R fallback. |
| Final combined evidence | Earlier checkpoint-bound GTEx, source, sanitizer, timing/RSS and platform artifacts | Re-run required combined checks and rebuild the local source bundle only after the remaining changes are combined. Earlier artifacts do not certify the final replacement source. |

## Is rcompat-lbfgsb the current precision limiter?

The configured dependency remains `rcompat-lbfgsb 0.3.0`. No present isolated
library defect has been demonstrated. The current large-shape failures were
reproduced in probability kernels with no optimizer invocation. Six existing
hard GTEx R-optim replays and the bounded numerical-gradient check pass their
unchanged bounds. Those bounds do not imply bitwise agreement in every hard case.

Optimization remains a sensitivity source for sparse/boundary/poorly conditioned
rows. The earlier [apeglm experiment](native-apeglm-development.md) illustrates
why a smaller gradient can still produce a worse reference MAP comparison:
matching the original iteration/stopping policy is distinct from converging
further on the same mathematical objective. Default apeglm's unbounded prefit
and stopping policy are implemented in this repository; `rcompat-lbfgsb` supplies
bounded fitting/fallbacks.

Before changing the dependency, isolate one failing objective with identical
start, bounds, gradient, scaling and stopping controls. Compare objective values,
projected gradients, active bounds and traces against R's `optim`, then compare
the full quantification outputs. The present evidence prioritizes numerical
kernels and method/argument policy. It does not rule out rare library edge cases.

No current quantification speed/RSS or latest-registry claim is added by this
audit. Work remains local; pushing/tagging/publication require a new user request.
