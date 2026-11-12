# Review corrections after the focused refactor

Development source after `f6099e8` addresses the seven-item review. This is a
working validation record; combined acceptance and branch CI remain pending.

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
subset enumeration. Expansion is limited to 16,384 terms, input to 1 MiB and
parenthesis nesting to 128 levels; excessive requests return an error before
unbounded expansion. Formula parsing now has a real private module with seven
explicit parent-facing functions. CLI execution, adapters, readers and writers
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
logs. R, dependencies and Rust's standard library remain uninstrumented and
leak detection remains disabled, as recorded in the driver evidence.

## Validation status

Focused CLI/design/formula/posterior checks pass locally. Six formula
simplifications match independent R 4.6.2 matrices exactly. Combined validation,
actual R interruption, sanitizer execution, precision/interface identity and the
single paired speed/memory pass are still required before acceptance.
