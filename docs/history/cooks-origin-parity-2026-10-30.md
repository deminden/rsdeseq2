# Cook's-distance ties and formula-origin parity

This checkpoint corrects two result-filtering differences found by auditing
original DESeq2 1.52.0 under R 4.6.2 (`deseq2-r462-native`). The preceding shared
policy checkpoint is `7c867f9`. Package versions remain 0.3.1; the wider
[Rust-only migration](rust-only-wrapper-migration.md) continues.

## Corrected decisions

DESeq2 selects the **first** tied maximum Cook's distance. The previous shared
kernel selected the last. When the tied samples have different counts, this
changes whether the low-count exception preserves a p-value. The corrected
kernel also treats positive and negative zero as equal when selecting a maximum.
The standalone primitive's existing finite-only diagnostic contract is retained;
the fitted R bridge continues to reject infinite Cook's values.

The low-count exception applies when the original design is a formula whose
only variable is a factor with two declared levels. Matching a two-column numeric
partition alone does not establish eligibility. Conversely, the original rule
does not require an intercept or exactly two model-matrix columns.

| Original design/context | Low-count exception |
| --- | --- |
| `~ condition`, two-level factor | Eligible |
| Same factor with a different reference level | Eligible |
| `~ relevel(condition, "B")` | Eligible |
| `~ 0 + condition` | Eligible |
| `~ numeric` or `~ I(numeric)` | Ineligible |
| `~ factor(numeric)` with an originally numeric variable | Ineligible |
| `~ condition + offset(log(exposure))` | Ineligible: two original variables |
| `~ condition + nuisance - nuisance` | Ineligible: the removed variable still belongs to the expression |
| Directly supplied numeric matrix | Ineligible |
| One factor with three levels, or `~ batch + condition` | Ineligible |

All masking, tie selection, count comparisons and numerical eligibility remain
in Rust. R extracts formula identifiers using `all.vars()` and transports the
original factor names, declared levels and codes as interface metadata. The
standalone formula helper extracts original identifiers after its numerical
parser accepts the supported formula grammar. Function names, named arguments,
string literals and constants are excluded; offsets and removed terms are kept.
Ten identifier expectations are independently checked with base R 4.6.2.

## Shared interfaces and stored fits

`DesignOrigin` records formula, matrix or explicit legacy origin separately from
numeric `DesignMatrix` values. Formula/factor design helpers attach origin;
primitive `from_row_major` matrices default to matrix origin. Native importers
that know the source formula may attach `DesignOrigin::formula(...)` and supply
the original `FormulaModelFrame`. A model frame alone does not turn a matrix into
a formula. Existing explicit low-count heuristic overrides remain available.
Clones preserve origin; matrix equality continues to compare numeric values and
coefficient names, with origin available for separate comparison.

New R fits retain `schemaVersion = 2L` and validated `designOrigin` metadata.
Saved dataset views reject changes to original formula variables even when the
numeric model matrix remains identical. Save/load, row subsets and raw matrix
fits are covered. Older dataset fits can recover origin from their retained live
formula without mutating the saved object. Older raw fit lists without a formula
retain explicit legacy partition inference; that compatibility mode does not
claim strict original-formula provenance. Unknown schemas, missing new-schema
origin and malformed origin descriptors are rejected.

The synchronous C/Rust descriptor borrows factor names and formula identifiers
for one call. No foreign pointers are retained and no R numerical callbacks are
introduced. Boundary tests cover invalid kinds, duplicate/empty/missing
identifiers, illegal matrix variables and unnamed formula factors, alongside
valid formula/matrix/legacy controls.

## Independent reference and verification

`scripts/references/generate_cooks_policy_reference.R` requires the pinned original R and
DESeq2 versions and never loads rsdeseq2. It fits 40 genes in six samples using
original `nbinomWaldTest`, then supplies explicit diagnostic Cook's values and
maxima to isolate result policy from optimizer variation. Counts, fitted
coefficients and raw p-values stay unchanged. It records exact IEEE p-value and
adjusted-p-value bits, numeric designs, original factor metadata, observed
eligibility and missingness decisions for twelve cases.

The matrix case uses the identical fitted numerical model and changes its design
representation before calling original `results()`. Eligibility is observed
from a uniquely low-count row's retained result rather than reproduced by a
local predicate. Two additional rows distinguish opposite first/last-tie
outcomes; ordinary rows retain adjustment coverage.

The frozen oracle is checked into both the Rust and R tests, with identical
bytes. The generator refuses to overwrite existing output paths; supply fresh
paths to reproduce it. The initial uncommitted eleven-case panel was extended
with the numeric-to-factor case only after verifying every preceding case and
common input remained identical. Existing committed references and precision
bounds are unchanged.

Running the preceding implementation on this controlled panel gives **20
missingness disagreements** across the twelve cases. The corrected standalone
and wrapped paths match every reference mask and retained p-value/adjustment bit.
This diagnoses the edge cases; it does not estimate their frequency in GTEx.

The final checked R source reports **Status OK**, **1,584 assertions**, no
failures, warnings or skips. Its 266 packaged files match working sources,
allowing only generated DESCRIPTION fields/DCF whitespace. The direct identity
driver passes **113** cases: all preceding 101 plus twelve strict-origin cases,
with production R numerical functions guarded against execution. Standalone
Rust tests, backend tests, Clippy, formatting, minimum Rust 1.97.1 checks, vendor
integrity, Python support tests, documentation and crate packaging pass.
Instrumented C AddressSanitizer/UndefinedBehaviorSanitizer and Rust
AddressSanitizer boundary checks pass; their recorded scope excludes rebuilding
R, installed R dependencies and the Rust standard library, and disables leak
detection.

All **121** existing GTEx contrast precision gates and **9** global gates pass.
All **22** result/diagnostic exports are byte-identical to `7c867f9`; these runs
disable Cook's filtering and therefore complement, rather than replace, the new
independent enabled-filtering oracle.

## Performance evidence and comparison scope

One paired full eleven-contrast native pass and one paired eighteen-workflow R
pass compare this candidate with `7c867f9`. Compiler, package-check,
instrumentation and conformance jobs finish before timing. Existing measurements
and previously observed costs remain unchanged; an unchanged candidate is not
rerun to replace a timing result.

The R profiler preserves complete raw RDS payloads. Its default comparison
remains strict. The explicit `--allow-design-origin-migration` option permits
only schema 1→2 plus the added origin field, after validating that origin against
the registered workflow formula and its actual numeric model matrix. It then
requires zero-tolerance equality of all remaining fields, nested fitted states,
attributes and numerical outputs. Tests prove that wrong origin/schema, changed
p-values, changed outputs, extra fields and corrupted nested states still fail.
Raw RDS hashes may differ because they retain the new metadata.

Native summed elapsed time is **35.920 → 35.984 s (+0.18%)**, with a
**5.25×** DESeq2/Rust process-time ratio (DESeq2 totals 189.087 s). Maximum
native RSS is **302.016 → 302.563 MiB (+0.547 MiB)**. Throughput changes
−0.18%; **no 2% optimization round is counted**.

The eighteen R workflows total **78.355 → 78.591 s (+0.30%)**. Maximum RSS
is **1,117.578 → 1,117.266 MiB (−0.313 MiB)**, below the 2,048 MiB budget.
All eighteen complete RDS payloads become **byte-identical serialized objects**
after only the separately validated schema/origin metadata normalization; the
raw original payloads and their different hashes are retained. Twenty-two
fitted-state origin additions are validated across those eighteen workflows,
including nested states in rlog metadata.

The full-size two-group ten-result stage is **1.507 → 1.504 s**, its process
**4.925 → 4.923 s**, and RSS **954.859 → 954.297 MiB**. Controls with costs
remain explicit: full-size individual extraction is **7.473 → 7.510 s
(+0.49%)**; the smallest rlog process is **2.165 → 2.280 s (+5.33%)**, with its
rlog stage 0.106 → 0.110 s. The largest rlog process is **21.353 → 21.304 s
(−0.23%)**, with its stage 18.851 → 18.812 s. These single observations do
not establish causes or erase earlier checkpoint costs. No measurement was
repeated to replace an observed result.

See [native summary](../data/cooks_origin_2026-10-30_native_summary.tsv),
[native process measurements](../data/cooks_origin_2026-10-30_native_timings.tsv),
[R workflows](../data/cooks_origin_2026-10-30_r_workflows.tsv),
[R stages](../data/cooks_origin_2026-10-30_r_stages.tsv),
[validated metadata migrations](../data/cooks_origin_2026-10-30_r_validated_origin_migrations.tsv),
[serialized payload equality](../data/cooks_origin_2026-10-30_r_normalized_serialization_identity.tsv),
[original-implementation reference disagreements](../data/cooks_origin_2026-10-30_reference_baseline_mismatches.tsv)
and [validation/artifact identities](../data/cooks_origin_2026-10-30_validation.json).
The [GTEx gates](../data/cooks_origin_2026-10-30_gtex_precision.tsv),
[global gates](../data/cooks_origin_2026-10-30_gtex_global_precision.tsv) and
[unchanged exports](../data/cooks_origin_2026-10-30_native_output_identity.tsv)
are separate from timing. Local immutable artifacts use the `cooks-origin-*`
prefix under `results/codebase_improvements`. Main README figures stay unchanged.

## Remaining work

This resolves the specific tie/origin gaps recorded in the preceding
[Cook/refit report](rust-only-cooks-policy-2026-10-30.md). It does not complete
the remaining normal-shrinkage zero check, offset/normalization/formula
arithmetic or optional posterior/shrinkage migration. Enabled-but-no-refit
allocation, redundant initial native result adjustments and the final combined
performance audit also remain open. Older builder coefficient/numeric paths
still choose a policy helper from the requested contrast. Their intercept,
no-intercept, mixed/scaled numeric and transformed-factor routes need the same
contrast-independent origin classification now used by stored `ResultRequest`
and the fitted R policy. This checkpoint does not claim those older builders
already have complete policy parity.

The subsequent [builder-policy correction](builder-cooks-policy-2026-10-30.md)
resolves the contrast-dependent dispatch gap recorded above across all 21 native
builder masking sites, while preserving this checkpoint's artifacts and measured
costs.
