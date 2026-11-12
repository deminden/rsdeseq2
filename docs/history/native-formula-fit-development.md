# Complete formula fits use R term order in Rust

The standalone formula workflow previously grouped factor columns ahead of numeric
columns. For `~ exposure + condition`, that changed the default last coefficient
and the numerical fitting path relative to the R wrapper. The earlier parser
checks mapped columns explicitly and therefore could not detect this complete
workflow discrepancy.

Complete native formula workflows now retain first-occurrence variable order,
sort terms by interaction degree while preserving term order, and assemble
columns through the existing shared `ModelMatrixPlan`. An interaction factor
uses contrasts when its complementary variables are contained in a preceding
term; otherwise it uses full indicators. This follows original R's `TermCode`
rule, including cases where an earlier interaction supplies the complement.
No numerical calculation delegates to R.

`poly()` remains one matrix-valued variable. Its component order is retained
when multiplying by factor indicators or other numeric variables. Expression
identity uses syntax tokens rather than equality of evaluated values. Repeated
expressions differing only in whitespace deduplicate; genuinely distinct calls
retain separate variables. Colliding compatibility aliases gain unique suffixes.
Rank validation occurs after subtraction and interaction expansion, so removed
aliases do not cause a premature rank error.

The public `formula_design_for_fitting_from_model_frame()` helper exposes the
same preparation as complete formula fits. Existing `expanded_formula_design*`
compatibility helpers retain their explicit historical column layout. Callers
preparing positional contrasts for complete formula fits should use the fitting
helper or the actual fitted design. Existing coefficient aliases are retained
except where distinct expressions require disambiguation.

## Normalization and reduced models

Formula offsets are applied once to the effective gene/sample normalization
matrix. When that matrix is present, unused sample-only size-factor metadata is
canonicalized to ones, matching the wrapped workflow. Explicit matrix-based
builder input retains its existing metadata contract.

LRT formula workflows accept matching offsets in the reduced formula. A reduced
formula without offsets remains valid when the full formula supplies them; a
nonempty reduced offset vector that differs from the full vector is rejected.
All six default/named/contrast and replacement LRT entry paths share that rule.

## Evidence and reproduction

`generate_formula_fitting_reference.R` uses only original R 4.6.2 model-frame,
model-matrix and offset calculations. Its 32 cases check complete matrix and
offset bits without column mapping: numeric/factor order, removed/readded terms,
interaction complements, no intercept, polynomial grouping, nesting, powers,
ordered transforms, dot expansion, expression aliases and offsets. Previously
committed reference files are not modified. This local checkpoint does not
claim a new independently generated Linux oracle.

```sh
Rscript --vanilla scripts/references/generate_formula_fitting_reference.R "$FRESH_REFERENCE_JSON"
cargo test --locked -p rsdeseq2 --test formula_fitting --test native_formula --test native_lrt_pipeline
cargo build --release --locked -p rsdeseq2 --example evaluate_r_interface
Rscript --vanilla scripts/validation/check_formula_fit_identity.R \
  target/release/examples/evaluate_r_interface "$INSTALLED_R_LIBRARY" "$FRESH_OUTPUT_DIR"
```

The complete formula checker supplies primitive native frame columns, never an
R-prepared design matrix. It compares actual public native fits and later
numeric/character queries with actual wrapped fits while R numerical calls are
blocked. Independent original-R matrices and controlled rank/offset rejections
remain separate from cross-interface identity. Object construction with original
DESeq2 occurs before the firewall; candidate fitting and queries run inside it.
All 33 requests (29 successful fits and four controlled rejections), with 96
later queries, agree in every encoded numerical bit. Both R CI lanes now run
this additional checker.

The existing 31 matrix-fit requests and 113 later queries also remain exact
between rebuilt interfaces and unchanged in all previously encoded fields
against `567eec3`. The 189 primitive cases, 347 focused Rust tests, 210 connected
tests after alias handling, targeted R formula/LRT/normalization tests and both
Clippy scopes pass. All 213 previously committed core reference files remain
unchanged. The [local report](../data/native-formula-fit-local.json) binds runtime
sources and artifacts.

The shared exact encoder also adds effective normalization, full-design and
reduced-design arrays to the existing matrix-workflow comparison. Old-schema
payload comparisons exclude only these three previously absent fields; all new
fields must agree between current interfaces, including signed-zero and missing
value bits. No finite tolerance or existing oracle is relaxed.

Current-source speed/RSS, full package/platform/sanitizer checks and the final
source bundle remain C23 requirements. The broader supported warning, option and
state inventory remains open under C19/C21. Ordered inputs use the explicit
`ordered()` transform in this primitive frame API; arbitrary R classes/functions
and a complete R formula grammar are not newly claimed. No optimizer, dependency,
version, publication or optimization-round claim is part of this correction.
