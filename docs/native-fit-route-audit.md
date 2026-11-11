# Fit/result authority audit: local checkpoint

The [C19–C23 plan](codebase-improvement-plan-2026-10-22.md) requires shared
statistical decisions as well as preserved public APIs. This checkpoint removes
one remaining R statistical calculation and one duplicate fit-dispatch decision.
It does not close the full method/state/warning inventory or final combined gates.

## Corrections

The R wrapper previously calculated `ncol(design) - ncol(reduced)` for stored
LRT degrees of freedom. It now imports the actual `LrtOutput` value produced by
the standalone Rust fit. The private C/Rust fit transport gains one checked
integer output; Wald exposes NULL and LRT retains its existing integer type.
The registered `.Call` argument count, fitted-object schema and public fields
remain unchanged. C attaches the scalar without computing the statistic.

Ordinary fits now enter `DeseqBuilder::fit_with_results` for both Wald and LRT.
The adapter stores the reduced design once and lets the core choose the fit
route, removing its repeated test dispatch and reduced-design clone. Replacement
fits already use the core's test-aware replacement entry point. No core numerical
operation, optimizer setting or dependency changes.

## Inspected authority paths

| Production responsibility | Authoritative implementation | Adapter role |
| --- | --- | --- |
| Ordinary Wald/LRT fitting | Core `builder/workflows/entry.rs` and native fit stages | Argument/object translation, buffer validation and forwarding to the core dispatcher |
| Replacement/refit fitting | Core `builder/workflows/replacement.rs` and Cook/refit kernels | Copies selected refitted storage and attaches metadata; numerical decisions remain in core |
| LRT degrees of freedom | Core `glm/lrt.rs`, stored `LrtOutput` | Checked scalar transport; no R recomputation |
| Stored result requests | `results/request.rs`, `core/fitted_results.rs`, shared GLM projection/inference and Cook/filter kernels | Selector, row alignment and presentation metadata |
| Single/bulk R projection | Shared `PreparedContrastView` and `ContrastZeroPolicy` | Borrowed buffers and bounded request chunks |
| Distributions, adjustments and filtering | Core probability, inference and independent-filtering modules | No direct R statistical calls or R numerical imports in inspected bridge paths |

The remaining R arithmetic search finds indexing, enum translation, shape and
label handling, plus quoted likelihood syntax. The demonstrated statistical
dimension subtraction is now rejected by the lexical guard. Its eighth negative
control catches the former expression; a positive control keeps zero-based
index translation valid. All 39 production R files pass. This lexical scope
does not prove every possible future alias or expression safe; actual runtime
comparisons remain separately required.

## Validation

- All **189** actual result-interface cases agree exactly with a separate R-free
  driver while scoped R numerical calls are blocked.
- Existing fit integration and isolated malformed-boundary tests pass. LRT
  coverage now checks one and two degrees of freedom, exact integer type,
  independent chi-square probabilities and selected-coefficient behavior. Wald
  retains NULL degrees of freedom.
- One paired pass covers Wald, LRT, actual replacement and enabled replacement
  with no refitted rows. All four complete R workflow payloads, including fitted
  arrays, metadata and missingness, are byte-identical to `39021c1`.
- Those four candidate processes peak at **838.203 MiB**, below the 2 GiB budget.
  Summed time is 9.104 s versus 9.000 s, an observed **1.16% increase** in this
  small targeted pass. Maximum RSS changes by 0.078 MiB. These observations make
  no total-quantification speed claim and count as no optimization round.
- Backend Clippy and vendor integrity pass. All 210 existing oracle files and
  the standalone core runtime sources/executable remain unchanged.

```sh
Rscript scripts/audit_r_numerical_boundary.R r-pkg/rsdeseq2/R "$FRESH_AUDIT_JSON"
Rscript scripts/check_rust_r_identity.R "$STANDALONE_DRIVER" "$R_LIBRARY" "$FRESH_IDENTITY_OUT"
python3 scripts/profile_r_workflows.py --rscript "$RSCRIPT" \
  --baseline-library "$BASELINE_LIBRARY" --candidate-library "$R_LIBRARY" \
  --output-dir "$FRESH_PAIRED_OUT" --cases wald_1200x16 lrt_1200x16 \
  replacement_40x8 replacement_none_1200x16
```

The [source/artifact report](data/native-fit-route-local.json) records this
checkpoint. Remaining work includes the complete warning/option/supported-state
inventory, full fresh-fit cross-interface authority reconciliation, combined
GTEx/held-out precision and speed/RSS checks, sanitizer/MSRV/platform/package
gates and the final source bundle. Earlier checkpoint artifacts certify their
recorded source only. No release or publication occurs here.
