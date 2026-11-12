# Frozen local VST migration: source audit

The 0.4.0 candidate still evaluates `dispersionFunction(object)` in R for a
frozen local VST. It also constructs normalized counts and the sinh grid in R.
This path does not satisfy the shared Rust numerical-engine requirement.

In the pinned DESeq2 1.52.0 environment, `localDispersionFit` creates a closure
whose environment holds a fitted `locfit` object. Its function evaluates that
fit at log means and exponentiates the prediction. DESeq2's local VST evaluates
that closure on a 999-point sinh grid, then integrates, interpolates and anchors
the transform. The Rust core already has local trend fitting and
`vst_local_from_dispersion_grid`, but reconstructing a fit from training data
is not proven equivalent to using the exact saved R fit.

The wrapper also accepts arbitrary closures tagged `fitType = "local"`. The
existing integration tests intentionally exercise a custom analytic closure,
including a normalization-factor matrix. A tag alone therefore cannot identify
a genuine DESeq2 locfit closure. Calling it to discover its values would keep
the R calculation and could execute arbitrary code.

The next implementation needs to:

1. Generate independent fixtures from actual DESeq2 local fits, including the
   all-minimum fallback, serialization/reload, both normalization forms, and
   out-of-range grid predictions. Record the original fit representation and
   grid/output values without invoking rsdeseq2.
2. Determine whether the saved locfit evaluation grid and coefficients can be
   validated and imported into a Rust prediction representation without refitting.
   Check interpolation, extrapolation and error behavior against those fixtures.
3. Move maximum normalized count and grid construction into shared Rust, then
   pass only validated, owned fit primitives across the R boundary. Test direct
   standalone/wrapped bit identity independently from DESeq2 precision.
4. Reject unsupported or unrecognized closures without executing them. Preserve
   support only for representations that can be translated safely and verified;
   document any compatibility change explicitly. The earlier user decision to
   reject custom *formula* functions does not by itself establish what to do
   with custom frozen dispersion functions.
5. Re-run the full source-package, GTEx, cross-interface and paired workflow gates
   before making performance or release claims. Keep the 0.4.0 evidence intact.

This audit establishes the concrete representation problem; it does not complete
the frozen-local VST port, the broader Rust-only migration, or the release audit.

## Implementation follow-up after 0.4.0

The development implementation now imports the standard one-dimensional
quadratic Gaussian locfit profile. `SavedLocalDispersionFit` owns the original
log-mean knots, fitted residual values/slopes and centered global polynomial.
It sorts the original knots, uses cubic Hermite interpolation and boundary
tangents, then adds the saved polynomial. It never refits the training data.
Invalid shapes, duplicate knots, unsupported fit profiles, altered prediction
transforms and masked prediction functions fail explicitly without execution.

The user selected restricted expression translation for custom dispersion
functions. The wrapper inspects a single-expression body and primitive captured
constants, then sends its AST to `FrozenLocalDispersion`. Supported built-ins,
including numeric `pmax`/`pmin`, use the shared standalone expression engine.
Masked/custom functions, active bindings and delayed expressions fail before
execution. Unresolved globals require a local captured constant or a literal;
this restriction avoids forcing a hidden promise during translation. The
wrapper constructs no numerical grid and never calls the dispersion closure.

`scripts/references/generate_frozen_local_vst_reference.R` uses original R 4.6.2,
DESeq2 1.52.0 and locfit 1.5.9.12, without loading rsdeseq2. The committed fixture
contains two actual saved fits (including uneven/tied training means),
serialization equality, both normalization
forms, all 999 grid/prediction points, and complete transform outputs. It also
records the original all-near-minimum fallback. Rust prediction/grid checks
use a relative/absolute `2e-12` bound; full VST retains the existing `2e-6` bound.
These are independent accuracy checks, separate from interface identity.

`scripts/validation/check_frozen_local_identity.R` compares eight complete cases from an
R-free standalone executable with actual production wrapper calls. Saved and
reloaded fits, restricted expressions and constant minimum curves cover both
normalization forms. Every grid, dispersion and transform output bit agrees,
while R `sinh`, `asinh`, `exp`, `sqrt`, `mean`, `sweep`, division and `predict` are blocked during
wrapper execution. The required R CI lane now runs this script. Adversarial
translation tests verify that unsupported bodies, active bindings and delayed
expressions remain unexecuted.

This closes the implemented numerical path for these frozen-local
representations. It does not close optional posterior/s-value ports, general
warning/method coverage, the broader Rust-only migration or its final audit.

## Performance observations

All measurements below are one paired pass per implementation, using R 4.6.2
in `deseq2-r462-native`, one Rust/BLAS worker and whole-process peak RSS. No
contrast or slower case was discarded. The 0.4.0 artifacts remain unchanged.

The first implementation serialized a full normalized matrix through the
primitive request. Its [target measurements](../data/frozen_local_vst_initial_2026-11-10.tsv)
showed a saved-fit VST stage increasing from 0.110 to 0.248 s with sample size
factors and from 0.245 to 0.389 s with matrix factors. That implementation was
replaced by a borrowed count/factor scan in Rust. Only the resulting 999-point
grid and small fit/expression primitives cross the request boundary.

The [borrowed-grid measurements](../data/frozen_local_vst_borrowed_2026-11-10.tsv)
cover 56,200 genes and 16 samples:

| Curve / normalization | Release VST stage | Borrowed-grid VST stage | Release RSS | Borrowed-grid RSS |
| --- | ---: | ---: | ---: | ---: |
| Saved fit / sample factors | 0.110 s | 0.099 s | 731.000 MiB | 717.422 MiB |
| Saved fit / matrix factors | 0.246 s | 0.249 s | 768.141 MiB | 768.125 MiB |
| Expression / sample factors | 0.101 s | 0.099 s | 730.641 MiB | 723.938 MiB |
| Expression / matrix factors | 0.246 s | 0.248 s | 771.609 MiB | 764.922 MiB |

Summed target process time is 9.374 s versus 9.360 s. All four paired inputs
are identical, and transformed values remain within the existing VST bound.
The source validation additionally requires exact standalone/wrapped output
bits. These observations do not estimate run-to-run variation.

The initial implementation's [twenty ordinary R workflows](../data/frozen_local_vst_r_workflows_initial_2026-11-10.tsv)
retain exact complete payloads: 86.475 s versus 86.957 s, with maximum RSS
1,150.109 versus 1,130.828 MiB, below the 2 GiB budget. Its
[native eleven-contrast pass](../data/frozen_local_vst_native_initial_2026-11-10.tsv)
measures 37.132 versus 36.026 s (**+3.07%**), predominantly the unblocked second
split and two blocked splits. Native maximum RSS is 301.969 versus 302.141 MiB;
original DESeq2 takes 190.584 s and 1,764.500 MiB, a 5.13× time ratio. The timing
cost is retained, and no 2% quantification optimization round is counted.
These initial-candidate measurements do not certify the final combined
performance gate. Later compact transport and validation changes require the
plan's final combined reconciliation after the outstanding ports/adjustment work.

## Reproduction

Use a native R 4.6.2 environment with the pinned reference dependencies:

```sh
Rscript scripts/references/generate_frozen_local_vst_reference.R /tmp/frozen-local-reference.json
cmp /tmp/frozen-local-reference.json crates/rsdeseq2/tests/data/frozen_local_vst_r462.json
cargo test --locked -p rsdeseq2 --test frozen_local_vst
cargo build --release --locked -p rsdeseq2 --example evaluate_r_interface
Rscript scripts/validation/check_frozen_local_identity.R \
  target/release/examples/evaluate_r_interface R_CHECK_LIBRARY /tmp/frozen-local-identity
```

`scripts/benchmarks/profile_frozen_local_vst.R` accepts a baseline or candidate installed
library, `saved|expression`, `size|matrix`, dimensions and two output paths.
For example, on macOS (use `/usr/bin/time -v` on Linux):

```sh
/usr/bin/time -l Rscript scripts/benchmarks/profile_frozen_local_vst.R \
  BASELINE_LIBRARY saved matrix 56200 16 /tmp/frozen-before.rds /tmp/frozen-before.tsv
/usr/bin/time -l Rscript scripts/benchmarks/profile_frozen_local_vst.R \
  CANDIDATE_LIBRARY saved matrix 56200 16 /tmp/frozen-after.rds /tmp/frozen-after.tsv
Rscript scripts/benchmarks/profile_frozen_local_vst.R compare /tmp/frozen-before.rds /tmp/frozen-after.rds
```

Alternate baseline/candidate order across the four combinations. Fixture
preparation uses independent R calculations in the benchmark process; the timed
stage is the production wrapper. RSS includes fixture preparation and the whole
process. The comparison checks complete input identity, dimensions/missingness
and a maximum absolute transform difference of `2e-6`; ordinary workflow
comparisons remain exact.
