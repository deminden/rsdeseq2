# Platform numerical parity

The standalone engine and R wrapper use the same Rust arithmetic. Independent
references reflect the original R environment on each affected platform; R is
never a production numerical backend. [Current CI](testing-parity.md#ci-status)
passes Linux, both macOS architectures and Windows at unchanged numerical bounds.

## Linux formulas

Linux x86_64 uses R 4.6.2 and OpenBLAS 0.3.26 with one thread and
`OPENBLAS_CORETYPE=Haswell`. Pinning the kernel prevents hosted CPU changes from
silently selecting a different original-R oracle. The Rust scalar implementation
of the reduction schedule requires no AVX CPU.

R's Linux means, finite per-argument/column sums and BLAS norms use extended
accumulation. Rust implements these reductions directly and retains the extended
result until the corresponding binary64 conversion. Polynomial QR uses the
reference kernel's reduction order and fused operations. General fitting/rank
policies remain separate from this formula arithmetic.

On the matching Linux x86_64 reference environment:

```sh
python3 scripts/validation/check_linux_formula_references.py \
  --rscript Rscript --output-dir results/checks/linux-formula-reference
```

This independently regenerates all seven formula panels with original R and
compares every non-provenance field with the frozen Linux snapshots: inputs,
output bits, attributes, order and rejection decisions. The complete fitting
panel checks R term/matrix order without column remapping. macOS snapshots stay
unchanged. Metadata displays library basenames rather than installation paths.
The elementary-function boundary fixture is portable and uses its common name;
required platform-specific snapshots do not silently fall back when absent.

## Bounded posterior density

Original Linux x86_64 R rounds the far-tail normal-density square and addition
separately; the macOS reference uses a fused operation. A one-ulp log-density
difference can grow when posterior variance subtracts close moments. The shared
Rust helper uses the target's original arithmetic schedule. Tests check the
independent density probe's exact bits and existing posterior bounds; no broader
CDF or SD tolerance substitutes for the correct expression.

## Windows and Intel macOS CDFs

Ultra-narrow uniform intervals can differ across original R arithmetic schedules.
For one nearly coincident normal-probability case, original Windows R returns
component CDF `0.5`, while original ARM R returns `0.50113895216400917` with the
same mixture weights. A single platform's saved expectation cannot represent
both results at the existing `2e-8` bound.

`generate_uniform_cdf_reference.R` independently evaluates the model using only
base R. It never loads rsdeseq2 or selects expected values from candidate output.
Windows and Intel macOS CI regenerate all 6,500 CDFs under R 4.6.2. Rust validates
every input/value and coverage at the unchanged `2e-8` bound; all other posterior
fields retain their original fixtures. On ARM the independent values match the
frozen original-ashr CDF panel within `1.12e-16`.

```sh
Rscript --vanilla scripts/references/generate_uniform_cdf_reference.R \
  "$PWD/results/uniform-cdf.tsv"
RSDESEQ2_UNIFORM_CDF_REFERENCE="$PWD/results/uniform-cdf.tsv" \
  cargo test --locked -p rsdeseq2 --test bounded_posterior
```

Use an absolute reference path because Cargo executes tests from the crate.
CI saves independently generated inputs, outputs and environment provenance.
