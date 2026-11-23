# Platform numerical parity

The standalone engine and R wrapper use the same Rust arithmetic. Independent
references reflect the original R environment on each affected platform; R is
never a production numerical backend. `r-host` selects the R arithmetic for
the executing host. Linux x86_64, Linux arm64 and macOS have specialized
profiles. Windows uses `r-binary64` with the system math library.
The portable profile uses the same Rust math library on every target.
See the [validation record](testing-parity.md#ci-status) for checked revisions.

## Linux formulas

Linux CI and release builds that run exact numerical checks use Ubuntu 26.04.
The Linux references use R 4.6.2 and one OpenBLAS 0.3.32 thread. Linux
x86_64 pins `OPENBLAS_CORETYPE=Haswell` and uses 80-bit x87 scalar reductions.
Linux arm64 pins `OPENBLAS_CORETYPE=ARMV8` and uses 128-bit scalar reductions
with a 113-bit significand. These settings keep the original R reference
independent of the machine's automatic kernel selection.

The Rust kernels follow the reference reduction order through the binary64
conversion. Polynomial QR follows the target dot product, norm and fused
operation schedules, including the fused scalar tail of the Ubuntu 26.04
Haswell axpy kernel. The Rust Haswell implementation also runs on x86_64
CPUs without AVX. Snapshots were regenerated with original base R for both
architectures. [Reference capture](https://github.com/deminden/rsdeseq2/actions/runs/45277310935)

Linux R-host golden manifests use glibc 2.43. The portable profile uses a
single manifest across platforms.
[Golden capture](https://github.com/deminden/rsdeseq2/actions/runs/45276866316)

Run the regeneration check on a matching Linux R installation. For arm64:

```sh
python3 scripts/validation/check_linux_formula_references.py \
  --arch aarch64 --rscript Rscript --output-dir results/checks/linux-arm64-formula-reference
```

Use `--arch x86_64` with the x86_64 R installation. The check regenerates
all seven formula panels with original R and compares every non-provenance
field with the snapshots, including input and output bits, attributes, column
order and rejection decisions. The elementary-function boundary fixture
uses its common name on all targets.

## Bounded posterior density

Original Linux x86_64 R rounds the far-tail normal-density square and addition
separately; the macOS reference uses a fused operation. A one-ulp log-density
difference can grow when posterior variance subtracts close moments. The shared
Rust helper uses the schedule selected by the active R profile. Tests check the
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
  cargo test --locked -p rsdeseq2 --test integration bounded_posterior::
```

Use an absolute reference path because Cargo executes tests from the crate.
CI saves independently generated inputs, outputs and environment provenance.
