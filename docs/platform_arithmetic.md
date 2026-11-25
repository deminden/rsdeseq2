# Platform arithmetic

The standalone engine and R wrapper use the same Rust arithmetic. Independent
references reflect the original R environment they were recorded in; R is
never a production numerical backend. `r-host` selects the R arithmetic for
the executing host: `r-macos` on macOS and `r-binary64` everywhere else,
which uses binary64 reductions and the system maths library. The portable
profile uses the same Rust maths library on every target.
See the [validation record](testing-parity.md#ci-status) for checked revisions.

## Formula references

Exact formula references (`formula_expression`, `model_matrix`,
`native_formula`, `native_factor`, `numeric_transform`, `formula_fitting` and
the uniform posterior panel) were recorded with R 4.6.2 on macOS arm64. Tests
compare them under the `r-macos` and `portable` profiles. Other hosts have no
recorded reference, and those tests return early there.

Release 0.7.0 removed the Linux x86_64 and Linux arm64 R profiles. Their
reduction kernels reproduced summation schedules that had been taken from
OpenBLAS source, so they could not stay in an MIT crate. Linux now uses
`r-binary64`: results agree with R at the tolerances of the reference tests,
but Linux R's last bits are no longer reproduced.

## Bounded posterior density

The far-tail normal log-density uses one fused multiply-add on every profile.
A one-ulp log-density difference can grow when posterior variance subtracts
close moments, so the rounding is fixed rather than left to the compiler.

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
