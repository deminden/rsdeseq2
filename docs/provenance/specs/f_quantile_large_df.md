# Specification: F quantile for very large degrees of freedom

Status: written from black-box observation of R 4.6.2 only.

## Observation

`scripts/provenance/probe_qf_large_df.R` compares `qf(0.99, df1, df2)` with
the two chi-square limits of the F distribution (standard results, e.g.
Johnson, Kotz and Balakrishnan, Continuous Univariate Distributions, vol. 2,
chapter 27):

- df2 -> infinity: F(df1, df2) -> chi-square(df1) / df1;
- df1 -> infinity: F(df1, df2) -> df2 / chi-square(df2).

Output for macOS arm64: `docs/provenance/observations/qf_large_df_r462_macos.tsv`.

Observed facts:

1. For df1 in {1, 2, 5, 20, 200}, the smallest integer df2 at which
   `qf(0.99, df1, df2)` equals `qchisq(0.99, df1) / df1` bit for bit is
   400001.
2. For df2 in {1, 2, 5, 20, 200}, the smallest integer df1 at which
   `qf(0.99, df1, df2)` equals `df2 / qchisq(0.99, df2, lower.tail = FALSE)`
   bit for bit is 400001.
3. When both exceed 400000, the first limit applies when df1 <= df2
   (pairs (5e5, 6e5), (5e5, 5e5), (4e5, 400001)) and the second when
   df1 > df2 (pair (6e5, 5e5)).

## Required behaviour

`f_99_quantile(df1, df2)`:

- if df1 <= df2 and df2 > 400000: return the 0.99 lower-tail quantile of the
  chi-square distribution with df1 degrees of freedom, divided by df1;
- else if df1 > 400000: return df2 divided by the 0.99 upper-tail quantile
  (equivalently the 0.01 lower-tail quantile) of the chi-square distribution
  with df2 degrees of freedom;
- otherwise use the existing general path.

The chi-square quantile is computed in Rust from this crate's own chi-square
tail functions (for example by bisection on the upper tail), not by calling
R.

## Acceptance

`crates/rsdeseq2/tests/inference_primitives.rs`
(`standalone_probability_primitives_match_pinned_r462_references`) passes at
its existing tolerance, including the case m = 400002, p = 1. Add a unit test
with the probe's pairs against `qf` values recorded by a reference script.
