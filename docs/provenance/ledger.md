# Provenance ledger

This ledger records how each part of the workspace crates and the R bridge
came to be. It is checked by `scripts/validation/check_provenance.py`: every
tracked source file must match a row, and every provenance statement in a
source file (phrases such as "adapted from", names of upstream projects,
licence file names) must fall under a row whose status is resolved.

Rows are matched by path pattern (`*` within a path segment, `**` across
segments). When several rows match a file, the row with the longest pattern
wins.

Knowledge sources: `paper`, `manual` (help pages and reference manuals, for
behaviour only), `black-box` (outputs of running R or a host library on chosen
inputs), `public-domain code`, `MIT code`, `upstream source`, `own design`,
`unknown`.

Statuses: `clean`, `notice kept`, `to replace`, `to remove`,
`owner to confirm`. Only `clean` and `notice kept` are resolved.

"Owner statement" means the owner states that the module was written from
published papers, reference manuals and help pages, and black-box runs of R,
without upstream source code open.

## Rows

| Path | Computes | Knowledge source | Evidence | Status |
|---|---|---|---|---|
| `crates/**` | Workspace crates not listed more specifically below: statistical core, formula and model-matrix code, shrinkage engines, I/O and CLI | paper, manual, black-box | Owner statement (2026-11-25); R reference fixtures under `crates/rsdeseq2/tests/data` and `scripts/references` | clean |
| `r-pkg/rsdeseq2/**` | R bridge, R wrapper functions and their tests | manual, black-box | Owner statement (2026-11-25) | clean |
| `crates/rcompat-numeric/src/math/rank.rs` | `matrix_rank_r_qr`: rank by Householder QR with R's column-deferral rule | upstream source | Doc comment states it follows R's `src/appl/dqrdc2.f` (GPL) | to replace |
| `crates/rcompat-formula/src/design/implementation/basics.rs` | Calls the R-style rank rule; doc comment names LINPACK `dqrdc2` | paper, manual, black-box | Owner statement; the rank rule itself is the `rank.rs` row | clean |
| `crates/rcompat-formula/src/design/numeric_transform/extended.rs` | x87 extended reductions and norm for the `RLinuxX86` profile | upstream source | OpenBLAS 0.3.26 Haswell schedules; BSD notice removed in 0.6.2 while the code was reshaped | to remove |
| `crates/rcompat-formula/src/design/numeric_transform/extended_unavailable.rs` | Binary64 stand-in for the x87 kernels | own design | Exists only to serve the `extended.rs` row | to remove |
| `crates/rcompat-formula/src/design/numeric_transform/binary128.rs` | Binary128 reductions for the `RLinuxArm64` profile | upstream source | OpenBLAS ARMV8 schedules; BSD notice removed in 0.6.2 | to remove |
| `crates/rcompat-formula/src/design/numeric_transform/binary128/wide.rs` | Binary128 round-to-nearest arithmetic replacing `rustc_apfloat` | unknown | Owner to confirm it was written from the IEEE 754 definition; removed with the `RLinuxArm64` profile | to remove |
| `crates/rcompat-formula/src/design/numeric_transform/reduction.rs` | Dot product, norm and axpy order used by polynomial QR on every profile | upstream source | 0.6.1 header: "Arithmetic schedules adapted from OpenBLAS 0.3.26"; reshaped in 0.6.2 | to replace |
| `crates/rcompat-formula/src/design/numeric_transform/polynomial_qr.rs` | Householder QR with column deferral for `poly()` and `contr.poly()` | unknown | Mirrors the rank rule of the `rank.rs` row and uses `reduction.rs` | to replace |
| `crates/rcompat-numeric/src/numeric/profile.rs` | Numeric profile selection; describes the OpenBLAS Haswell and ARMV8 schedules | own design | The two Linux R profiles are being removed | to replace |
| `crates/rcompat-numeric/src/probability/f_quantile.rs` | 0.99 F quantile for the Cook's cutoff, by inverting the upper tail of the beta-logistic density; F from the beta relation | paper, black-box | The large degrees-of-freedom branch that cited R's `qf.c` is removed; `f_quantile_r462.tsv` from `qf()` checks df1 up to 200 and df2 up to 100,000 | clean |
| `crates/rcompat-numeric/src/numeric/special.rs` | Gamma, digamma, incomplete gamma, beta, normal and gamma quantiles | MIT code | Port of `statrs` 0.19.1 (MIT); `erf_inv` uses Boost rational approximations by way of statrs | to replace |
| `crates/rcompat-shrinkage/src/shrinkage/heavy_tail/unbounded_lbfgs.rs` | Unbounded L-BFGS with strong-Wolfe line search | MIT code | Adaptation of LBFGSpp (MIT); notice in `LBFGSPP-LICENSE` | to replace |
| `crates/rcompat-shrinkage/src/shrinkage/heavy_tail/LBFGSPP-LICENSE` | Licence notice for the LBFGSpp adaptation | MIT code | Covers `unbounded_lbfgs.rs` | to remove |
| `crates/rcompat-shrinkage/src/shrinkage/mixture_fit/sqp/compatibility.rs` | Active-set and floored EM policy of the mixture SQP solver | MIT code | Implements mixsqp 0.3-54 (MIT); notice in `MIXSQP-LICENSE` | to replace |
| `crates/rcompat-shrinkage/src/shrinkage/mixture_fit/sqp/MIXSQP-LICENSE` | Licence notice for the mixsqp algorithm | MIT code | Covers `sqp/compatibility.rs` | to remove |
| `crates/rcompat-shrinkage/src/shrinkage/gamma_normalization.rs` | Gamma function series for the R profiles | public-domain code | Fullerton's SLATEC `dgamma.f` tables, cited; SLATEC guide section 4 | clean |
| `crates/rcompat-numeric/src/probability/normal_log.rs` | Log normal tail probabilities | paper | Cody (1969), Math. Comp. | clean |
| `crates/rcompat-numeric/src/probability/digamma.rs` | Digamma for positive arguments | paper | DLMF 5.5 and 5.11 | clean |
| `crates/rcompat-shrinkage/src/shrinkage/heavy_tail/random.rs` | MT19937 and CMRG random starts | paper | Matsumoto and Nishimura (1998); L'Ecuyer (1999) | clean |
| `crates/rcompat-shrinkage/src/shrinkage/mixture_fit/squarem.rs` | SQUAREM extrapolation | paper | Varadhan and Roland (2008) | clean |
| `crates/rcompat-shrinkage/src/shrinkage/log_f_quadrature.rs` | 15-point Gauss-Kronrod integration on infinite intervals | public-domain code | QUADPACK `dqk15i.f` (netlib) nodes and weights, cited | clean |
| `crates/rcompat-shrinkage/src/shrinkage/location_fit.rs` | Location preprocessing; reproduces the mixture solver's warning text | black-box | Warning text as emitted by R; owner statement | clean |
| `crates/rcompat-numeric/src/lib.rs` | Crate documentation; mentions the statrs ports | own design | Follows the `special.rs` row | to replace |
| `r-pkg/rsdeseq2/tests/testthat/test-size-factors.R` | Size-factor tests; one block is a port of an upstream DESeq2 test file | upstream source | Comment names `external/DESeq2/tests/testthat/test_size_factor.R` (LGPL) | to replace |
| `THIRD_PARTY_NOTICES.md` | Third-party notices of the workspace | MIT code | statrs, LBFGSpp and mixsqp sections | to replace |
| `crates/*/THIRD_PARTY_NOTICES.md` | Per-crate third-party notices | MIT code | Copies of the workspace sections | to replace |
| `r-pkg/rsdeseq2/inst/COPYRIGHTS` | Third-party notices of the R package | MIT code | Copies of the workspace sections | to replace |
| `rcompat-locfit 0.2.0` (crates.io) | Local likelihood dispersion trend | unknown | Separate repository, declared MIT; upstream locfit is GPL | owner to confirm |
| `rcompat-lbfgsb 0.3.0` (crates.io) | Bounded L-BFGS-B for `optim` | unknown | Separate repository, declared MIT | owner to confirm |
