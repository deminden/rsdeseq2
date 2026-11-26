# Specification: Householder QR with limited column pivoting

Status: written from allowed inputs only: R's help page for `qr` (behaviour),
the LINPACK Users' Guide (Dongarra, Bunch, Moler and Stewart, 1979, SIAM,
chapter 9, QR decomposition, describing the public-domain routines DQRDC and
DQRSL), and black-box runs of R 4.6.2. No R, LINPACK, LAPACK or BLAS source
was consulted.

## Purpose

Three callers need R's default `qr()` (not the LAPACK variant):

1. `matrix_rank_r_qr` in `rcompat-numeric/src/math/rank.rs`: the rank that
   `qr(x)$rank` reports, which DESeq2's full-rank check compares with the
   number of columns.
2. The orthogonal polynomial basis of `poly()`, and
3. the polynomial contrasts of `contr.poly()`, both in `rcompat-formula`.
   Both need, for each column j of the pivoted matrix, the vector
   `Q (d_j e_j)`, where `d_j` is the j-th diagonal entry of the compact QR
   matrix: the columns of `qr.qy(QR, diag(diag(QR$qr), n, p))`.

## Behaviour from the help page

R's `?qr`: the default method uses a modification of LINPACK's DQRDC that
uses "a limited column pivoting strategy": columns with near-zero 2-norm are
moved to the right-hand edge of `x`, and `tol` (default 1e-7) is the
tolerance for detecting linear dependence. `rank` is the computed rank and
`pivot` the column order. The owner's statement of the rule: a column whose
remaining norm (the part not explained by the columns before it) falls below
`tol` times its original norm is moved to the end and does not count towards
the rank.

## Factorization (LINPACK Users' Guide, chapter 9)

Matrix x is n by p, stored by columns; positions are 1-based here. All
arithmetic is binary64; `dot`, `norm` and `axpy` are the reductions of
`docs/provenance/specs/reductions.md`, applied to contiguous column segments.

Before the loop, record `orig_j = norm(x(1:n, j))` for every column.
`active` starts at p. For l = 1, 2, ..., min(n, p):

1. Deferral (help page rule): while l <= active and the remaining norm
   `rn_l` of column l is below `tol * orig_l`, move column l (with its
   `orig` and `rn`) to the end, shifting columns l+1..p one place left, and
   decrease `active`. Columns with orig = 0 are moved (observed).
2. If l = n, stop: the last row has no reflector.
3. `nrmxl = norm(x(l:n, l))`. If it is 0, set qraux_l = 0 and go to the next l.
4. If x(l,l) != 0, give nrmxl the sign of x(l,l).
5. Scale x(l:n, l) by multiplying each entry by `1 / nrmxl` (the reciprocal
   is computed once).
6. x(l,l) = 1 + x(l,l).
7. For each later column j = l+1..p (including moved columns):
   t = -dot(x(l:n, l), x(l:n, j)) / x(l,l);
   x(l:n, j) = axpy(t, x(l:n, l), x(l:n, j)).
8. qraux_l = x(l,l); x(l,l) = -nrmxl.

### Remaining norms

`rn_j` starts as `orig_j`. After step l has transformed column j (j > l),
the part of column j below row l loses the entry x(l,j), so its norm is
downdated as in QR with column pivoting (Golub and Van Loan, Matrix
Computations, 4th ed., section 5.4.2, norm downdating):
`rn_j <- rn_j * sqrt(max(0, 1 - (|x(l,j)| / rn_j)^2))`, skipped when
`rn_j = 0`. Because the update loses accuracy under heavy cancellation, it
is replaced by a direct recomputation `rn_j = norm(x(l+1:n, j))` when the
update becomes unreliable.

The observations fix the details that this description leaves open: the
condition that triggers recomputation (candidates include a threshold on the
factor `1 - ratio^2` and a threshold on `rn_j` relative to `orig_j`), and
the rounding of the update (fused or separate operations). R stores `rn_j`
in `qraux_j` for positions without a reflector (l >= n and the extra columns
of wide matrices), so those qraux values expose the updated norms directly.
Choose the variant that reproduces rank, pivot and qraux of every case in
both observation files, and document it. If none does, report the cases.


The rank is the number of columns not moved, but never more than n. Rank
and pivot of the cases in the observation file decide the details (when the
loop stops for wide matrices, zero columns, one-row matrices).

## Applying Q (LINPACK DQRSL, "qy")

Given the compact matrix, qraux, k = rank, and a vector y of length n: let
ju = min(k, n - 1). For l = ju, ju-1, ..., 1: if qraux_l != 0, temporarily
use qraux_l as x(l,l), compute t = -dot(x(l:n, l), y(l:n)) / x(l,l) and
y(l:n) = axpy(t, x(l:n, l), y(l:n)), then restore x(l,l).

The basis vector for column j (1-based, pivoted order, j = 1..p) is this
product applied to y = d_j e_j with d_j = x(j,j) of the compact matrix
(y = 0 when j > n).

## Observation and acceptance

`--downdate` mode of the same script writes
`crates/rcompat-numeric/tests/data/qr_downdate_cases_r462_macos.json`: 170 larger
cases (`contr.poly` inputs for 11 to 120 levels, centred polynomial bases up
to degree 20, near-dependent random matrices) with rank, pivot and qraux.
The `contr_poly_n` inputs are not stored; they are
`outer(seq_len(n) - (n + 1) / 2, 0:(n - 1), "^")` evaluated in R, so compute
them with R or check that Rust `powi`/`powf` reproduces those bits.

`scripts/provenance/probe_qr_rules.R` writes
`crates/rcompat-numeric/tests/data/qr_cases_r462_macos.json`: for 65 matrices
(identity, zero columns, all-zero, one row, one column, wide, dependent,
near-collinear at 1e-5 to 1e-11 in two scalings, other tolerances, mixed
magnitudes, centred Vandermonde matrices like `poly()` builds, `contr.poly()`
inputs and random matrices) it records the input, `rank`, `pivot`, the
compact `qr`, `qraux` and the `qr.qy` basis above, as big-endian binary64
hex (matrices column-major; scalars wrapped in one-element arrays).

Acceptance:

- `rank`, `pivot` and `qraux` match every case of both files exactly.
- `qr` and the basis match bit for bit for every case of the first file.
  Report any case that does not; do not consult any source.
- The existing tests pass unchanged: rank fixtures
  (`crates/rsdeseq2/tests/data/rank_r462.tsv`, `rank_reference_r462.json`,
  integration modules `design_rank::` and `model_validity::`), the exact
  formula references (`numeric_transform`, `model_matrix`,
  `formula_expression`, `native_formula`, `formula_fitting`) and the golden
  manifests (`cargo golden --check` with the `r-macos` and `portable`
  profiles).
- Property tests: rank <= min(n, p); rank unchanged when a column is scaled
  by a power of two; full rank for orthogonal designs.
