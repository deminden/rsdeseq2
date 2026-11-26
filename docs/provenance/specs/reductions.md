# Specification: BLAS-order reductions for Householder QR

Status: written from black-box observation only. No BLAS, LAPACK, LINPACK
or R source was consulted.

## Purpose

R's `qr()` reduces columns with the BLAS routines `ddot`, `dnrm2` and
`daxpy` of the BLAS that R is linked against. To reproduce R's `poly()` and
`contr.poly()` bases bit for bit, the Rust code must perform these three
operations in the same order and with the same rounding as that BLAS.

## Observation

`scripts/provenance/probe_reduction_order.R` compiles a small C shim against
the BLAS of the running R (`$(BLAS_LIBS)`, as an R package would) and calls
the three routines on chosen inputs. Its output for R 4.6.2 on macOS arm64
with the conda-forge OpenBLAS is in
`crates/rcompat-numeric/tests/data/reductions_r462_macos/`.

- `environment.tsv`: R version, platform and BLAS library.
- `ddot_tree.tsv`: the summation tree of `ddot(x, 1)` for n = 1..80, found
  with the method of Xie et al. (2024), "FPRev: Revealing the order of
  floating-point summation" (x[i] = 2^60, x[j] = -2^60, other terms 1; the
  result is n minus the size of the subtree joining i and j).
- `ddot_fma.tsv`, `ddot_fma_pairs.tsv`: whether a product is rounded before
  it is added (a = 1 + 2^-30, b = 1 - 2^-30, so a*b = 1 - 2^-60 rounds to 1;
  a fused update leaves -2^-60 after a -1 in the same accumulator).
- `daxpy_fma.tsv`: the same test for `y + alpha*x`.
- `dnrm2.tsv`: norms of chosen vectors, including 1e-200 and 1e200 entries.
- `qr_norm_crosscheck.tsv`: `qr(cbind(x))$qr[1,1]` equals
  `-sign(x[1]) * dnrm2(x)` bit for bit, so `qr()` uses this `dnrm2`.
- `ddot_random.tsv`, `dnrm2_random.tsv`: 300 random cases each (lengths
  1..100, 200, 257 and, for `ddot`, 1000) with the BLAS result, as acceptance
  data. Values are hexadecimal big-endian binary64 bit patterns separated by
  commas.

## Required behaviour

All arithmetic is IEEE binary64 with round-to-nearest-even. "fma(a, b, c)"
is a fused multiply-add (one rounding), `f64::mul_add` in Rust.

### dot(x, y), equal lengths n

Observed in `ddot_tree.tsv` and the FMA files:

1. Let m = 32 * floor(n / 32) (the largest multiple of 32 not above n).
2. Sixteen accumulators a_1..a_16 start at 0. For each position p = 1..m
   in increasing order, accumulator k = ((p - 1) mod 16) + 1 is updated as
   a_k = fma(x_p, y_p, a_k). So a_1 receives terms 1, 17, 33, 49, ...;
   a_16 receives terms 16, 32, 48, ...
3. The accumulators are combined by plain additions in this exact tree:
   odd = ((a_1 + a_3) + (a_5 + a_7)) + ((a_9 + a_11) + (a_13 + a_15));
   even = ((a_2 + a_4) + (a_6 + a_8)) + ((a_10 + a_12) + (a_14 + a_16));
   s = odd + even.
   When m = 0 the accumulators are all 0 and s = 0.
4. The remaining positions p = m+1..n are added in increasing order with
   s = fma(x_p, y_p, s).
5. The result is s.

For n < 32 this is a left-to-right fused accumulation from 0.

### axpy(alpha, x, y)

Observed in `daxpy_fma.tsv`: every element becomes y_i = fma(alpha, x_i, y_i).

### norm(x)

Observed in `dnrm2.tsv` and `dnrm2_random.tsv`:

- Vectors of 1e-200 and of 1e200 return the correct norm without underflow
  or overflow, so the routine scales.
- The result depends on the order of the entries (`1:10` and `10:1` differ
  in the last bit), so it is a sequential scan.
- Zero entries do not change the result (`c(0, 0, 5, 0)` gives exactly 5).

The published one-pass scaled sum-of-squares method fits these facts (Blue,
J. L., 1978, "A portable Fortran program to find the Euclidean norm of a
vector", ACM TOMS 4(1):15-23; Higham, N. J., 2002, Accuracy and Stability of
Numerical Algorithms, 2nd ed., SIAM, discussion of computing the 2-norm
without overflow): keep `scale`, the largest magnitude seen so far, and
`ssq`, the sum of squares in units of `scale`, starting from scale = 0 and
ssq = 1; for each nonzero entry v, if |v| > scale set
ssq = 1 + ssq * (scale/|v|)^2 and scale = |v|, otherwise add (|v|/scale)^2
to ssq; return scale * sqrt(ssq). The facts above do not fix whether each
update rounds once (fused multiply-add) or more than once, nor the exact
grouping of the update. The implementation must choose the variant that
reproduces every row of `dnrm2.tsv` and `dnrm2_random.tsv` bit for bit, and
document the choice. If no variant reproduces all rows, report the
mismatching rows rather than consulting any source.

## Acceptance

- `dot` reproduces every row of `ddot_random.tsv` bit for bit.
- `norm` reproduces every row of `dnrm2.tsv` and `dnrm2_random.tsv`.
- `axpy` reproduces `daxpy_fma.tsv`.
- These are unit tests that read the observation files.

## Scope

These kernels describe one host BLAS. They are used for every numeric
profile, so that `poly()` and `contr.poly()` give the same basis everywhere;
the exact formula references were recorded on this host.
