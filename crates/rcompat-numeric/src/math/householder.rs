//! Householder QR with limited column pivoting, as R's default `qr()`.
//!
//! The factorization follows the Householder method of the LINPACK Users'
//! Guide (chapter 9) with the deferral rule of R's `?qr` help page: a column
//! whose remaining norm falls below `tolerance` times its original norm is
//! moved to the right-hand edge and does not count towards the rank. All
//! column reductions use [`super::reduction`], so results are reproducible
//! bit for bit. The behaviour is specified in
//! `docs/provenance/specs/householder_qr.md`.

use super::reduction::{axpy, dot, norm};
use crate::errors::{Error, ErrorKind, Requirement, invalid_parameter};
use crate::matrix::try_alloc_matrix;

/// Compact Householder QR factorization of an `n` by `p` matrix.
///
/// Positions and pivots are 0-based. Column `j` of the compact matrix holds
/// `R` on and above the diagonal and the scaled Householder vector below it;
/// `qraux[l]` holds the leading entry of reflector `l`.
#[doc(hidden)]
#[derive(Debug, Clone)]
pub struct HouseholderQr {
    rows: usize,
    columns: Vec<Vec<f64>>,
    qraux: Vec<f64>,
    pivot: Vec<usize>,
    rank: usize,
}

impl HouseholderQr {
    /// Factorize the matrix whose columns are `columns`, each of length `rows`.
    ///
    /// For each position `l` (up to `min(rows, p)`), columns whose remaining
    /// norm `‖x[l.., l]‖` is below `tolerance` times their original norm, and
    /// columns that are entirely zero, are moved to the end while columns
    /// that have not been moved remain at `l`. A reflector is then formed from
    /// `x[l.., l]` and applied to every later column (moved ones included).
    /// The last row has no reflector. The rank is the number of columns not
    /// moved, at most `rows`.
    ///
    /// The remaining norm of each column starts as its original norm and is
    /// downdated after every reflector, with a direct recomputation when the
    /// downdate cancels heavily (described in the module source); positions
    /// without a reflector (`l >= rows - 1`) record it in `qraux`.
    pub fn factorize(columns: Vec<Vec<f64>>, rows: usize, tolerance: f64) -> Result<Self, Error> {
        if !tolerance.is_finite() || tolerance < 0.0 {
            return Err(invalid_parameter(
                ErrorKind::Options,
                "QR tolerance",
                Requirement::FiniteAndNonNegative,
            ));
        }
        if let Some(bad) = columns.iter().find(|column| column.len() != rows) {
            return Err(Error::InvalidDimensions {
                context: "QR column length".to_string(),
                expected: rows,
                actual: bad.len(),
            });
        }
        let mut columns = columns;
        let p = columns.len();
        let mut original = try_alloc_matrix(p, 1, "QR column norms")?;
        original.extend(columns.iter().map(|column| norm(column)));
        let mut remaining = try_alloc_matrix(p, 1, "QR remaining norms")?;
        remaining.extend_from_slice(&original);
        let mut pivot = try_alloc_matrix(p, 1, "QR pivot")?;
        pivot.extend(0..p);
        let mut qraux = try_alloc_matrix(p, 1, "QR auxiliary values")?;
        qraux.resize(p, 0.0);
        let mut active = p;

        for l in 0..rows.min(p) {
            while l < active {
                let moved = original[l] == 0.0 || remaining[l] < tolerance * original[l];
                if !moved {
                    break;
                }
                columns[l..].rotate_left(1);
                original[l..].rotate_left(1);
                remaining[l..].rotate_left(1);
                pivot[l..].rotate_left(1);
                active -= 1;
            }
            if l + 1 == rows {
                break;
            }
            let (done, rest) = columns.split_at_mut(l + 1);
            let reflector = &mut done[l][l..];
            let mut length = norm(reflector);
            if length == 0.0 {
                qraux[l] = 0.0;
                continue;
            }
            if reflector[0] != 0.0 {
                length = length.copysign(reflector[0]);
            }
            let inverse = 1.0 / length;
            for value in reflector.iter_mut() {
                *value *= inverse;
            }
            reflector[0] += 1.0;
            let lead = reflector[0];
            for (column, norm_left) in rest.iter_mut().zip(&mut remaining[l + 1..]) {
                let target = &mut column[l..];
                let t = -dot(reflector, target) / lead;
                axpy(t, reflector, target);
                *norm_left = downdate(*norm_left, target);
            }
            qraux[l] = lead;
            reflector[0] = -length;
        }
        if rows > 0 {
            let first = (rows - 1).min(p);
            qraux[first..].copy_from_slice(&remaining[first..]);
        }
        Ok(Self {
            rows,
            columns,
            qraux,
            pivot,
            rank: active.min(rows),
        })
    }

    /// Number of rows `n`.
    pub fn rows(&self) -> usize {
        self.rows
    }

    /// Number of columns `p`.
    pub fn cols(&self) -> usize {
        self.columns.len()
    }

    /// Numerical rank: the number of columns not moved, at most `n`.
    pub fn rank(&self) -> usize {
        self.rank
    }

    /// Original (0-based) index of the column at each pivoted position.
    pub fn pivot(&self) -> &[usize] {
        &self.pivot
    }

    /// Leading entries of the reflectors (see [`HouseholderQr::factorize`]).
    pub fn qraux(&self) -> &[f64] {
        &self.qraux
    }

    /// Column `j` of the compact factorization, in pivoted order.
    pub fn column(&self, j: usize) -> &[f64] {
        &self.columns[j]
    }

    /// Diagonal entry `j` of the compact matrix, or 0 when `j >= n`.
    pub fn diagonal(&self, j: usize) -> f64 {
        if j < self.rows {
            self.columns[j][j]
        } else {
            0.0
        }
    }

    /// Overwrite `y` (length `n`) with `Q y`, using the first `rank` reflectors.
    ///
    /// Reflectors are applied from the last used one back to the first; a
    /// reflector whose `qraux` entry is zero is skipped.
    pub fn apply_q(&mut self, y: &mut [f64]) -> Result<(), Error> {
        if y.len() != self.rows {
            return Err(Error::InvalidDimensions {
                context: "QR product vector".to_string(),
                expected: self.rows,
                actual: y.len(),
            });
        }
        let used = self.rank.min(self.rows.saturating_sub(1));
        for l in (0..used).rev() {
            let lead = self.qraux[l];
            if lead == 0.0 {
                continue;
            }
            let reflector = &mut self.columns[l][l..];
            let saved = reflector[0];
            reflector[0] = lead;
            let target = &mut y[l..];
            let t = -dot(reflector, target) / lead;
            axpy(t, reflector, target);
            reflector[0] = saved;
        }
        Ok(())
    }

    /// Visit `Q (d_j e_j)` for every pivoted column `j`, in order.
    ///
    /// `d_j` is the `j`-th diagonal entry of the compact matrix (zero beyond
    /// the last row). These are the columns of R's
    /// `qr.qy(QR, diag(diag(QR$qr), n, p))`, the basis behind `poly()` and
    /// `contr.poly()`.
    pub fn visit_basis(&mut self, mut visit: impl FnMut(usize, &[f64])) -> Result<(), Error> {
        let mut y = try_alloc_matrix(self.rows, 1, "QR basis vector")?;
        y.resize(self.rows, 0.0);
        for j in 0..self.cols() {
            y.fill(0.0);
            if j < self.rows {
                y[j] = self.diagonal(j);
            }
            self.apply_q(&mut y)?;
            visit(j, &y);
        }
        Ok(())
    }
}

/// Factor `1 - (|x_l| / rn)^2` below which the downdate is replaced by a
/// direct recomputation.
///
/// Fitted to the observations: R recomputed at a factor of
/// 9.999999997177891e-7 and downdated at 1.0000000001618781e-6 (the two
/// attainable factors either side of 1e-6 for a two-column probe), and the
/// switch does not depend on how far the remaining norm has already fallen
/// relative to the original norm (probed at ratios 1, 0.71 and 0.01).
const RECOMPUTE_BELOW: f64 = 1e-6;

/// Remaining norm of a column after one reflector has been applied.
///
/// `target` is the transformed column from the reflector's row down; its
/// first entry leaves the remaining part. Following the norm downdating of
/// QR with column pivoting (Golub and Van Loan, Matrix Computations, 4th ed.,
/// section 5.4.2), `rn <- rn * sqrt(max(0, 1 - (|x_l| / rn)^2))`, with the
/// ratio rounded, `1 - ratio^2` evaluated as one fused multiply-add, the
/// square root and the product rounded. When that factor is below
/// [`RECOMPUTE_BELOW`] the norm is recomputed from the entries below. A zero
/// norm stays zero.
///
/// This variant was selected because it reproduces rank, pivot and every
/// `qraux` value of both observation files
/// (`crates/rcompat-numeric/tests/data/qr_cases_r462_macos.json` and
/// `qr_downdate_cases_r462_macos.json`). Separately rounded `1 - r * r` and
/// the factored `(1 - r)(1 + r)` change pivots; triggers on the norm relative
/// to the original (or last recomputed) norm leave stored norms that differ.
fn downdate(remaining: f64, target: &[f64]) -> f64 {
    if remaining == 0.0 {
        return 0.0;
    }
    let ratio = target[0].abs() / remaining;
    let factor = (-ratio).mul_add(ratio, 1.0).max(0.0);
    if factor < RECOMPUTE_BELOW {
        norm(&target[1..])
    } else {
        remaining * factor.sqrt()
    }
}

#[cfg(test)]
mod tests {
    use super::HouseholderQr;
    use crate::math::rank::matrix_rank_r_qr;

    const CASES: &str = include_str!("../../tests/data/qr_cases_r462_macos.json");
    const DOWNDATE_CASES: &str = include_str!("../../tests/data/qr_downdate_cases_r462_macos.json");

    /// One recorded case; matrices are column-major. The larger file stores
    /// neither `qr` nor the basis (left empty here), nor the `contr.poly()`
    /// inputs, which are rebuilt by [`contr_poly_input`].
    struct Case {
        name: String,
        rows: usize,
        cols: usize,
        tolerance: f64,
        x: Vec<f64>,
        rank: usize,
        pivot: Vec<usize>,
        qr: Vec<f64>,
        qraux: Vec<f64>,
        qy: Vec<f64>,
    }

    // The observation file holds flat arrays of strings or integers only, so
    // a field is the text between `"key":[` and the next `]`.
    // A missing key reads as an empty array.
    fn field<'a>(object: &'a str, key: &str) -> Vec<&'a str> {
        let marker = format!("\"{key}\":[");
        let Some(found) = object.find(&marker) else {
            return Vec::new();
        };
        let start = found + marker.len();
        let end = start + object[start..].find(']').unwrap();
        object[start..end]
            .split(',')
            .filter(|item| !item.is_empty())
            .map(|item| item.trim_matches('"'))
            .collect()
    }

    /// `outer(seq_len(n) - (n + 1) / 2, 0:(n - 1), "^")`, column-major. The
    /// powers are of integers or half-integers. R computed them with the macOS
    /// C library's `pow`, which is not always correctly rounded, so `powf`
    /// gives R's inputs only on macOS (checked against R for n = 11..120);
    /// other hosts skip these cases.
    fn contr_poly_input(n: usize) -> Vec<f64> {
        let center = (n as f64 + 1.0) / 2.0;
        (0..n)
            .flat_map(|degree| {
                (0..n).map(move |row| (row as f64 + 1.0 - center).powf(degree as f64))
            })
            .collect()
    }

    fn doubles(object: &str, key: &str) -> Vec<f64> {
        field(object, key)
            .into_iter()
            .map(|hex| f64::from_bits(u64::from_str_radix(hex, 16).unwrap()))
            .collect()
    }

    fn integers(object: &str, key: &str) -> Vec<usize> {
        field(object, key)
            .into_iter()
            .map(|text| text.parse().unwrap())
            .collect()
    }

    fn cases(text: &str) -> Vec<Case> {
        let body = &text[text.find("\"cases\":[").unwrap()..];
        body.split("{\"name\":")
            .skip(1)
            .map(|chunk| {
                let object = format!("{{\"name\":{chunk}");
                let name = field(&object, "name")[0].to_string();
                let rows = integers(&object, "rows")[0];
                let mut x = doubles(&object, "x");
                if x.is_empty() && name.starts_with("contr_poly_") {
                    x = contr_poly_input(rows);
                }
                Case {
                    name,
                    rows,
                    cols: integers(&object, "cols")[0],
                    tolerance: doubles(&object, "tolerance")[0],
                    x,
                    rank: integers(&object, "rank")[0],
                    pivot: integers(&object, "pivot"),
                    qr: doubles(&object, "qr"),
                    qraux: doubles(&object, "qraux"),
                    qy: doubles(&object, "qy"),
                }
            })
            .collect()
    }

    fn split_columns(values: &[f64], rows: usize, cols: usize) -> Vec<Vec<f64>> {
        (0..cols)
            .map(|j| values[j * rows..(j + 1) * rows].to_vec())
            .collect()
    }

    fn row_major(columns: &[Vec<f64>], rows: usize) -> Vec<f64> {
        (0..rows)
            .flat_map(|i| columns.iter().map(move |column| column[i]))
            .collect()
    }

    fn assert_bits(got: &[f64], want: &[f64], what: &str, name: &str) {
        assert_eq!(got.len(), want.len(), "{name}: {what} length");
        for (index, (g, w)) in got.iter().zip(want).enumerate() {
            assert_eq!(
                g.to_bits(),
                w.to_bits(),
                "{name}: {what}[{index}] = {g:e}, recorded {w:e}"
            );
        }
    }

    /// Factorize a recorded case and check rank, pivot and `qraux` exactly.
    fn check_factorization(case: &Case) -> HouseholderQr {
        let (n, p, name) = (case.rows, case.cols, case.name.as_str());
        assert_eq!(case.x.len(), n * p, "{name}: input size");
        let columns = split_columns(&case.x, n, p);
        let qr = HouseholderQr::factorize(columns.clone(), n, case.tolerance).unwrap();
        assert_eq!(qr.rank(), case.rank, "{name}: rank");
        let pivot: Vec<usize> = qr.pivot().iter().map(|j| j + 1).collect();
        assert_eq!(pivot, case.pivot, "{name}: pivot");
        let rank = matrix_rank_r_qr(&row_major(&columns, n), n, p, case.tolerance).unwrap();
        assert_eq!(rank, case.rank, "{name}: rank from row-major input");
        // Includes the downdated remaining norms stored at positions
        // without a reflector.
        assert_bits(qr.qraux(), &case.qraux, "qraux", name);
        qr
    }

    #[test]
    fn reproduces_recorded_downdate_cases() {
        let cases = cases(DOWNDATE_CASES);
        assert_eq!(cases.len(), 170);
        for case in &cases {
            if case.name.starts_with("contr_poly_") && !cfg!(target_os = "macos") {
                continue;
            }
            check_factorization(case);
        }
    }

    #[test]
    fn reproduces_recorded_factorizations() {
        let cases = cases(CASES);
        // The observation file currently records 65 cases.
        assert_eq!(cases.len(), 65);
        for case in &cases {
            let (n, p, name) = (case.rows, case.cols, case.name.as_str());
            let mut qr = check_factorization(case);

            let compact: Vec<f64> = (0..p).flat_map(|j| qr.column(j).to_vec()).collect();
            assert_bits(&compact, &case.qr, "qr", name);

            let mut basis = Vec::with_capacity(n * p);
            qr.visit_basis(|j, column| {
                assert_eq!(j * n, basis.len(), "{name}: visiting order");
                basis.extend_from_slice(column);
            })
            .unwrap();
            assert_bits(&basis, &case.qy, "basis", name);
        }
    }

    #[test]
    fn recomputation_threshold_matches_r_probe() {
        // qr(cbind(c(1, 0), c(1, s)))$qraux[2] in R 4.6.2: the factor
        // 1 - (1/rn)^2 is 9.999999997177891e-7 for the first s (R stores the
        // recomputed |s|) and 1.0000000001618781e-6 for the second (R stores
        // the downdated norm, which differs from |s| in the last bits).
        for (s, recomputed) in [(0.001000000499947672, true), (0.0010000004999476916, false)] {
            let qr =
                HouseholderQr::factorize(vec![vec![1.0, 0.0], vec![1.0, s]], 2, 1e-300).unwrap();
            assert_eq!(qr.qraux()[1] == s, recomputed, "s = {s:e}");
        }
    }

    /// Small deterministic generator for the property tests.
    struct Lcg(u64);

    impl Lcg {
        fn next(&mut self) -> f64 {
            self.0 = self
                .0
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            ((self.0 >> 11) as f64) / ((1_u64 << 53) as f64) * 2.0 - 1.0
        }

        fn below(&mut self, bound: usize) -> usize {
            ((self.next() + 1.0) / 2.0 * bound as f64) as usize % bound
        }
    }

    fn random_matrix(rng: &mut Lcg) -> (usize, usize, Vec<Vec<f64>>) {
        let n = 1 + rng.below(9);
        let p = 1 + rng.below(9);
        let mut columns: Vec<Vec<f64>> = (0..p)
            .map(|_| (0..n).map(|_| rng.next()).collect())
            .collect();
        // Make some columns zero or copies of earlier ones.
        for j in 0..p {
            match rng.below(5) {
                0 => columns[j].fill(0.0),
                1 if j > 0 => {
                    let source = rng.below(j);
                    let factor = rng.next() * 3.0;
                    columns[j] = columns[source].iter().map(|v| v * factor).collect();
                }
                _ => {}
            }
        }
        (n, p, columns)
    }

    #[test]
    fn rank_never_exceeds_either_dimension() {
        let mut rng = Lcg(17);
        for _ in 0..500 {
            let (n, p, columns) = random_matrix(&mut rng);
            let qr = HouseholderQr::factorize(columns, n, 1e-7).unwrap();
            assert!(qr.rank() <= n.min(p), "{n} x {p}: rank {}", qr.rank());
            let mut pivot = qr.pivot().to_vec();
            pivot.sort_unstable();
            assert_eq!(pivot, (0..p).collect::<Vec<_>>());
        }
    }

    #[test]
    fn rank_is_invariant_to_power_of_two_column_scaling() {
        let mut rng = Lcg(29);
        for case in 0..500 {
            let (n, p, columns) = random_matrix(&mut rng);
            let base = HouseholderQr::factorize(columns.clone(), n, 1e-7)
                .unwrap()
                .rank();
            let mut scaled = columns;
            let column = rng.below(p);
            let exponent = rng.below(41) as i32 - 20;
            for value in &mut scaled[column] {
                *value *= 2.0_f64.powi(exponent);
            }
            let rank = HouseholderQr::factorize(scaled, n, 1e-7).unwrap().rank();
            assert_eq!(
                rank, base,
                "case {case}: column {column} scaled by 2^{exponent}"
            );
        }
    }

    #[test]
    fn orthogonal_designs_have_full_rank() {
        // Columns of a Hadamard-type design (+1/-1 entries, mutually
        // orthogonal), with any subset of columns and any column scaling.
        let n = 16;
        let hadamard = |i: usize, j: usize| {
            if (i & j).count_ones().is_multiple_of(2) {
                1.0
            } else {
                -1.0
            }
        };
        for p in 1..=n {
            let columns: Vec<Vec<f64>> = (0..p)
                .map(|j| (0..n).map(|i| hadamard(i, j) * (j as f64 + 0.5)).collect())
                .collect();
            let qr = HouseholderQr::factorize(columns, n, 1e-7).unwrap();
            assert_eq!(qr.rank(), p);
            assert_eq!(qr.pivot(), (0..p).collect::<Vec<_>>().as_slice());
        }
        let identity: Vec<f64> = (0..5 * 5)
            .map(|k| if k / 5 == k % 5 { 1.0 } else { 0.0 })
            .collect();
        assert_eq!(matrix_rank_r_qr(&identity, 5, 5, 1e-7).unwrap(), 5);
    }

    #[test]
    fn rejects_invalid_input() {
        assert!(HouseholderQr::factorize(vec![vec![1.0]], 1, -1.0).is_err());
        assert!(HouseholderQr::factorize(vec![vec![1.0]], 1, f64::NAN).is_err());
        assert!(HouseholderQr::factorize(vec![vec![1.0, 2.0]], 1, 1e-7).is_err());
        assert!(matrix_rank_r_qr(&[1.0, f64::INFINITY], 1, 2, 1e-7).is_err());
        assert!(matrix_rank_r_qr(&[1.0, 2.0], 1, 3, 1e-7).is_err());
        assert!(matrix_rank_r_qr(&[1.0], 1, 1, f64::INFINITY).is_err());
        assert_eq!(matrix_rank_r_qr(&[], 0, 3, 1e-7).unwrap(), 0);
        assert_eq!(matrix_rank_r_qr(&[], 3, 0, 1e-7).unwrap(), 0);
    }
}
