//! Special functions and the gamma distribution used by the pipeline.
//!
//! Written from published methods (see `docs/provenance/specs/special_functions.md`):
//!
//! - log-gamma and gamma reduce the argument to `2 + z` with `|z| <= 1/2`
//!   by the recurrence `Gamma(x + 1) = x Gamma(x)` (DLMF 5.5.1), carrying the
//!   product of the recurrence factors in double-length arithmetic (Dekker,
//!   1971, Numer. Math. 18:224-242), and evaluate `ln Gamma(2 + z)` from its
//!   Taylor series, whose coefficients are `(-1)^k (zeta(k) - 1) / k`
//!   (DLMF 5.7.3 combined with `ln(1 + z)`). Large arguments use Stirling's
//!   series (DLMF 5.11.1) and very negative ones the reflection formula
//!   (DLMF 5.5.3).
//! - digamma uses its Taylor series about the positive zero
//!   `x0 = 1.46163...`, the recurrence (DLMF 5.5.2), the asymptotic
//!   expansion (DLMF 5.11.2) and the reflection formula (DLMF 5.5.4).
//! - the regularized incomplete gamma ratios use the power series
//!   (DLMF 8.7.1), the even part of the continued fraction (DLMF 8.9.2)
//!   evaluated with the modified Lentz method (Thompson and Barnett, 1986,
//!   J. Comput. Phys. 64:490-509), and for large shapes near the transition
//!   Temme's uniform asymptotic expansion (Temme, 1979, SIAM J. Math. Anal.
//!   10:757-766; DLMF section 8.12).
//! - the normal quantile is Wichura's Algorithm AS 241 (Applied Statistics
//!   37:477-484, 1988).
//!
//! All series coefficients were computed for this file in 50-digit
//! arithmetic from their defining formulas and rounded to the nearest double.
//! Elementary functions go through [`ProfileMath`], so the `portable`
//! profile gives the same bits on every host.

use crate::numeric::math::ProfileMath;
use crate::probability::normal_upper_tail;

const PI: f64 = std::f64::consts::PI;
const HALF_LN_2PI: f64 = 0.9189385332046728;
const LN_PI: f64 = 1.1447298858494002;
const FRAC_1_SQRT_2PI: f64 = 0.3989422804014327;

/// Positive arguments below this make R's `gamma` report overflow (the
/// thresholds here were found by bisection on R 4.6.2).
const GAMMA_TINY: f64 = 2.2474362225598545e-308;
/// Smallest argument for which R's `gamma` reports overflow.
const GAMMA_OVERFLOW: f64 = 171.614478871823;
/// Arguments below this give 0 from R's `gamma` (the value would be
/// close to or below the smallest normal double).
const GAMMA_UNDERFLOW: f64 = -170.5674972726612;

/// Evaluate a polynomial with ascending coefficients by Horner's rule.
fn polynomial(x: f64, ascending: &[f64]) -> f64 {
    ascending.iter().rev().fold(0.0, |acc, &c| acc * x + c)
}

/// A product carried as an unevaluated sum `hi + lo` (Dekker's
/// double-length multiplication via fused multiply-add).
#[derive(Clone, Copy)]
struct Product {
    hi: f64,
    lo: f64,
}

impl Product {
    fn new(start: f64) -> Self {
        Self { hi: start, lo: 0.0 }
    }

    fn times(self, factor: f64) -> Self {
        let hi = self.hi * factor;
        let lo = self.hi.mul_add(factor, -hi) + self.lo * factor;
        let sum = hi + lo;
        Self {
            hi: sum,
            lo: lo - (sum - hi),
        }
    }

    fn value(self) -> f64 {
        self.hi + self.lo
    }

    fn ln_abs(self) -> f64 {
        self.hi.abs().ln_p() + self.lo / self.hi
    }
}

/// Writes `x = n + z` with `n` the nearest integer, so `|z| <= 1/2`, and
/// returns `z` with the product of recurrence factors linking `Gamma(x)` to
/// `Gamma(2 + z)`: `Gamma(x) = Gamma(2 + z) * product` when the flag is false
/// and `Gamma(2 + z) / product` when it is true. Every factor `x - j`
/// (`x >= 1.5`) or `x + m` (`x < 1.5`) is exact unless its magnitude exceeds
/// `|x|`, which happens only for `|x| < 1.5`. `start` scales the product to
/// avoid overflow.
fn recurrence(x: f64, start: f64) -> (f64, Product, bool) {
    let n = x.round();
    let z = x - n;
    let mut product = Product::new(start);
    if n >= 2.0 {
        let mut j = 1.0;
        while j <= n - 2.0 {
            product = product.times(x - j);
            j += 1.0;
        }
        (z, product, false)
    } else {
        let mut m = 0.0;
        while m <= 1.0 - n {
            product = product.times(x + m);
            m += 1.0;
        }
        (z, product, true)
    }
}

/// Taylor coefficients of `ln Gamma(2 + z) / z` about `z = 0`: `1 - gamma`, then
/// `(-1)^k (zeta(k) - 1) / k` for `k >= 2` (DLMF 5.7.3).
const LN_GAMMA_2P: [f64; 30] = [
    0.42278433509846713,
    0.3224670334241132,
    -0.0673523010531981,
    0.020580808427784546,
    -0.007385551028673986,
    0.0028905103307415234,
    -0.001192753911703261,
    0.0005096695247430425,
    -0.00022315475845357939,
    9.945751278180853e-05,
    -4.492623673813314e-05,
    2.050721277567069e-05,
    -9.439488275268397e-06,
    4.374866789907488e-06,
    -2.039215753801366e-06,
    9.55141213040742e-07,
    -4.492469198764566e-07,
    2.1207184805554665e-07,
    -1.0043224823968099e-07,
    4.7698101693639804e-08,
    -2.2711094608943164e-08,
    1.0838659214896955e-08,
    -5.183475041970047e-09,
    2.4836745438024785e-09,
    -1.1921401405860912e-09,
    5.731367241678862e-10,
    -2.7595228851242334e-10,
    1.330476437424449e-10,
    -6.4229645638381e-11,
    3.1044247747322276e-11,
];
/// The positive zero of digamma, `1.4616321449683623412626595423257...`,
/// as an unevaluated sum of two doubles.
const PSI_ROOT_HI: f64 = 1.4616321449683622;
const PSI_ROOT_LO: f64 = 9.549995429965697e-17;
/// Taylor coefficients `psi^(k)(x0) / k!`, `k >= 1`, of digamma about its
/// positive zero `x0`.
const PSI_ROOT_TAYLOR: [f64; 43] = [
    0.9676722454476212,
    -0.4427631689835921,
    0.258499760955651,
    -0.16394270544240652,
    0.10782405069126237,
    -0.07219956125645471,
    0.04880428816414311,
    -0.03316112647484736,
    0.022597648232218104,
    -0.01542476590494896,
    0.010538791616612175,
    -0.007204534386356869,
    0.004926781395729853,
    -0.003369801655439328,
    0.002305126326734928,
    -0.0015769367714301972,
    0.0010788252019162967,
    -0.0007380709389960052,
    0.000504953265834602,
    -0.0003454680251063077,
    0.00023635601564027053,
    -0.00016170622091974803,
    0.0001106337276874741,
    -7.569179582195066e-05,
    5.178575795222081e-05,
    -3.5430070947659604e-05,
    2.424006611860132e-05,
    -1.6584242271854135e-05,
    1.134638458466385e-05,
    -7.762817668462094e-06,
    5.3110609208898636e-06,
    -3.6336507898010456e-06,
    2.486022733129538e-06,
    -1.7008538854332607e-06,
    1.1636675363548843e-06,
    -7.96142543124197e-07,
    5.446941930669446e-07,
    -3.7266161283438227e-07,
    2.549626552021554e-07,
    -1.7443695117727745e-07,
    1.1934394829830244e-07,
    -8.165115189488409e-08,
    5.586299683532171e-08,
];
/// Taylor coefficients in `eta` of Temme's `c_0(eta), ..., c_10(eta)`
/// (DLMF section 8.12), truncated where they no longer matter for
/// `a >= 20` and `|eta| <= 1`.
const TEMME: [&[f64]; 11] = [
    &[
        -0.3333333333333333,
        0.08333333333333333,
        -0.014814814814814815,
        0.0011574074074074073,
        0.0003527336860670194,
        -0.0001787551440329218,
        3.919263178522438e-05,
        -2.185448510679992e-06,
        -1.85406221071516e-06,
        8.296711340953087e-07,
        -1.7665952736826078e-07,
        6.707853543401498e-09,
        1.0261809784240309e-08,
        -4.382036018453353e-09,
        9.14769958223679e-10,
        -2.5514193994946248e-11,
        -5.830772132550426e-11,
        2.4361948020667415e-11,
        -5.0276692801141755e-12,
        1.1004392031956135e-13,
        3.371763262400985e-13,
        -1.392388722418162e-13,
        2.8534893807047445e-14,
        -5.139111834242572e-16,
        -1.9752288294349442e-15,
        8.099521156704561e-16,
        -1.6522531216398162e-16,
        2.5305430097478883e-18,
        1.1686939738559576e-17,
        -4.770037049820485e-18,
        9.699126059056237e-19,
    ],
    &[
        -0.001851851851851852,
        -0.003472222222222222,
        0.0026455026455026454,
        -0.0009902263374485596,
        0.00020576131687242798,
        -4.018775720164609e-07,
        -1.8098550334489977e-05,
        7.64916091608111e-06,
        -1.6120900894563446e-06,
        4.647127802807434e-09,
        1.378633446915721e-07,
        -5.752545603517705e-08,
        1.1951628599778148e-08,
        -1.7543241719747647e-11,
        -1.0091543710600413e-09,
        4.162792991842583e-10,
        -8.56390702649298e-11,
        6.067215101604758e-14,
        7.1624989648114856e-12,
        -2.933186643771437e-12,
        5.996696365683689e-13,
        -2.1671786527323313e-16,
        -4.978339972369262e-14,
        2.0291628823713425e-14,
        -4.13125571381061e-15,
        8.286516239883097e-19,
        3.4100308869333327e-16,
        -1.3854195302893971e-16,
        2.812346653228875e-17,
        -3.406444194143029e-21,
        -2.3109797315115572e-18,
    ],
    &[
        0.004133597883597883,
        -0.0026813271604938273,
        0.0007716049382716049,
        2.0093878600823047e-06,
        -0.0001073665322636516,
        5.2923448829120125e-05,
        -1.2760635188618728e-05,
        3.423578734096138e-08,
        1.3721957309062934e-06,
        -6.298992138380055e-07,
        1.4280614206064242e-07,
        -2.0477098421990866e-10,
        -1.409252991086752e-08,
        6.228974084922022e-09,
        -1.3670488396617114e-09,
        9.428356159014678e-13,
        1.2872252400089318e-10,
        -5.5645956134363323e-11,
        1.197593554636698e-11,
        -4.1689782251838634e-15,
        -1.0940640427884595e-12,
        4.662239946390136e-13,
        -9.905105763906907e-14,
        1.8931876768373515e-17,
        8.859221872591127e-15,
        -3.737820398046405e-15,
        7.868833639035156e-16,
        -9.000027395741211e-20,
        -6.928881229347671e-17,
    ],
    &[
        0.0006494341563786008,
        0.00022947209362139917,
        -0.0004691894943952557,
        0.00026772063206283885,
        -7.561801671883977e-05,
        -2.396505113867297e-07,
        1.1082654115347302e-05,
        -5.6749528269915965e-06,
        1.4230900732435883e-06,
        -2.7861080291528143e-11,
        -1.6958404091930278e-07,
        8.099464905388083e-08,
        -1.9111168485973655e-08,
        2.3928620439808118e-12,
        2.0620131815488797e-09,
        -9.460496661855133e-10,
        2.1541049775774907e-10,
        -1.388823336813903e-14,
        -2.1894761681963938e-11,
        9.790998951171684e-12,
        -2.178219188018096e-12,
        6.208819573407901e-17,
        2.126978363279737e-13,
        -9.344688791517433e-14,
        2.045367122678285e-14,
        -2.58260790403495e-19,
        -1.9405297673344544e-15,
        8.415979290484816e-16,
    ],
    &[
        -0.0008618882909167117,
        0.0007840392217200666,
        -0.0002990724803031902,
        -1.4638452578843418e-06,
        6.641498215465122e-05,
        -3.968365047179435e-05,
        1.1375726970678419e-05,
        2.507497226237533e-10,
        -1.6954149536558305e-06,
        8.907507532205309e-07,
        -2.292934834000805e-07,
        2.956794137544049e-11,
        2.8865829742708783e-08,
        -1.4189739437803219e-08,
        3.4463580499464896e-09,
        -2.3024517174528067e-13,
        -3.9409233028046403e-10,
        1.86023389685045e-10,
        -4.356323005056618e-11,
        1.278600101629623e-15,
        4.67927502665792e-12,
        -2.149246470613483e-12,
        4.908815614809652e-13,
        -6.33859148489156e-18,
        -5.045332069080094e-14,
        2.2722958222901286e-14,
    ],
    &[
        -0.00033679855336635813,
        -6.972813758365857e-05,
        0.0002772753244959392,
        -0.00019932570516188847,
        6.797780477937208e-05,
        1.419062920643967e-07,
        -1.3594048189768693e-05,
        8.018470256334202e-06,
        -2.291481176508095e-06,
        -3.252473551298454e-10,
        3.4652846491085265e-07,
        -1.8447187191171344e-07,
        4.8240967037894184e-08,
        -1.7989466721743514e-14,
        -6.306194500013523e-09,
        3.162417628774568e-09,
        -7.840924253697429e-10,
        5.192679165254041e-15,
        9.358944242306784e-11,
        -4.513426216163278e-11,
        1.0799129993116828e-11,
        -3.661886712685252e-17,
        -1.210902069055155e-12,
        5.680743584990564e-13,
    ],
    &[
        0.0005313079364639922,
        -0.0005921664373536939,
        0.0002708782096718045,
        7.902353232660328e-07,
        -8.153969367561969e-05,
        5.61168275310625e-05,
        -1.8329116582843375e-05,
        -3.0796134506033047e-09,
        3.465155368803609e-06,
        -2.0291327396058603e-06,
        5.788792863149004e-07,
        2.338630673826657e-13,
        -8.828600746330484e-08,
        4.7435958880408125e-08,
        -1.2545415020710383e-08,
        8.649648858010293e-14,
        1.6846058979264062e-09,
        -8.575492823577594e-10,
        2.1598224929232125e-10,
        -7.613230520476153e-16,
        -2.6639822008536144e-11,
        1.3065700536611057e-11,
    ],
    &[
        0.00034436760689237765,
        5.171790908260592e-05,
        -0.00033493161081142234,
        0.0002812695154763237,
        -0.00010976582244684731,
        -1.2741009095484485e-07,
        2.7744451511563645e-05,
        -1.8263488805711332e-05,
        5.7876949497350525e-06,
        4.93875893393627e-10,
        -1.0595367014026043e-06,
        6.166714376110408e-07,
        -1.7562973359060463e-07,
        -1.297447328701544e-12,
        2.695423606288966e-08,
        -1.4578352908731272e-08,
        3.887645959386175e-09,
        -3.881002251019412e-17,
        -5.327994173877286e-10,
        2.7437977643314844e-10,
    ],
    &[
        -0.0006526239185953094,
        0.0008394987206720873,
        -0.000438297098541721,
        -6.969091458420552e-07,
        0.00016644846642067547,
        -0.00012783517679769218,
        4.629953263691304e-05,
        4.557909867922708e-09,
        -1.0595271125805195e-05,
        6.783342904865167e-06,
        -2.1075476666258803e-06,
        -1.7213731432817144e-11,
        3.773587741611098e-07,
        -2.1867506700122867e-07,
        6.220228804018927e-08,
        6.597703826733e-16,
        -9.590386497425686e-09,
        5.213214492280807e-09,
    ],
    &[
        -0.0005967612901927463,
        -7.204895416020011e-05,
        0.0006782308837667328,
        -0.0006401475260262758,
        0.00027750107634328704,
        1.819700838046515e-07,
        -8.479507117068503e-05,
        6.105192082501531e-05,
        -2.1073920183404862e-05,
        -8.858589014125599e-10,
        4.5284535953805374e-06,
        -2.8427815022504407e-06,
        8.708234177864641e-07,
        3.6886101871706966e-12,
        -1.534469519070206e-07,
        8.862466778790695e-08,
    ],
    &[
        0.0013324454494800656,
        -0.0019144384985654776,
        0.0011089369134596636,
        9.9324041226423e-07,
        -0.0005087450129309319,
        0.00042735056665392886,
        -0.00016858853767910798,
        -8.1301893922785e-09,
        4.5284402370562144e-05,
        -3.127053674781734e-05,
        1.044986828530338e-05,
        4.8435226265680926e-11,
        -2.148256587345626e-06,
        1.329369701097492e-06,
    ],
];

/// `ln Gamma(2 + z)` for `|z| <= 1/2` from its Taylor series about 2.
fn ln_gamma_2p(z: f64) -> f64 {
    z * polynomial(z, &LN_GAMMA_2P)
}

/// `ln Gamma(1 + a)` for `a >= 0`, accurate near `a = 0` and `a = 1`.
fn ln_gamma_1p(a: f64) -> f64 {
    if a < 0.5 {
        ln_gamma_2p(a) - a.ln_1p_p()
    } else if a <= 1.5 {
        ln_gamma_2p(a - 1.0)
    } else {
        ln_gamma(a) + a.ln_p()
    }
}

/// Remainder of Stirling's series, `ln Gamma(x) - [(x - 1/2) ln x - x +
/// ln(2 pi) / 2]`, for `x >= 10` (DLMF 5.11.1, coefficients
/// `B_2k / (2k (2k - 1))`).
fn stirling_series(x: f64) -> f64 {
    let inverse = 1.0 / x;
    let square = inverse * inverse;
    inverse
        * polynomial(
            square,
            &[
                1.0 / 12.0,
                -1.0 / 360.0,
                1.0 / 1260.0,
                -1.0 / 1680.0,
                1.0 / 1188.0,
                -691.0 / 360360.0,
                1.0 / 156.0,
                -3617.0 / 122400.0,
                43867.0 / 244188.0,
            ],
        )
}

/// `sin(pi r)` for `|r| <= 1/2`.
fn sin_pi_reduced(r: f64) -> f64 {
    (PI * r).sin_p()
}

/// `cos(pi r)` for `|r| <= 1/2`, accurate near `|r| = 1/2`.
fn cos_pi_reduced(r: f64) -> f64 {
    if r.abs() < 0.25 {
        (PI * r).cos_p()
    } else {
        sin_pi_reduced(0.5 - r.abs())
    }
}

/// `sin(pi x)` with exact argument reduction.
fn sin_pi(x: f64) -> f64 {
    // x - 2 round(x / 2) is exact and lies in [-1, 1].
    let mut r = x - 2.0 * (0.5 * x).round();
    if r > 0.5 {
        r = 1.0 - r;
    } else if r < -0.5 {
        r = -1.0 - r;
    }
    sin_pi_reduced(r)
}

/// `pi / tan(pi x)` with exact argument reduction (period 1).
fn pi_cot_pi(x: f64) -> f64 {
    let r = x - x.round();
    PI * cos_pi_reduced(r) / sin_pi_reduced(r)
}

/// Natural log of the absolute value of the gamma function.
///
/// `+Inf` at zero, at the negative integers and at both infinities; NaN
/// propagates. The relative error is a few units in the last place except
/// near the zeros at 1 and 2 (where the absolute error is about 1e-16) and
/// near the zeros on the negative axis.
pub fn ln_gamma(x: f64) -> f64 {
    if x.is_nan() {
        return x;
    }
    if x.is_infinite() || (x <= 0.0 && x == x.floor()) {
        return f64::INFINITY;
    }
    if x >= 10.0 {
        return (x - 0.5) * x.ln_p() - x + HALF_LN_2PI + stirling_series(x);
    }
    if x > GAMMA_UNDERFLOW {
        let (z, product, divide) = recurrence(x, 1.0);
        let log_product = product.ln_abs();
        let reduced = ln_gamma_2p(z);
        return if divide {
            reduced - log_product
        } else {
            reduced + log_product
        };
    }
    // Reflection, DLMF 5.5.3: Gamma(x) Gamma(1 - x) = pi / sin(pi x).
    LN_PI - sin_pi(x).abs().ln_p() - ln_gamma(1.0 - x)
}

/// The gamma function for `GAMMA_UNDERFLOW <= x < GAMMA_OVERFLOW`, `x` not
/// a non-positive integer, without R's range limits.
fn gamma_in_range(x: f64) -> f64 {
    // Scale the product of up to 170 factors so it cannot overflow before the
    // final multiplication.
    const SCALE: f64 = 1.0 / 4294967296.0;
    let start = if x > 160.0 { SCALE } else { 1.0 };
    let (z, product, divide) = recurrence(x, start);
    let reduced = ln_gamma_2p(z).exp_p();
    if divide {
        reduced / product.value()
    } else {
        reduced * product.value() / start
    }
}

/// The gamma function.
///
/// NaN at zero, at the negative integers, at `-Inf` and for NaN; `+Inf` at
/// `+Inf`. As in R, positive arguments below about `2.2e-308` and arguments
/// from about 171.6145 report overflow (`+Inf`), and non-integer arguments
/// below about -170.567 give 0. The relative error is a few units in the last
/// place.
pub fn gamma(x: f64) -> f64 {
    if x.is_nan() {
        return x;
    }
    if x == f64::INFINITY {
        return x;
    }
    if x <= 0.0 && x == x.floor() {
        return f64::NAN;
    }
    if x > 0.0 && x < GAMMA_TINY {
        return f64::INFINITY;
    }
    if x >= GAMMA_OVERFLOW {
        return f64::INFINITY;
    }
    if x < GAMMA_UNDERFLOW {
        return 0.0;
    }
    gamma_in_range(x)
}

/// Digamma for `x >= 10` from the asymptotic expansion (DLMF 5.11.2).
fn digamma_asymptotic(x: f64) -> f64 {
    let inverse = 1.0 / x;
    let square = inverse * inverse;
    x.ln_p()
        - 0.5 * inverse
        - square
            * polynomial(
                square,
                &[
                    1.0 / 12.0,
                    -1.0 / 120.0,
                    1.0 / 252.0,
                    -1.0 / 240.0,
                    1.0 / 132.0,
                    -691.0 / 32760.0,
                    1.0 / 12.0,
                    -3617.0 / 8160.0,
                ],
            )
}

/// Digamma from its Taylor series in `t = x - x0` about the positive zero
/// `x0`, for `|t| <= 0.54`.
fn digamma_near_root(t: f64) -> f64 {
    t * polynomial(t, &PSI_ROOT_TAYLOR)
}

/// Digamma for finite `x > 0`.
fn digamma_positive_argument(x: f64) -> f64 {
    if x >= 10.0 {
        return digamma_asymptotic(x);
    }
    if x >= 2.0 {
        // Downward recurrence psi(y + 1) = psi(y) + 1 / y with exact y.
        let mut y = x;
        let mut sum = 0.0;
        while y >= 2.0 {
            y -= 1.0;
            sum += 1.0 / y;
        }
        return digamma_near_root((y - PSI_ROOT_HI) - PSI_ROOT_LO) + sum;
    }
    if x >= 1.0 {
        return digamma_near_root((x - PSI_ROOT_HI) - PSI_ROOT_LO);
    }
    // psi(x) = psi(x + 1) - 1 / x, with x + 1 - x0 formed without rounding
    // x + 1 first.
    digamma_near_root((x - (PSI_ROOT_HI - 1.0)) - PSI_ROOT_LO) - 1.0 / x
}

/// The digamma function `psi(x)`.
///
/// NaN at zero, the negative integers, `-Inf` and NaN; `+Inf` at `+Inf`.
/// Relative error a few units in the last place for positive arguments,
/// including near the zero at 1.4616.
pub fn digamma(x: f64) -> f64 {
    if x.is_nan() || x == f64::INFINITY {
        return x;
    }
    if x <= 0.0 && x == x.floor() {
        return f64::NAN;
    }
    if x < 0.0 {
        // Reflection, DLMF 5.5.4: psi(1 - x) - psi(x) = pi cot(pi x).
        return digamma_positive_argument(1.0 - x) - pi_cot_pi(x);
    }
    digamma_positive_argument(x)
}

/// `ln(1 + u) - u` for `u > -1`. For `-1/2 <= u <= 1` it uses the series
/// of `ln(1 + u) = 2 artanh(u / (2 + u))`, so the leading `-u^2/2`
/// is formed without cancellation.
fn log1pmx(u: f64) -> f64 {
    if !(-0.5..=1.0).contains(&u) {
        return u.ln_1p_p() - u;
    }
    let t = u / (2.0 + u);
    let square = t * t;
    let mut power = t * square;
    let mut odd = 3.0;
    let mut sum = 0.0;
    for _ in 0..64 {
        let term = power / odd;
        sum += term;
        if term.abs() <= sum.abs() * 1e-17 {
            break;
        }
        power *= square;
        odd += 2.0;
    }
    (2.0 * sum - u * t).min(0.0)
}

/// `a ln(x / a) - (x - a)`, the exponent of `x^a e^-x` relative to its value
/// at `x = a`.
fn scaled_log1pmx(a: f64, x: f64) -> f64 {
    let u = (x - a) / a;
    if (-0.5..=1.0).contains(&u) {
        a * log1pmx(u)
    } else {
        a * (x / a).ln_p() - (x - a)
    }
}

/// `ln(x^a e^-x / Gamma(a))` for `a > 0`, `x > 0`. For `a >= 10` Stirling's
/// series (DLMF 5.11.1) cancels the large terms analytically.
fn ln_power_exp_over_gamma(a: f64, x: f64) -> f64 {
    if a < 10.0 {
        a * x.ln_p() - x - ln_gamma(a)
    } else {
        scaled_log1pmx(a, x) + 0.5 * (a / (2.0 * PI)).ln_p() - stirling_series(a)
    }
}

/// `P(a, x)` from the power series of DLMF 8.7.1,
/// `P = x^a e^-x / Gamma(a + 1) * sum_k x^k / ((a + 1) ... (a + k))`.
fn lower_series(a: f64, x: f64) -> f64 {
    let log_front = if a < 10.0 {
        a * x.ln_p() - x - ln_gamma_1p(a)
    } else {
        ln_power_exp_over_gamma(a, x) - a.ln_p()
    };
    let mut term = 1.0;
    let mut sum = 1.0;
    let mut k = 1.0;
    while k < 100_000.0 {
        term *= x / (a + k);
        sum += term;
        if term <= sum * 1e-17 {
            break;
        }
        k += 1.0;
    }
    log_front.exp_p() * sum
}

/// `Q(a, x)` for `a < 1` and `x < 2` from the series of DLMF 8.7.1 for
/// `gamma(a, x)`: with `t = a ln x - ln Gamma(1 + a)` and
/// `S = sum_{k>=1} (-x)^k / (k! (a + k))`,
/// `Q = -expm1(t) - a e^t S`, which keeps relative accuracy when `Q` is small.
fn upper_small_shape(a: f64, x: f64) -> f64 {
    let t = a * x.ln_p() - ln_gamma_1p(a);
    let mut power = 1.0;
    let mut sum = 0.0;
    let mut k = 1.0;
    while k < 200.0 {
        power *= -x / k;
        let term = power / (a + k);
        sum += term;
        if term.abs() <= sum.abs() * 1e-17 {
            break;
        }
        k += 1.0;
    }
    -t.exp_m1_p() - a * t.exp_p() * sum
}

/// `Q(a, x)` for `x >= a + 1` from the even part of the continued fraction
/// of DLMF 8.9.2,
/// `Q = x^a e^-x / Gamma(a) / (x + 1 - a - 1 (1 - a) / (x + 3 - a - 2 (2 - a) / ...))`,
/// evaluated with the modified Lentz method.
fn upper_fraction(a: f64, x: f64) -> f64 {
    const TINY: f64 = 1e-300;
    let mut value = x + 1.0 - a;
    if value.abs() < TINY {
        value = TINY;
    }
    let mut c = value;
    let mut d = 0.0;
    let mut k = 1.0;
    while k < 100_000.0 {
        let numerator = -k * (k - a);
        let denominator = x + 2.0 * k + 1.0 - a;
        d = denominator + numerator * d;
        if d.abs() < TINY {
            d = TINY;
        }
        c = denominator + numerator / c;
        if c.abs() < TINY {
            c = TINY;
        }
        d = 1.0 / d;
        let delta = c * d;
        value *= delta;
        if (delta - 1.0).abs() <= 1e-16 {
            break;
        }
        k += 1.0;
    }
    ln_power_exp_over_gamma(a, x).exp_p() / value
}

/// `(P, Q)` from Temme's uniform expansion (DLMF section 8.12) for
/// `a >= 20` and `|eta| <= 1`, where `log_ratio = ln(x/a) - (x/a - 1)` and
/// `u = x/a - 1`, so `eta = sign(u) sqrt(-2 log_ratio)`:
/// `Q = erfc(eta sqrt(a/2)) / 2 + R` and
/// `P = erfc(-eta sqrt(a/2)) / 2 - R`, with
/// `R ~ exp(-a eta^2 / 2) / sqrt(2 pi a) sum_k c_k(eta) a^-k`,
/// `c_0 = 1 / (x/a - 1) - 1 / eta` and
/// `c_k = c_{k-1}'(eta) / eta + (-1)^k gamma_k / (x/a - 1)` with `gamma_k` the
/// Stirling coefficients. The Taylor coefficients of `c_k` in `eta` were derived in
/// exact rational arithmetic.
fn temme(a: f64, u: f64, log_ratio: f64) -> (f64, f64) {
    let eta = (-2.0 * log_ratio).sqrt().copysign(u);
    let z = (-2.0 * a * log_ratio).sqrt().copysign(u);
    let inverse = 1.0 / a;
    let series = TEMME.iter().rev().fold(0.0, |acc, coefficients| {
        acc * inverse + polynomial(eta, coefficients)
    });
    let remainder = (a * log_ratio).exp_p() / (2.0 * PI * a).sqrt() * series;
    (
        normal_upper_tail(-z) - remainder,
        normal_upper_tail(z) + remainder,
    )
}

/// `(P(a, x), Q(a, x))` for finite `a > 0` and finite `x >= 0`. The smaller
/// of the two is always computed directly; the larger may be one minus the
/// smaller, which loses nothing in relative terms.
fn gamma_tails(a: f64, x: f64) -> (f64, f64) {
    if x == 0.0 {
        return (0.0, 1.0);
    }
    if a >= 20.0 {
        let u = (x - a) / a;
        if (-0.7..=1.2).contains(&u) {
            let log_ratio = log1pmx(u);
            if log_ratio >= -0.5 {
                return temme(a, u, log_ratio);
            }
        }
    }
    if x < a + 1.0 {
        let lower = lower_series(a, x);
        let upper = if a < 1.0 {
            upper_small_shape(a, x)
        } else {
            1.0 - lower
        };
        (lower, upper)
    } else {
        let upper = upper_fraction(a, x);
        (1.0 - upper, upper)
    }
}

fn check_incomplete_gamma_arguments(a: f64, x: f64) {
    assert!(
        a > 0.0 && a.is_finite() && x >= 0.0 && x.is_finite(),
        "incomplete gamma needs finite a > 0 and finite x >= 0, got a = {a}, x = {x}"
    );
}

/// Upper regularized incomplete gamma function `Q(a, x) = Gamma(a, x) /
/// Gamma(a)`.
///
/// NaN if either argument is NaN; `Q(a, 0) = 1`.
///
/// # Panics
///
/// Unless `a > 0` and finite and `x >= 0` and finite.
pub fn gamma_ur(a: f64, x: f64) -> f64 {
    if a.is_nan() || x.is_nan() {
        return f64::NAN;
    }
    check_incomplete_gamma_arguments(a, x);
    gamma_tails(a, x).1
}

/// Lower regularized incomplete gamma function `P(a, x) = gamma(a, x) /
/// Gamma(a)`.
///
/// NaN if either argument is NaN; `P(a, 0) = 0`.
///
/// # Panics
///
/// Unless `a > 0` and finite and `x >= 0` and finite.
pub fn gamma_lr(a: f64, x: f64) -> f64 {
    if a.is_nan() || x.is_nan() {
        return f64::NAN;
    }
    check_incomplete_gamma_arguments(a, x);
    gamma_tails(a, x).0
}

fn check_beta_arguments(a: f64, b: f64) {
    assert!(
        a > 0.0 && b > 0.0,
        "beta needs a > 0 and b > 0, got a = {a}, b = {b}"
    );
}

/// `Gamma(p) Gamma(q) / Gamma(p + q)` for `1e-300 <= p <= q` and
/// `p + q < GAMMA_OVERFLOW`.
fn beta_from_gammas(p: f64, q: f64) -> f64 {
    gamma_in_range(p) / gamma_in_range(p + q) * gamma_in_range(q)
}

/// Natural log of the beta function.
///
/// For a large argument `q` the log-gamma difference `ln Gamma(q) -
/// ln Gamma(p + q)` is formed from Stirling's series (DLMF 5.11.1) with the
/// leading terms combined analytically (`ln(1 + p/q)` and
/// `ln(1 + p/q) - p/q`), so nothing large is subtracted. `-Inf` when either
/// argument is infinite.
///
/// # Panics
///
/// Unless `a > 0` and `b > 0`.
pub fn ln_beta(a: f64, b: f64) -> f64 {
    check_beta_arguments(a, b);
    let (p, q) = if a <= b { (a, b) } else { (b, a) };
    if q.is_infinite() {
        return f64::NEG_INFINITY;
    }
    if q < 10.0 {
        if p < 1e-300 {
            return ln_gamma(p) + ln_gamma(q) - ln_gamma(p + q);
        }
        return beta_from_gammas(p, q).ln_p();
    }
    if p < 10.0 {
        let r = p / q;
        return ln_gamma(p) - p * q.ln_p() - p * r.ln_1p_p() + 0.5 * r - (q - 0.5) * log1pmx(r)
            + (stirling_series(q) - stirling_series(q + p));
    }
    HALF_LN_2PI
        - 0.5 * (p + q).ln_p()
        - (p - 0.5) * (q / p).ln_1p_p()
        - (q - 0.5) * (p / q).ln_1p_p()
        + (stirling_series(p) + stirling_series(q) - stirling_series(p + q))
}

/// The beta function `B(a, b) = Gamma(a) Gamma(b) / Gamma(a + b)`.
///
/// With `p <= q` the smaller and larger argument: for `q < 10`, or both at
/// least 10 with `p + q` below the gamma overflow threshold, the ratio of
/// gamma functions. For `p < 10 <= q`, `Gamma(p) q^-p` times the exponential
/// of the small remaining Stirling terms (DLMF 5.11.1), which never forms
/// `p + q`. Otherwise Stirling's series with the two large factors raised to
/// their powers separately. 0 when either argument is infinite.
///
/// # Panics
///
/// Unless `a > 0` and `b > 0`.
pub fn beta(a: f64, b: f64) -> f64 {
    check_beta_arguments(a, b);
    let (p, q) = if a <= b { (a, b) } else { (b, a) };
    if q.is_infinite() {
        return 0.0;
    }
    if p < 1e-307 {
        return ln_beta(p, q).exp_p();
    }
    if q < 10.0 || (p >= 10.0 && p + q < GAMMA_OVERFLOW) {
        return beta_from_gammas(p, q);
    }
    let r = p / q;
    if p < 10.0 {
        let small = -p * r.ln_1p_p() + 0.5 * r - (q - 0.5) * log1pmx(r)
            + (stirling_series(q) - stirling_series(q + p));
        return gamma_in_range(p) * q.powf_p(-p) * small.exp_p();
    }
    // sqrt(2 pi / s) (1 + q/p)^-(p - 1/2) (1 + p/q)^-(q - 1/2) exp(S(p) + S(q) - S(s)).
    let s = p + q;
    let corrections = stirling_series(p) + stirling_series(q) - stirling_series(s);
    (2.0 * PI / s).sqrt()
        * (1.0 + q / p).powf_p(-(p - 0.5))
        * (-(q - 0.5) * r.ln_1p_p() + corrections).exp_p()
}

/// Density of the standard normal distribution, `exp(-x^2/2) / sqrt(2 pi)`.
///
/// `x^2` is split exactly into `hi + lo` with a fused multiply-add, so the
/// exponent carries no rounding error even in the far tails.
pub fn standard_normal_pdf(x: f64) -> f64 {
    if x.is_nan() {
        return x;
    }
    if x.is_infinite() {
        return 0.0;
    }
    let hi = x * x;
    let lo = x.mul_add(x, -hi);
    FRAC_1_SQRT_2PI * (-0.5 * hi).exp_p() * (-0.5 * lo).exp_p()
}

/// Quantile of the standard normal distribution, Wichura's Algorithm AS 241
/// (`PPND16`, about 16 significant digits), with the coefficients as printed
/// in Applied Statistics 37:477-484 (1988).
///
/// `-Inf` at 0, `+Inf` at 1, NaN for NaN.
///
/// # Panics
///
/// If `p` lies outside `[0, 1]`.
pub fn standard_normal_inverse_cdf(p: f64) -> f64 {
    if p.is_nan() {
        return p;
    }
    assert!(
        (0.0..=1.0).contains(&p),
        "normal quantile needs 0 <= p <= 1, got {p}"
    );
    if p == 0.0 {
        return f64::NEG_INFINITY;
    }
    if p == 1.0 {
        return f64::INFINITY;
    }
    let q = p - 0.5;
    if q.abs() <= 0.425 {
        let r = 0.180625 - q * q;
        return q * polynomial(r, &AS241_A) / polynomial(r, &AS241_B);
    }
    let tail = if q < 0.0 { p } else { 1.0 - p };
    let mut r = (-tail.ln_p()).sqrt();
    let value = if r <= 5.0 {
        r -= 1.6;
        polynomial(r, &AS241_C) / polynomial(r, &AS241_D)
    } else {
        r -= 5.0;
        polynomial(r, &AS241_E) / polynomial(r, &AS241_F)
    };
    if q < 0.0 { -value } else { value }
}

// Algorithm AS 241, PPND16 coefficients (Wichura, 1988), ascending powers.
#[allow(clippy::excessive_precision)] // printed with 20 digits in the paper
const AS241_A: [f64; 8] = [
    3.3871328727963666080e0,
    1.3314166789178437745e2,
    1.9715909503065514427e3,
    1.3731693765509461125e4,
    4.5921953931549871457e4,
    6.7265770927008700853e4,
    3.3430575583588128105e4,
    2.5090809287301226727e3,
];
#[allow(clippy::excessive_precision)] // printed with 20 digits in the paper
const AS241_B: [f64; 8] = [
    1.0,
    4.2313330701600911252e1,
    6.8718700749205790830e2,
    5.3941960214247511077e3,
    2.1213794301586595867e4,
    3.9307895800092710610e4,
    2.8729085735721942674e4,
    5.2264952788528545610e3,
];
#[allow(clippy::excessive_precision)] // printed with 20 digits in the paper
const AS241_C: [f64; 8] = [
    1.42343711074968357734e0,
    4.63033784615654529590e0,
    5.76949722146069140550e0,
    3.64784832476320460504e0,
    1.27045825245236838258e0,
    2.41780725177450611770e-1,
    2.27238449892691845833e-2,
    7.74545014278341407640e-4,
];
#[allow(clippy::excessive_precision)] // printed with 20 digits in the paper
const AS241_D: [f64; 8] = [
    1.0,
    2.05319162663775882187e0,
    1.67638483018380384940e0,
    6.89767334985100004550e-1,
    1.48103976427480074590e-1,
    1.51986665636164571966e-2,
    5.47593808499534494600e-4,
    1.05075007164441684324e-9,
];
#[allow(clippy::excessive_precision)] // printed with 20 digits in the paper
const AS241_E: [f64; 8] = [
    6.65790464350110377720e0,
    5.46378491116411436990e0,
    1.78482653991729133580e0,
    2.96560571828504891230e-1,
    2.65321895265761230930e-2,
    1.24266094738807843860e-3,
    2.71155556874348757815e-5,
    2.01033439929228813265e-7,
];
#[allow(clippy::excessive_precision)] // printed with 20 digits in the paper
const AS241_F: [f64; 8] = [
    1.0,
    5.99832206555887937690e-1,
    1.36929880922735805310e-1,
    1.48753612908506148525e-2,
    7.86869131145613259100e-4,
    1.84631831751005468180e-5,
    1.42151175831644588870e-7,
    2.04426310338993978564e-15,
];

/// A gamma distribution with shape and rate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GammaDistribution {
    shape: f64,
    rate: f64,
}

impl GammaDistribution {
    /// A gamma distribution, or `None` unless `shape` and `rate` are finite
    /// and positive.
    pub fn new(shape: f64, rate: f64) -> Option<Self> {
        (shape.is_finite() && shape > 0.0 && rate.is_finite() && rate > 0.0)
            .then_some(Self { shape, rate })
    }

    /// The chi-squared distribution with `freedom` degrees of freedom
    /// (shape `freedom / 2`, rate `1/2`).
    pub fn chi_squared(freedom: f64) -> Option<Self> {
        Self::new(freedom / 2.0, 0.5)
    }

    /// Degrees of freedom when this is a chi-squared distribution
    /// (twice the shape).
    pub fn freedom(&self) -> f64 {
        2.0 * self.shape
    }

    /// Lower cumulative probability `P(X <= x)`: 0 for `x <= 0`, 1 at
    /// `+Inf`, NaN for NaN.
    pub fn cdf(&self, x: f64) -> f64 {
        if x.is_nan() {
            return x;
        }
        if x <= 0.0 {
            return 0.0;
        }
        let scaled = self.rate * x;
        if scaled.is_infinite() {
            return 1.0;
        }
        gamma_tails(self.shape, scaled).0
    }

    /// Upper tail probability `P(X > x)`: 1 for `x <= 0`, 0 at `+Inf`, NaN
    /// for NaN.
    pub fn sf(&self, x: f64) -> f64 {
        if x.is_nan() {
            return x;
        }
        if x <= 0.0 {
            return 1.0;
        }
        let scaled = self.rate * x;
        if scaled.is_infinite() {
            return 0.0;
        }
        gamma_tails(self.shape, scaled).1
    }

    /// Quantile function: the `x` with `P(X <= x) = p`.
    ///
    /// 0 at `p = 0`, `+Inf` at `p = 1`, NaN for NaN. Otherwise Newton's
    /// method on the log of whichever tail is smaller, started from the
    /// better of the Wilson-Hilferty approximation (Wilson and Hilferty, 1931,
    /// PNAS 17:684-688) and the leading terms of the tail expansions, and
    /// safeguarded by a bracket that falls back to bisection.
    ///
    /// # Panics
    ///
    /// If `p` lies outside `[0, 1]`.
    pub fn inverse_cdf(&self, p: f64) -> f64 {
        if p.is_nan() {
            return p;
        }
        assert!(
            (0.0..=1.0).contains(&p),
            "gamma quantile needs 0 <= p <= 1, got {p}"
        );
        if p == 0.0 {
            return 0.0;
        }
        if p == 1.0 {
            return f64::INFINITY;
        }
        standard_gamma_quantile(self.shape, p) / self.rate
    }
}

/// Quantile of the gamma distribution with shape `a` and rate 1, for
/// `0 < p < 1`.
fn standard_gamma_quantile(a: f64, p: f64) -> f64 {
    let lower = p <= 0.5;
    let target = if lower { p } else { 1.0 - p };
    let tail = |y: f64| {
        let (lower_tail, upper_tail) = gamma_tails(a, y);
        if lower { lower_tail } else { upper_tail }
    };
    if lower && tail(f64::MIN_POSITIVE) >= target {
        // The quantile is below the smallest normal double.
        return 0.0;
    }
    let mut y = quantile_start(a, p, lower, target, &tail);
    let (mut low, mut high) = (0.0_f64, f64::INFINITY);
    for _ in 0..400 {
        let current = tail(y);
        if current == target {
            return y;
        }
        // The lower tail increases with y and the upper tail decreases.
        if (current > target) == lower {
            high = y;
        } else {
            low = y;
        }
        // Newton on ln(tail): d ln(tail) / dy = -+ density / tail.
        let log_slope = ln_power_exp_over_gamma(a, y) - y.ln_p() - current.ln_p();
        let step = (current / target).ln_p() / log_slope.exp_p();
        let mut next = if lower { y - step } else { y + step };
        if !(next > low && next < high) {
            next = if high.is_infinite() {
                4.0 * low
            } else if low == 0.0 {
                0.0625 * high
            } else if high > 2.0 * low {
                (low * high).sqrt()
            } else {
                low + 0.5 * (high - low)
            };
            if next <= low || next >= high {
                return y;
            }
        } else if (next - y).abs() <= 2.0 * f64::EPSILON * y {
            return next;
        }
        y = next;
    }
    y
}

/// Starting value for [`standard_gamma_quantile`]: whichever candidate has
/// the tail closest to the target on the log scale.
fn quantile_start(a: f64, p: f64, lower: bool, target: f64, tail: &impl Fn(f64) -> f64) -> f64 {
    let mut candidates = [f64::NAN; 3];
    // Wilson-Hilferty: (X / a)^(1/3) is close to normal.
    let c = 1.0 / (9.0 * a);
    let base = 1.0 - c + standard_normal_inverse_cdf(p) * c.sqrt();
    if base > 0.0 {
        candidates[0] = a * base * base * base;
    }
    if lower {
        // P(a, y) <= y^a / Gamma(a + 1), so this lies just below the quantile.
        candidates[1] = ((target.ln_p() + ln_gamma_1p(a)) / a).exp_p();
    } else {
        // Q(a, y) ~ y^(a-1) e^-y / Gamma(a) for large y.
        let first = -target.ln_p() - ln_gamma(a);
        let first = first.max(1.0);
        candidates[1] = -target.ln_p() - ln_gamma(a) + (a - 1.0) * first.ln_p();
    }
    candidates[2] = a;
    let mut best = a;
    let mut best_distance = f64::INFINITY;
    for &y in &candidates {
        if !(y.is_finite() && y > 0.0) {
            continue;
        }
        let distance = (tail(y) / target).ln_p().abs();
        if distance < best_distance {
            best = y;
            best_distance = distance;
        }
    }
    best
}

#[cfg(test)]
mod tests {
    use super::*;

    // R 4.6.2 values (scripts/references/generate_special_functions_reference.R).
    // Columns: function, first argument, second argument ("-" if none), and the
    // bits of R's result.
    const REFERENCE: &str = include_str!("../../tests/data/special_functions_r462.tsv");

    struct Row {
        a: f64,
        b: f64,
        expected: f64,
    }

    fn rows(function: &str) -> Vec<Row> {
        let rows: Vec<Row> = REFERENCE
            .lines()
            .filter(|line| !line.starts_with('#'))
            .skip(1)
            .map(|line| line.split('\t').collect::<Vec<_>>())
            .filter(|fields| fields[0] == function)
            .map(|fields| Row {
                a: fields[1].parse().unwrap(),
                b: if fields[2] == "-" {
                    f64::NAN
                } else {
                    fields[2].parse().unwrap()
                },
                expected: f64::from_bits(u64::from_str_radix(fields[3], 16).unwrap()),
            })
            .collect();
        assert!(!rows.is_empty(), "no reference rows for {function}");
        rows
    }

    /// Compares every reference row of `function`: NaN, infinite and zero
    /// values must match exactly; otherwise `|actual - expected|` may be at
    /// most `relative * |expected| + absolute`, with `(relative, absolute)`
    /// from `tolerance`, plus a few subnormal units for subnormal values.
    /// Reports the largest relative difference among the normal values that
    /// are held to a purely relative target.
    fn compare(
        function: &str,
        evaluate: impl Fn(f64, f64) -> f64,
        tolerance: impl Fn(&Row) -> (f64, f64),
    ) {
        let mut worst = (0.0_f64, f64::NAN, f64::NAN);
        for row in rows(function) {
            let actual = evaluate(row.a, row.b);
            let expected = row.expected;
            if expected.is_nan() || expected.is_infinite() || expected == 0.0 {
                assert!(
                    actual.to_bits() == expected.to_bits()
                        || (actual.is_nan() && expected.is_nan()),
                    "{function}({}, {}) = {actual:e}, R gives {expected:e}",
                    row.a,
                    row.b
                );
                continue;
            }
            let (relative, mut absolute) = tolerance(&row);
            if expected.abs() < f64::MIN_POSITIVE {
                absolute += 1e-322;
            } else if absolute == 0.0 {
                let difference = ((actual - expected) / expected).abs();
                if difference > worst.0 {
                    worst = (difference, row.a, row.b);
                }
            }
            assert!(
                (actual - expected).abs() <= relative * expected.abs() + absolute,
                "{function}({}, {}) = {actual:e}, R gives {expected:e}",
                row.a,
                row.b
            );
        }
        eprintln!(
            "{function}: largest relative difference from R {:.2e} at ({:e}, {:e})",
            worst.0, worst.1, worst.2
        );
    }

    #[test]
    fn ln_gamma_matches_r() {
        // Within 0.1 of the zeros at 1 and 2, and on the negative axis where
        // the value is small, the target is absolute.
        compare(
            "lgamma",
            |x, _| ln_gamma(x),
            |row| {
                let near_zero = (row.a - 1.0).abs() <= 0.1
                    || (row.a - 2.0).abs() <= 0.1
                    || (row.a < 0.0 && row.expected.abs() < 0.1);
                (4e-15, if near_zero { 1e-15 } else { 0.0 })
            },
        );
    }

    #[test]
    fn gamma_matches_r() {
        // Checked against 40-digit values, R 4.6.2's gamma itself is off by
        // up to about 1.1e-15 * |x| in relative terms for |x| > 10 (1.8e-13
        // at -169.05), while this implementation stays within 2.3e-16 on the
        // whole grid; the comparison allows for R's error there.
        compare(
            "gamma",
            |x, _| gamma(x),
            |row| {
                let relative = if row.a.abs() <= 10.0 {
                    4e-15
                } else {
                    1.5e-15 * row.a.abs()
                };
                (relative, 0.0)
            },
        );
    }

    #[test]
    fn digamma_matches_r() {
        // Near the zero at 1.4616 the target is absolute. For negative x R's
        // reflection loses about |x| (1 + cot^2(pi x)) units of 1e-16 in
        // absolute terms (6e-11 relative at -1e6 - 0.5 against 40-digit
        // values); this implementation reduces pi x exactly and stays within
        // 8e-15 relative there, so the comparison allows for R's error.
        compare(
            "digamma",
            |x, _| digamma(x),
            |row| {
                if row.a < 0.0 {
                    let cot = pi_cot_pi(row.a) / PI;
                    (4e-15, 2e-15 * row.a.abs() * (1.0 + cot * cot))
                } else {
                    let near_root = (row.a - PSI_ROOT_HI).abs() < 0.25;
                    (4e-15, if near_root { 4e-16 } else { 0.0 })
                }
            },
        );
    }

    fn incomplete_gamma_tolerance(row: &Row) -> (f64, f64) {
        if row.expected >= 1e-6 {
            (1e-13, 0.0)
        } else if row.expected >= 1e-300 {
            (1e-12, 0.0)
        } else {
            (1e-12, 1e-310)
        }
    }

    #[test]
    fn lower_incomplete_gamma_matches_r() {
        compare("pgamma_lower", gamma_lr, incomplete_gamma_tolerance);
    }

    #[test]
    fn upper_incomplete_gamma_matches_r() {
        compare("pgamma_upper", gamma_ur, incomplete_gamma_tolerance);
    }

    #[test]
    fn gamma_quantile_matches_r() {
        for rate in ["1", "0.5", "3.75"] {
            let rate_value: f64 = rate.parse().unwrap();
            compare(
                &format!("qgamma_rate_{rate}"),
                |shape, p| {
                    GammaDistribution::new(shape, rate_value)
                        .unwrap()
                        .inverse_cdf(p)
                },
                |_| (1e-12, 0.0),
            );
        }
    }

    #[test]
    fn ln_beta_matches_r() {
        compare("lbeta", ln_beta, |row| {
            (4e-15, if row.expected.abs() < 0.1 { 1e-15 } else { 0.0 })
        });
    }

    #[test]
    fn beta_matches_r() {
        // For a larger argument q >= 10 R's beta inherits the error of its
        // gamma (about 1e-15 * q relative while q < 171.6) or of exp(lbeta)
        // (about 1.3e-16 * |ln B|), up to 8.6e-14 on this grid against
        // 40-digit values. This implementation is within 9e-15 there; the
        // comparison allows for R's error.
        compare("beta", beta, |row| {
            let q = row.a.max(row.b);
            if q < 10.0 {
                return (4e-15, 0.0);
            }
            let gamma_error = if q < GAMMA_OVERFLOW { 1.5e-15 * q } else { 0.0 };
            (4e-15 + gamma_error + 2e-16 * row.expected.ln_p().abs(), 0.0)
        });
    }

    #[test]
    fn normal_density_matches_r() {
        compare("dnorm", |x, _| standard_normal_pdf(x), |_| (4e-15, 0.0));
    }

    #[test]
    fn normal_quantile_matches_r() {
        compare(
            "qnorm",
            |p, _| standard_normal_inverse_cdf(p),
            |_| (4e-15, 0.0),
        );
    }

    #[test]
    fn edge_values_follow_r() {
        assert_eq!(ln_gamma(0.0), f64::INFINITY);
        assert_eq!(ln_gamma(-3.0), f64::INFINITY);
        assert_eq!(ln_gamma(f64::NEG_INFINITY), f64::INFINITY);
        assert!(ln_gamma(f64::NAN).is_nan());
        assert!(gamma(0.0).is_nan() && gamma(-0.0).is_nan() && gamma(-4.0).is_nan());
        assert!(gamma(f64::NEG_INFINITY).is_nan());
        assert_eq!(gamma(f64::INFINITY), f64::INFINITY);
        assert_eq!(gamma(200.0), f64::INFINITY);
        assert!(digamma(0.0).is_nan() && digamma(-2.0).is_nan());
        assert!(digamma(f64::NEG_INFINITY).is_nan());
        assert_eq!(digamma(f64::INFINITY), f64::INFINITY);
        assert_eq!(standard_normal_pdf(f64::INFINITY), 0.0);
        assert!(standard_normal_pdf(f64::NAN).is_nan());
        assert_eq!(standard_normal_inverse_cdf(0.0), f64::NEG_INFINITY);
        assert_eq!(standard_normal_inverse_cdf(1.0), f64::INFINITY);
        assert!(standard_normal_inverse_cdf(f64::NAN).is_nan());
        assert_eq!(ln_beta(f64::INFINITY, 2.0), f64::NEG_INFINITY);
        assert_eq!(beta(3.0, f64::INFINITY), 0.0);
        // Factorials up to 22! are exact doubles.
        let mut factorial = 1.0;
        for n in 1..=23 {
            assert_eq!(gamma(f64::from(n)), factorial, "gamma({n})");
            factorial *= f64::from(n);
        }
    }

    #[test]
    fn incomplete_gamma_boundaries() {
        assert_eq!(gamma_lr(2.5, 0.0), 0.0);
        assert_eq!(gamma_ur(2.5, 0.0), 1.0);
        assert!(gamma_lr(f64::NAN, 1.0).is_nan());
        assert!(gamma_ur(1.0, f64::NAN).is_nan());
        // P + Q = 1 across the methods.
        for &(a, x) in &[
            (0.3, 0.2),
            (0.3, 5.0),
            (4.0, 3.0),
            (4.0, 9.0),
            (50.0, 49.0),
            (50.0, 140.0),
        ] {
            assert!(
                (gamma_lr(a, x) + gamma_ur(a, x) - 1.0).abs() < 1e-15,
                "{a} {x}"
            );
        }
    }

    #[test]
    #[should_panic(expected = "incomplete gamma")]
    fn incomplete_gamma_rejects_zero_shape() {
        gamma_lr(0.0, 1.0);
    }

    #[test]
    #[should_panic(expected = "incomplete gamma")]
    fn incomplete_gamma_rejects_infinite_argument() {
        gamma_ur(1.0, f64::INFINITY);
    }

    #[test]
    #[should_panic(expected = "beta")]
    fn ln_beta_rejects_zero() {
        ln_beta(0.0, 1.0);
    }

    #[test]
    #[should_panic(expected = "beta")]
    fn beta_rejects_negative() {
        beta(2.0, -1.0);
    }

    #[test]
    #[should_panic(expected = "normal quantile")]
    fn normal_quantile_rejects_probability_above_one() {
        standard_normal_inverse_cdf(1.5);
    }

    #[test]
    #[should_panic(expected = "gamma quantile")]
    fn gamma_quantile_rejects_negative_probability() {
        GammaDistribution::new(2.0, 1.0).unwrap().inverse_cdf(-0.1);
    }

    #[test]
    fn gamma_distribution_contract() {
        assert!(GammaDistribution::new(0.0, 1.0).is_none());
        assert!(GammaDistribution::new(1.0, 0.0).is_none());
        assert!(GammaDistribution::new(f64::INFINITY, 1.0).is_none());
        assert!(GammaDistribution::new(1.0, f64::NAN).is_none());
        let chi = GammaDistribution::chi_squared(3.0).unwrap();
        assert_eq!(chi, GammaDistribution::new(1.5, 0.5).unwrap());
        assert_eq!(chi.freedom(), 3.0);
        assert_eq!(chi.cdf(0.0), 0.0);
        assert_eq!(chi.cdf(-1.0), 0.0);
        assert_eq!(chi.sf(-1.0), 1.0);
        assert_eq!(chi.cdf(f64::INFINITY), 1.0);
        assert_eq!(chi.sf(f64::INFINITY), 0.0);
        assert!(chi.cdf(f64::NAN).is_nan());
        assert_eq!(chi.inverse_cdf(0.0), 0.0);
        assert_eq!(chi.inverse_cdf(1.0), f64::INFINITY);
        // cdf(x) = P(shape, rate x) and sf(x) = Q(shape, rate x).
        let d = GammaDistribution::new(2.5, 3.75).unwrap();
        assert_eq!(d.cdf(0.8), gamma_lr(2.5, 3.75 * 0.8));
        assert_eq!(d.sf(0.8), gamma_ur(2.5, 3.75 * 0.8));
        let q = d.inverse_cdf(0.3);
        assert!((d.cdf(q) - 0.3).abs() < 1e-15);
    }
}
