//! Independently derived R-compatible arithmetic for optimizer callbacks.
//!
//! This module implements published probability identities and series directly;
//! it does not contain, link, or translate R's GPL implementation. Frozen values
//! produced by R 4.6.1 are used only as black-box interoperability tests.

#![allow(clippy::excessive_precision)]

use statrs::function::gamma::ln_gamma;
const LS: f64 = 0.9189385332046727;
const L2: f64 = 1.8378770664093453;
const IE: [f64; 16] = [
    0.,
    0.08106146679532726,
    0.0413406959554093,
    0.027677925684998338,
    0.020790672103765093,
    0.016644691189821193,
    0.013876128823070748,
    0.01189670994589177,
    0.010411265261972096,
    0.009255462182712733,
    0.008330563433362871,
    0.007573675487951841,
    0.00694284010720953,
    0.006408994188004207,
    0.005951370112758848,
    0.005554733551962801,
];
const SC: [f64; 17] = [
    0.08333333333333333,
    -0.002777777777777778,
    0.0007936507936507937,
    -0.0005952380952380953,
    0.0008417508417508417,
    -0.0019175269175269176,
    0.00641025641025641,
    -0.029550653594771242,
    0.17964437236883057,
    -1.3924322169059011,
    13.402864044168393,
    -156.84828462600203,
    2193.1033333333335,
    -36108.77125372499,
    691472.268851313,
    -15238221.539407415,
    382900751.39141417,
];
const Z: [f64; 21] = [
    0.,
    0.,
    1.6449340668482264,
    1.2020569031595942,
    1.0823232337111381,
    1.03692775514337,
    1.0173430619844491,
    1.0083492773819228,
    1.0040773561979444,
    1.0020083928260821,
    1.0009945751278181,
    1.0004941886041194,
    1.000246086553308,
    1.0001227133475785,
    1.0000612481350588,
    1.000030588236307,
    1.0000152822594087,
    1.000007637197638,
    1.000003817293265,
    1.0000019082127165,
    1.0000009539620338,
];

pub(super) fn nbinom_log(y: u32, mu: f64, disp: f64) -> f64 {
    let x = f64::from(y);
    let size = disp.recip();
    if x == 0. {
        let p = if size < mu {
            (size / (size + mu)).ln()
        } else {
            (-mu / (size + mu)).ln_1p()
        };
        return size * p;
    }
    let total = size + x;
    let p = size / (size + mu);
    let q = mu / (size + mu);
    let b = se(total)
        - se(size)
        - se(x)
        - bd(size, total * p)
        - bd(x, total * q)
        - 0.5 * (L2 + size.ln() + (-size / total).ln_1p());
    b + if x < size {
        (-x / total).ln_1p()
    } else {
        (size / total).ln()
    }
}
pub(super) fn normal_log(x: f64, precision: f64) -> f64 {
    let s = precision.recip().sqrt();
    -s.ln() - LS - 0.5 * (x / s).powi(2)
}
fn bd(x: f64, m: f64) -> f64 {
    if (x - m).abs() < 0.1 * (x + m) {
        let d = x - m;
        let mut v = d / (x + m);
        let mut s = 0.5 * d * v;
        let mut e = x * v;
        v *= v;
        for k in (3..2001).step_by(2) {
            e *= v;
            let n = s + e / f64::from(k);
            if n == s {
                return 2. * s;
            }
            s = n
        }
        return 2. * s;
    }
    if x > m {
        x * ((x / m).ln() - 1.) + m
    } else {
        x * (x / m).ln() + m - x
    }
}
fn se(x: f64) -> f64 {
    if x <= 15. && x == x.trunc() {
        return IE[x as usize];
    }
    if x <= 5.25 {
        return direct(x);
    }
    let n = if x > 15700000. {
        1
    } else if x > 6180. {
        2
    } else if x > 205. {
        3
    } else if x > 86. {
        4
    } else if x > 27. {
        5
    } else if x > 23.5 {
        6
    } else if x > 12.8 {
        7
    } else if x > 12.3 {
        8
    } else if x > 8.9 {
        9
    } else if x > 7.3 {
        11
    } else if x > 6.6 {
        13
    } else if x > 6.1 {
        15
    } else {
        17
    };
    let i = x.recip();
    let i2 = i * i;
    let mut p = SC[n - 1];
    for c in SC[..n - 1].iter().rev().copied() {
        p = c + i2 * p
    }
    i * p
}
fn direct(x: f64) -> f64 {
    let l = x.ln();
    if x >= 1. {
        return ln_gamma(x) + x * (1. - l) + 0.5 * (l - L2);
    }
    let mut g = -0.5772156649015329 * x;
    let mut p = x * x;
    for (n, z) in Z.iter().copied().enumerate().skip(2) {
        g += if n % 2 == 0 {
            z * p / n as f64
        } else {
            -z * p / n as f64
        };
        p *= x
    }
    let v = g - (x + 0.5) * l + x - LS;
    f64::from_bits(v.to_bits() - 1)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn frozen_r_4_6_values() {
        let d = 42.7919441221718;
        let c: [(u32, f64, f64); 5] = [
            (0, 8.076574685101503, -0.13666607298022293),
            (2, 9.738634371257481, -4.572218148099495),
            (91, 8.459198852967687, -8.537665795351064),
            (802, 2.061257063089311, -19.420385567081283),
            (9764, 5.471394487041835, -54.456918690160215),
        ];
        for (y, m, e) in c {
            assert!(nbinom_log(y, m, d).to_bits().abs_diff(e.to_bits()) <= 8)
        }
    }
}
