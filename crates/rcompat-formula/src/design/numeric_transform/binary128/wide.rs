//! Round-to-nearest, ties-to-even arithmetic at binary128 precision.
//!
//! The inputs are binary64 values and slice lengths. Their sums, differences,
//! and means stay inside the binary128 exponent range, so finite nonzero
//! values can always use a normalized integer significand.

const PRECISION: u32 = 113;
const GUARD: u32 = 3;
const WORK_BITS: u32 = PRECISION + GUARD;
const FRACTION_MASK: u64 = (1 << 52) - 1;
const QUIET_BIT: u64 = 1 << 51;
const SIGN_BIT: u64 = 1 << 63;

#[derive(Clone, Copy)]
pub(super) enum Wide {
    Finite {
        significand: u128,
        exponent: i32,
        negative: bool,
    },
    Special(f64),
}

impl Wide {
    pub(super) const ZERO: Self = Self::Special(0.0);

    pub(super) fn from_double(value: f64) -> Self {
        let bits = value.to_bits();
        if !value.is_finite() {
            return Self::Special(if value.is_nan() {
                f64::from_bits(bits | QUIET_BIT)
            } else {
                value
            });
        }
        if value == 0.0 {
            return Self::Special(value);
        }
        let encoded_exponent = ((bits >> 52) & 0x7ff) as i32;
        let fraction = bits & FRACTION_MASK;
        let (integer, exponent) = if encoded_exponent == 0 {
            (fraction, -1074)
        } else {
            (fraction | (1 << 52), encoded_exponent - 1075)
        };
        Self::round(u128::from(integer), exponent, bits & SIGN_BIT != 0)
    }

    pub(super) fn to_double(self) -> f64 {
        let (significand, exponent, negative) = match self {
            Self::Finite {
                significand,
                exponent,
                negative,
            } => (significand, exponent, negative),
            Self::Special(value) => return value,
        };
        let sign = if negative { SIGN_BIT } else { 0 };
        let leading_exponent = exponent + PRECISION as i32 - 1;
        if leading_exponent < -1022 {
            let mantissa = rounded_shift(significand, (-1074 - exponent) as u32);
            return f64::from_bits(sign | mantissa as u64);
        }
        let mantissa = rounded_shift(significand, PRECISION - 53) as u64;
        // Add the full mantissa so a rounding carry increments the exponent.
        let encoded = if leading_exponent > 1023 {
            0x7ff << 52
        } else {
            (((leading_exponent + 1022) as u64) << 52) + mantissa
        };
        f64::from_bits(sign | encoded)
    }

    pub(super) fn add(self, rhs: Self) -> Self {
        match (self, rhs) {
            (Self::Special(a), Self::Special(b)) => {
                let value = if a.is_nan() {
                    a
                } else if b.is_nan() {
                    b
                } else if a.is_infinite() && b == -a {
                    f64::from_bits(0x7fff_ffff_ffff_ffff)
                } else {
                    a + b
                };
                Self::Special(value)
            }
            (Self::Special(x), finite) | (finite, Self::Special(x)) => {
                if x == 0.0 {
                    finite
                } else {
                    Self::Special(x)
                }
            }
            (
                Self::Finite {
                    significand: a,
                    exponent: ea,
                    negative: na,
                },
                Self::Finite {
                    significand: b,
                    exponent: eb,
                    negative: nb,
                },
            ) => {
                let exponent = ea.max(eb);
                let a = sticky_shift(a << GUARD, (exponent - ea) as u32);
                let b = sticky_shift(b << GUARD, (exponent - eb) as u32);
                let (magnitude, negative) = if na == nb {
                    (a + b, na)
                } else if a >= b {
                    (a - b, na)
                } else {
                    (b - a, nb)
                };
                Self::round(magnitude, exponent - GUARD as i32, negative)
            }
        }
    }

    pub(super) fn subtract(self, rhs: Self) -> Self {
        self.add(match rhs {
            Self::Special(value) => Self::Special(-value),
            Self::Finite {
                significand,
                exponent,
                negative,
            } => Self::Finite {
                significand,
                exponent,
                negative: !negative,
            },
        })
    }

    pub(super) fn divide(self, count: usize) -> Self {
        let (significand, exponent, negative) = match self {
            // Preserve NaN payloads and use a positive quiet NaN for 0 / 0.
            Self::Special(value) if value.is_nan() => return self,
            Self::Special(value) if value == 0.0 && count == 0 => return Self::Special(f64::NAN),
            Self::Special(value) => return Self::Special(value / count as f64),
            Self::Finite {
                significand,
                exponent,
                negative,
            } => (significand, exponent, negative),
        };
        if count == 0 {
            return Self::Special(if negative {
                f64::NEG_INFINITY
            } else {
                f64::INFINITY
            });
        }
        let divisor = count as u128;
        let mut quotient = significand / divisor;
        let mut remainder = significand % divisor;
        let shift = WORK_BITS - bit_length(quotient);
        // Each remainder is smaller than a slice length. Shifting it at most
        // 64 bits leaves enough space in u128 even on 64-bit targets.
        for bits in [shift.min(64), shift.saturating_sub(64)] {
            remainder <<= bits;
            quotient = (quotient << bits) | (remainder / divisor);
            remainder %= divisor;
        }
        quotient |= u128::from(remainder != 0);
        Self::round(quotient, exponent - shift as i32, negative)
    }

    pub(super) fn exceeds_double(self) -> bool {
        let Self::Finite {
            significand,
            exponent,
            ..
        } = self
        else {
            return false;
        };
        // The largest binary64 value is a 53-bit all-ones significand
        // multiplied by 2^(1023 - 52), then normalized to this precision.
        let max_significand = ((1 << 53) - 1) << (PRECISION - 53);
        let max_exponent = 1024 - PRECISION as i32;
        (exponent, significand) > (max_exponent, max_significand)
    }

    fn round(magnitude: u128, mut exponent: i32, negative: bool) -> Self {
        if magnitude == 0 {
            return Self::ZERO;
        }
        let adjustment = bit_length(magnitude) as i32 - PRECISION as i32;
        exponent += adjustment;
        let mut significand = if adjustment < 0 {
            magnitude << -adjustment
        } else {
            rounded_shift(magnitude, adjustment as u32)
        };
        if significand == 1 << PRECISION {
            significand >>= 1;
            exponent += 1;
        }
        Self::Finite {
            significand,
            exponent,
            negative,
        }
    }
}

fn bit_length(value: u128) -> u32 {
    u128::BITS - value.leading_zeros()
}

fn sticky_shift(value: u128, shift: u32) -> u128 {
    if shift >= u128::BITS {
        return u128::from(value != 0);
    }
    let retained = value >> shift;
    retained | u128::from(value != retained << shift)
}

fn rounded_shift(value: u128, shift: u32) -> u128 {
    match shift {
        0 => value,
        1..u128::BITS => {
            let retained = value >> shift;
            let remainder = value - (retained << shift);
            let halfway = 1 << (shift - 1);
            retained + u128::from(remainder > halfway || remainder == halfway && retained & 1 != 0)
        }
        u128::BITS => u128::from(value > 1 << (u128::BITS - 1)),
        _ => 0,
    }
}
