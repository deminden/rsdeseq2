//! Portable random-start state, using published MT19937/CMRG equations.
//! R seed scrambling/state layout and two-uniform normal construction are
//! interoperability data. No R/GPL generator implementation is copied or linked.
use super::prior::invalid;
use crate::{DeseqError, probability::normal_upper_quantile};

use std::{
    sync::atomic::{AtomicU64, Ordering},
    time::{SystemTime, UNIX_EPOCH},
};
const MT_N: usize = 624;
const M1: i64 = 4_294_967_087;
const M2: i64 = 4_294_944_443;
#[derive(Clone, Debug)]
enum Uniform {
    Mt {
        words: Box<[u32; MT_N]>,
        position: usize,
    },
    Cmrg([i64; 6]),
}
/// Portable uniform families available with Inversion normal starts.
#[derive(Clone, Copy, Debug)]
pub enum PosteriorUniformKind {
    MersenneTwister,
    LecuyerCmrg,
}
/// Serializable native RNG for the supported R Inversion-normal states.
/// Uniforms and exported states are exact; normal quantiles use the same
/// standalone probability kernel as other inference, with no R math calls.
#[derive(Clone, Debug)]
pub struct PosteriorRng {
    header: i32,
    uniform: Uniform,
    draws: usize,
}
impl PosteriorRng {
    /// Reproduce the documented R default set.seed() initialization.
    pub fn from_seed(seed: i32) -> Result<Self, DeseqError> {
        Self::from_seed_with_kind(seed, PosteriorUniformKind::MersenneTwister)
    }
    /// Seed either portable family without requiring R or an imported state.
    pub fn from_seed_with_kind(seed: i32, kind: PosteriorUniformKind) -> Result<Self, DeseqError> {
        if seed == i32::MIN {
            return Err(invalid("random seed must not be NA_integer_"));
        }
        let mut state = seed as u32;
        let step = |s: &mut u32| {
            *s = s.wrapping_mul(69069).wrapping_add(1);
            *s
        };
        for _ in 0..50 {
            step(&mut state);
        }
        match kind {
            PosteriorUniformKind::MersenneTwister => {
                step(&mut state); // The position slot is then fixed to 624.
                let mut words = Box::new([0; MT_N]);
                for word in words.iter_mut() {
                    *word = step(&mut state);
                }
                Ok(Self {
                    header: 10403,
                    uniform: Uniform::Mt {
                        words,
                        position: MT_N,
                    },
                    draws: 0,
                })
            }
            PosteriorUniformKind::LecuyerCmrg => {
                let mut words = [0; 6];
                for word in &mut words {
                    loop {
                        let candidate = step(&mut state);
                        if i64::from(candidate) < M2 {
                            *word = i64::from(candidate);
                            break;
                        }
                    }
                }
                Ok(Self {
                    header: 10407,
                    uniform: Uniform::Cmrg(words),
                    draws: 0,
                })
            }
        }
    }
    /// Restore MT19937 or L'Ecuyer-CMRG with Inversion normals. All signed
    /// state-word bits, including the integer-NA bit pattern, are preserved.
    pub fn from_r_seed(seed: &[i32]) -> Result<Self, DeseqError> {
        let header = *seed
            .first()
            .ok_or_else(|| invalid("missing random-state header"))?;
        if header < 0 || header / 10000 > 1 || header / 100 % 100 != 4 {
            return Err(invalid(
                "native random starts require an Inversion-normal RNG state",
            ));
        }
        let uniform = match header % 100 {
            3 => {
                if seed.len() != MT_N + 2 || !(0..=MT_N as i32).contains(&seed[1]) {
                    return Err(invalid(
                        "invalid Mersenne-Twister state dimensions/position",
                    ));
                }
                let mut words = Box::new([0; MT_N]);
                for (word, &value) in words.iter_mut().zip(&seed[2..]) {
                    *word = value as u32;
                }
                if words.iter().all(|&v| v == 0) {
                    return Err(invalid("Mersenne-Twister state must not be all zero"));
                }
                Uniform::Mt {
                    words,
                    position: if seed[1] == 0 { MT_N } else { seed[1] as usize },
                }
            }
            7 => {
                if seed.len() != 7 {
                    return Err(invalid("invalid L'Ecuyer-CMRG state dimensions"));
                }
                let values = std::array::from_fn(|i| i64::from(seed[i + 1] as u32));
                if values[..3].iter().any(|&v| v >= M1)
                    || values[3..].iter().any(|&v| v >= M2)
                    || values[..3].iter().all(|&v| v == 0)
                    || values[3..].iter().all(|&v| v == 0)
                {
                    return Err(invalid("invalid L'Ecuyer-CMRG state words"));
                }
                Uniform::Cmrg(values)
            }
            _ => {
                return Err(invalid(
                    "native random starts support Mersenne-Twister and L'Ecuyer-CMRG",
                ));
            }
        };
        Ok(Self {
            header,
            uniform,
            draws: 0,
        })
    }
    /// Initialize an unseeded request entirely in Rust. Reproducible callers
    /// should supply a seed/state; an unseeded R session receives the result.
    pub fn from_entropy() -> Self {
        static COUNTER: AtomicU64 = AtomicU64::new(0);
        let time = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .map(|v| v.as_nanos())
            .unwrap_or(0);
        let mix = time as u64
            ^ (time >> 64) as u64
            ^ u64::from(std::process::id())
            ^ COUNTER
                .fetch_add(1, Ordering::Relaxed)
                .wrapping_mul(0x9e3779b97f4a7c15);
        let seed = (mix ^ mix.rotate_left(29)) as i32;
        Self::from_seed(if seed == i32::MIN { 0 } else { seed }).unwrap()
    }
    pub fn r_seed(&self) -> Vec<i32> {
        let mut out = vec![self.header];
        match &self.uniform {
            Uniform::Mt { words, position } => {
                out.push(*position as i32);
                out.extend(words.iter().map(|&v| v as i32));
            }
            Uniform::Cmrg(v) => out.extend(v.iter().map(|&v| v as u32 as i32)),
        }
        out
    }
    pub fn has_drawn(&self) -> bool {
        self.draws != 0
    }
    pub fn uniform(&mut self) -> f64 {
        self.draws = self.draws.saturating_add(1);
        match &mut self.uniform {
            Uniform::Mt { words, position } => {
                if *position == MT_N {
                    for i in 0..MT_N {
                        let y = (words[i] & 0x80000000) | (words[(i + 1) % MT_N] & 0x7fffffff);
                        words[i] = words[(i + 397) % MT_N]
                            ^ (y >> 1)
                            ^ if y & 1 == 0 { 0 } else { 0x9908b0df };
                    }
                    *position = 0;
                }
                let mut y = words[*position];
                *position += 1;
                y ^= y >> 11;
                y ^= (y << 7) & 0x9d2c5680;
                y ^= (y << 15) & 0xefc60000;
                y ^= y >> 18;
                let p = f64::from(y) * (1.0 / 4294967296.0);
                if p == 0.0 {
                    0.5 * 2.328306437080797e-10
                } else {
                    p
                }
            }
            Uniform::Cmrg(v) => {
                let first = (1403580 * v[1] - 810728 * v[0]).rem_euclid(M1);
                let second = (527612 * v[5] - 1370589 * v[3]).rem_euclid(M2);
                *v = [v[1], v[2], first, v[4], v[5], second];
                let delta = first - second;
                let positive = if delta <= 0 { delta + M1 } else { delta };
                positive as f64 * (1.0 / 4294967088.0)
            }
        }
    }
    /// One N(0, 0.5^2) initialization coordinate, consuming two uniforms.
    pub fn start_coordinate(&mut self) -> f64 {
        let high = (134217728.0 * self.uniform()).floor();
        let p = (high + self.uniform()) / 134217728.0;
        if p <= 0.0 {
            return f64::NEG_INFINITY;
        }
        if p >= 1.0 {
            return f64::INFINITY;
        }
        let tail = p.min(1.0 - p);
        let q = normal_upper_quantile(tail);
        (if p < 0.5 { -q } else { q }) * 0.5
    }
    pub fn start_vector(&mut self, coefficients: usize) -> Vec<f64> {
        (0..coefficients).map(|_| self.start_coordinate()).collect()
    }
}
