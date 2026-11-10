//! Integer state transport only; generation and probability calculations are core.
use super::{Data, Value};
use rsdeseq2::shrinkage::PosteriorRng;
pub(super) fn restore(value: Option<&Value>) -> Result<PosteriorRng, String> {
    let Some(value) = value.filter(|v| !matches!(v.data, Data::Null)) else {
        return Ok(PosteriorRng::from_entropy());
    };
    let Data::Number {
        values,
        integer: true,
    } = &value.data
    else {
        return Err("random state requires integer storage".into());
    };
    if value.dimensions.is_some() {
        return Err("random state must be a vector".into());
    }
    let words = values
        .iter()
        .map(|&v| {
            // The signed-INT_MIN bit pattern is a legitimate RNG word. R's
            // general integer wire represents that bit pattern as numeric NA.
            if v.is_nan() {
                Ok(i32::MIN)
            } else if v.is_finite()
                && v > f64::from(i32::MIN)
                && v <= f64::from(i32::MAX)
                && v.fract() == 0.0
            {
                Ok(v as i32)
            } else {
                Err("invalid integer random-state word")
            }
        })
        .collect::<Result<Vec<_>, _>>()?;
    PosteriorRng::from_r_seed(&words).map_err(|e| e.to_string())
}
pub(super) fn encode(state: Option<&[i32]>) -> Value {
    state
        .map(|state| {
            Value::new(Data::Number {
                integer: true,
                values: state
                    .iter()
                    .map(|&v| {
                        if v == i32::MIN {
                            f64::from_bits(0x7ff0_0000_0000_07a2)
                        } else {
                            f64::from(v)
                        }
                    })
                    .collect(),
            })
        })
        .unwrap_or_else(|| Value::new(Data::Null))
}
