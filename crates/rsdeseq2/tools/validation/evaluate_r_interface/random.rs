//! Portable state/draw conformance transport, with no R dependency.
use super::bits;
use rsdeseq2::shrinkage::{PosteriorRng, PosteriorUniformKind};
use serde::Deserialize;
use serde_json::{Value, json};
use std::error::Error;
#[derive(Deserialize)]
pub struct Input {
    pub random_cases: Vec<Case>,
}
#[derive(Deserialize)]
pub struct Case {
    pub name: String,
    pub seed: Option<i32>,
    pub kind: Option<String>,
    pub state: Vec<String>,
}
pub fn state(v: &[String]) -> Result<Vec<i32>, Box<dyn Error>> {
    v.iter()
        .map(|v| Ok(u32::from_str_radix(v, 16)? as i32))
        .collect()
}
pub fn encode_state(v: &[i32]) -> Vec<String> {
    v.iter().map(|v| format!("{:08x}", *v as u32)).collect()
}
pub fn generator(c: &Case) -> Result<PosteriorRng, Box<dyn Error>> {
    Ok(if let Some(seed) = c.seed {
        PosteriorRng::from_seed_with_kind(
            seed,
            match c.kind.as_deref() {
                None | Some("Mersenne-Twister") => PosteriorUniformKind::MersenneTwister,
                Some("L'Ecuyer-CMRG") => PosteriorUniformKind::LecuyerCmrg,
                _ => return Err("unsupported seeded RNG kind".into()),
            },
        )?
    } else {
        PosteriorRng::from_r_seed(&state(&c.state)?)?
    })
}
pub fn evaluate(input: Input) -> Result<Value, Box<dyn Error>> {
    let cases = input
        .random_cases
        .into_iter()
        .map(|c| {
            let original = generator(&c)?;
            let initial = encode_state(&original.r_seed());
            let mut u = original.clone();
            let mut n = original;
            let uniforms = (0..64).map(|_| bits(u.uniform())).collect::<Vec<_>>();
            let normals = (0..64)
                .map(|_| bits(n.start_coordinate()))
                .collect::<Vec<_>>();
            Ok(
                json!({"name":c.name,"initial":initial,"uniforms":uniforms,"normals":normals,
            "uniform_state":encode_state(&u.r_seed()),"normal_state":encode_state(&n.r_seed())}),
            )
        })
        .collect::<Result<Vec<Value>, Box<dyn Error>>>()?;
    Ok(json!({"cases":cases}))
}
