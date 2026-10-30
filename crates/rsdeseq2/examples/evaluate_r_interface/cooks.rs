//! Standalone Cook/refit policy conformance, independent of R objects/runtime.
use rsdeseq2::prelude::*;
use serde::Deserialize;
use serde_json::{Value, json};
use std::error::Error;

#[derive(Deserialize)]
pub(super) struct Input {
    cooks_policy: Policy,
}
#[derive(Deserialize)]
struct Factor {
    codes: Vec<i32>,
    levels: usize,
}
#[derive(Deserialize)]
struct Policy {
    mode: String,
    samples: usize,
    pvalue: Vec<Option<String>>,
    maximum: Vec<Option<String>>,
    effect: Vec<String>,
    se: Vec<String>,
    stat: Vec<String>,
    cutoff: Option<String>,
    counts: Vec<u32>,
    cooks: Vec<String>,
    design: Option<Vec<String>>,
    coefficients: usize,
    factors: Vec<Factor>,
    low_count: Option<bool>,
    refitted: Vec<bool>,
    new_zero: Vec<bool>,
}
fn parse(values: Vec<String>) -> Result<Vec<f64>, std::num::ParseFloatError> {
    values.into_iter().map(|v| v.parse()).collect()
}
pub(super) fn evaluate(input: Input) -> Result<Value, Box<dyn Error>> {
    let p = input.cooks_policy;
    let genes = p.pvalue.len();
    let mut effect = parse(p.effect)?;
    let mut se = parse(p.se)?;
    let mut stat = parse(p.stat)?;
    let mut pvalue = p
        .pvalue
        .into_iter()
        .map(|v| {
            v.map(|v| v.parse::<f64>())
                .transpose()
                .map(|v| v.unwrap_or(f64::NAN))
        })
        .collect::<Result<Vec<_>, _>>()?;
    let maximum = p
        .maximum
        .into_iter()
        .map(|v| v.map(|v| v.parse::<f64>()).transpose())
        .collect::<Result<Vec<_>, _>>()?;
    let counts = CountView::from_unsigned(MatrixView::new(
        genes,
        p.samples,
        &p.counts,
        MatrixOrder::RowMajor,
    )?);
    let cooks = parse(p.cooks)?;
    let cooks = MatrixView::new(genes, p.samples, &cooks, MatrixOrder::RowMajor)?;
    let design = p.design.map(parse).transpose()?;
    let mut eligible = false;
    if let Some(values) = &design {
        let design = MatrixView::new(p.samples, p.coefficients, values, MatrixOrder::RowMajor)?;
        for factor in &p.factors {
            if factor.codes.contains(&i32::MIN) {
                continue;
            }
            eligible |= two_level_design_matches(design, &factor.codes, factor.levels)?;
        }
    }
    let heuristic = if p.low_count.unwrap_or(eligible) {
        Some(LowCountHeuristic::new(counts, cooks)?)
    } else {
        None
    };
    let rows = ReplacementRows::new(
        genes,
        Some(RowFlags::from_bools(&p.refitted)),
        Some(RowFlags::from_bools(&p.new_zero)),
    )?;
    rows.apply_columns(&mut effect, &mut se, &mut stat, &mut pvalue)?;
    let cutoff = p.cutoff.map(|v| v.parse::<f64>()).transpose()?;
    let mask = CooksPolicy::new(cutoff, heuristic)?
        .with_replacement(rows)
        .evaluate(
            pvalue
                .iter()
                .map(|&v| if v.is_nan() { None } else { Some(v) }),
            maximum,
        )?;
    let hex = |values: Vec<f64>| {
        values
            .into_iter()
            .map(|v| {
                if v.is_nan() {
                    None
                } else {
                    Some(super::bits(v))
                }
            })
            .collect::<Vec<_>>()
    };
    let mut output = json!({"pvalue":mask.pvalue.into_iter().map(|v|v.map(super::bits)).collect::<Vec<_>>(),"cooksOutlier":mask.outlier});
    if p.mode == "present" {
        output["log2FoldChange"] = json!(hex(effect));
        output["lfcSE"] = json!(hex(se));
        output["stat"] = json!(hex(stat));
    } else if p.mode != "mask" {
        return Err("unknown Cook's conformance mode".into());
    }
    Ok(output)
}
