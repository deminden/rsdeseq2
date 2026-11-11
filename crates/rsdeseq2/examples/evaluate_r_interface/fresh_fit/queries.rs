//! Later requests reuse one freshly fitted model; all arithmetic is in the core.
use super::{bits, floats, numbers};
use rsdeseq2::{cooks::CooksRefitPlan, prelude::*};
use serde::Deserialize;
use serde_json::{Value, json};
use std::error::Error;

#[derive(Deserialize)]
pub struct Query {
    name: String,
    #[serde(default)]
    contrast: Vec<String>,
    character: Option<[String; 3]>,
    threshold: String,
    alternative: String,
    df: Option<String>,
    cooks: String,
    filtering: bool,
    alpha: String,
    missing: String,
}

pub fn evaluate(
    fit: &DeseqFit,
    counts: &CountMatrix,
    plan: Option<&CooksRefitPlan>,
    queries: &[Query],
) -> Result<Value, Box<dyn Error>> {
    let mut output = Vec::with_capacity(queries.len());
    for q in queries {
        let missing = super::super::parse_numeric_input(&q.missing)?;
        if !missing.is_nan() {
            return Err("query transport missing marker must be NaN".into());
        }
        let alternative = match q.alternative.as_str() {
            "greaterAbs" => WaldAlternative::GreaterAbs,
            "lessAbs" => WaldAlternative::LessAbs,
            "greater" => WaldAlternative::Greater,
            "less" => WaldAlternative::Less,
            _ => return Err("unknown query alternative".into()),
        };
        let wald = if let Some(df) = &q.df {
            WaldTestOptions::t_degrees_of_freedom(super::super::parse_numeric_input(df)?)
        } else {
            WaldTestOptions::normal()
        }
        .with_lfc_threshold(
            super::super::parse_numeric_input(&q.threshold)?,
            alternative,
        );
        let cooks = match q.cooks.as_str() {
            "auto" => CooksCutoff::Default,
            "off" => CooksCutoff::Disabled,
            v => CooksCutoff::Threshold(super::super::parse_numeric_input(v)?),
        };
        let contrast = if let Some([factor, numerator, denominator]) = &q.character {
            ResultsContrast::Character {
                factor: factor.clone(),
                numerator: numerator.clone(),
                denominator: denominator.clone(),
                reference: None,
            }
        } else {
            ResultsContrast::Numeric(numbers(&q.contrast)?)
        };
        let request = ResultRequest::new(
            contrast,
            wald,
            cooks,
            IndependentFilteringOptions {
                enabled: q.filtering,
                alpha: super::super::parse_numeric_input(&q.alpha)?,
                theta: None,
            },
        )?;
        let result = if let Some(plan) = plan {
            fit.results_for_request_with_replacement(counts, &request, None, plan)?
        } else {
            fit.results_for_request(counts, &request, None)?
        };
        let column = |field: fn(&rsdeseq2::results::DeseqResultRow) -> Option<f64>, absent: f64| {
            result
                .rows
                .iter()
                .map(|r| bits(field(r).unwrap_or(absent)))
                .collect::<Vec<_>>()
        };
        let m = result
            .independent_filtering
            .as_ref()
            .ok_or("missing filtering metadata")?;
        // Match the public ABI's explicit missing marker. Stored LRT tests and
        // active Cook's masking already use the core's canonical NaN marker.
        let absent_stat = if fit.lrt.is_some() { f64::NAN } else { missing };
        let absent_p = if q.cooks != "off" || fit.lrt.is_some() {
            f64::NAN
        } else {
            missing
        };
        let padj = result
            .rows
            .iter()
            .map(|r| {
                bits(rsdeseq2::independent_filtering::format_adjusted_pvalue(
                    r.padj,
                    r.pvalue.unwrap_or(absent_p),
                    r.base_mean,
                    m.filter_threshold,
                    missing,
                ))
            })
            .collect::<Vec<_>>();
        output.push(json!({"name":q.name,"value":{
            "baseMean":result.rows.iter().map(|r|bits(r.base_mean)).collect::<Vec<_>>(),
            "log2FoldChange":column(|r|r.log2_fold_change,missing),"lfcSE":column(|r|r.lfc_se,missing),
            "stat":column(|r|r.stat,absent_stat),"pvalue":column(|r|r.pvalue,absent_p),"padj":padj,
            "filtered":result.rows.iter().map(|r|r.filtered).collect::<Vec<_>>(),
            "filtering":{"enabled":m.enabled,"theta":floats(&m.theta),"numRej":m.num_rejections,
                "selected":m.selected_index.map(|v|v+1),"filterTheta":m.filter_theta.map(bits),
                "filterThreshold":m.filter_threshold.map(bits),"lowessX":m.lowess_fit.as_ref().map(|_|floats(&m.theta)),
                "lowessY":m.lowess_fit.as_ref().map(|v|floats(v)),"alpha":bits(m.alpha)}
        }}));
    }
    Ok(Value::Array(output))
}
