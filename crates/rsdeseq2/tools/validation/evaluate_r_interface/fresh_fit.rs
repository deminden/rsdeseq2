//! Full fresh-fit transport for independent native/wrapped bit comparisons.
use super::{bits, heavy_tail::Matrix, parse_numeric_inputs as numbers};
use rsdeseq2::{cooks::CooksRefitPlan, prelude::*};
use serde::Deserialize;
use serde_json::{Value, json};
use std::error::Error;
#[path = "fresh_fit/formula_case.rs"]
mod formula_case;
#[path = "fresh_fit/queries.rs"]
mod queries;

#[derive(Deserialize)]
pub struct Input {
    pub fresh_fit_cases: Vec<Case>,
}
#[derive(Deserialize)]
pub struct Case {
    pub name: String,
    pub counts: Matrix,
    pub design: Option<Matrix>,
    pub formula: Option<formula_case::FormulaCase>,
    pub reduced: Option<Matrix>,
    pub size_factors: Option<Vec<String>>,
    pub normalization_factors: Option<Matrix>,
    pub weights: Option<Matrix>,
    pub fit_type: String,
    pub test: String,
    pub threads: usize,
    pub min_replicates: usize,
    pub replacement_cutoff: String,
    #[serde(default)]
    pub queries: Vec<queries::Query>,
}
fn matrix(m: &Matrix) -> Result<RowMajorMatrix<f64>, Box<dyn Error>> {
    let values = numbers(&m.values)?;
    let view = m.view(&values)?;
    Ok(RowMajorMatrix::from_row_major(
        m.rows,
        m.columns,
        (0..m.rows)
            .flat_map(|i| (0..m.columns).map(move |j| *view.get(i, j).unwrap()))
            .collect(),
    )?)
}
fn floats(v: &[f64]) -> Vec<String> {
    v.iter().copied().map(bits).collect()
}
fn optional(v: Option<&[f64]>, n: usize) -> Vec<String> {
    v.map(floats).unwrap_or_else(|| vec![bits(f64::NAN); n])
}
pub fn run(c: &Case) -> Result<Value, Box<dyn Error>> {
    let raw = matrix(&c.counts)?;
    rsdeseq2::matrix::CountView::from_doubles((&raw).into())?;
    let counts = CountMatrix::from_row_major_u32(
        raw.n_rows(),
        raw.n_cols(),
        raw.as_slice().iter().map(|v| *v as u32).collect(),
    )?;
    let mut builder = DeseqBuilder::new()
        .automatic_dispersion_means()
        .threads(c.threads)
        .fit_type(match c.fit_type.as_str() {
            "mean" => FitType::Mean,
            "local" => FitType::Local,
            "parametric" => FitType::Parametric,
            _ => return Err("unknown fit type".into()),
        })
        .test(match c.test.as_str() {
            "Wald" => TestType::Wald,
            "LRT" => TestType::Lrt,
            _ => return Err("unknown test".into()),
        });
    // Match the public wrapper's argument precedence before constructing the
    // lower-level builder, which can retain unused sample-scale metadata.
    if c.normalization_factors.is_none()
        && let Some(v) = &c.size_factors
    {
        builder = builder.size_factors(numbers(v)?);
    }
    if let Some(v) = &c.normalization_factors {
        builder = builder.normalization_factors(matrix(v)?);
    }
    if let Some(v) = &c.weights {
        builder = builder.observation_weights(matrix(v)?);
    }
    if let Some(v) = &c.reduced {
        let m = matrix(v)?;
        builder = builder.reduced_design(DesignMatrix::from_row_major(
            m.n_rows(),
            m.n_cols(),
            m.into_values(),
            None,
        )?);
    }
    let (fit, results, plan) = if let Some(formula) = &c.formula {
        formula.fit(builder, &counts, c)?
    } else {
        let dm = matrix(
            c.design
                .as_ref()
                .ok_or("matrix design or formula is required")?,
        )?;
        let design =
            DesignMatrix::from_row_major(dm.n_rows(), dm.n_cols(), dm.into_values(), None)?;
        if c.min_replicates == 0 {
            let (f, r) = builder.fit_with_results(&counts, &design)?;
            (f, r, None)
        } else {
            let mut options =
                CooksReplacementOptions::new(super::parse_numeric_input(&c.replacement_cutoff)?);
            options.min_replicates = c.min_replicates;
            let (f, r, p) = builder
                .fit_with_test_results_with_cooks_replacement(&counts, &design, &options)?
                .into_merged_fit()?;
            (f, r, Some(p))
        }
    };
    let queries = queries::evaluate(&fit, &counts, plan.as_ref(), &c.queries)?;
    let mut output = encode(fit, results, plan)?;
    if !c.queries.is_empty() {
        output["queries"] = queries;
    }
    Ok(output)
}
fn encode(
    fit: DeseqFit,
    results: DeseqResults,
    plan: Option<CooksRefitPlan>,
) -> Result<Value, Box<dyn Error>> {
    let g = fit.base_mean.len();
    let rows = &results.rows;
    let numeric = |field: fn(&rsdeseq2::results::DeseqResultRow) -> Option<f64>| {
        rows.iter()
            .map(|r| bits(field(r).unwrap_or(f64::NAN)))
            .collect::<Vec<_>>()
    };
    let flags = |field: fn(&rsdeseq2::results::DeseqResultRow) -> Option<bool>| {
        rows.iter().map(field).collect::<Vec<_>>()
    };
    let raw = (0..g)
        .map(|i| {
            if plan
                .as_ref()
                .is_some_and(|p| p.new_all_zero_rows.contains(&i))
            {
                return (0.0, 1.0);
            }
            if fit.all_zero[i] {
                return (f64::NAN, f64::NAN);
            }
            if let Some(lrt) = &fit.lrt {
                (
                    lrt.deviance[i].unwrap_or(f64::NAN),
                    lrt.pvalue[i].unwrap_or(f64::NAN),
                )
            } else {
                let wald = fit.wald.as_ref().unwrap();
                (
                    wald.stat[i].unwrap_or(f64::NAN),
                    wald.pvalue[i].unwrap_or(f64::NAN),
                )
            }
        })
        .collect::<Vec<_>>();
    let eligible = plan
        .as_ref()
        .is_some_and(CooksRefitPlan::has_replaceable_samples);
    let n_refit = plan.as_ref().filter(|_| eligible).map_or(0, |p| p.n_refit);
    let replacement = plan.as_ref().filter(|_| n_refit > 0).map(|p| {
        p.replacement
            .replaced_counts
            .as_slice()
            .iter()
            .map(|v| bits(f64::from(*v)))
            .collect::<Vec<_>>()
    });
    let matrix = |m: Option<&RowMajorMatrix<f64>>| m.map(|m| floats(m.as_slice()));
    Ok(json!({
        "results":{"baseMean":rows.iter().map(|r|bits(r.base_mean)).collect::<Vec<_>>(),
            "log2FoldChange":numeric(|r|r.log2_fold_change),"lfcSE":numeric(|r|r.lfc_se),
            "stat":numeric(|r|r.stat),"pvalue":numeric(|r|r.pvalue),"padj":numeric(|r|r.padj),
            "dispersion":rows.iter().enumerate().map(|(i,r)|bits(r.dispersion.unwrap_or_else(||fit.dispersion.as_ref().map_or(f64::NAN,|v|v[i])))).collect::<Vec<_>>(),
            "converged":flags(|r|r.converged),"maxCooks":numeric(|r|r.max_cooks),"cooksOutlier":flags(|r|r.cooks_outlier),"filtered":flags(|r|r.filtered)},
        "rawStat":raw.iter().map(|r|bits(r.0)).collect::<Vec<_>>(),"rawPvalue":raw.iter().map(|r|bits(r.1)).collect::<Vec<_>>(),
        "baseMean":floats(&fit.base_mean),"baseVar":floats(&fit.base_var),"allZero":fit.all_zero,
        "sizeFactors":floats(&fit.size_factors),"dispersions":optional(fit.dispersion.as_deref(),g),"dispFit":optional(fit.disp_fit.as_deref(),g),
        "beta":matrix(fit.beta.as_ref()),"betaSE":matrix(fit.beta_se.as_ref()),"betaCovariance":matrix(fit.beta_covariance.as_ref()),
        "normalizationFactors":matrix(fit.normalization_factors.as_ref()),
        "design":fit.design.as_ref().map(|d|floats(d.matrix().as_slice())),
        "reduced":fit.reduced_design.as_ref().map(|d|floats(d.matrix().as_slice())),
        "mu":matrix(fit.mu.as_ref()),"cooks":matrix(fit.cooks.as_ref()),"lrtDf":fit.lrt.as_ref().map(|v|v.degrees_of_freedom),
        "nRefit":n_refit,"replaceCounts":replacement,
        "replace":(0..g).map(|i|plan.as_ref().filter(|_|eligible).and_then(|p|p.replacement.replace[i])).collect::<Vec<_>>(),
        "replaceable":plan.as_ref().map(|p|p.replacement.replaceable_samples.clone()).unwrap_or_else(||vec![false;fit.counts_summary.n_samples]),
        "refitRows":(0..g).map(|i|eligible && plan.as_ref().is_some_and(|p|p.refit_rows.contains(&i))).collect::<Vec<_>>(),
        "newAllZeroRows":(0..g).map(|i|eligible && plan.as_ref().is_some_and(|p|p.new_all_zero_rows.contains(&i))).collect::<Vec<_>>()
    }))
}
pub fn evaluate(input: Input) -> Value {
    json!({"cases":input.fresh_fit_cases.iter().map(|c|match run(c){Ok(v)=>json!({"name":c.name,"value":v}),Err(e)=>json!({"name":c.name,"error":e.to_string()})}).collect::<Vec<_>>()})
}
