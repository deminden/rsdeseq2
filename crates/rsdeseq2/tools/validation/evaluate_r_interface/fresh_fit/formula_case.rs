//! Public formula workflows, using only primitive model-frame input.
use super::{Case, numbers};
use rsdeseq2::{cooks::CooksRefitPlan, prelude::*};
use serde::Deserialize;
use std::error::Error;

#[derive(Deserialize)]
pub struct FormulaCase {
    full: String,
    reduced: Option<String>,
    factors: Vec<Factor>,
    numeric: Vec<Numeric>,
}
#[derive(Deserialize)]
struct Factor {
    name: String,
    sample_levels: Vec<String>,
    levels: Option<Vec<String>>,
    reference: Option<String>,
}
#[derive(Deserialize)]
struct Numeric {
    name: String,
    values: Vec<String>,
}
type Output = (DeseqFit, DeseqResults, Option<CooksRefitPlan>);
impl FormulaCase {
    pub fn fit(
        &self,
        builder: DeseqBuilder,
        counts: &CountMatrix,
        case: &Case,
    ) -> Result<Output, Box<dyn Error>> {
        let frame = FormulaModelFrame {
            factors: self
                .factors
                .iter()
                .map(|f| FormulaFactorColumn {
                    name: f.name.clone(),
                    sample_levels: f.sample_levels.clone(),
                    levels: f.levels.clone(),
                    reference: f.reference.clone(),
                })
                .collect(),
            numeric_covariates: self
                .numeric
                .iter()
                .map(|v| {
                    Ok(FormulaNumericColumn {
                        name: v.name.clone(),
                        values: numbers(&v.values)?,
                    })
                })
                .collect::<Result<Vec<_>, Box<dyn Error>>>()?,
        };
        let builder = builder.try_model_frame(frame)?;
        let reduced = if case.test == "LRT" {
            Some(
                self.reduced
                    .as_deref()
                    .ok_or("LRT requires a reduced formula")?,
            )
        } else {
            None
        };
        let mut request = Analysis::formula(counts, &self.full);
        if let Some(reduced) = reduced {
            request = request.lrt_formula(reduced);
        }
        let mut options = CooksReplacementOptions::new(super::super::parse_numeric_input(
            &case.replacement_cutoff,
        )?);
        options.min_replicates = case.min_replicates;
        if case.min_replicates != 0 {
            request = request.cooks_replacement(&options);
        }
        match builder.run(&request)? {
            RunOutput::Results { fit, results } => Ok((fit, results, None)),
            RunOutput::Replacement(output) => {
                let (f, r, p) = output.into_merged_fit()?;
                Ok((f, r, Some(p)))
            }
            _ => Err("unexpected run output".into()),
        }
    }
}
