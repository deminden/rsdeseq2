//! Standalone numerical conformance driver; it neither links nor launches R.
#[path = "evaluate_r_interface/adaptive.rs"]
mod adaptive;
#[path = "evaluate_r_interface/batch.rs"]
mod batch;
#[path = "evaluate_r_interface/cooks.rs"]
mod cooks;
#[path = "evaluate_r_interface/factor_expansion.rs"]
mod factor_expansion;
#[path = "evaluate_r_interface/formula.rs"]
mod formula;
#[path = "evaluate_r_interface/frozen_local.rs"]
mod frozen_local;
#[path = "evaluate_r_interface/heavy_tail.rs"]
mod heavy_tail;
#[path = "evaluate_r_interface/likelihood.rs"]
mod likelihood;
#[path = "evaluate_r_interface/location_shrinkage.rs"]
mod location_shrinkage;
#[path = "evaluate_r_interface/mixture_fit.rs"]
mod mixture_fit;
#[path = "evaluate_r_interface/mixture_grid.rs"]
mod mixture_grid;
#[path = "evaluate_r_interface/model_matrix.rs"]
mod model_matrix;
#[path = "evaluate_r_interface/native_formula.rs"]
mod native_formula;
#[path = "evaluate_r_interface/normal_log.rs"]
mod normal_log;
#[path = "evaluate_r_interface/posterior.rs"]
mod posterior;
#[path = "evaluate_r_interface/preparation.rs"]
mod preparation;
#[path = "evaluate_r_interface/random.rs"]
mod random;
#[path = "evaluate_r_interface/regression.rs"]
mod regression;
#[path = "evaluate_r_interface/uniform.rs"]
mod uniform;
use rsdeseq2::independent_filtering::{
    IndependentFilteringOptions, evaluate_independent_filtering,
};
use serde::Deserialize;
use serde_json::json;
use std::error::Error;

#[derive(Deserialize)]
struct FilteringInput {
    filter: Vec<String>,
    pvalues: Vec<Option<String>>,
    theta: Option<Vec<String>>,
    alpha: String,
    enabled: bool,
}

fn bits(value: f64) -> String {
    format!("{:016x}", value.to_bits())
}

fn parse_numeric_input(value: &str) -> Result<f64, Box<dyn Error>> {
    Ok(if let Some(bits) = value.strip_prefix("bits:") {
        f64::from_bits(u64::from_str_radix(bits, 16)?)
    } else if matches!(value, "NA" | "NaN") {
        f64::NAN
    } else {
        value.parse()?
    })
}

fn parse_numeric_inputs(values: &[String]) -> Result<Vec<f64>, Box<dyn Error>> {
    values.iter().map(|v| parse_numeric_input(v)).collect()
}

fn filtering(input: FilteringInput) -> Result<serde_json::Value, Box<dyn Error>> {
    let parse = |values: Vec<String>| -> Result<Vec<f64>, std::num::ParseFloatError> {
        values.into_iter().map(|v| v.parse()).collect()
    };
    let filter = parse(input.filter)?;
    let pvalues = input
        .pvalues
        .into_iter()
        .map(|v| v.map(|v| v.parse()).transpose())
        .collect::<Result<Vec<_>, _>>()?;
    let options = IndependentFilteringOptions {
        enabled: input.enabled,
        alpha: input.alpha.parse()?,
        theta: input.theta.map(parse).transpose()?,
    };
    let evaluated = evaluate_independent_filtering(&filter, &pvalues, &options)?;
    let metadata = evaluated.metadata;
    Ok(json!({
            "padj": evaluated.padj.into_iter().map(|v| v.map(bits)).collect::<Vec<_>>(),
            "filtered": evaluated.filtered,
            "enabled": metadata.enabled,
            "theta": metadata.theta.iter().copied().map(bits).collect::<Vec<_>>(),
            "numRej": metadata.num_rejections,
            "selected": metadata.selected_index.map(|v| v + 1),
            "filterTheta": metadata.filter_theta.map(bits),
            "filterThreshold": metadata.filter_threshold.map(bits),
            "lowessX": metadata.lowess_fit.as_ref().map(|_| metadata.theta.iter().copied().map(bits).collect::<Vec<_>>()),
            "lowessY": metadata.lowess_fit.map(|v| v.into_iter().map(bits).collect::<Vec<_>>()),
            "alpha": bits(metadata.alpha),
    }))
}

#[derive(Deserialize)]
struct NormalizationInput {
    genes: usize,
    samples: usize,
    counts: Vec<u32>,
    method: String,
    geo_means: Option<Vec<String>>,
    control_genes: Option<Vec<usize>>,
    size_factors: Option<Vec<String>>,
    normalization_factors: Option<Vec<String>>,
    weights: Option<Vec<String>>,
}

fn normalization(input: NormalizationInput) -> Result<serde_json::Value, Box<dyn Error>> {
    use rsdeseq2::normalization::{
        base_mean, base_mean_with_weights, base_variance, base_variance_with_weights,
        estimate_size_factors_with_options, normalized_counts, normalized_counts_with_factors,
    };
    use rsdeseq2::prelude::{CountMatrix, RowMajorMatrix, SizeFactorMethod};
    let parse = |values: Vec<String>| {
        values
            .into_iter()
            .map(|v| v.parse::<f64>())
            .collect::<Result<Vec<_>, _>>()
    };
    let counts = CountMatrix::from_row_major_u32(input.genes, input.samples, input.counts)?;
    let geo = input.geo_means.map(parse).transpose()?;
    let method = match input.method.as_str() {
        "ratio" => SizeFactorMethod::Ratio,
        "poscounts" => SizeFactorMethod::PosCounts,
        _ => return Err("unknown size-factor method".into()),
    };
    let factors = estimate_size_factors_with_options(
        &counts,
        method,
        geo.as_deref(),
        input.control_genes.as_deref(),
    )?;
    let matrix = |values| -> Result<RowMajorMatrix<f64>, Box<dyn Error>> {
        Ok(RowMajorMatrix::from_row_major(
            input.genes,
            input.samples,
            parse(values)?,
        )?)
    };
    let normalization_factors = input.normalization_factors.map(matrix).transpose()?;
    let weights = input.weights.map(matrix).transpose()?;
    let supplied = input.size_factors.map(parse).transpose()?;
    let normalized = match normalization_factors.as_ref() {
        Some(matrix) => normalized_counts_with_factors(&counts, matrix)?,
        None => normalized_counts(&counts, supplied.as_deref().unwrap_or(&factors))?,
    };
    let (means, variances) = match weights.as_ref() {
        Some(weights) => (
            base_mean_with_weights(&normalized, weights)?,
            base_variance_with_weights(&normalized, weights)?,
        ),
        None => (base_mean(&normalized)?, base_variance(&normalized)?),
    };
    Ok(json!({
        "sizeFactors": factors.into_iter().map(bits).collect::<Vec<_>>(),
        "normalized": normalized.as_slice().iter().copied().map(bits).collect::<Vec<_>>(),
        "baseMean": means.into_iter().map(bits).collect::<Vec<_>>(),
        "baseVar": variances.into_iter().map(bits).collect::<Vec<_>>(),
        "allZero": counts.all_zero_flags(),
    }))
}

#[derive(Deserialize)]
struct ZeroInput {
    samples: usize,
    counts: Vec<u32>,
    design: Option<Vec<String>>,
    levels: Option<Vec<String>>,
    numerator: Option<String>,
    denominator: Option<String>,
    rule: String,
}
impl ZeroInput {
    fn apply(
        self,
        genes: usize,
        coefficients: usize,
        contrast: &[f64],
        effect: &mut [Option<f64>],
        inference: Option<&mut rsdeseq2::glm::WaldOutput>,
    ) -> Result<(), Box<dyn Error>> {
        use rsdeseq2::prelude::{
            ContrastCounts, ContrastZeroPolicy, CountMatrix, MatrixOrder, MatrixView,
        };
        let counts = CountMatrix::from_row_major_u32(genes, self.samples, self.counts)?;
        let counts = ContrastCounts::from(&counts);
        let design = self
            .design
            .map(|v| {
                v.into_iter()
                    .map(|v| v.parse::<f64>())
                    .collect::<Result<Vec<_>, _>>()
            })
            .transpose()?;
        let flags = match self.rule.as_str() {
            "numeric" => counts.numeric_flags(
                design
                    .as_deref()
                    .map(|d| MatrixView::new(self.samples, coefficients, d, MatrixOrder::RowMajor))
                    .transpose()?,
                contrast,
            )?,
            "character" => counts.factor_flags(
                self.levels.as_deref().ok_or("missing sample levels")?,
                self.numerator.as_deref().ok_or("missing numerator")?,
                self.denominator.as_deref().ok_or("missing denominator")?,
            )?,
            _ => return Err("unknown zero rule".into()),
        };
        let policy = ContrastZeroPolicy::new(&flags, counts.all_zero())?;
        if let Some(inference) = inference {
            policy.apply_wald(effect, &mut inference.stat, &mut inference.pvalue)?;
        } else {
            policy.apply_effects(effect.iter_mut())?;
        }
        Ok(())
    }
}

#[derive(Deserialize)]
#[serde(tag = "operation", rename_all = "snake_case")]
enum StatisticsInput {
    Wald {
        beta: Vec<String>,
        se: Vec<String>,
        threshold: String,
        alternative: String,
        df: Option<Vec<String>>,
    },
    Projection {
        genes: usize,
        coefficients: usize,
        beta: Vec<String>,
        se: Option<Vec<String>>,
        covariance: Option<Vec<String>>,
        contrast: Vec<String>,
        threshold: Option<String>,
        alternative: Option<String>,
        df: Option<Vec<String>>,
        zero: Option<Box<ZeroInput>>,
    },
    Lrt {
        statistics: Vec<String>,
        df: Vec<String>,
    },
    Cooks {
        samples: usize,
        coefficients: usize,
    },
}

fn statistics(input: StatisticsInput) -> Result<serde_json::Value, Box<dyn Error>> {
    use rsdeseq2::glm::{
        WaldAlternative, WaldTestOptions, lrt_pvalues_from_statistics, wald_test_from_estimates,
    };
    let parse = |values: Vec<String>| {
        values
            .into_iter()
            .map(|v| v.parse::<f64>())
            .collect::<Result<Vec<_>, _>>()
    };
    match input {
        StatisticsInput::Wald {
            beta,
            se,
            threshold,
            alternative,
            df,
        } => {
            let alternative = match alternative.as_str() {
                "greaterAbs" => WaldAlternative::GreaterAbs,
                "greaterAbsUPSHOT" => WaldAlternative::GreaterAbsUpshot,
                "greaterAbs2014" => WaldAlternative::GreaterAbs2014,
                "lessAbs" => WaldAlternative::LessAbs,
                "greater" => WaldAlternative::Greater,
                "less" => WaldAlternative::Less,
                _ => return Err("unknown Wald alternative".into()),
            };
            let df = df.map(parse).transpose()?;
            let options = match df {
                None => WaldTestOptions::normal(),
                Some(df) if df.len() == 1 => WaldTestOptions::t_degrees_of_freedom(df[0]),
                Some(df) => WaldTestOptions::t_per_gene_degrees_of_freedom(df),
            }
            .with_lfc_threshold(threshold.parse()?, alternative);
            let result = wald_test_from_estimates(&parse(beta)?, &parse(se)?, &options)?;
            Ok(
                json!({"stat": result.stat.into_iter().map(|v| v.map(bits)).collect::<Vec<_>>(),
                "pvalue": result.pvalue.into_iter().map(|v| v.map(bits)).collect::<Vec<_>>() }),
            )
        }
        StatisticsInput::Projection {
            genes,
            coefficients,
            beta,
            se,
            covariance,
            contrast,
            threshold,
            alternative,
            df,
            zero,
        } => {
            use rsdeseq2::glm::{
                CovarianceView, InferenceEstimates, wald_test_contrast_from_estimates,
            };
            use rsdeseq2::matrix::{MatrixOrder, MatrixView};
            let beta = parse(beta)?;
            let se = se.map(parse).transpose()?;
            let covariance = covariance.map(parse).transpose()?;
            let matrix =
                |values| MatrixView::new(genes, coefficients, values, MatrixOrder::RowMajor);
            let covariance = covariance
                .as_deref()
                .map(|values| {
                    CovarianceView::new(
                        coefficients,
                        MatrixView::new(
                            genes,
                            coefficients
                                .checked_mul(coefficients)
                                .ok_or("covariance shape overflow")?,
                            values,
                            MatrixOrder::RowMajor,
                        )?,
                    )
                    .map_err(|e| -> Box<dyn Error> { Box::new(e) })
                })
                .transpose()?;
            let estimates = InferenceEstimates::new(
                matrix(&beta)?,
                se.as_deref().map(matrix).transpose()?,
                covariance,
            )?;
            let contrast = parse(contrast)?;
            if let Some(threshold) = threshold {
                let alternative = match alternative.as_deref().unwrap_or("greaterAbs") {
                    "greaterAbs" => WaldAlternative::GreaterAbs,
                    "greaterAbsUPSHOT" => WaldAlternative::GreaterAbsUpshot,
                    "greaterAbs2014" => WaldAlternative::GreaterAbs2014,
                    "lessAbs" => WaldAlternative::LessAbs,
                    "greater" => WaldAlternative::Greater,
                    "less" => WaldAlternative::Less,
                    _ => return Err("unknown Wald alternative".into()),
                };
                let options = match df.map(parse).transpose()? {
                    None => WaldTestOptions::normal(),
                    Some(df) if df.len() == 1 => WaldTestOptions::t_degrees_of_freedom(df[0]),
                    Some(df) => WaldTestOptions::t_per_gene_degrees_of_freedom(df),
                }
                .with_lfc_threshold(threshold.parse()?, alternative);
                let mut output =
                    wald_test_contrast_from_estimates(&estimates, &contrast, &options)?;
                if let Some(zero) = zero {
                    (*zero).apply(
                        genes,
                        coefficients,
                        &contrast,
                        &mut output.log2_fold_change,
                        Some(&mut output.wald),
                    )?;
                }
                Ok(
                    json!({"log2FoldChange":output.log2_fold_change.into_iter().map(|v|v.map(bits)).collect::<Vec<_>>(),
                    "lfcSE":output.lfc_se.into_iter().map(|v|v.map(bits)).collect::<Vec<_>>(),
                    "stat":output.wald.stat.into_iter().map(|v|v.map(bits)).collect::<Vec<_>>(),
                    "pvalue":output.wald.pvalue.into_iter().map(|v|v.map(bits)).collect::<Vec<_>>() }),
                )
            } else {
                let mut output = estimates.project_contrast(&contrast)?;
                if let Some(zero) = zero {
                    (*zero).apply(
                        genes,
                        coefficients,
                        &contrast,
                        &mut output.log2_fold_change,
                        None,
                    )?;
                }
                Ok(
                    json!({"log2FoldChange":output.log2_fold_change.into_iter().map(|v|v.map(bits)).collect::<Vec<_>>(),
                    "lfcSE":output.lfc_se.into_iter().map(|v|v.map(bits)).collect::<Vec<_>>() }),
                )
            }
        }
        StatisticsInput::Lrt { statistics, df } => {
            let result = lrt_pvalues_from_statistics(&parse(statistics)?, &parse(df)?)?;
            Ok(json!({"pvalue": result.into_iter().map(|v| v.map(bits)).collect::<Vec<_>>() }))
        }
        StatisticsInput::Cooks {
            samples,
            coefficients,
        } => Ok(
            json!({"cutoff": rsdeseq2::results::default_cooks_cutoff(samples, coefficients)?.map(bits)}),
        ),
    }
}

#[derive(Deserialize)]
#[serde(untagged)]
enum Input {
    Random(random::Input),
    Regression(regression::Input),
    HeavyTail(heavy_tail::Input),
    Adaptive(adaptive::Input),
    Uniform(uniform::Input),
    MixtureGrid(mixture_grid::Input),
    LocationShrinkage(location_shrinkage::Input),
    NormalLog(normal_log::Input),
    MixtureFit(mixture_fit::Input),
    Posterior(posterior::Input),
    FrozenLocal(frozen_local::Input),
    NativeFormula(native_formula::Input),
    Formula(formula::Input),
    ModelMatrix(model_matrix::Case),
    Preparation(preparation::Input),
    CooksPolicy(cooks::Input),
    Batch(batch::BatchInput),
    Statistics(StatisticsInput),
    Filtering(FilteringInput),
    Normalization(NormalizationInput),
}

fn main() -> Result<(), Box<dyn Error>> {
    let path = std::env::args()
        .nth(1)
        .ok_or("expected an input JSON path")?;
    let input: Input = serde_json::from_slice(&std::fs::read(path)?)?;
    let result = match input {
        Input::Random(input) => random::evaluate(input)?,
        Input::Regression(input) => regression::evaluate(input)?,
        Input::HeavyTail(input) => heavy_tail::evaluate(input)?,
        Input::Adaptive(input) => adaptive::evaluate(input)?,
        Input::Uniform(input) => uniform::evaluate(input)?,
        Input::MixtureGrid(input) => mixture_grid::evaluate(input)?,
        Input::LocationShrinkage(input) => location_shrinkage::evaluate(input)?,
        Input::NormalLog(input) => normal_log::evaluate(input)?,
        Input::MixtureFit(input) => mixture_fit::evaluate(input)?,
        Input::Posterior(input) => posterior::evaluate(input)?,
        Input::FrozenLocal(input) => frozen_local::evaluate(input)?,
        Input::NativeFormula(input) => native_formula::evaluate(&input),
        Input::Formula(input) => formula::panel(input),
        Input::ModelMatrix(input) => {
            model_matrix::evaluate(&input, rsdeseq2::matrix::MatrixOrder::RowMajor)
        }
        Input::Preparation(input) => preparation::evaluate(input)?,
        Input::CooksPolicy(input) => cooks::evaluate(input)?,
        Input::Batch(input) => batch::evaluate(input)?,
        Input::Statistics(input) => statistics(input)?,
        Input::Filtering(input) => filtering(input)?,
        Input::Normalization(input) => normalization(input)?,
    };
    println!("{result}");
    Ok(())
}
