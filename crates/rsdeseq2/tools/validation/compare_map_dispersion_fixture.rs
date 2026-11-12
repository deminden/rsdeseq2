use std::collections::HashMap;
use std::env;
use std::fs;
use std::path::{Path, PathBuf};

use rsdeseq2::prelude::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let fixture = env::args().nth(1).map(PathBuf::from).unwrap_or_else(|| {
        PathBuf::from("../../results/fixtures/se_covariance_hard_real_2026-06-05")
    });

    let manifest = read_key_value_table(&fixture.join("manifest.tsv"))?;
    let disp_prior_var = required_value(&manifest, "disp_prior_var")?;
    let var_log_disp_estimates = required_value(&manifest, "var_log_disp_estimates")?;

    let (genes, samples, counts_values) = read_count_matrix(&fixture.join("selected_counts.tsv"))?;
    let counts = CountMatrix::from_row_major_u32_with_names(
        genes.len(),
        samples.len(),
        counts_values,
        Some(genes.clone()),
        Some(samples),
    )?;
    let (coefficient_names, design_values) =
        read_design_matrix(&fixture.join("design_matrix.tsv"))?;
    let design = DesignMatrix::from_row_major(
        counts.n_samples(),
        coefficient_names.len(),
        design_values,
        Some(coefficient_names),
    )?;
    let map_mu = read_named_matrix(&fixture.join("map_input_mu.tsv"))?;
    let selected = read_gene_table(&fixture.join("selected_dispersions.tsv"))?;
    let line_search_reference = read_gene_table(&fixture.join("map_line_search.tsv"))?;
    let mu = RowMajorMatrix::from_row_major(
        counts.n_genes(),
        counts.n_samples(),
        genes
            .iter()
            .flat_map(|gene| map_mu.rows[gene].iter().copied())
            .collect(),
    )?;
    let disp_gene_est = genes
        .iter()
        .map(|gene| selected[gene]["dispGeneEst"])
        .collect::<Vec<_>>();
    let disp_fit = genes
        .iter()
        .map(|gene| selected[gene]["dispFit"])
        .collect::<Vec<_>>();
    let all_zero = vec![false; counts.n_genes()];

    let map = estimate_map_dispersions(
        MapDispersionInput {
            counts: &counts,
            design: &design,
            mu: &mu,
            disp_gene_est: &disp_gene_est,
            disp_fit: &disp_fit,
            all_zero: &all_zero,
            observation_weights: None,
            disp_prior_var,
            var_log_disp_estimates,
        },
        MapDispersionOptions::default(),
    )?;

    summarize_gene_values("dispMAP", &genes, &map.disp_map, &selected, "dispMAP");
    summarize_gene_values(
        "dispersion",
        &genes,
        &map.dispersion,
        &selected,
        "dispersion",
    );
    summarize_gene_values("dispInit", &genes, &map.disp_init, &selected, "dispInit");
    summarize_gene_iters("dispIter", &genes, &map.disp_iter, &selected, "dispIter");
    let line_search =
        run_line_search_trace(&genes, &counts, &design, &mu, &selected, disp_prior_var)?;
    summarize_gene_values(
        "lineLogAlpha",
        &genes,
        &line_search.log_alpha,
        &line_search_reference,
        "logAlpha",
    );
    summarize_gene_values(
        "initialLp",
        &genes,
        &line_search.initial_lp,
        &line_search_reference,
        "initialLp",
    );
    summarize_gene_values(
        "initialDlp",
        &genes,
        &line_search.initial_dlp,
        &line_search_reference,
        "initialDlp",
    );
    summarize_gene_values(
        "lastLp",
        &genes,
        &line_search.last_lp,
        &line_search_reference,
        "lastLp",
    );
    summarize_gene_values(
        "lastDlp",
        &genes,
        &line_search.last_dlp,
        &line_search_reference,
        "lastDlp",
    );
    summarize_gene_iters(
        "lineIter",
        &genes,
        &line_search.iter,
        &line_search_reference,
        "iter",
    );
    summarize_gene_iters(
        "lineIterAccept",
        &genes,
        &line_search.iter_accept,
        &line_search_reference,
        "iterAccept",
    );
    if env::var_os("RSDESEQ2_PRINT_MISMATCHES").is_some() {
        print_line_search_mismatches(&genes, &line_search, &line_search_reference);
    }
    print_component_breakdown(
        "MT-ND4",
        &genes,
        &counts,
        &design,
        &mu,
        &selected,
        disp_prior_var,
    )?;
    print_component_breakdown(
        "ANKLE2",
        &genes,
        &counts,
        &design,
        &mu,
        &selected,
        disp_prior_var,
    )?;
    if let Ok(gene) = env::var("RSDESEQ2_TRACE_GENE") {
        print_wrapper_vs_line_search(
            &gene,
            &genes,
            &map,
            &line_search,
            &selected,
            &line_search_reference,
        );
        print_line_search_trace(
            &gene,
            &genes,
            &counts,
            &design,
            &mu,
            &selected,
            disp_prior_var,
        )?;
    }
    if let Ok(points) = env::var("RSDESEQ2_COMPONENT_POINTS") {
        let view = MapFixtureView {
            genes: &genes,
            counts: &counts,
            design: &design,
            mu: &mu,
            selected: &selected,
            disp_prior_var,
        };
        print_component_points(&points, view)?;
    }
    Ok(())
}

fn print_line_search_mismatches(
    genes: &[String],
    line_search: &LineSearchTrace,
    reference: &HashMap<String, HashMap<String, f64>>,
) {
    for (gene_idx, gene) in genes.iter().enumerate() {
        let expected_iter = reference[gene]["iter"] as usize;
        let iter_diff = line_search.iter[gene_idx].abs_diff(expected_iter);
        let log_diff = line_search.log_alpha[gene_idx] - reference[gene]["logAlpha"];
        let dlp_diff = line_search.last_dlp[gene_idx] - reference[gene]["lastDlp"];
        if iter_diff > 0 || log_diff.abs() > 1.0e-12 || dlp_diff.abs() > 1.0e-12 {
            println!(
                "lineMismatch\tgene={gene}\titer={}\texpectedIter={expected_iter}\titerDiff={iter_diff}\tlogAlpha={:.16e}\texpectedLogAlpha={:.16e}\tlogDiff={log_diff:.16e}\tlastDlp={:.16e}\texpectedLastDlp={:.16e}\tdlpDiff={dlp_diff:.16e}",
                line_search.iter[gene_idx],
                line_search.log_alpha[gene_idx],
                reference[gene]["logAlpha"],
                line_search.last_dlp[gene_idx],
                reference[gene]["lastDlp"],
            );
        }
    }
}

fn print_wrapper_vs_line_search(
    gene: &str,
    genes: &[String],
    map: &MapDispersionOutput,
    line_search: &LineSearchTrace,
    selected: &HashMap<String, HashMap<String, f64>>,
    line_search_reference: &HashMap<String, HashMap<String, f64>>,
) {
    let Some(gene_idx) = genes.iter().position(|candidate| candidate == gene) else {
        return;
    };
    println!(
        "wrapperVsLine\tgene={gene}\twrapperDispMap={:.16e}\twrapperLogAlpha={:.16e}\tdirectLogAlpha={:.16e}\tdirectDispMap={:.16e}\texpectedLogAlpha={:.16e}\texpectedDispMap={:.16e}\twrapperIter={}\tdirectIter={}\texpectedIter={:.0}\tdirectLastLp={:.16e}\tdirectLastDlp={:.16e}\texpectedLastLp={:.16e}\texpectedLastDlp={:.16e}",
        map.disp_map[gene_idx],
        map.disp_map[gene_idx].ln(),
        line_search.log_alpha[gene_idx],
        line_search.log_alpha[gene_idx].exp(),
        line_search_reference[gene]["logAlpha"],
        selected[gene]["dispMAP"],
        map.disp_iter[gene_idx],
        line_search.iter[gene_idx],
        line_search_reference[gene]["iter"],
        line_search.last_lp[gene_idx],
        line_search.last_dlp[gene_idx],
        line_search_reference[gene]["lastLp"],
        line_search_reference[gene]["lastDlp"],
    );
}

#[derive(Clone, Copy)]
struct MapFixtureView<'a> {
    genes: &'a [String],
    counts: &'a CountMatrix,
    design: &'a DesignMatrix,
    mu: &'a RowMajorMatrix<f64>,
    selected: &'a HashMap<String, HashMap<String, f64>>,
    disp_prior_var: f64,
}

struct LineSearchTrace {
    log_alpha: Vec<f64>,
    iter: Vec<usize>,
    iter_accept: Vec<usize>,
    initial_lp: Vec<f64>,
    initial_dlp: Vec<f64>,
    last_lp: Vec<f64>,
    last_dlp: Vec<f64>,
}

fn run_line_search_trace(
    genes: &[String],
    counts: &CountMatrix,
    design: &DesignMatrix,
    mu: &RowMajorMatrix<f64>,
    selected: &HashMap<String, HashMap<String, f64>>,
    disp_prior_var: f64,
) -> Result<LineSearchTrace, Box<dyn std::error::Error>> {
    let mut trace = LineSearchTrace {
        log_alpha: Vec::with_capacity(genes.len()),
        iter: Vec::with_capacity(genes.len()),
        iter_accept: Vec::with_capacity(genes.len()),
        initial_lp: Vec::with_capacity(genes.len()),
        initial_dlp: Vec::with_capacity(genes.len()),
        last_lp: Vec::with_capacity(genes.len()),
        last_dlp: Vec::with_capacity(genes.len()),
    };
    let options = GeneWiseDispersionOptions::default();
    for (gene_idx, gene) in genes.iter().enumerate() {
        let prior = DispersionPrior::new(selected[gene]["dispFit"].ln(), disp_prior_var)?;
        let fit = fit_dispersion_line_search_with_prior(
            counts.row(gene_idx)?,
            mu.row(gene_idx)?,
            design,
            selected[gene]["dispInit"],
            options,
            counts.n_samples(),
            prior,
        )?;
        trace.log_alpha.push(fit.log_alpha);
        trace.iter.push(fit.iter);
        trace.iter_accept.push(fit.iter_accept);
        trace.initial_lp.push(fit.initial_lp);
        trace.initial_dlp.push(fit.initial_dlp);
        trace.last_lp.push(fit.last_lp);
        trace.last_dlp.push(fit.last_dlp);
    }
    Ok(trace)
}

fn print_component_breakdown(
    gene: &str,
    genes: &[String],
    counts: &CountMatrix,
    design: &DesignMatrix,
    mu: &RowMajorMatrix<f64>,
    selected: &HashMap<String, HashMap<String, f64>>,
    disp_prior_var: f64,
) -> Result<(), Box<dyn std::error::Error>> {
    let Some(gene_idx) = genes.iter().position(|candidate| candidate == gene) else {
        return Ok(());
    };
    let log_alpha = selected[gene]["dispInit"].ln();
    let log_mean = selected[gene]["dispFit"].ln();
    let prior = DispersionPrior::new(log_mean, disp_prior_var)?;
    let likelihood =
        dispersion_nb_log_likelihood_kernel(counts.row(gene_idx)?, mu.row(gene_idx)?, log_alpha)?;
    let likelihood_derivative = dispersion_nb_log_likelihood_kernel_derivative(
        counts.row(gene_idx)?,
        mu.row(gene_idx)?,
        log_alpha,
    )?;
    let cox_reid = cox_reid_adjustment(design, mu.row(gene_idx)?, log_alpha)?;
    let cox_reid_derivative = cox_reid_adjustment_derivative(design, mu.row(gene_idx)?, log_alpha)?;
    let prior_density = dispersion_prior_log_density(log_alpha, prior)?;
    let prior_derivative = dispersion_prior_derivative(log_alpha, prior)?;
    println!(
        "components\tgene={gene}\tlog_alpha={log_alpha:.16e}\tlikelihood={likelihood:.16e}\tcox_reid={cox_reid:.16e}\tprior={prior_density:.16e}\ttotal={:.16e}\tlikelihood_d={likelihood_derivative:.16e}\tcox_reid_d={cox_reid_derivative:.16e}\tprior_d={prior_derivative:.16e}\ttotal_d={:.16e}",
        likelihood + cox_reid + prior_density,
        likelihood_derivative + cox_reid_derivative + prior_derivative
    );
    Ok(())
}

fn print_component_points(
    points: &str,
    view: MapFixtureView<'_>,
) -> Result<(), Box<dyn std::error::Error>> {
    for spec in points.split(';').filter(|spec| !spec.is_empty()) {
        let Some((gene, values)) = spec.split_once(':') else {
            continue;
        };
        for value in values.split(',').filter(|value| !value.is_empty()) {
            let log_alpha = value.parse::<f64>()?;
            print_component_at_log_alpha(gene, log_alpha, view)?;
        }
    }
    Ok(())
}

fn print_component_at_log_alpha(
    gene: &str,
    log_alpha: f64,
    view: MapFixtureView<'_>,
) -> Result<(), Box<dyn std::error::Error>> {
    let Some(gene_idx) = view.genes.iter().position(|candidate| candidate == gene) else {
        return Ok(());
    };
    let log_mean = view.selected[gene]["dispFit"].ln();
    let prior = DispersionPrior::new(log_mean, view.disp_prior_var)?;
    let likelihood = dispersion_nb_log_likelihood_kernel(
        view.counts.row(gene_idx)?,
        view.mu.row(gene_idx)?,
        log_alpha,
    )?;
    let likelihood_derivative = dispersion_nb_log_likelihood_kernel_derivative(
        view.counts.row(gene_idx)?,
        view.mu.row(gene_idx)?,
        log_alpha,
    )?;
    let cox_reid = cox_reid_adjustment(view.design, view.mu.row(gene_idx)?, log_alpha)?;
    let cox_reid_derivative =
        cox_reid_adjustment_derivative(view.design, view.mu.row(gene_idx)?, log_alpha)?;
    let prior_density = dispersion_prior_log_density(log_alpha, prior)?;
    let prior_derivative = dispersion_prior_derivative(log_alpha, prior)?;
    let total = likelihood + cox_reid + prior_density;
    let total_derivative = likelihood_derivative + cox_reid_derivative + prior_derivative;
    println!(
        "componentPoint\tgene={gene}\tlog_alpha={log_alpha:.16e}\tlikelihood={likelihood:.16e}\tcox_reid={cox_reid:.16e}\tprior={prior_density:.16e}\ttotal={total:.16e}\ttheta={:.16e}\tlikelihood_d={likelihood_derivative:.16e}\tcox_reid_d={cox_reid_derivative:.16e}\tprior_d={prior_derivative:.16e}\ttotal_d={total_derivative:.16e}",
        -total
    );
    Ok(())
}

fn print_line_search_trace(
    gene: &str,
    genes: &[String],
    counts: &CountMatrix,
    design: &DesignMatrix,
    mu: &RowMajorMatrix<f64>,
    selected: &HashMap<String, HashMap<String, f64>>,
    disp_prior_var: f64,
) -> Result<(), Box<dyn std::error::Error>> {
    let Some(gene_idx) = genes.iter().position(|candidate| candidate == gene) else {
        return Ok(());
    };
    let row_counts = counts.row(gene_idx)?;
    let row_mu = mu.row(gene_idx)?;
    let prior = DispersionPrior::new(selected[gene]["dispFit"].ln(), disp_prior_var)?;
    let mut log_alpha = selected[gene]["dispInit"].ln();
    let mut lp = map_objective(row_counts, row_mu, design, log_alpha, prior)?;
    let mut dlp = map_derivative(row_counts, row_mu, design, log_alpha, prior)?;
    let mut kappa = 1.0;
    let epsilon = 1.0e-4;
    let mut iter_accept = 0_usize;
    println!("traceStart\tgene={gene}\ta={log_alpha:.16e}\tlp={lp:.16e}\tdlp={dlp:.16e}");
    for iter in 1..=100 {
        let mut proposed = log_alpha + kappa * dlp;
        if proposed < -30.0 {
            kappa = (-30.0 - log_alpha) / dlp;
            proposed = log_alpha + kappa * dlp;
        }
        if proposed > 10.0 {
            kappa = (10.0 - log_alpha) / dlp;
            proposed = log_alpha + kappa * dlp;
        }
        let theta = -map_objective(row_counts, row_mu, design, proposed, prior)?;
        let theta_hat = -lp - kappa * epsilon * dlp * dlp;
        let accept = theta <= theta_hat;
        println!(
            "traceStep\titer={iter}\taccept={accept}\ta={log_alpha:.16e}\tproposed={proposed:.16e}\tkappa={kappa:.16e}\tlp={lp:.16e}\tdlp={dlp:.16e}\ttheta={theta:.16e}\tthetaHat={theta_hat:.16e}\tdelta={:.16e}",
            theta - theta_hat
        );
        if accept {
            iter_accept += 1;
            log_alpha = proposed;
            let lp_new = map_objective(row_counts, row_mu, design, log_alpha, prior)?;
            let change = lp_new - lp;
            println!(
                "traceAccept\titer={iter}\tacceptCount={iter_accept}\ta={log_alpha:.16e}\tlpNew={lp_new:.16e}\tchange={change:.16e}"
            );
            lp = lp_new;
            if change < 1.0e-6 {
                break;
            }
            if log_alpha < (1.0e-8_f64 / 10.0).ln() {
                break;
            }
            dlp = map_derivative(row_counts, row_mu, design, log_alpha, prior)?;
            kappa = (kappa * 1.1).min(1.0);
            if iter_accept.is_multiple_of(5) {
                kappa /= 2.0;
            }
        } else {
            kappa /= 2.0;
        }
    }
    Ok(())
}

fn map_objective(
    counts: &[u32],
    mu: &[f64],
    design: &DesignMatrix,
    log_alpha: f64,
    prior: DispersionPrior,
) -> Result<f64, Box<dyn std::error::Error>> {
    Ok(dispersion_nb_log_likelihood_kernel(counts, mu, log_alpha)?
        + cox_reid_adjustment(design, mu, log_alpha)?
        + dispersion_prior_log_density(log_alpha, prior)?)
}

fn map_derivative(
    counts: &[u32],
    mu: &[f64],
    design: &DesignMatrix,
    log_alpha: f64,
    prior: DispersionPrior,
) -> Result<f64, Box<dyn std::error::Error>> {
    Ok(
        dispersion_nb_log_likelihood_kernel_derivative(counts, mu, log_alpha)?
            + cox_reid_adjustment_derivative(design, mu, log_alpha)?
            + dispersion_prior_derivative(log_alpha, prior)?,
    )
}

fn summarize_gene_values(
    label: &str,
    genes: &[String],
    actual: &[f64],
    expected: &HashMap<String, HashMap<String, f64>>,
    expected_column: &str,
) {
    let mut diffs = Vec::with_capacity(genes.len());
    let mut worst = ("", 0.0, 0.0, 0.0);
    for (gene_idx, gene) in genes.iter().enumerate() {
        let expected_value = expected[gene][expected_column];
        let actual_value = actual[gene_idx];
        if !actual_value.is_finite() || !expected_value.is_finite() {
            continue;
        }
        let diff = (actual_value - expected_value).abs();
        diffs.push(diff);
        if diff > worst.1 {
            worst = (gene.as_str(), diff, actual_value, expected_value);
        }
    }
    summarize_diffs(label, &mut diffs, worst);
}

fn summarize_gene_iters(
    label: &str,
    genes: &[String],
    actual: &[usize],
    expected: &HashMap<String, HashMap<String, f64>>,
    expected_column: &str,
) {
    let mut worst = ("", 0_usize, 0_usize, 0_usize);
    let diffs = genes
        .iter()
        .enumerate()
        .map(|(gene_idx, gene)| {
            let expected_iter = expected[gene][expected_column] as usize;
            let diff = actual[gene_idx].abs_diff(expected_iter);
            if diff > worst.1 {
                worst = (gene.as_str(), diff, actual[gene_idx], expected_iter);
            }
            diff
        })
        .collect::<Vec<_>>();
    let mismatches = diffs.iter().filter(|diff| **diff > 0).count();
    let max = diffs.iter().copied().max().unwrap_or(0);
    println!(
        "{label}\tmismatches={mismatches}\tmax_abs_iter_diff={max}\tworst_gene={}\tactual={}\texpected={}",
        worst.0, worst.2, worst.3
    );
}

fn summarize_diffs(label: &str, diffs: &mut [f64], worst: (&str, f64, f64, f64)) {
    diffs.sort_by(|left, right| left.total_cmp(right));
    let mean = diffs.iter().sum::<f64>() / diffs.len() as f64;
    let median = diffs[diffs.len() / 2];
    let p95 = diffs[((diffs.len() as f64 * 0.95).ceil() as usize).saturating_sub(1)];
    println!(
        "{label}\tmean={mean:.6e}\tmedian={median:.6e}\tp95={p95:.6e}\tmax={:.6e}\tworst_gene={}\tactual={:.16e}\texpected={:.16e}",
        worst.1, worst.0, worst.2, worst.3
    );
}

fn required_value(
    values: &HashMap<String, String>,
    key: &str,
) -> Result<f64, Box<dyn std::error::Error>> {
    values
        .get(key)
        .ok_or_else(|| format!("missing manifest key {key}"))?
        .parse()
        .map_err(Into::into)
}

fn read_key_value_table(
    path: &Path,
) -> Result<HashMap<String, String>, Box<dyn std::error::Error>> {
    let (_, rows) = read_tsv(path)?;
    let mut values = HashMap::with_capacity(rows.len());
    for row in rows {
        values.insert(row[0].clone(), row[1].clone());
    }
    Ok(values)
}

fn read_gene_table(
    path: &Path,
) -> Result<HashMap<String, HashMap<String, f64>>, Box<dyn std::error::Error>> {
    let (header, rows) = read_tsv(path)?;
    let mut values = HashMap::with_capacity(rows.len());
    for row in rows {
        let mut columns = HashMap::with_capacity(header.len().saturating_sub(1));
        for (name, value) in header[1..].iter().zip(&row[1..]) {
            if let Ok(parsed) = parse_f64(value) {
                columns.insert(name.clone(), parsed);
            }
        }
        values.insert(row[0].clone(), columns);
    }
    Ok(values)
}

type TsvTable = (Vec<String>, Vec<Vec<String>>);

fn read_tsv(path: &Path) -> Result<TsvTable, Box<dyn std::error::Error>> {
    let text = fs::read_to_string(path)?;
    let mut lines = text.lines();
    let header = lines
        .next()
        .ok_or_else(|| format!("empty TSV: {}", path.display()))?
        .split('\t')
        .map(str::to_string)
        .collect::<Vec<_>>();
    let rows = lines
        .filter(|line| !line.is_empty())
        .map(|line| line.split('\t').map(str::to_string).collect::<Vec<_>>())
        .collect::<Vec<_>>();
    Ok((header, rows))
}

type CountMatrixParts = (Vec<String>, Vec<String>, Vec<u32>);

fn read_count_matrix(path: &Path) -> Result<CountMatrixParts, Box<dyn std::error::Error>> {
    let (header, rows) = read_tsv(path)?;
    let samples = header[1..].to_vec();
    let mut genes = Vec::with_capacity(rows.len());
    let mut values = Vec::with_capacity(rows.len() * samples.len());
    for row in rows {
        genes.push(row[0].clone());
        for value in &row[1..] {
            values.push(value.parse()?);
        }
    }
    Ok((genes, samples, values))
}

fn read_design_matrix(path: &Path) -> Result<(Vec<String>, Vec<f64>), Box<dyn std::error::Error>> {
    let (header, rows) = read_tsv(path)?;
    let coefficient_names = header[1..].to_vec();
    let mut values = Vec::with_capacity(rows.len() * coefficient_names.len());
    for row in rows {
        for value in &row[1..] {
            values.push(value.parse()?);
        }
    }
    Ok((coefficient_names, values))
}

struct NamedMatrix {
    rows: HashMap<String, Vec<f64>>,
}

fn read_named_matrix(path: &Path) -> Result<NamedMatrix, Box<dyn std::error::Error>> {
    let (_, rows) = read_tsv(path)?;
    let mut values = HashMap::with_capacity(rows.len());
    for row in rows {
        values.insert(
            row[0].clone(),
            row[1..]
                .iter()
                .map(|value| parse_f64(value))
                .collect::<Result<Vec<f64>, _>>()?,
        );
    }
    Ok(NamedMatrix { rows: values })
}

fn parse_f64(value: &str) -> Result<f64, Box<dyn std::error::Error>> {
    if value == "NA" || value.is_empty() {
        Ok(f64::NAN)
    } else {
        Ok(value.parse()?)
    }
}
