use rsdeseq2::matrix::{MatrixOrder, MatrixView};
use rsdeseq2::shrinkage::{MixtureFitOptions, NormalMixture, fit_mixture_proportions};
use serde::Deserialize;
#[derive(Deserialize)]
struct References {
    cases: Vec<Case>,
}
#[derive(Deserialize)]
struct Case {
    name: String,
    observations: usize,
    components: usize,
    means: Vec<String>,
    sds: Vec<String>,
    estimates: Vec<String>,
    standard_errors: Vec<String>,
    dirichlet: Vec<String>,
    observation_weights: Vec<String>,
    log_likelihoods: Vec<String>,
    weighted_log_likelihood: Vec<String>,
    posterior_mean: Vec<String>,
    posterior_sd: Vec<String>,
    converged: bool,
}
fn numbers(x: &[String]) -> Vec<f64> {
    x.iter().map(|s| s.parse().unwrap()).collect()
}
#[test]
fn convex_mixture_fits_match_independent_ashr_likelihoods_and_posteriors() {
    let references: References =
        serde_json::from_str(include_str!("data/mixture_fit_r462.json")).unwrap();
    assert_eq!(references.cases.len(), 8);
    for case in references.cases {
        assert!(case.converged);
        let log_likelihoods = numbers(&case.log_likelihoods);
        let fit = fit_mixture_proportions(
            MatrixView::new(
                case.observations,
                case.components,
                &log_likelihoods,
                MatrixOrder::RowMajor,
            )
            .unwrap(),
            Some(&numbers(&case.observation_weights)),
            &numbers(&case.dirichlet),
            None,
            MixtureFitOptions {
                likelihood_offset: 1e-6,
                ..Default::default()
            },
        )
        .unwrap();
        if case.name == "shifted" {
            let column_major = (0..case.components)
                .flat_map(|j| {
                    (0..case.observations).map({
                        let data = &log_likelihoods;
                        move |i| data[i * case.components + j]
                    })
                })
                .collect::<Vec<_>>();
            let other = fit_mixture_proportions(
                MatrixView::new(
                    case.observations,
                    case.components,
                    &column_major,
                    MatrixOrder::ColumnMajor,
                )
                .unwrap(),
                Some(&numbers(&case.observation_weights)),
                &numbers(&case.dirichlet),
                None,
                MixtureFitOptions {
                    likelihood_offset: 1e-6,
                    ..Default::default()
                },
            )
            .unwrap();
            assert_eq!(fit.weights, other.weights);
            assert_eq!(fit.log_likelihood.to_bits(), other.log_likelihood.to_bits());
            assert_eq!(fit.dual_residual.to_bits(), other.dual_residual.to_bits());
        }
        let reference = numbers(&case.weighted_log_likelihood)[0];
        let prior = NormalMixture::new(
            fit.weights.clone(),
            numbers(&case.means),
            numbers(&case.sds),
        )
        .unwrap();
        let expected_means = numbers(&case.posterior_mean);
        let expected_sds = numbers(&case.posterior_sd);
        let mut mean_error = 0.0_f64;
        let mut sd_error = 0.0_f64;
        for (i, (&estimate, se)) in numbers(&case.estimates)
            .iter()
            .zip(numbers(&case.standard_errors))
            .enumerate()
        {
            let posterior = prior.posterior(estimate, se).unwrap();
            mean_error = mean_error.max((posterior.mean() - expected_means[i]).abs());
            sd_error = sd_error.max((posterior.standard_deviation() - expected_sds[i]).abs());
        }
        eprintln!(
            "{}: converged {}, grad {:.3e}, ll_error {:.3e}, mean {:.3e}, sd {:.3e}, calls {}: {}",
            case.name,
            fit.converged,
            fit.dual_residual,
            (fit.log_likelihood - reference).abs(),
            mean_error,
            sd_error,
            fit.function_evaluations,
            fit.message
        );
        assert!(fit.converged, "{}: {}", case.name, fit.message);
        assert!((fit.weights.iter().sum::<f64>() - 1.0).abs() < 1e-14);
        assert!(
            (fit.log_likelihood - reference).abs() < 1e-5,
            "{} likelihood",
            case.name
        );
        assert!(
            mean_error <= 2e-6,
            "{} posterior mean {mean_error}",
            case.name
        );
        assert!(sd_error <= 2e-6, "{} posterior SD {sd_error}", case.name);
    }
}
#[test]
fn categorical_weighted_dirichlet_and_prior_only_modes_have_known_solutions() {
    let logs = [0.0, f64::NEG_INFINITY, f64::NEG_INFINITY, 0.0];
    let fit = fit_mixture_proportions(
        MatrixView::new(2, 2, &logs, MatrixOrder::RowMajor).unwrap(),
        Some(&[3., 1.]),
        &[2., 2.],
        None,
        MixtureFitOptions::default(),
    )
    .unwrap();
    assert!(fit.converged);
    assert!((fit.weights[0] - 2. / 3.).abs() < 1e-7);
    assert!((fit.weights[1] - 1. / 3.).abs() < 1e-7);
    let expected = 4. * (2.0_f64 / 3.).ln() + 2. * (1.0_f64 / 3.).ln();
    assert!((fit.penalized_log_likelihood - expected).abs() < 1e-12);
    let fit = fit_mixture_proportions(
        MatrixView::new(0, 2, &[], MatrixOrder::ColumnMajor).unwrap(),
        None,
        &[2., 4.],
        None,
        MixtureFitOptions::default(),
    )
    .unwrap();
    assert!(fit.converged);
    assert!((fit.weights[0] - 0.25).abs() < 1e-7);
    assert_eq!(fit.log_likelihood, 0.);
}
#[test]
fn invalid_mixture_problems_fail_and_iteration_limits_remain_visible() {
    fn view(v: &[f64]) -> MatrixView<'_, f64> {
        MatrixView::new(2, 2, v, MatrixOrder::RowMajor).unwrap()
    }
    let good = [0., -3., -4., 0.];
    for (a, p, w) in [
        (Some(vec![1.]), vec![1., 1.], None),
        (Some(vec![-1., 1.]), vec![1., 1.], None),
        (Some(vec![0., 0.]), vec![1., 1.], None),
        (None, vec![1.], None),
        (None, vec![0.5, 1.], None),
        (None, vec![1., f64::NAN], None),
        (None, vec![1., 1.], Some(vec![0., 0.])),
        (None, vec![1., 1.], Some(vec![1.])),
    ] {
        assert!(
            fit_mixture_proportions(
                view(&good),
                a.as_deref(),
                &p,
                w.as_deref(),
                MixtureFitOptions::default()
            )
            .is_err()
        );
    }
    for v in [
        [0., f64::NAN, 0., 0.],
        [0., f64::INFINITY, 0., 0.],
        [f64::NEG_INFINITY, f64::NEG_INFINITY, 0., 0.],
    ] {
        assert!(
            fit_mixture_proportions(
                view(&v),
                None,
                &[1., 1.],
                None,
                MixtureFitOptions::default()
            )
            .is_err()
        );
    }
    let fit = fit_mixture_proportions(
        view(&good),
        Some(&[7., 1.]),
        &[1., 1.],
        None,
        MixtureFitOptions {
            max_iterations: 1,
            em_iterations: 0,
            ..Default::default()
        },
    )
    .unwrap();
    assert!(!fit.converged);
    let skipped = [f64::NEG_INFINITY, f64::NEG_INFINITY, 0., -2.];
    let fit = fit_mixture_proportions(
        view(&skipped),
        Some(&[0., 1.]),
        &[1., 1.],
        None,
        MixtureFitOptions::default(),
    )
    .unwrap();
    assert!(fit.converged);
    assert_eq!(fit.weights, vec![1., 0.]);
}
