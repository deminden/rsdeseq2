//! `line_search_dispersion` / `grid_search_dispersion` with a
//! `DispersionFitInput` against the twelve deprecated wrappers they replace.

use rsdeseq2::prelude::*;

fn shown<T: std::fmt::Debug>(value: Result<T, DeseqError>) -> String {
    format!("{value:?}")
}

#[test]
fn dispersion_fit_input_matches_the_deprecated_wrappers() {
    let design = DesignMatrix::from_row_major(
        6,
        2,
        vec![1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        None,
    )
    .unwrap();
    let genes: [[u32; 6]; 3] = [
        [3, 9, 4, 20, 31, 18],
        [100, 140, 90, 60, 70, 55],
        [0, 1, 0, 5, 2, 7],
    ];
    let mus: [[f64; 6]; 3] = [
        [5.2, 5.2, 5.2, 23.0, 23.0, 23.0],
        [110.0, 110.0, 110.0, 62.0, 62.0, 62.0],
        [0.4, 0.4, 0.4, 4.6, 4.6, 4.6],
    ];
    let weights = [1.0, 0.8, 1.0, 0.6, 1.0, 0.9];
    let prior = DispersionPrior::new(0.1f64.ln(), 0.5).unwrap();
    let options = GeneWiseDispersionOptions::default();
    for (y, mu) in genes.iter().zip(&mus) {
        let base = DispersionFitInput::new(y, mu, 0.3, options, 6);
        let weighted = WeightedDispersionFitInput {
            counts: y,
            mu,
            design: &design,
            initial_dispersion: 0.3,
            options,
            n_samples: 6,
            prior,
            weights: &weights,
        };
        let line = |input| shown(line_search_dispersion(input));
        let grid = |input| shown(grid_search_dispersion(input));
        let pairs = [
            (
                line(base.design(&design)),
                shown(fit_dispersion_line_search(y, mu, &design, 0.3, options, 6)),
            ),
            (
                line(base.design(&design).prior(prior)),
                shown(fit_dispersion_line_search_with_prior(
                    y, mu, &design, 0.3, options, 6, prior,
                )),
            ),
            (
                line(base.design(&design).prior(prior).weights(&weights)),
                shown(fit_dispersion_line_search_with_prior_and_weights(weighted)),
            ),
            (
                line(base),
                shown(fit_dispersion_line_search_no_cr(y, mu, 0.3, options, 6)),
            ),
            (
                line(base.prior(prior)),
                shown(fit_dispersion_line_search_no_cr_with_prior(
                    y, mu, 0.3, options, 6, prior,
                )),
            ),
            (
                line(base.prior(prior).weights(&weights)),
                shown(fit_dispersion_line_search_no_cr_with_prior_and_weights(
                    y, mu, 0.3, options, 6, prior, &weights,
                )),
            ),
            (
                grid(base.design(&design)),
                shown(fit_dispersion_grid(y, mu, &design, 0.3, options, 6)),
            ),
            (
                grid(base.design(&design).prior(prior)),
                shown(fit_dispersion_grid_with_prior(
                    y, mu, &design, 0.3, options, 6, prior,
                )),
            ),
            (
                grid(base.design(&design).prior(prior).weights(&weights)),
                shown(fit_dispersion_grid_with_prior_and_weights(weighted)),
            ),
            (
                grid(base),
                shown(fit_dispersion_grid_no_cr(y, mu, 0.3, options, 6)),
            ),
            (
                grid(base.prior(prior)),
                shown(fit_dispersion_grid_no_cr_with_prior(
                    y, mu, 0.3, options, 6, prior,
                )),
            ),
            (
                grid(base.prior(prior).weights(&weights)),
                shown(fit_dispersion_grid_no_cr_with_prior_and_weights(
                    y, mu, 0.3, options, 6, prior, &weights,
                )),
            ),
        ];
        for (index, (new, old)) in pairs.iter().enumerate() {
            assert_eq!(new, old, "pair {index}");
            assert!(new.starts_with("Ok("), "pair {index}: {new}");
        }
    }
}
