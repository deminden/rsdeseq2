//! `line_search_dispersion` / `grid_search_dispersion` with every
//! `DispersionFitInput` combination of design, prior and weights.

use rsdeseq2::prelude::*;

fn shown<T: std::fmt::Debug>(value: Result<T, DeseqError>) -> String {
    format!("{value:?}")
}

#[test]
fn dispersion_fit_input_combinations_all_fit() {
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
        let inputs = [
            base.design(&design),
            base.design(&design).prior(prior),
            base.design(&design).prior(prior).weights(&weights),
            base,
            base.prior(prior),
            base.prior(prior).weights(&weights),
        ];
        let fits = inputs
            .iter()
            .map(|&input| shown(line_search_dispersion(input)))
            .chain(
                inputs
                    .iter()
                    .map(|&input| shown(grid_search_dispersion(input))),
            );
        for (index, fit) in fits.enumerate() {
            assert!(fit.starts_with("Ok("), "fit {index}: {fit}");
        }
    }
}
