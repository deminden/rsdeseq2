use super::*;
use crate::support::run::run_results;
use rsdeseq2::prelude::{Analysis, ResultColumn};

#[test]
fn supplied_fixed_dispersions_from_optional_reference_are_preserved() {
    let dispersions = read_fixed_dispersions();
    let size_factors = read_size_factors("size_factors_ratio.tsv");

    let (fit, results) = run_results(
        DeseqBuilder::new()
            .size_factors(size_factors)
            .disable_cooks_cutoff()
            .disable_independent_filtering()
            .run(
                &Analysis::new(&reference_counts(), &reference_full_design())
                    .fixed_dispersion(&dispersions)
                    .result(ResultColumn::Coefficient(1)),
            ),
    )
    .unwrap();

    let fit_dispersions = fit.dispersion.as_ref().unwrap();
    assert_eq!(fit_dispersions.len(), dispersions.len());
    for (gene, (actual, expected)) in fit_dispersions
        .iter()
        .copied()
        .zip(dispersions.iter().copied())
        .enumerate()
    {
        assert_float_close(
            actual,
            expected,
            1e-15,
            1e-15,
            &format!("supplied dispersion gene {gene}"),
        );
        assert_eq!(results.rows[gene].dispersion, Some(expected));
    }
}
